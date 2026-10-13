"""
foldbot: Telegram protein-structure agent.

- Selectable prediction method per user (/method): ESMFold remote API, ESMFold local, AlphaFold DB lookup
- pLDDT-colored PNG render, mutation effect report (structure comparison + residue/burial analysis
  + plain-language explanation from your own LLM: Ollama or Gemini)
- LLM-guided design loop, per-method job queues, disk cache, user allow-list
- Conversational mode: plain messages ("fold this", "what if L43 became D?") are interpreted by the LLM
  and routed to the same tools as the slash commands (which still work)
"""
import asyncio
import csv
import functools
import hashlib
import importlib.util
import io
import json
import logging
import os
import re
import shlex
import shutil
import tempfile
import time
from pathlib import Path

import httpx
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from collections import Counter
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBIO, PDBParser, Select, Superimposer
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from telegram import InlineKeyboardButton, InlineKeyboardMarkup, Update
from telegram.constants import ChatAction
from telegram.ext import (
    ApplicationBuilder, CallbackQueryHandler, CommandHandler, ContextTypes,
    MessageHandler, filters,
)
from telegram.request import HTTPXRequest

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("foldbot")

# ---------------------------------------------------------------- .env loader
def _load_dotenv():
    """Read KEY=VALUE lines from .env next to bot.py (or in the current folder). Real environment variables win."""
    vals = {}
    for f in (Path(__file__).resolve().parent / ".env", Path.cwd() / ".env"):
        if not f.exists():
            continue
        for line in f.read_text(encoding="utf-8-sig").splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            k, v = k.strip(), v.strip().strip('"').strip("'")
            if v:
                vals[k] = v  # a later line overrides an earlier one
    for k, v in vals.items():
        os.environ.setdefault(k, v)


_load_dotenv()

# ---------------------------------------------------------------- config
TOKEN = os.environ["TELEGRAM_BOT_TOKEN"]
ALLOWED = {int(x) for x in os.getenv("ALLOWED_USER_IDS", "").split(",") if x.strip()}
DEFAULT_METHOD = "esmfold_local" if os.getenv("FOLD_BACKEND", "api").lower() == "local" else "esmfold_api"
MAX_LEN_API = int(os.getenv("MAX_LEN_API", "400"))
MAX_LEN_LOCAL = int(os.getenv("MAX_LEN_LOCAL", "700"))
LLM_PROVIDER = os.getenv("LLM_PROVIDER", "ollama").lower()  # "ollama", "gemini" or "openai" (any OpenAI-compatible API)
OLLAMA_URL = os.getenv("OLLAMA_URL", "http://localhost:11434")
LLM_MODEL = os.getenv("LLM_MODEL", "qwen2.5:7b-instruct")
GEMINI_KEY = os.getenv("GEMINI_API_KEY", "")
GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-flash-latest")  # alias for the newest stable Flash model
OPENAI_BASE_URL = os.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1").rstrip("/")
OPENAI_KEY = os.getenv("OPENAI_API_KEY", "")
OPENAI_MODEL = os.getenv("OPENAI_MODEL", "")
GROQ_KEY = os.getenv("GROQ_API_KEY", "")
GROQ_MODEL = os.getenv("GROQ_MODEL", "qwen/qwen3.8-27b")
RENDER_STYLE = os.getenv("RENDER_STYLE", "cartoon").lower()  # "cartoon" or "trace" (simple CA line)
DESIGN_ROUNDS = int(os.getenv("DESIGN_ROUNDS", "2"))
DESIGN_VARIANTS = int(os.getenv("DESIGN_VARIANTS", "4"))
CACHE = Path(os.getenv("CACHE_DIR", ".cache"))
CACHE.mkdir(exist_ok=True)

# Extra prediction methods (see .env.example). CLI tools are found on PATH or via *_BIN.
BOLTZ_BIN = os.getenv("BOLTZ_BIN", "boltz")
BOLTZ_MSA_SERVER = os.getenv("BOLTZ_USE_MSA_SERVER", "1") == "1"  # 1 = sends sequence to api.colabfold.com
COLABFOLD_BIN = os.getenv("COLABFOLD_BIN", "colabfold_batch")
COLABFOLD_MSA_MODE = os.getenv("COLABFOLD_MSA_MODE", "mmseqs2_uniref_env")  # or single_sequence (offline, worse)
OMEGAFOLD_BIN = os.getenv("OMEGAFOLD_BIN", "omegafold")
NVIDIA_KEY = os.getenv("NVIDIA_API_KEY", "")
NVIDIA_ESMFOLD_URL = os.getenv("NVIDIA_ESMFOLD_URL", "https://health.api.nvidia.com/v1/biology/nvidia/esmfold")
CLI_TIMEOUT = int(os.getenv("CLI_TIMEOUT", "1800"))

# AlphaFold 3: you must request the weights from Google DeepMind yourself and install the databases (~630 GB).
# AF3_CMD is a command template; {json} and {out} are filled in by the bot. Docker is fine, e.g.
#   python /opt/alphafold3/run_alphafold.py --json_path={json} --model_dir=/data/af3/models --db_dir=/data/af3/db --output_dir={out}
AF3_CMD = os.getenv("AF3_CMD", "")
AF3_TIMEOUT = int(os.getenv("AF3_TIMEOUT", "7200"))

# Stability (ddG) tools. FoldX is licensed software you download yourself; ThermoMPNN is an open-source repo.
FOLDX_BIN = os.getenv("FOLDX_BIN", "foldx")
FOLDX_DIR = os.getenv("FOLDX_DIR", "")            # folder containing rotabase.txt (defaults to the binary's folder)
FOLDX_RUNS = int(os.getenv("FOLDX_RUNS", "3"))
THERMOMPNN_DIR = os.getenv("THERMOMPNN_DIR", "")  # checkout of github.com/Kuhlman-Lab/ThermoMPNN
THERMOMPNN_CMD = os.getenv("THERMOMPNN_CMD", "")  # e.g. "/path/env/bin/python custom_inference.py --pdb {pdb} --chain {chain}"

ESM_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{}"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
UNIPROT_RE = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")


# ---------------------------------------------------------------- access control
def restricted(fn):
    @functools.wraps(fn)
    async def wrapper(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
        uid = update.effective_user.id
        if uid not in ALLOWED:
            await update.message.reply_text(
                f"Not authorised. Your Telegram ID is {uid}. Add it to ALLOWED_USER_IDS to use this bot."
            )
            return
        return await fn(update, ctx)
    return wrapper


# ---------------------------------------------------------------- sequence helpers
def clean_sequence(text: str) -> str:
    lines = [l.strip() for l in text.splitlines() if l.strip() and not l.startswith(">")]
    return "".join(lines).upper().replace(" ", "")


def validate(seq: str, method: str) -> str | None:
    if not seq:
        return "Empty sequence."
    bad = set(seq) - VALID_AA
    if bad:
        return f"Invalid characters: {', '.join(sorted(bad))}"
    limit = METHODS[method]["max_len"]
    if len(seq) > limit:
        return f"Sequence too long ({len(seq)} aa). Max for {method} is {limit}."
    return None


MUT_RE = re.compile(r"([A-Z])(\d+)([A-Z])")


def parse_mutations(args: list[str]) -> list[str]:
    """Accept 'L43D,K11R', 'L43D, K11R' or 'L43D K11R'."""
    return [m.upper() for m in re.split(r"[,\s;]+", " ".join(args)) if m]


def apply_mutations(seq: str, muts: list[str]) -> str:
    """Mutations like A45G (wild-type residue, 1-based position, new residue)."""
    s = list(seq)
    for m in muts:
        g = MUT_RE.fullmatch(m.strip().upper())
        if not g:
            raise ValueError(f"Bad mutation format: {m}")
        wt, pos, new = g.group(1), int(g.group(2)), g.group(3)
        if not 1 <= pos <= len(s):
            raise ValueError(f"{m}: position out of range (1-{len(s)})")
        if s[pos - 1] != wt:
            raise ValueError(f"{m}: residue {pos} is {s[pos - 1]}, not {wt}")
        if new not in VALID_AA:
            raise ValueError(f"{m}: invalid target residue")
        s[pos - 1] = new
    return "".join(s)


# ---------------------------------------------------------------- structure helpers
def residues(pdb: str):
    """Return CA atoms and per-residue pLDDT (0-100) from the B-factor column."""
    s = PDBParser(QUIET=True).get_structure("s", io.StringIO(pdb))
    cas = [r["CA"] for r in s.get_residues() if "CA" in r]
    pl = np.array([a.get_bfactor() for a in cas], dtype=float)
    if pl.size and pl.max() <= 1.0:
        pl = pl * 100
    return cas, pl


def mean_plddt(pdb: str) -> float:
    _, pl = residues(pdb)
    return float(pl.mean()) if pl.size else float("nan")


def confidence_label(p: float) -> str:
    return "very high" if p >= 90 else "confident" if p >= 70 else "low" if p >= 50 else "very low"


def plddt_color(p: float) -> str:  # AlphaFold DB palette
    return "#0053D6" if p >= 90 else "#65CBF3" if p >= 70 else "#FFDB13" if p >= 50 else "#FF7D45"


def render_trace_png(pdb: str, title: str) -> io.BytesIO:
    """Simple CA-trace render colored by pLDDT."""
    cas, pl = residues(pdb)
    xyz = np.array([a.coord for a in cas])
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")
    for i in range(len(xyz) - 1):
        ax.plot(*xyz[i:i + 2].T, color=plddt_color(pl[i]), lw=2.5)
    ax.set_axis_off()
    ax.set_title(title)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    buf.name = "render.png"
    return buf


# ---------------------------------------------------------------- cartoon renderer (ribbons: helices, strands, loops)
def backbone(pdb: str):
    """Residues with full backbone (N, CA, C, O) from the first model: coords (L,4,3), pLDDT (L,), connected (L,) flags."""
    s = PDBParser(QUIET=True).get_structure("s", io.StringIO(pdb))
    coords, pl, chain_of, resname = [], [], [], []
    for ch in s[0]:
        for r in ch:
            if r.id[0] != " " or not all(a in r for a in ("N", "CA", "C", "O")):
                continue
            coords.append([r[a].coord for a in ("N", "CA", "C", "O")])
            pl.append(r["CA"].get_bfactor())
            chain_of.append(ch.id)
            resname.append(r.get_resname())
    coords, pl = np.array(coords, dtype=float), np.array(pl, dtype=float)
    if pl.size and pl.max() <= 1.0:
        pl = pl * 100
    n = len(coords)
    connected = np.zeros(n, bool)  # residue i is bonded to residue i-1
    for i in range(1, n):
        connected[i] = chain_of[i] == chain_of[i - 1] and np.linalg.norm(coords[i - 1, 2] - coords[i, 0]) < 2.0
    return coords, pl, connected, resname


def assign_ss(coords, connected, resname):
    """Simplified DSSP (Kabsch-Sander hydrogen-bond energies): returns 'H' (helix), 'E' (strand) or '-' per residue."""
    n = len(coords)
    N, C, O = coords[:, 0], coords[:, 2], coords[:, 3]
    H = np.full((n, 3), np.nan)
    for i in range(1, n):
        if connected[i] and resname[i] != "PRO":
            v = C[i - 1] - O[i - 1]
            H[i] = N[i] + v / np.linalg.norm(v)
    # E[i, j]: energy of the bond between C=O of residue i (acceptor) and N-H of residue j (donor)
    d = lambda a, b: np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        E = 27.888 * (1 / d(O, N) + 1 / d(C, H) - 1 / d(O, H) - 1 / d(C, N))
    A = np.nan_to_num(E, nan=0.0) < -0.5
    idx = np.arange(n)
    A &= np.abs(idx[:, None] - idx[None, :]) >= 2
    ss = np.array(["-"] * n)

    def shift(M, di, dj):
        out = np.zeros_like(M)
        i0, i1, j0, j1 = max(0, -di), min(n, n - di), max(0, -dj), min(n, n - dj)
        out[i0:i1, j0:j1] = M[i0 + di:i1 + di, j0 + dj:j1 + dj]
        return out

    # helices: consecutive n-turns (n=4 alpha, n=3 3-10)
    for turn, length in ((3, 3), (4, 4)):  # alpha (4) is applied last so it takes priority
        t = np.array([A[i, i + turn] if i + turn < n else False for i in range(n)])
        for i in range(1, n):
            if t[i - 1] and t[i]:
                ss[i:i + length] = "H"
    # strands: residues in at least two neighbouring beta bridges
    AT = A.T
    par = (shift(A, -1, 0) & shift(AT, 1, 0)) | (shift(AT, 0, -1) & shift(A, 0, 1))
    anti = (A & AT) | (shift(A, -1, 1) & shift(AT, 1, -1))
    br = (par | anti) & (np.abs(idx[:, None] - idx[None, :]) >= 3)
    inb = br.any(axis=1)
    for i in range(n):
        if inb[i] and ((i > 0 and inb[i - 1]) or (i < n - 1 and inb[i + 1])) and ss[i] == "-":
            ss[i] = "E"
    return ss


def _catmull(P, sub):
    """Catmull-Rom spline through the rows of P with `sub` samples per segment (plus the last point)."""
    L = len(P)
    ext = np.vstack([P[0], P, P[-1]])
    out = []
    for i in range(L - 1):
        p0, p1, p2, p3 = ext[i], ext[i + 1], ext[i + 2], ext[i + 3]
        for s in range(sub):
            t = s / sub
            out.append(0.5 * (2 * p1 + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t ** 2
                              + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    out.append(P[-1])
    return np.array(out)


def _color(p):
    return to_rgb("#0053D6" if p >= 90 else "#65CBF3" if p >= 70 else "#FFDB13" if p >= 50 else "#FF7D45")


def render_cartoon_png(pdb, title, elev=55, azim=-70):
    coords, pl, connected, resname = backbone(pdb)
    n = len(coords)
    if n < 5:
        raise ValueError("too few residues with a full backbone")
    ss = assign_ss(coords, connected, resname)
    sub = 10 if n < 300 else 6 if n < 700 else 4

    # PCA orientation so the protein's longest dimension is horizontal
    ca = coords[:, 1]
    mean = ca.mean(axis=0)
    _, _, Vt = np.linalg.svd(ca - mean, full_matrices=False)
    if np.linalg.det(Vt) < 0:
        Vt[2] *= -1
    xyz = (coords - mean) @ Vt.T  # (n, 4, 3)

    # light direction (camera-relative: upper left, towards the viewer)
    e, a = np.radians(elev), np.radians(azim)
    cam = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    up = np.array([-np.sin(e) * np.cos(a), -np.sin(e) * np.sin(a), np.cos(e)])
    left = np.cross(up, cam)
    light = cam * 0.7 + up * 0.5 + left * 0.35
    light /= np.linalg.norm(light)

    polys, colors = [], []
    starts = [0] + [i for i in range(1, n) if not connected[i]] + [n]
    for f0, f1 in zip(starts[:-1], starts[1:]):
        m = f1 - f0
        if m < 2:
            continue
        P = xyz[f0:f1, 1]
        V = xyz[f0:f1, 3] - xyz[f0:f1, 2]
        V /= np.linalg.norm(V, axis=1, keepdims=True)
        for i in range(1, m):                       # keep carbonyl directions consistent along the chain
            if np.dot(V[i], V[i - 1]) < 0:
                V[i] = -V[i]
        S, W = _catmull(P, sub), _catmull(V, sub)
        T = np.gradient(S, axis=0)
        T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-9
        W = W - (W * T).sum(axis=1, keepdims=True) * T
        W /= np.linalg.norm(W, axis=1, keepdims=True) + 1e-9
        Nrm = np.cross(T, W)

        res_of = np.array([min(i + (1 if s >= sub // 2 else 0), m - 1) for i in range(m - 1) for s in range(sub)] + [m - 1])
        kind = np.array([ss[f0 + r] for r in res_of])
        hw = np.where(kind == "H", 1.15, np.where(kind == "E", 0.95, 0.2))
        ht = np.where(kind == "H", 0.3, np.where(kind == "E", 0.25, 0.2))
        k = max(3, sub)                              # smooth width changes between secondary-structure types
        pad = np.pad(hw, (k // 2, k // 2), mode="edge")
        hw = np.convolve(pad, np.ones(k) / k, mode="valid")[:len(hw)]
        for j in range(m):                           # arrow head on the last residue of each strand
            if ss[f0 + j] == "E" and (j == m - 1 or ss[f0 + j + 1] != "E") and j >= 1:
                a0 = (j - 1) * sub
                base = 0.95
                for s in range(sub):
                    hw[a0 + s] = base * (1.75 - 1.55 * s / sub)
        col = np.array([_color(pl[f0 + r]) for r in res_of])

        corners = np.stack([S + W * hw[:, None] + Nrm * ht[:, None],
                            S + W * hw[:, None] - Nrm * ht[:, None],
                            S - W * hw[:, None] - Nrm * ht[:, None],
                            S - W * hw[:, None] + Nrm * ht[:, None]], axis=1)  # (K, 4, 3)
        for kk in range(len(S) - 1):
            for a_, b_ in ((0, 1), (1, 2), (2, 3), (3, 0)):
                quad = np.array([corners[kk, a_], corners[kk, b_], corners[kk + 1, b_], corners[kk + 1, a_]])
                nv = np.cross(quad[1] - quad[0], quad[3] - quad[0])
                nv /= np.linalg.norm(nv) + 1e-9
                shade = 0.5 + 0.5 * abs(float(np.dot(nv, light)))
                polys.append(quad)
                colors.append(np.clip(np.array(col[kk]) * shade, 0, 1))

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Poly3DCollection(polys, facecolors=colors, edgecolors=colors, linewidths=0.3))
    allp = np.concatenate(polys)
    c, r = (allp.max(0) + allp.min(0)) / 2, (allp.max(0) - allp.min(0)).max() / 2 * 0.74
    ax.set_xlim(c[0] - r, c[0] + r); ax.set_ylim(c[1] - r, c[1] + r); ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    ax.set_title(title)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def render_png(pdb: str, title: str) -> io.BytesIO:
    """Cartoon picture colored by pLDDT (falls back to a simple CA trace if the cartoon cannot be drawn)."""
    if RENDER_STYLE == "cartoon":
        try:
            buf = render_cartoon_png(pdb, title)
            buf.name, buf.style = "render.png", "cartoon"
            return buf
        except ValueError as e:  # e.g. a file without full backbone atoms
            log.warning("cartoon not possible (%s), using CA trace", e)
        except Exception:  # noqa: BLE001
            log.exception("cartoon render failed, falling back to CA trace")
    buf = render_trace_png(pdb, title)
    buf.name, buf.style = "render.png", "trace"
    return buf


def compare(wt_pdb: str, mut_pdb: str, site_positions: list[int]):
    a, pa = residues(wt_pdb)
    b, pb = residues(mut_pdb)
    sup = Superimposer()
    sup.set_atoms(a, b)
    sup.apply(b)
    disp = np.array([np.linalg.norm(x.coord - y.coord) for x, y in zip(a, b)])
    delta = pb - pa

    fig, axs = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
    x = np.arange(1, len(pa) + 1)
    axs[0].plot(x, pa, label="wild type")
    axs[0].plot(x, pb, label="mutant")
    axs[0].set_ylabel("pLDDT")
    axs[0].legend()
    axs[1].plot(x, disp, color="tab:red")
    axs[1].set_ylabel("CA shift (Å)")
    axs[1].set_xlabel("residue")
    for ax in axs:
        for p in site_positions:
            ax.axvline(p, color="gray", ls="--", lw=0.8)
    buf = io.BytesIO()
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    buf.seek(0)
    buf.name = "mutation_compare.png"

    lines = [
        "• Method: Kabsch 3D coordinate superposition (BioPython Superimposer)",
        f"• Cα RMSD (after superposition): {sup.rms:.2f} Å",
        f"• Mean pLDDT (model confidence): WT {pa.mean():.1f} → mutant {pb.mean():.1f} ({pb.mean() - pa.mean():+.1f})",
    ]
    for p in site_positions:
        lo, hi = max(0, p - 6), min(len(pa), p + 5)
        lines.append(f"• Local pLDDT around pos {p} (±5 aa): {delta[lo:hi].mean():+.1f}")
        lines.append(f"• Cα backbone shift at pos {p}: {disp[p - 1]:.2f} Å")
    top = np.argsort(disp)[::-1][:5]
    lines.append("• Largest Cα backbone shifts: " + ", ".join(f"pos {i + 1} ({disp[i]:.1f} Å)" for i in top))
    return "\n".join(lines), buf


# ---------------------------------------------------------------- mutation analysis
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)   # Kyte-Doolittle
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1, H=153.2, I=166.7,
           L=166.7, K=168.6, M=162.9, F=189.9, P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)  # Å^3
CHG = {**{a: 0 for a in VALID_AA}, "D": -1, "E": -1, "K": 1, "R": 1}


def site_analysis(seq: str, muts: list[str], wt_pdb: str) -> str:
    """Rule-based description of each substitution using the wild-type predicted structure."""
    cas, pl = residues(wt_pdb)
    xyz = np.array([a.coord for a in cas])
    dist = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    nb = (dist < 10).sum(axis=1) - 1  # CA neighbours within 10 Å = crude burial proxy
    out = [
        "Method: 10 Å Cα coordination sphere (burial) | Kyte-Doolittle scale (hydropathy) | Chothia scale (van der Waals volume)",
    ]
    for m in muts:
        g = MUT_RE.fullmatch(m.strip().upper())
        wt, pos, new = g.group(1), int(g.group(2)), g.group(3)
        i = pos - 1
        pct = float((nb < nb[i]).mean())
        burial = "buried core" if pct >= 0.67 else "exposed surface" if pct <= 0.33 else "partially buried"
        dh, dv, dq = KD[new] - KD[wt], VOL[new] - VOL[wt], CHG[new] - CHG[wt]
        notes, risk = [], 0
        if pl[i] < 50:
            notes.append("site is in a low-confidence (possibly disordered) region")
        if new == "P":
            notes.append("proline restricts backbone angles & breaks secondary structure"); risk += 1
        if wt == "P":
            notes.append("loses proline's backbone rigidity")
        if wt == "G":
            notes.append("glycine may be needed for tight turns/flexibility"); risk += 1
        if new == "G":
            notes.append("glycine adds backbone flexibility")
        if "C" in (wt, new):
            notes.append("cysteine change affects disulfide bonds/redox sensitivity")
        if burial == "buried core":
            if dh <= -3:
                notes.append("polar/charged residue replacing hydrophobic core is strongly destabilising"); risk += 2
            if dv >= 30:
                notes.append("bulkier side chain in packed core causes steric clashes"); risk += 1
            if dv <= -30:
                notes.append("smaller side chain leaves a cavity in the core"); risk += 1
            if dq != 0:
                notes.append("uncompensated charge in buried core is energetically costly"); risk += 2
        elif burial == "exposed surface":
            if dq != 0:
                notes.append("surface charge change alters local electrostatics/solubility")
            if dh >= 3:
                notes.append("new surface hydrophobic patch may promote aggregation")
            if not notes:
                notes.append("conservative surface change, typically well tolerated")
        if not notes:
            notes.append("no specific structural risk flags raised")
        verdict = "HIGH" if risk >= 3 else "MODERATE" if risk >= 1 else "LOW"
        out.append(
            f"• {m}: {burial} (Cα neighbors: {nb[i]}, percentile: {pct*100:.0f}%, local pLDDT: {pl[i]:.0f}/100)\n"
            f"  - Hydropathy Δ: {dh:+.1f} | Volume Δ: {dv:+.0f} Å³ | Net charge Δ: {dq:+d}\n"
            f"  - Destabilisation risk: {verdict} ({'; '.join(notes)})."
        )
    return "\n".join(out)


async def llm_describe(facts: str) -> str | None:
    prompt = (
        "You are a structural biology assistant. Below are computed facts about point mutations in a protein "
        "(from PREDICTED structures, not experiments). Write a plain-language explanation in at most 150 words: "
        "what changes chemically, whether the site is buried or exposed, what predicted structural change was seen, "
        "likely stability/function implications, and how uncertain this is. Do not invent experimental data, "
        "database entries or citations.\n\nFacts:\n" + facts
    )
    try:
        return (await llm_text(prompt, as_json=False)).strip()
    except Exception:  # noqa: BLE001
        log.exception("LLM description failed")
        return None


# ---------------------------------------------------------------- compute backends
async def _predict_api(seq: str) -> str:
    async with httpx.AsyncClient(timeout=300) as client:
        last_exc = None
        for attempt in range(3):
            try:
                r = await client.post(ESM_URL, content=seq, headers={"Content-Type": "text/plain"})
                if r.status_code in (502, 503, 504) and attempt < 2:
                    wait_s = (attempt + 1) * 3
                    log.warning("ESM Atlas returned %d, retrying in %ds (attempt %d/3)...", r.status_code, wait_s, attempt + 1)
                    await asyncio.sleep(wait_s)
                    continue
                r.raise_for_status()
                return r.text
            except httpx.HTTPStatusError as e:
                last_exc = e
                if e.response.status_code in (502, 503, 504) and attempt < 2:
                    await asyncio.sleep((attempt + 1) * 3)
                    continue
                raise
            except (httpx.TimeoutException, httpx.NetworkError) as e:
                last_exc = e
                if attempt < 2:
                    await asyncio.sleep((attempt + 1) * 3)
                    continue
                raise
        if last_exc:
            raise last_exc
        raise RuntimeError("ESM Atlas request failed after retries.")


_local: dict = {}


def _load_local():
    if "model" in _local:
        return
    import torch
    from transformers import EsmForProteinFolding
    log.info("Loading ESMFold weights (first run downloads several GB)...")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(dev).eval()
    if dev == "cuda":
        model.esm = model.esm.half()
    model.trunk.set_chunk_size(64)
    _local.update(model=model, device=dev)
    log.info("ESMFold ready on %s", dev)


def _predict_local(seq: str) -> str:
    import torch
    _load_local()
    with torch.no_grad():
        return _local["model"].infer_pdb(seq)


async def fetch_afdb(acc: str):
    async with httpx.AsyncClient(timeout=60, follow_redirects=True) as client:
        meta = await client.get(AFDB_API.format(acc))
        meta.raise_for_status()
        entry = meta.json()[0]
        pdb = (await client.get(entry["pdbUrl"])).text
    return pdb, entry.get("uniprotDescription", acc)


def split_cmd(cmd: str, windows: bool | None = None) -> list[str]:
    """Split a command template into arguments; keeps Windows backslashes intact."""
    win = (os.name == "nt") if windows is None else windows
    if not win:
        return shlex.split(cmd)
    return [t[1:-1] if len(t) > 1 and t[0] == t[-1] and t[0] in "\"'" else t
            for t in shlex.split(cmd, posix=False)]


async def _run(cmd: list[str], cwd=None, timeout: int | None = None) -> str:
    """Run a command-line tool without blocking the bot."""
    try:
        proc = await asyncio.create_subprocess_exec(
            *cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.STDOUT, cwd=cwd
        )
    except (FileNotFoundError, NotADirectoryError, PermissionError) as e:
        raise RuntimeError(f"Could not start '{cmd[0]}' ({e.__class__.__name__}). Check the path in your .env file.")
    try:
        out, _ = await asyncio.wait_for(proc.communicate(), timeout or CLI_TIMEOUT)
    except asyncio.TimeoutError:
        proc.kill()
        raise RuntimeError(f"{Path(cmd[0]).name} timed out after {timeout or CLI_TIMEOUT}s")
    text = out.decode(errors="replace")
    if proc.returncode != 0:
        raise RuntimeError(f"{Path(cmd[0]).name} failed:\n{text[-600:]}")
    return text


async def _predict_boltz(seq: str) -> str:
    with tempfile.TemporaryDirectory() as d:
        inp = Path(d) / "query.yaml"
        yaml = f"version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: {seq}\n"
        if not BOLTZ_MSA_SERVER:
            yaml += "      msa: empty\n"  # single-sequence mode, lower accuracy
        inp.write_text(yaml)
        cmd = [BOLTZ_BIN, "predict", str(inp), "--out_dir", d, "--output_format", "pdb"]
        if BOLTZ_MSA_SERVER:
            cmd.append("--use_msa_server")
        cmd += split_cmd(os.getenv("BOLTZ_EXTRA_ARGS", ""))  # e.g. "--accelerator cpu --no_kernels"
        await _run(cmd)
        hits = sorted(Path(d).rglob("*_model_0.pdb"))
        if not hits:
            raise RuntimeError("Boltz finished but no PDB was found")
        return hits[0].read_text()


async def _predict_colabfold(seq: str) -> str:
    with tempfile.TemporaryDirectory() as d:
        fa, out = Path(d) / "query.fasta", Path(d) / "out"
        fa.write_text(f">query\n{seq}\n")
        await _run([COLABFOLD_BIN, str(fa), str(out), "--num-models", "1", "--num-recycle", "3",
                    "--msa-mode", COLABFOLD_MSA_MODE])
        hits = sorted(out.glob("*rank_001*.pdb"))
        if not hits:
            raise RuntimeError("ColabFold finished but no ranked PDB was found")
        return hits[0].read_text()


async def _predict_omegafold(seq: str) -> str:
    with tempfile.TemporaryDirectory() as d:
        fa, out = Path(d) / "query.fasta", Path(d) / "out"
        fa.write_text(f">query\n{seq}\n")
        out.mkdir()
        await _run([OMEGAFOLD_BIN, str(fa), str(out)])
        hits = sorted(out.glob("*.pdb"))
        if not hits:
            raise RuntimeError("OmegaFold finished but no PDB was found")
        return hits[0].read_text()


def cif_to_pdb(cif_text: str) -> str:
    """AlphaFold 3 writes mmCIF; convert to PDB so the rest of the bot (and Telegram users) can use it."""
    from Bio.PDB import MMCIFParser, PDBIO
    structure = MMCIFParser(QUIET=True).get_structure("s", io.StringIO(cif_text))
    out = io.StringIO()
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(out)
    return out.getvalue()


async def _predict_af3(seq: str) -> str:
    with tempfile.TemporaryDirectory() as d:
        inp, out = Path(d) / "fold_input.json", Path(d) / "out"
        out.mkdir()
        inp.write_text(json.dumps({
            "name": "query", "modelSeeds": [1], "dialect": "alphafold3", "version": 1,
            "sequences": [{"protein": {"id": "A", "sequence": seq}}],
        }))
        await _run(split_cmd(AF3_CMD.format(json=str(inp), out=str(out))), timeout=AF3_TIMEOUT)
        hits = sorted(out.rglob("*_model.cif"), key=lambda f: len(f.parts))  # top-ranked model sits at the job's root
        if not hits:
            raise RuntimeError("AlphaFold 3 finished but no *_model.cif was found")
        return cif_to_pdb(hits[0].read_text())


async def _predict_nvidia(seq: str) -> str:
    headers = {"Authorization": f"Bearer {NVIDIA_KEY}", "Accept": "application/json",
               "NVCF-POLL-SECONDS": "300"}
    async with httpx.AsyncClient(timeout=600) as c:
        r = await c.post(NVIDIA_ESMFOLD_URL, headers=headers, json={"sequence": seq})
        if r.status_code == 202:  # long job: poll until done
            poll = f"https://health.api.nvidia.com/v1/status/{r.headers['nvcf-reqid']}"
            for _ in range(120):
                await asyncio.sleep(5)
                r = await c.get(poll, headers=headers)
                if r.status_code != 202:
                    break
        r.raise_for_status()
        pdbs = r.json().get("pdbs")
    if not pdbs:
        raise RuntimeError("NVIDIA API returned no structure")
    return pdbs[0] if isinstance(pdbs, list) else pdbs


# Registry of prediction methods. To add one, write a predict function and add an entry here.
METHODS = {
    "esmfold_api": dict(family="ESMFold", label="ESMFold – remote API (no GPU needed)", fn=_predict_api,
                        max_len=MAX_LEN_API, workers=2),
    "esmfold_nvidia": dict(family="ESMFold", label="ESMFold – NVIDIA hosted (needs free API key)", fn=_predict_nvidia,
                           max_len=int(os.getenv("MAX_LEN_NVIDIA", "1000")), workers=1),
    "esmfold_local": dict(family="ESMFold", label="ESMFold – local GPU/CPU", fn=_predict_local,
                          max_len=MAX_LEN_LOCAL, workers=1),
    "boltz2_local": dict(family="Boltz-2", label="Boltz-2 – local GPU (open AlphaFold3-class)", fn=_predict_boltz,
                         max_len=int(os.getenv("MAX_LEN_BOLTZ", "800")), workers=1),
    "colabfold_local": dict(family="AlphaFold2", label="ColabFold / AlphaFold2 – local GPU", fn=_predict_colabfold,
                            max_len=int(os.getenv("MAX_LEN_COLABFOLD", "800")), workers=1),
    "omegafold_local": dict(family="OmegaFold", label="OmegaFold – local GPU (no MSA)", fn=_predict_omegafold,
                            max_len=int(os.getenv("MAX_LEN_OMEGAFOLD", "600")), workers=1),
    "alphafold3_local": dict(family="AlphaFold3", label="AlphaFold 3 – local install (weights from Google)", fn=_predict_af3,
                             max_len=int(os.getenv("MAX_LEN_AF3", "1000")), workers=1),
    "afdb": dict(family="AlphaFold DB", label="AlphaFold DB – lookup by UniProt ID", fn=None, max_len=0, workers=0),
}

CLI_BINS = {"boltz2_local": BOLTZ_BIN, "colabfold_local": COLABFOLD_BIN, "omegafold_local": OMEGAFOLD_BIN}


def method_available(name: str):
    if name == "esmfold_local":
        if not (importlib.util.find_spec("torch") and importlib.util.find_spec("transformers")):
            return False, "Local ESMFold needs torch and transformers installed on the machine running the bot."
    if name in CLI_BINS and not shutil.which(CLI_BINS[name]):
        return False, f"'{CLI_BINS[name]}' was not found on PATH. Install it (see .env.example) or set its *_BIN path."
    if name == "alphafold3_local" and not AF3_CMD:
        return False, "AlphaFold 3 isn't set up. It needs weights requested from Google DeepMind, ~630 GB of databases and an NVIDIA GPU; then set AF3_CMD in .env."
    if name == "esmfold_nvidia" and not NVIDIA_KEY:
        return False, "Set NVIDIA_API_KEY (free key from build.nvidia.com)."
    return True, ""


class JobQueue:
    """Serialises heavy jobs so a GPU (or the public API) is not overwhelmed."""

    def __init__(self):
        self.q: asyncio.Queue = asyncio.Queue()

    def start(self, workers: int):
        for _ in range(workers):
            asyncio.create_task(self._worker())

    async def _worker(self):
        while True:
            fn, arg, fut = await self.q.get()
            try:
                res = await fn(arg) if asyncio.iscoroutinefunction(fn) else await asyncio.to_thread(fn, arg)
                fut.set_result(res)
            except Exception as e:  # noqa: BLE001
                fut.set_exception(e)
            finally:
                self.q.task_done()

    async def submit(self, fn, arg):
        fut = asyncio.get_running_loop().create_future()
        await self.q.put((fn, arg, fut))
        return await fut


queues = {n: JobQueue() for n, m in METHODS.items() if m["fn"]}


async def predict(seq: str, method: str) -> str:
    key = hashlib.sha1(f"{method}:{seq}".encode()).hexdigest()
    f = CACHE / f"{key}.pdb"
    if f.exists():
        return f.read_text()
    pdb = await queues[method].submit(METHODS[method]["fn"], seq)
    f.write_text(pdb)
    return pdb


async def resolve_method(update: Update, ctx: ContextTypes.DEFAULT_TYPE, need_sequence: bool):
    """Return the user's chosen method (or a sequence-capable fallback), or None if unusable."""
    m = ctx.user_data.get("method", DEFAULT_METHOD)
    if need_sequence and m == "afdb":
        m = DEFAULT_METHOD
        await update.message.reply_text(
            f"AlphaFold DB can't fold custom or mutated sequences, so I'm using {METHODS[m]['label']} here."
        )
    ok, why = method_available(m)
    if not ok:
        await update.message.reply_text(f"❌ {why} Pick another method with /method.")
        return None
    return m


# ---------------------------------------------------------------- LLM
async def _call_openai_compat(
    client: httpx.AsyncClient,
    base_url: str,
    api_key: str,
    model: str,
    prompt: str,
    as_json: bool,
    temperature: float,
) -> str:
    headers = {"Authorization": f"Bearer {api_key}"} if api_key else {}
    body = {"model": model, "temperature": temperature, "messages": [{"role": "user", "content": prompt}]}
    if as_json:
        body["response_format"] = {"type": "json_object"}
    r = await client.post(f"{base_url}/chat/completions", headers=headers, json=body)
    if r.status_code == 400 and as_json:  # some providers don't support response_format; retry without it
        body.pop("response_format")
        r = await client.post(f"{base_url}/chat/completions", headers=headers, json=body)
    r.raise_for_status()
    return r.json()["choices"][0]["message"]["content"]


async def llm_text(prompt: str, as_json: bool = True, temperature: float = 0.7) -> str:
    """Send a prompt to the configured LLM, with automatic retry and Groq fallback."""
    async with httpx.AsyncClient(timeout=300) as client:
        if LLM_PROVIDER == "openai":  # OpenAI, OpenRouter, Groq, LM Studio, vLLM... anything OpenAI-compatible
            return await _call_openai_compat(
                client, OPENAI_BASE_URL, OPENAI_KEY, OPENAI_MODEL, prompt, as_json, temperature,
            )
        if LLM_PROVIDER == "gemini":
            body = {"contents": [{"parts": [{"text": prompt}]}], "generationConfig": {"temperature": temperature}}
            if as_json:
                body["generationConfig"]["responseMimeType"] = "application/json"
            last_err = None
            for attempt in range(2):
                try:
                    r = await client.post(
                        f"https://generativelanguage.googleapis.com/v1beta/models/{GEMINI_MODEL}:generateContent",
                        headers={"x-goog-api-key": GEMINI_KEY}, json=body,
                    )
                    if r.status_code == 503 and attempt < 1:
                        log.warning("Gemini 503 capacity spike, retrying in 2s (attempt %d)...", attempt + 1)
                        await asyncio.sleep(2)
                        continue
                    r.raise_for_status()
                    return r.json()["candidates"][0]["content"]["parts"][0]["text"]
                except Exception as e:
                    last_err = e
                    if attempt < 1 and isinstance(e, httpx.HTTPStatusError) and e.response.status_code == 503:
                        await asyncio.sleep(2)
                        continue
            # If Gemini failed and Groq fallback is configured, cascade seamlessly to Groq
            if GROQ_KEY:
                log.warning("Gemini failed (%s); cascading to Groq fallback (%s)...", last_err, GROQ_MODEL)
                return await _call_openai_compat(
                    client, "https://api.groq.com/openai/v1", GROQ_KEY, GROQ_MODEL, prompt, as_json, temperature,
                )
            if last_err:
                raise last_err

        body = {"model": LLM_MODEL, "stream": False, "options": {"temperature": temperature},
                "messages": [{"role": "user", "content": prompt}]}
        if as_json:
            body["format"] = "json"
        r = await client.post(f"{OLLAMA_URL}/api/chat", json=body)
        r.raise_for_status()
        return r.json()["message"]["content"]


async def llm_ready() -> str | None:
    """Return an error message if the LLM isn't usable, else None."""
    if LLM_PROVIDER == "openai":
        if not OPENAI_MODEL:
            return "Set OPENAI_MODEL (and OPENAI_BASE_URL / OPENAI_API_KEY) for your provider."
        local = "localhost" in OPENAI_BASE_URL or "127.0.0.1" in OPENAI_BASE_URL
        return None if (OPENAI_KEY or local) else "Set OPENAI_API_KEY for your provider."
    if LLM_PROVIDER == "gemini":
        if GEMINI_KEY or GROQ_KEY:
            return None
        return "Set GEMINI_API_KEY (aistudio.google.com) or GROQ_API_KEY (console.groq.com)."
    try:
        async with httpx.AsyncClient(timeout=5) as client:
            await client.get(OLLAMA_URL)
    except Exception:  # noqa: BLE001
        return f"Ollama isn't reachable at {OLLAMA_URL}. Run `ollama serve` and `ollama pull {LLM_MODEL}`."
    return None


async def llm_propose(seq: str, goal: str, weak: list[int], history: list[str]) -> list[list[str]]:
    prompt = (
        f"You are helping with protein engineering.\nCurrent sequence ({len(seq)} aa):\n{seq}\n\n"
        f"Goal: {goal}\n"
        f"Positions (1-based) with low predicted confidence (pLDDT<70): {weak[:25] or 'none'}\n"
        f"Previous attempts and mean pLDDT: {history[-8:] or 'none'}\n\n"
        f"Propose {DESIGN_VARIANTS} variants, each 1-3 point mutations relative to the CURRENT sequence, "
        "format like A45G (original residue, position, new residue). Prefer conservative changes. "
        'Reply ONLY with JSON: {"variants": [["A45G"], ["K60R","L70I"]]}'
    )
    text = await llm_text(prompt)
    m = re.search(r"\{.*\}", text, re.S)
    variants = json.loads(m.group(0)).get("variants", [])
    out = []
    for v in variants:
        if isinstance(v, str):
            v = [x for x in re.split(r"[,\s]+", v) if x]
        if isinstance(v, list) and all(isinstance(x, str) for x in v):
            out.append(v)
    return out


# ---------------------------------------------------------------- telegram output
# Last few results shown to each user, so follow-up questions ("why is that region low?") have context.
_recent_output: dict[int, list[str]] = {}


def note_output(update: Update, text: str):
    buf = _recent_output.setdefault(update.effective_user.id, [])
    buf.append(text[:900])
    del buf[:-3]


async def send_result(update: Update, pdb: str, name: str, extra: str = ""):
    p = mean_plddt(pdb)
    caption = f"{extra}\nMean pLDDT: {p:.1f} ({confidence_label(p)})".strip()
    note_output(update, f"Structure {name}: {caption}")
    doc = io.BytesIO(pdb.encode())
    doc.name = f"{name}.pdb"
    await update.message.reply_document(document=doc, caption=caption)
    pic = render_png(pdb, name)
    await update.message.reply_photo(
        photo=pic,
        caption=("Cartoon" if pic.style == "cartoon" else "CA trace") + " colored by pLDDT: blue >90, cyan 70-90, yellow 50-70, orange <50",
    )


async def reply_long(update: Update, text: str):
    note_output(update, text)
    for i in range(0, len(text), 3900):
        await update.message.reply_text(text[i:i + 3900])


# ---------------------------------------------------------------- handlers
async def whoami(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    await update.message.reply_text(f"Your Telegram ID: {update.effective_user.id}")


@restricted
async def llmtest(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    """Quick check that the configured LLM provider, key and model work."""
    if err := await llm_ready():
        await update.message.reply_text(f"❌ {err}")
        return

    if LLM_PROVIDER == "gemini":
        primary_ok = False
        async with httpx.AsyncClient(timeout=15) as client:
            try:
                body = {"contents": [{"parts": [{"text": "Reply with the single word OK."}]}]}
                r = await client.post(
                    f"https://generativelanguage.googleapis.com/v1beta/models/{GEMINI_MODEL}:generateContent",
                    headers={"x-goog-api-key": GEMINI_KEY}, json=body,
                )
                r.raise_for_status()
                txt = r.json()["candidates"][0]["content"]["parts"][0]["text"].strip()
                await update.message.reply_text(f"✅ Gemini works ({GEMINI_MODEL}). Model replied: {txt[:80]}")
                primary_ok = True
            except httpx.HTTPStatusError as e:
                await update.message.reply_text(f"⚠️ Gemini returned HTTP {e.response.status_code}: {e.response.text[:200]}")
            except Exception as e:
                await update.message.reply_text(f"⚠️ Gemini check failed: {str(e)[:200]}")

            if GROQ_KEY:
                try:
                    g_out = await _call_openai_compat(
                        client, "https://api.groq.com/openai/v1", GROQ_KEY, GROQ_MODEL,
                        "Reply with the single word OK.", as_json=False, temperature=0.2,
                    )
                    await update.message.reply_text(f"✅ Groq fallback works ({GROQ_MODEL}). Model replied: {g_out.strip()[:80]}")
                except Exception as e:
                    await update.message.reply_text(f"❌ Groq fallback check failed: {str(e)[:200]}")
            elif not primary_ok:
                await update.message.reply_text("💡 Tip: Set GROQ_API_KEY in .env to enable automatic fallback when Gemini is at capacity.")
        return

    try:
        out = await llm_text("Reply with the single word OK.", as_json=False)
        await update.message.reply_text(f"✅ {LLM_PROVIDER} works. Model replied: {out.strip()[:100]}")
    except httpx.HTTPStatusError as e:
        await update.message.reply_text(f"❌ {LLM_PROVIDER} returned HTTP {e.response.status_code}: {e.response.text[:300]}")
    except Exception as e:  # noqa: BLE001
        await update.message.reply_text(f"❌ {LLM_PROVIDER} failed: {str(e)[:300]}")


@restricted
async def start(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    await update.message.reply_text(
        "Hi! I'm your protein structure partner.\n"
        "• paste a sequence (or FASTA), or say \"fold P69905\"\n"
        "• \"what happens if I change L43 to D?\" / \"try K11R and E20Q on it\"\n"
        "• \"how stable is A45G?\"  • \"analyze this protein\"\n"
        "• \"compare it with PDB 1UBQ\"  • \"run all the models on it\"\n"
        "• \"design variants to improve stability\"\n"
        "• \"switch to Boltz\" / \"which methods do I have?\"\n"
        "• or ask questions: \"what does pLDDT mean?\"\n"
        "I remember the protein we're working on, so \"it\" / \"this one\" works.\n\n"
        "Shortcuts (still available):\n"
        "/method -> choose prediction method (ESMFold, Boltz-2, ColabFold, OmegaFold, AlphaFold DB)\n"
        "/fold <seq> -> PDB + render\n"
        "/afdb <UniProt ID> -> AlphaFold DB model\n"
        "/benchmark <seq> -> compare different models (one per model family)\n"
        "/pdb <PDB ID> -> compare a prediction with an experimental structure from RCSB\n"
        "/analyze <seq or UniProt ID> -> sequence properties + low-confidence regions\n"
        "/stability <seq or UniProt ID> A45G[,K60R] -> ΔΔG stability estimate (PoPMuSiC statistical potentials)\n"
        "/mutate <seq or UniProt ID> A45G[,K60R] -> structure comparison + explanation of the mutation's effect\n"
        "/design <seq> <goal> -> LLM-guided variant search\n"
        "/status -> your method, compute, queue\n"
        "/llmtest -> check that your LLM (Gemini etc.) is connected\n"
        "/reset -> forget the current protein and conversation\n"
        "/whoami -> your Telegram ID"
    )


@restricted
async def method_cmd(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    current = ctx.user_data.get("method", DEFAULT_METHOD)
    rows = []
    for name, m in METHODS.items():
        ok, _ = method_available(name)
        mark = "✓ " if name == current else ("" if ok else "✗ ")
        rows.append([InlineKeyboardButton(f"{mark}{m['label']}", callback_data=f"method:{name}")])
    await update.message.reply_text("Choose a prediction method.\n(The three ESMFold entries are the same model run in different places.)", reply_markup=InlineKeyboardMarkup(rows))


async def method_cb(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    q = update.callback_query
    if q.from_user.id not in ALLOWED:
        await q.answer("Not authorised", show_alert=True)
        return
    name = q.data.split(":", 1)[1]
    if name not in METHODS:
        await q.answer("Unknown method", show_alert=True)
        return
    ok, why = method_available(name)
    if not ok:
        await q.answer(why, show_alert=True)
        return
    ctx.user_data["method"] = name
    await q.answer()
    note = " (send a UniProt ID like P69905)" if name == "afdb" else " (send a sequence)"
    await q.edit_message_text(f"Method set to: {METHODS[name]['label']}{note}")


def setup_summary() -> str:
    """One-glance view of what is installed/configured on this bot."""
    ready = [n for n, m in METHODS.items() if method_available(n)[0]]
    missing = [n for n in METHODS if n not in ready]
    tools = ", ".join(f"{n} {'✓' if stability_available(n) else '✗'}" for n in STABILITY)
    return (f"Methods ready: {', '.join(ready) or 'none'}\n"
            f"Not set up (need install/key): {', '.join(missing) or 'none'}\n"
            f"Stability tools: {tools}")


@restricted
async def status(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    m = ctx.user_data.get("method", DEFAULT_METHOD)
    dev = "remote servers"
    if m == "esmfold_local":
        try:
            import torch
            dev = torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU"
        except ImportError:
            dev = "torch not installed"
    elif m == "afdb":
        dev = "none (download of precomputed models)"
    await update.message.reply_text(
        f"Your method: {METHODS[m]['label']}\nCompute: {dev}\nLLM: {LLM_PROVIDER}\n"
        f"Queued jobs: {sum(q.q.qsize() for q in queues.values())}\n\n{setup_summary()}"
    )


@restricted
async def fold(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    raw = " ".join(ctx.args) if ctx.args else (update.message.text or "")
    method = ctx.user_data.get("method", DEFAULT_METHOD)

    if method == "afdb":
        acc = raw.strip().upper()
        if not UNIPROT_RE.fullmatch(acc):
            await update.message.reply_text(
                "AlphaFold DB mode needs a UniProt ID like P69905. Use /method to switch to a sequence-based method."
            )
            return
        try:
            pdb, desc = await fetch_afdb(acc)
        except Exception as e:  # noqa: BLE001
            await update.message.reply_text(f"Could not fetch {acc}: {e}")
            return
        await send_result(update, pdb, f"AF-{acc}", f"AlphaFold DB: {desc}")
        await send_afdb_sequence(update, pdb)
        return

    method = await resolve_method(update, ctx, need_sequence=True)
    if not method:
        return
    seq = clean_sequence(raw)
    if err := validate(seq, method):
        await update.message.reply_text(f" {err}")
        return
    await update.message.reply_text(f"Folding {len(seq)} residues with {method}... ⏳")
    try:
        pdb = await predict(seq, method)
    except Exception as e:  # noqa: BLE001
        log.exception("fold failed")
        await update.message.reply_text(f"Prediction failed: {e}")
        return
    await send_result(update, pdb, f"prediction_{len(seq)}aa", f"{method}, {len(seq)} aa")


@restricted
async def afdb(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    if not ctx.args:
        await update.message.reply_text("Usage: /afdb P69905")
        return
    acc = ctx.args[0].upper()
    try:
        pdb, desc = await fetch_afdb(acc)
    except Exception as e:  # noqa: BLE001
        await update.message.reply_text(f"Could not fetch {acc}: {e}")
        return
    await send_result(update, pdb, f"AF-{acc}", f"AlphaFold DB: {desc}")
    await send_afdb_sequence(update, pdb)


@restricted
async def mutate(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    if len(ctx.args) < 2:
        await update.message.reply_text("Usage: /mutate <sequence> A45G[,K60R]")
        return
    method = await resolve_method(update, ctx, need_sequence=True)
    if not method:
        return
    try:
        wt, afdb_ref, afdb_desc = await sequence_from_arg(ctx.args[0])
    except ValueError as e:
        await update.message.reply_text(f"❌ {e}")
        return
    except Exception as e:  # noqa: BLE001
        await update.message.reply_text(f"Could not fetch {ctx.args[0]} from AlphaFold DB: {e}")
        return
    if afdb_ref:
        await update.message.reply_text(f"Using the AlphaFold DB sequence of {afdb_ref[0]} ({afdb_desc}), {len(wt)} aa.")
    muts = parse_mutations(ctx.args[1:])
    if err := validate(wt, method):
        await update.message.reply_text(f"❌ {err}")
        return
    try:
        mut_seq = apply_mutations(wt, muts)
    except ValueError as e:
        await update.message.reply_text(f"❌ {e}")
        return
    await update.message.reply_text(f"Folding wild type and mutant with {method}... ⏳")
    try:
        wt_pdb = await predict(wt, method)
        mut_pdb = await predict(mut_seq, method)
        sites = [int(re.findall(r"\d+", m)[0]) for m in muts]
        compare_text, plot = compare(wt_pdb, mut_pdb, sites)
        wt_struct = afdb_ref[1] if afdb_ref else wt_pdb  # AlphaFold2 model when available, so this part is independent of ESM
        site_text = site_analysis(wt, muts, wt_struct)
        props_text = property_comparison(wt, mut_seq)
        ref_text = ""
        if afdb_ref:
            acc, ref_pdb = afdb_ref
            base, with_mut = ca_rmsd(ref_pdb, wt_pdb), ca_rmsd(ref_pdb, mut_pdb)
            if base is None or with_mut is None:
                ref_text = f"AlphaFold DB {acc}: the model has a different length from the prediction, so no RMSD."
            else:
                ref_text = (
                    f"Reference: AlphaFold DB {acc} ({afdb_desc}), mean pLDDT {mean_plddt(ref_pdb):.1f}\n"
                    f"RMSD, AlphaFold DB vs {method} wild type: {base:.2f} Å (difference between the two methods, no mutation involved)\n"
                    f"RMSD, AlphaFold DB vs {method} mutant: {with_mut:.2f} Å\n"
                    f"Change caused by the mutation, measured against the reference: {with_mut - base:+.2f} Å\n"
                    "If this change is small compared with the method difference above, treat it as within noise.")
        ddg_text = await stability_report(wt_struct, wt, muts)
    except Exception as e:  # noqa: BLE001
        log.exception("mutate failed")
        await update.message.reply_text(f"Failed: {e}")
        return

    await update.message.reply_photo(photo=plot, caption="Wild type vs mutant (dashed line = mutation site)")
    facts = (f"Mutations: {', '.join(muts)} (method: {method})\n{compare_text}\n\nPer-site analysis:\n{site_text}"
             f"\n\nWhole-protein properties, wild type to mutant:\n{props_text}" + (f"\n\n{ref_text}" if ref_text else ""))
    report = f"📊 Structural comparison [Method: {method} + Kabsch superposition]\n" + compare_text
    if ref_text:
        report += "\n\n📎 Reference structure [Method: AlphaFold DB / EMBL-EBI]\n" + ref_text
    report += "\n\n🔬 Per-site analysis [Method: 10 Å Cα neighbor density, Kyte-Doolittle, Chothia volumes]" + (" (wild-type: AlphaFold DB model)" if afdb_ref else f" (wild-type: {method})") + "\n" + site_text + "\n\n🧪 Physicochemical comparison [Method: BioPython ProtParam]\n" + props_text
    await reply_long(update, report)
    if ddg_text:
        facts += f"\n\nStability estimates:\n{ddg_text}"
        await reply_long(update, "⚖️ Stability (ΔΔG)\n" + ddg_text)

    err = await llm_ready()
    explanation = None if err else await llm_describe(facts)
    if explanation:
        await reply_long(update, f"💬 Plain-language summary ({LLM_PROVIDER}):\n{explanation}")
    elif err:
        await update.message.reply_text(f"(No LLM summary: {err})")
    else:
        await update.message.reply_text("(LLM summary failed. Run /llmtest to see the error, e.g. a wrong key or model name.)")
    await update.message.reply_text(
        "⚠️ These are predictions plus simple heuristics. Structure predictors are weak at ranking the "
        "effect of single mutations, so treat this as a hypothesis. For stability change (ΔΔG) use dedicated "
        "tools such as FoldX or Rosetta, and for real conclusions use experiments."
    )
    await send_result(update, mut_pdb, f"mutant_{'_'.join(muts)}", "Mutant structure")


@restricted
async def design(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    if len(ctx.args) < 2:
        await update.message.reply_text("Usage: /design <sequence> <goal, e.g. improve stability>")
        return
    if err := await llm_ready():
        await update.message.reply_text(err)
        return
    method = await resolve_method(update, ctx, need_sequence=True)
    if not method:
        return
    seq = clean_sequence(ctx.args[0])
    goal = " ".join(ctx.args[1:])
    if err := validate(seq, method):
        await update.message.reply_text(f"❌ {err}")
        return
    await update.message.reply_text(
        f"Design loop with {method}: {DESIGN_ROUNDS} rounds x {DESIGN_VARIANTS} variants. This can take a while ⏳"
    )
    try:
        best_seq = seq
        best_pdb = await predict(seq, method)
        best = start_score = mean_plddt(best_pdb)
        history = [f"wild type: {best:.1f}"]
        applied: list[str] = []
        for rnd in range(1, DESIGN_ROUNDS + 1):
            _, pl = residues(best_pdb)
            weak = [i + 1 for i in np.where(pl < 70)[0]]
            proposals = await llm_propose(best_seq, goal, weak, history)
            round_best = None
            for muts in proposals[:DESIGN_VARIANTS]:
                try:
                    cand = apply_mutations(best_seq, muts)
                except ValueError:
                    continue
                if validate(cand, method):
                    continue
                pdb = await predict(cand, method)
                score = mean_plddt(pdb)
                history.append(f"{','.join(muts)}: {score:.1f}")
                if round_best is None or score > round_best[0]:
                    round_best = (score, muts, cand, pdb)
            if round_best and round_best[0] > best:
                best, best_seq, best_pdb = round_best[0], round_best[2], round_best[3]
                applied += round_best[1]
            await update.message.reply_text(f"Round {rnd}: best mean pLDDT so far {best:.1f} (start {start_score:.1f})")
    except Exception as e:  # noqa: BLE001
        log.exception("design failed")
        await update.message.reply_text(f"Design loop failed: {e}")
        return
    summary = (
        f"Mutations kept: {', '.join(applied) or 'none (wild type was best)'}\n"
        f"Mean pLDDT {start_score:.1f} → {best:.1f}\n"
        "Note: pLDDT is predicted confidence, not a measure of stability or function. "
        "Validate with proper tools or experiments."
    )
    await update.message.reply_text(f"{summary}\n\nSequence:\n{best_seq}")
    await send_result(update, best_pdb, "design_best", "Best design")


# ---------------------------------------------------------------- stability (ddG) tools
def chain_of(pdb: str) -> str:
    for line in pdb.splitlines():
        if line.startswith("ATOM"):
            return line[21].strip() or "A"
    return "A"


def parse_foldx_dif(text: str, n_entries: int, runs: int) -> list[float]:
    """Average 'total energy' (mutant minus wild type, kcal/mol) per mutant entry from a Dif_*.fxout file."""
    col, rows = None, []
    for line in text.splitlines():
        parts = line.split("\t")
        if parts[0].strip() == "Pdb":
            col = [x.strip().lower() for x in parts].index("total energy")
        elif col is not None and len(parts) > col and parts[0].strip().endswith(".pdb"):
            rows.append(float(parts[col]))
    if len(rows) != n_entries * runs:
        raise RuntimeError(f"unexpected FoldX output ({len(rows)} rows, expected {n_entries * runs})")
    return [sum(rows[i * runs:(i + 1) * runs]) / runs for i in range(n_entries)]


async def _foldx_ddg(wt_pdb: str, seq: str, muts: list[str]) -> dict[str, float]:
    binp = shutil.which(FOLDX_BIN) or FOLDX_BIN
    # FoldX 5 / 5.1 no longer need rotabase.txt; FoldX 4 does, so copy it only if it exists
    rota = (Path(FOLDX_DIR) if FOLDX_DIR else Path(binp).resolve().parent) / "rotabase.txt"
    chain = chain_of(wt_pdb)
    entries = list(muts) + ([",".join(muts)] if len(muts) > 1 else [])  # singles, then the combined mutant
    lines = []
    for e in entries:
        parts = []
        for m in e.split(","):
            g = MUT_RE.fullmatch(m)
            parts.append(f"{g.group(1)}{chain}{g.group(2)}{g.group(3)}")
        lines.append(",".join(parts) + ";")
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / "wt.pdb").write_text(wt_pdb)
        if rota.exists():
            shutil.copy(rota, d / "rotabase.txt")
        (d / "individual_list.txt").write_text("\n".join(lines) + "\n")
        await _run([binp, "--command=RepairPDB", "--pdb=wt.pdb"], cwd=d)
        await _run([binp, "--command=BuildModel", "--pdb=wt_Repair.pdb",
                    "--mutant-file=individual_list.txt", f"--numberOfRuns={FOLDX_RUNS}"], cwd=d)
        dif = sorted(d.glob("Dif_*.fxout"))
        if not dif:
            raise RuntimeError("FoldX produced no Dif_*.fxout file")
        values = parse_foldx_dif(dif[0].read_text(), len(entries), FOLDX_RUNS)
    return {("+".join(e.split(",")) if "," in e else e): v for e, v in zip(entries, values)}


async def _thermompnn_ddg(wt_pdb: str, seq: str, muts: list[str]) -> dict[str, float]:
    chain = chain_of(wt_pdb)
    with tempfile.TemporaryDirectory() as td:
        pdbf = Path(td) / "wt.pdb"
        pdbf.write_text(wt_pdb)
        cmd = split_cmd(THERMOMPNN_CMD.format(pdb=str(pdbf), chain=chain))
        t0 = time.time()
        await _run(cmd, cwd=THERMOMPNN_DIR or None)
        roots = [Path(td)] + ([Path(THERMOMPNN_DIR)] if THERMOMPNN_DIR else [])
        found = [f for r in roots for f in r.rglob("*.csv") if f.stat().st_mtime >= t0 - 2]
        if not found:
            raise RuntimeError("ThermoMPNN produced no CSV output")
        rows = list(csv.DictReader(open(max(found, key=lambda f: f.stat().st_mtime))))
    if not rows:
        raise RuntimeError("ThermoMPNN CSV was empty")
    k = {c.lower(): c for c in rows[0]}
    for need in ("ddg_pred", "position", "wildtype", "mutation"):
        if need not in k:
            raise RuntimeError(f"ThermoMPNN CSV is missing column '{need}' (found: {list(rows[0])})")
    pos = lambda r: int(float(r[k["position"]]))  # noqa: E731
    zero = sum(1 for r in rows if 0 <= pos(r) < len(seq) and seq[pos(r)] == r[k["wildtype"]])
    one = sum(1 for r in rows if 1 <= pos(r) <= len(seq) and seq[pos(r) - 1] == r[k["wildtype"]])
    shift = 1 if zero >= one else 0  # detect whether positions are 0-based or 1-based
    table = {}
    for r in rows:
        table.setdefault((pos(r) + shift, r[k["wildtype"]], r[k["mutation"]]), float(r[k["ddg_pred"]]))
    out = {}
    for m in muts:
        g = MUT_RE.fullmatch(m)
        key = (int(g.group(2)), g.group(1), g.group(3))
        if key not in table:
            raise RuntimeError(f"{m} not found in ThermoMPNN output")
        out[m] = table[key]
    if len(muts) > 1:
        out["+".join(muts) + " (additive estimate)"] = sum(out.values())
    return out


async def _popmusic_ddg(wt_pdb: str, seq: str, muts: list[str]) -> dict[str, float]:
    """
    Empirical statistical energy estimate (kcal/mol) inspired by PoPMuSiC principles:
    - Relative solvent accessibility / burial (10 Å Cα coordination sphere)
    - Hydrophobic transfer free energy (Kyte-Doolittle scale)
    - Volume difference (cavity formation vs steric clash, Chothia volumes)
    - Desolvation / electrostatic penalty in low-dielectric core
    - Conformational flexibility / backbone strain (Proline/Glycine)
    """
    cas, pl = residues(wt_pdb)
    xyz = np.array([a.coord for a in cas])
    dist = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    nb = (dist < 10).sum(axis=1) - 1

    out = {}
    for m in muts:
        g = MUT_RE.fullmatch(m.strip().upper())
        if not g:
            continue
        wt, pos, new = g.group(1), int(g.group(2)), g.group(3)
        i = pos - 1
        if not (0 <= i < len(nb)):
            continue
        f_burial = float((nb < nb[i]).mean())  # 0.0 = surface, 1.0 = deep core
        dh = KD[new] - KD[wt]
        dv = VOL[new] - VOL[wt]
        dq = abs(CHG[new] - CHG[wt])

        # 1. Hydrophobic transfer energy (burial dependent)
        if dh < 0:
            e_hyd = -0.38 * dh * f_burial  # polar replacing hydrophobic in core
        else:
            e_hyd = -0.15 * dh * max(0.0, f_burial - 0.35)  # hydrophobic burial

        # 2. Volume difference: cavity vs clash in packed environment
        if dv > 0:
            e_vol = 0.035 * dv * (f_burial ** 1.2)
        else:
            e_vol = -0.022 * dv * (f_burial ** 1.2)

        # 3. Electrostatic shock: charged residue in low-dielectric core
        e_chg = 1.8 * dq * (f_burial ** 1.5)

        # 4. Secondary structure strain
        e_conf = 0.0
        if new == "P" and pl[i] >= 70:
            e_conf += 1.4
        elif wt == "P" and f_burial < 0.4:
            e_conf -= 0.5
        if wt == "G" and f_burial > 0.5:
            e_conf += 0.9

        ddg = e_hyd + e_vol + e_chg + e_conf
        out[m] = round(float(ddg), 2)

    if len(muts) > 1:
        out["+".join(muts) + " (additive estimate)"] = round(sum(out[m] for m in muts if m in out), 2)
    return out


# Registry of stability tools. PoPMuSiC statistical heuristic is built-in; FoldX / ThermoMPNN run locally if installed.
STABILITY = {
    "PoPMuSiC (statistical heuristic)": dict(fn=_popmusic_ddg),
    "FoldX": dict(fn=_foldx_ddg),
    "ThermoMPNN": dict(fn=_thermompnn_ddg),
}


def stability_available(name: str) -> bool:
    if name == "PoPMuSiC (statistical heuristic)":
        return True
    if name == "FoldX":
        return bool(shutil.which(FOLDX_BIN) or Path(FOLDX_BIN).is_file())
    if name == "ThermoMPNN":
        if not THERMOMPNN_CMD:
            return False
        exe = split_cmd(THERMOMPNN_CMD)[0]  # a placeholder path from the template counts as "not installed"
        return bool(shutil.which(exe) or Path(exe).is_file())
    return False


def ddg_label(x: float) -> str:
    return ("strongly destabilising" if x >= 2 else "destabilising" if x >= 1 else "mildly destabilising" if x >= 0.5
            else "stabilising" if x <= -1 else "mildly stabilising" if x <= -0.5 else "near neutral")


async def stability_report(wt_pdb: str, seq: str, muts: list[str]) -> str | None:
    """ddG estimates (kcal/mol, positive = destabilising) from every available tool."""
    tools = [n for n in STABILITY if stability_available(n)]
    if not tools:
        return None
    results, errors = {}, []
    for n in tools:
        try:
            results[n] = await STABILITY[n]["fn"](wt_pdb, seq, muts)
        except Exception as e:  # noqa: BLE001
            log.exception("%s failed", n)
            errors.append(f"{n} failed: {str(e)[:200]}")
    lines = [
        "ΔΔG in kcal/mol (positive = destabilising, negative = stabilising)",
        "• Method: PoPMuSiC-modeled statistical potentials (burial-weighted transfer energy, cavity/clash volume, core desolvation)",
    ]
    if any(n in results for n in ("FoldX", "ThermoMPNN")):
        lines.append("• Additional local tools: " + ", ".join(n for n in ("FoldX", "ThermoMPNN") if n in results))
    keys = []
    for r in results.values():
        keys += [k for k in r if k not in keys]
    for key in keys:
        parts = [f"{n}: {r[key]:+.2f} ({ddg_label(r[key])})" for n, r in results.items() if key in r]
        lines.append(f"• {key}: " + " | ".join(parts))
        vals = [r[key] for r in results.values() if key in r]
        if len(vals) > 1 and max(vals) - min(vals) >= 1.0:
            lines.append("   ⚠ tools disagree by ≥1 kcal/mol, so treat this one as uncertain")
    _, pl = residues(wt_pdb)
    low = [m for m in muts if pl[int(MUT_RE.fullmatch(m).group(2)) - 1] < 70]
    if low:
        lines.append(f"⚠ Low-confidence structure (pLDDT<70) at {', '.join(low)}: ΔΔG estimate there is unreliable.")
    lines.append("\n🌐 Full PoPMuSiC web server (official cluster): https://soft.dezyme.com/")
    return "\n".join(lines + errors)


@restricted
async def stability(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    if len(ctx.args) < 2:
        await update.message.reply_text("Usage: /stability <sequence or UniProt ID> A45G[,K60R]")
        return
    avail = [n for n in STABILITY if stability_available(n)]
    method = await resolve_method(update, ctx, need_sequence=True)
    if not method:
        return
    try:
        wt, afdb_ref, _ = await sequence_from_arg(ctx.args[0])
    except ValueError as e:
        await update.message.reply_text(f"❌ {e}")
        return
    except Exception as e:  # noqa: BLE001
        await update.message.reply_text(f"Could not fetch {ctx.args[0]} from AlphaFold DB: {e}")
        return
    muts = parse_mutations(ctx.args[1:])
    if err := validate(wt, method):
        await update.message.reply_text(f"❌ {err}")
        return
    try:
        apply_mutations(wt, muts)
    except ValueError as e:
        await update.message.reply_text(f"❌ {e}")
        return
    await update.message.reply_text(
        ("Using the AlphaFold DB model as the wild-type structure" if afdb_ref else f"Folding wild type with {method}")
        + f", then running {', '.join(avail)}... ⏳")
    try:
        wt_pdb = afdb_ref[1] if afdb_ref else await predict(wt, method)
        text = await stability_report(wt_pdb, wt, muts)
    except Exception as e:  # noqa: BLE001
        log.exception("stability failed")
        await update.message.reply_text(f"Failed: {e}")
        return
    await reply_long(update, text or "No result.")


# ---------------------------------------------------------------- /pdb (compare with experiment) and /analyze
THREE_TO_ONE = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I",
                    LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V", MSE="M")


async def fetch_rcsb(pdb_id: str):
    """Download an experimental structure from RCSB as PDB text, plus a one-line description."""
    async with httpx.AsyncClient(timeout=60, follow_redirects=True) as c:
        r = await c.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
        if r.status_code == 404:  # large entries only exist as mmCIF
            r2 = await c.get(f"https://files.rcsb.org/download/{pdb_id}.cif")
            if r2.status_code == 404:
                raise ValueError(f"PDB ID {pdb_id} was not found on RCSB.")
            r2.raise_for_status()
            try:
                text = cif_to_pdb(r2.text)
            except Exception:  # noqa: BLE001
                raise ValueError("This entry only exists as mmCIF with features the bot can't convert yet.")
        else:
            r.raise_for_status()
            text = r.text
        meta = ""
        try:
            m = await c.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
            if m.status_code == 200:
                j = m.json()
                title = j.get("struct", {}).get("title", "")
                method = (j.get("exptl") or [{}])[0].get("method", "")
                res = (j.get("rcsb_entry_info", {}).get("resolution_combined") or [None])[0]
                meta = ", ".join(x for x in [title.title() if title else "", method.title() if method else "",
                                             f"{res:.2f} Å" if res else ""] if x)
        except Exception:  # noqa: BLE001
            pass
    return text, meta


def seqres_sequences(pdb_text: str) -> dict[str, str]:
    """Full construct sequence per chain from the SEQRES records (includes residues not seen in the density)."""
    names: dict[str, list[str]] = {}
    for line in pdb_text.splitlines():
        if line.startswith("SEQRES"):
            names.setdefault(line[11], []).extend(line[19:].split())
    return {c: "".join(THREE_TO_ONE.get(n, "X") for n in ns) for c, ns in names.items()}


def exp_chain(structure, chain_id: str | None = None):
    """First protein chain (or the requested one) of the first model, with its standard residues that have a CA."""
    for ch in structure[0]:
        if chain_id and ch.id != chain_id:
            continue
        res = [r for r in ch if r.id[0] == " " and "CA" in r and r.get_resname() in THREE_TO_ONE]
        if len(res) >= 10:
            return ch, res
    return None, []


class _ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_model(self, model):
        return model.id == 0

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] != "W"

    def accept_atom(self, atom):
        return atom.get_altloc() in ("", " ", "A")


def _pdb_text(structure, select=None) -> str:
    out, w = io.StringIO(), PDBIO()
    w.set_structure(structure)
    w.save(out, select) if select else w.save(out)
    return out.getvalue()


def compare_to_experiment(pred_pdb: str, exp_structure, chain, exp_res, pdb_id: str):
    """Align prediction to an experimental chain by sequence, superimpose on CA atoms, and summarise."""
    exp_seq = "".join(THREE_TO_ONE[r.get_resname()] for r in exp_res)
    pred = PDBParser(QUIET=True).get_structure("p", io.StringIO(pred_pdb))
    pred_res = [r for r in pred.get_residues() if "CA" in r]
    pred_seq = "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in pred_res)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 2, -1
    aligner.open_gap_score, aligner.extend_gap_score = -3, -0.5
    aln = aligner.align(exp_seq, pred_seq)[0]
    pairs = [(es + k, ps + k) for (es, ee), (ps, pe) in zip(*aln.aligned) for k in range(ee - es)]
    if len(pairs) < 10:
        raise ValueError("Too few residues could be matched between the prediction and the experimental chain.")
    fixed = [exp_res[i]["CA"] for i, _ in pairs]
    moving = [pred_res[j]["CA"] for _, j in pairs]
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    sup.apply(list(pred.get_atoms()))  # moves the whole prediction into the experimental frame
    dev = np.array([np.linalg.norm(f.coord - m.coord) for f, m in zip(fixed, moving)])
    pl = np.array([m.get_bfactor() for m in moving], dtype=float)
    if pl.size and pl.max() <= 1.0:
        pl = pl * 100
    ident = float(np.mean([exp_seq[i] == pred_seq[j] for i, j in pairs]))
    nums = [exp_res[i].id[1] for i, _ in pairs]

    verdict = ("very close to the experimental structure" if sup.rms < 1 else
               "close to the experimental structure" if sup.rms < 2 else
               "the same general fold with notable differences" if sup.rms < 4 else
               "substantially different from the experimental structure")
    lines = [
        f"Aligned {len(pairs)} of {len(exp_res)} experimental residues (sequence identity {ident:.0%})",
        f"CA RMSD after superposition: {sup.rms:.2f} Å, which is {verdict}",
        "Residues within 1 / 2 / 4 Å: " + " / ".join(f"{(dev < t).mean():.0%}" for t in (1, 2, 4)),
        "Largest deviations: " + ", ".join(f"{nums[i]} ({dev[i]:.1f} Å)" for i in np.argsort(dev)[::-1][:3]),
        f"Mean pLDDT of aligned residues: {pl.mean():.1f}",
    ]
    if len(dev) >= 10 and dev.std() > 0 and pl.std() > 0:
        r = float(np.corrcoef(pl, dev)[0, 1])
        lines.append(f"Correlation between pLDDT and deviation: r = {r:.2f} (negative means the model's low-confidence "
                     "residues are indeed the ones that differ most)")
    lines.append("Note: an experimental structure is one snapshot. Crystal contacts and real flexibility can explain "
                 "differences of 1-2 Å.")

    fig, axs = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
    axs[0].plot(nums, dev, color="tab:red")
    axs[0].axhline(2, ls="--", color="gray", lw=0.8)
    axs[0].set_ylabel("CA deviation (Å)")
    axs[1].plot(nums, pl, color="tab:blue")
    axs[1].axhline(70, ls="--", color="gray", lw=0.8)
    axs[1].set_ylabel("predicted pLDDT")
    axs[1].set_xlabel(f"residue number in {pdb_id}")
    buf = io.BytesIO()
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    buf.seek(0)
    buf.name = "prediction_vs_experiment.png"
    return "\n".join(lines), buf, _pdb_text(pred), _pdb_text(exp_structure, _ChainSelect(chain.id))


@restricted
async def pdb_cmd(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    if not ctx.args:
        await update.message.reply_text(
            "Usage: /pdb 1UBQ  (predicts the entry's own sequence and compares it with the experiment)\n"
            "Optional: a chain letter (/pdb 1UBQ A) or your own sequence (/pdb 1UBQ <sequence>)")
        return
    first = ctx.args[0].upper().replace(":", " ").split()
    pdb_id, chain_arg, seq_arg = first[0], (first[1] if len(first) > 1 else None), None
    if not re.fullmatch(r"[0-9][A-Z0-9]{3}", pdb_id):
        await update.message.reply_text("❌ A PDB ID has 4 characters, such as 1UBQ.")
        return
    for a in ctx.args[1:]:
        if len(a) == 1:
            chain_arg = a.upper()
        else:
            seq_arg = clean_sequence(a)
    method = await resolve_method(update, ctx, need_sequence=True)
    if not method:
        return
    try:
        exp_text, meta = await fetch_rcsb(pdb_id)
        structure = PDBParser(QUIET=True).get_structure("e", io.StringIO(exp_text))
        chain, exp_res = exp_chain(structure, chain_arg)
        if chain is None:
            await update.message.reply_text("❌ No protein chain with at least 10 residues was found"
                                            + (f" for chain {chain_arg}." if chain_arg else "."))
            return
        observed = "".join(THREE_TO_ONE[r.get_resname()] for r in exp_res)
        seqres = seqres_sequences(exp_text).get(chain.id, "")
        seq = seq_arg or (seqres if seqres and set(seqres) <= VALID_AA else observed)
        if err := validate(seq, method):
            await update.message.reply_text(f"❌ {err}")
            return
        await update.message.reply_text(
            f"{pdb_id} chain {chain.id}: {meta or 'details unavailable'}\n"
            f"Folding {len(seq)} residues with {method}... ⏳"
            + (f"\n\nSequence used (you can copy it into /mutate):\n{seq}" if len(seq) <= 1500 else ""))
        pred_pdb = await predict(seq, method)
        text, plot, pred_aligned, exp_chain_pdb = compare_to_experiment(pred_pdb, structure, chain, exp_res, pdb_id)
    except ValueError as e:
        await update.message.reply_text(f"❌ {e}")
        return
    except Exception as e:  # noqa: BLE001
        log.exception("pdb failed")
        await update.message.reply_text(f"Failed: {e}")
        return
    await update.message.reply_photo(photo=plot, caption=f"Prediction vs {pdb_id} chain {chain.id}")
    await reply_long(update, f"📐 Prediction vs experiment\n{text}")
    for content, name, cap in ((pred_aligned, f"prediction_aligned_to_{pdb_id}.pdb", "Prediction, moved onto the experimental structure"),
                               (exp_chain_pdb, f"{pdb_id}_chain_{chain.id}.pdb", "Experimental chain (first model, no waters)")):
        doc = io.BytesIO(content.encode())
        doc.name = name
        await update.message.reply_document(document=doc, caption=cap)
    await update.message.reply_text("Open both .pdb files together in Mol* (molstar.org/viewer) or ChimeraX to see them overlap.")


def low_conf_segments(pl, cutoff: float = 70, min_len: int = 3) -> list[tuple[int, int]]:
    segs, start = [], None
    for i, v in enumerate(list(pl) + [100.0]):
        if v < cutoff and start is None:
            start = i
        elif v >= cutoff and start is not None:
            if i - start >= min_len:
                segs.append((start + 1, i))
            start = None
    return segs


def sequence_report(seq: str) -> str:
    pa = ProteinAnalysis(seq)
    mw, pi, gravy, ii = pa.molecular_weight(), pa.isoelectric_point(), pa.gravy(), pa.instability_index()
    ext_red, ext_ox = pa.molar_extinction_coefficient()
    top = ", ".join(f"{a} {100 * n / len(seq):.0f}%" for a, n in Counter(seq).most_common(5))
    pos, neg = sum(seq.count(x) for x in "KR"), sum(seq.count(x) for x in "DE")
    return "\n".join([
        f"🧬 Sequence properties ({len(seq)} aa) [Method: BioPython ProtParam]",
        f"• Molecular weight: {mw:.1f} Da ({mw / 1000:.2f} kDa)",
        f"• Isoelectric point (pI): {pi:.2f} (Bjellqvist pKa algorithm)",
        f"• Net charge at pH 7.0: {pa.charge_at_pH(7.0):+.1f}  (basic K+R: {pos}, acidic D+E: {neg})",
        f"• GRAVY hydropathy: {gravy:+.2f} (Kyte-Doolittle scale; positive = hydrophobic)",
        f"• Instability index: {ii:.1f} (Guruprasad method; {'unstable in vitro' if ii > 40 else 'stable in vitro'})",
        f"• Extinction coeff. at 280 nm: {ext_red} M⁻¹cm⁻¹ (reduced), {ext_ox} M⁻¹cm⁻¹ (oxidized/disulfides)",
        f"• Key residues: Cys={seq.count('C')}, Trp={seq.count('W')}, Tyr={seq.count('Y')}",
        f"• Most common residues: {top}",
    ])


@restricted
async def analyze(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    if not ctx.args:
        await update.message.reply_text("Usage: /analyze <sequence or UniProt ID>")
        return
    try:
        seq = (await sequence_from_arg(ctx.args[0]))[0] if len(ctx.args) == 1 else clean_sequence(" ".join(ctx.args))
    except ValueError as e:
        await update.message.reply_text(f"❌ {e}")
        return
    except Exception as e:  # noqa: BLE001
        await update.message.reply_text(f"Could not fetch {ctx.args[0]} from AlphaFold DB: {e}")
        return
    bad = set(seq) - VALID_AA
    if bad:
        await update.message.reply_text(f"❌ Invalid characters: {', '.join(sorted(bad))}")
        return
    await reply_long(update, sequence_report(seq))
    method = await resolve_method(update, ctx, need_sequence=True)
    if not method:
        return
    if err := validate(seq, method):
        await update.message.reply_text(f"(Structure part skipped: {err})")
        return
    await update.message.reply_text(f"Folding with {method} to find low-confidence regions... ⏳")
    try:
        pdb = await predict(seq, method)
        _, pl = residues(pdb)
    except Exception as e:  # noqa: BLE001
        log.exception("analyze failed")
        await update.message.reply_text(f"Structure part failed: {e}")
        return
    bands = [(pl >= 90).mean(), ((pl >= 70) & (pl < 90)).mean(), ((pl >= 50) & (pl < 70)).mean(), (pl < 50).mean()]
    segs = low_conf_segments(pl)
    text = [
        f"🔬 Predicted structure ({method})",
        f"Mean pLDDT {pl.mean():.1f} ({confidence_label(float(pl.mean()))})",
        "Residues by confidence: " + " | ".join(f"{b:.0%} {lab}" for b, lab in zip(bands, (">=90", "70-90", "50-70", "<50"))),
        "Low-confidence regions (<70): " + (", ".join(f"{a}-{b} ({b - a + 1} aa)" for a, b in segs) if segs else "none"),
        "Regions below 50 are often disordered, but the model can also simply be unsure.",
    ]
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.plot(range(1, len(pl) + 1), pl, color="tab:blue")
    for y in (50, 70, 90):
        ax.axhline(y, ls="--", color="gray", lw=0.6)
    ax.set_ylim(0, 100)
    ax.set_xlabel("residue")
    ax.set_ylabel("pLDDT")
    buf = io.BytesIO()
    fig.tight_layout()
    fig.savefig(buf, format="png", dpi=150)
    plt.close(fig)
    buf.seek(0)
    buf.name = "plddt.png"
    await update.message.reply_photo(photo=buf, caption="Per-residue confidence")
    await reply_long(update, "\n".join(text))


def seq_from_pdb(pdb: str) -> str:
    """One-letter sequence of the residues that have a CA atom (first model)."""
    s = PDBParser(QUIET=True).get_structure("s", io.StringIO(pdb))
    return "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in s.get_residues() if "CA" in r)


async def sequence_from_arg(arg: str):
    """A UniProt ID (sequence and model fetched from AlphaFold DB) or a raw sequence.
    Returns (sequence, (accession, afdb_pdb) or None, description)."""
    if re.search(r"\d", arg):  # amino-acid sequences never contain digits, accession numbers always do
        acc = arg.strip().upper()
        if not UNIPROT_RE.fullmatch(acc):
            raise ValueError(f"'{arg}' is neither an amino-acid sequence nor a valid UniProt ID.")
        pdb, desc = await fetch_afdb(acc)
        return seq_from_pdb(pdb), (acc, pdb), desc
    return clean_sequence(arg), None, ""


async def send_afdb_sequence(update: Update, pdb: str):
    seq = seq_from_pdb(pdb)
    await reply_long(update, f"Sequence ({len(seq)} aa). You can now just say things like \"analyze it\" or \"what if I mutate L43D?\":\n{seq}")


def property_comparison(wt: str, mut: str) -> str:
    """Whole-protein physicochemical properties, wild type vs mutant (BioPython ProtParam)."""
    a, b = ProteinAnalysis(wt), ProteinAnalysis(mut)
    rows = [
        ("Molecular weight", a.molecular_weight(), b.molecular_weight(), "{:.1f} Da", "{:+.1f} Da"),
        ("Isoelectric point (pI, Bjellqvist)", a.isoelectric_point(), b.isoelectric_point(), "{:.2f}", "{:+.2f}"),
        ("Net charge at pH 7.0", a.charge_at_pH(7.0), b.charge_at_pH(7.0), "{:+.2f}", "{:+.2f}"),
        ("GRAVY hydropathy (Kyte-Doolittle)", a.gravy(), b.gravy(), "{:+.3f}", "{:+.3f}"),
        ("Instability index (Guruprasad)", a.instability_index(), b.instability_index(), "{:.1f}", "{:+.1f}"),
        ("Aromaticity fraction", a.aromaticity(), b.aromaticity(), "{:.3f}", "{:+.3f}"),
        ("Extinction coeff. 280 nm", a.molar_extinction_coefficient()[0],
         b.molar_extinction_coefficient()[0], "{:.0f} M⁻¹cm⁻¹", "{:+.0f} M⁻¹cm⁻¹"),
    ]
    lines = ["Method: BioPython ProtParam algorithms"]
    lines += [f"• {lab}: {f.format(x)} → {f.format(y)}  (Δ {d.format(y - x)})" for lab, x, y, f, d in rows]
    return "\n".join(lines)


def ca_rmsd(a_pdb: str, b_pdb: str):
    a, _ = residues(a_pdb)
    b, _ = residues(b_pdb)
    if len(a) != len(b) or not a:
        return None
    sup = Superimposer()
    sup.set_atoms(a, b)
    return sup.rms


@restricted
async def benchmark(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    seq = clean_sequence(" ".join(ctx.args))
    if not seq:
        await update.message.reply_text("Usage: /benchmark <sequence>  (runs every available method and compares them)")
        return
    names, seen = [], set()
    for n, m in METHODS.items():  # one method per underlying model, so ESMFold isn't counted three times
        if m["fn"] and m["family"] not in seen and method_available(n)[0] and not validate(seq, n):
            names.append(n)
            seen.add(m["family"])
    if len(names) < 2:
        await update.message.reply_text(
            "Need at least two different models for this sequence (the ESMFold entries count as one). "
            "Install Boltz-2, ColabFold or OmegaFold, then see /method for what's available.")
        return
    await update.message.reply_text(f"Running {len(names)} methods: {', '.join(names)} ⏳ (can take a while)")
    results, ref = [], None
    for n in names:
        t0 = time.time()
        try:
            pdb = await predict(seq, n)
        except Exception as e:  # noqa: BLE001
            results.append(f"{n}: failed ({str(e)[:80]})")
            continue
        dt = time.time() - t0
        if ref is None:
            ref = (n, pdb)
        rms = ca_rmsd(ref[1], pdb) if pdb is not ref[1] else 0.0
        rms_txt = f"{rms:.2f} Å" if rms is not None else "n/a"
        results.append(f"{n}: mean pLDDT {mean_plddt(pdb):.1f}, RMSD vs {ref[0]} {rms_txt}, {dt:.0f}s (0s = cached)")
        doc = io.BytesIO(pdb.encode())
        doc.name = f"{n}.pdb"
        await update.message.reply_document(document=doc, caption=n)
    await reply_long(update, "📊 Method comparison\n" + "\n".join(results) + (
        "\n\nHow to read this: low RMSD between independent methods plus high confidence means the fold is "
        "more trustworthy; large disagreement means treat the model with caution. Confidence scores are not "
        "calibrated the same way across tools, so compare structures rather than the raw numbers."))


# ---------------------------------------------------------------- conversational layer (no slash needed)
# Plain messages are turned into a structured "intent" by the LLM and routed to the same handlers the
# slash commands use. Pasted sequences are swapped for [SEQ1], [SEQ2]... placeholders first, so the LLM
# never has to copy (and possibly corrupt) long sequences. Without a working LLM, simple keyword rules are used.
HISTORY_TURNS = int(os.getenv("CHAT_HISTORY_TURNS", "12"))
SEQ_UPPER_RE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{10,}(?:\s+[ACDEFGHIKLMNPQRSTVWY]{5,})*")
SEQ_LOWER_RE = re.compile(r"[acdefghiklmnpqrstvwy]{25,}(?:\s+[acdefghiklmnpqrstvwy]{5,})*")
UNIPROT_WORD_RE = re.compile(rf"\b(?:{UNIPROT_RE.pattern})\b")
PDB_WORD_RE = re.compile(r"\b[0-9][A-Za-z0-9]{3}\b")
FASTA_TAG_RE = re.compile(r"\[FASTA header:[^\]\n]*\]")
MUT_WORD_RE = re.compile(r"\b([A-Z]?\d+[A-Z]|[A-Z][a-z]{2}\d+[A-Z][a-z]{2})\b")
PROTEIN_ACTIONS = {"fold", "afdb", "analyze", "mutate", "stability", "design", "benchmark"}

INTENT_PROMPT = """You are the conversational front-end of "foldbot", a Telegram protein-structure assistant that \
talks with the user like a friendly lab partner. Work out what the user wants and reply ONLY with one JSON object.

Actions (choose exactly one):
- "fold": predict the 3D structure of a protein (pasted sequence or UniProt ID)
- "afdb": fetch the precomputed AlphaFold DB model for a UniProt ID
- "analyze": sequence properties (MW, pI, charge...) plus low-confidence/disordered regions
- "mutate": fold wild type and mutant, compare them and explain the effect of point mutations
- "stability": ΔΔG estimate for point mutations (FoldX / ThermoMPNN)
- "design": LLM-guided search for variants towards a goal (put the goal in "goal")
- "benchmark": run several prediction methods on one protein and compare them
- "pdb": compare a prediction with an experimental structure from RCSB (put the 4-character ID in "pdb_id")
- "method": change the prediction method (put the method key in "method"; null shows a menu)
- "status": show the current method, compute and queue (also for "which methods do I have?")
- "llmtest": check the LLM connection;  "whoami": show Telegram ID;  "help": list what the bot can do;
  "reset": forget the current protein and conversation
- "chat": everything else - questions, explanations of earlier results, small talk, or when you need more
  information before running something. Put your full answer (max ~150 words) in "reply".

JSON fields:
{"action": "...", "protein": null, "mutations": [], "goal": null, "pdb_id": null, "chain": null, "method": null, "reply": "..."}
- protein: "SEQ1"/"SEQ2"... for a sequence the user pasted in THIS message, "CURRENT" for the protein we are
  already working on (the user may say "it", "this", "the same one"), or a UniProt accession such as P69905.
  If the user names a well-known protein and you are CERTAIN of its UniProt accession (human unless another
  organism is stated), you may use it and say which accession you chose in "reply". If unsure, use "chat" and ask.
  NEVER write out amino-acid sequences yourself.
- mutations: list like ["L43D", "K11R"] (original residue, 1-based position, new residue). If the user only
  gives a position and a new residue, write e.g. "43D" (the bot fills in the original residue).
- reply: for tool actions, one short friendly sentence saying what you are about to do. Never invent results;
  the tools produce them. For "chat", the answer itself.
- If something required is missing (no protein known, no mutation given, no PDB ID...), use "chat" and ask for it.
"""


def extract_sequences(text: str) -> tuple[str, dict[str, str]]:
    """Replace pasted sequences (raw or FASTA) by [SEQ1], [SEQ2]... Returns (masked text, {"SEQ1": seq}).."""
    seqs: dict[str, str] = {}
    text = re.sub(r"^\s*>(.*)$", lambda m: f"[FASTA header: {m.group(1).strip()[:80]}]", text, flags=re.M)

    def sub(m):
        key = f"SEQ{len(seqs) + 1}"
        seqs[key] = re.sub(r"\s+", "", m.group(0)).upper()
        return f"[{key}]"

    text = SEQ_UPPER_RE.sub(sub, text)
    text = SEQ_LOWER_RE.sub(sub, text)
    return text, seqs


def is_sequence(s: str | None) -> bool:
    return bool(s) and len(s) >= 10 and set(s) <= VALID_AA


def describe_subject(s: str | None) -> str:
    if not s:
        return "none yet"
    if is_sequence(s):
        return f"a pasted sequence of {len(s)} aa (starts {s[:12]}...)"
    return f"UniProt {s}"


def remember(ctx: ContextTypes.DEFAULT_TYPE, role: str, text: str):
    hist = ctx.user_data.setdefault("history", [])
    hist.append(f"{role}: {text[:500]}")
    del hist[:-HISTORY_TURNS]


async def llm_intent(update: Update, ctx: ContextTypes.DEFAULT_TYPE, masked: str, seqs: dict[str, str]) -> dict:
    methods = "\n".join(f"  {n}: {m['label']}" + ("" if method_available(n)[0] else " (NOT installed)")
                        for n, m in METHODS.items())
    pasted = ", ".join(f"{k} ({len(v)} aa)" for k, v in seqs.items()) or "none"
    recent = "\n---\n".join(_recent_output.get(update.effective_user.id, [])) or "none"
    prompt = (
        INTENT_PROMPT
        + f"\nAvailable method keys:\n{methods}\n"
        + f"User's current method: {ctx.user_data.get('method', DEFAULT_METHOD)}\n"
        + f"Current protein (CURRENT): {describe_subject(ctx.user_data.get('subject'))}\n"
        + f"Last PDB entry discussed: {ctx.user_data.get('last_pdb') or 'none'}\n"
        + f"Recent results shown to the user (truncated):\n{recent}\n\n"
        + "Conversation so far:\n" + ("\n".join(ctx.user_data.get("history", [])) or "(start of conversation)")
        + f"\n\nSequences pasted in this message: {pasted}\nUser: {masked}\n"
    )
    text = await llm_text(prompt, as_json=True, temperature=0.2)
    m = re.search(r"\{.*\}", text, re.S)
    if not m:
        raise ValueError(f"LLM did not return JSON: {text[:200]}")
    it = json.loads(m.group(0))
    if not isinstance(it, dict):
        raise ValueError("LLM JSON was not an object")
    return it


METHOD_ALIASES = [("nvidia", "esmfold_nvidia"), ("local esm", "esmfold_local"), ("esm", "esmfold_api"),
                  ("boltz", "boltz2_local"), ("colabfold", "colabfold_local"), ("alphafold2", "colabfold_local"),
                  ("omegafold", "omegafold_local"), ("alphafold 3", "alphafold3_local"), ("af3", "alphafold3_local"),
                  ("alphafold3", "alphafold3_local"), ("alphafold db", "afdb"), ("afdb", "afdb")]


def rule_intent(masked: str, seqs: dict[str, str], ctx: ContextTypes.DEFAULT_TYPE) -> dict:
    """Keyword fallback used when no LLM is configured or it fails."""
    t = masked.lower()
    acc = UNIPROT_WORD_RE.search(masked.upper())
    pdb_ids = [p for p in PDB_WORD_RE.findall(masked) if any(c.isalpha() for c in p)]
    plain = re.sub(r"\[SEQ\d+\]", "", masked)
    # "L43 to D" / "43 -> D" style first, then compact forms like L43D / 43D / Leu43Asp
    muts = [f"{a}{n}{b}".upper() for a, n, b in
            re.findall(r"\b([A-Za-z]?)(\d+)\s*(?:to|->|→|>)\s*([A-Za-z])\b", plain)]
    muts += [m for m in MUT_WORD_RE.findall(plain) if not m.isdigit()]
    protein = next(iter(seqs), None) or (acc.group(0) if acc else None)
    has = lambda *words: any(w in t for w in words)  # noqa: E731
    it = {"protein": protein or "CURRENT", "mutations": muts, "reply": ""}
    if has("reset", "start over", "forget"):
        return {"action": "reset"}
    if has("stab", "ddg", "δδg") and muts:
        return {**it, "action": "stability"}
    if muts or has("mutat", "mutant"):
        return {**it, "action": "mutate"}
    if has("design", "optimi", "improve"):
        goal = plain.strip() or "improve stability"
        return {**it, "action": "design", "goal": goal}
    if has("benchmark", "all models", "all the models", "every method", "compare method", "compare models"):
        return {**it, "action": "benchmark"}
    if pdb_ids and has("pdb", "experiment", "crystal", "rcsb", "compare"):
        return {**it, "action": "pdb", "pdb_id": pdb_ids[0], "protein": next(iter(seqs), None)}
    if has("analy", "propert", "disorder", "pi ", "molecular weight"):
        return {**it, "action": "analyze"}
    if has("method", "switch to", "use boltz", "use esm", "use colabfold", "use omegafold", "backend"):
        key = next((k for alias, k in METHOD_ALIASES if alias in t), None)
        return {"action": "method", "method": key}
    if has("status", "queue", "which methods"):
        return {"action": "status"}
    if has("whoami", "my id", "telegram id"):
        return {"action": "whoami"}
    if has("llm", "gemini", "ollama") and has("test", "check", "work"):
        return {"action": "llmtest"}
    if has("alphafold db", "afdb") and acc:
        return {**it, "action": "afdb"}
    if protein or has("fold", "structure", "predict", "model it"):
        return {**it, "action": "fold"}
    if has("help", "hello", "hi ", "hey", "what can you") or t.strip() in ("hi", "hey", "hello"):
        return {"action": "help"}
    llm_configured = bool(GEMINI_KEY or GROQ_KEY or OPENAI_KEY)
    status_note = (
        "The LLM service is temporarily experiencing high demand or unreachable, so I fell back to keyword mode."
        if llm_configured else
        "No LLM is configured right now, so I only understand simple keywords."
    )
    return {"action": "chat", "reply": (
        "I couldn't tell what you'd like me to do. Paste a sequence, give a UniProt ID, or say something like "
        f"\"mutate L43D\" or \"analyze it\". ({status_note} Send /llmtest to check.)")}


def normalize_mutations(raw, seq: str | None) -> list[str]:
    """Accept L43D, 43D (original residue filled in from the sequence) and Leu43Asp."""
    if isinstance(raw, str):
        raw = parse_mutations([raw])
    out = []
    for m in raw or []:
        m = str(m).strip()
        g3 = re.fullmatch(r"([A-Za-z]{3})(\d+)([A-Za-z]{3})", m)
        if g3 and g3.group(1).upper() in THREE_TO_ONE and g3.group(3).upper() in THREE_TO_ONE:
            m = THREE_TO_ONE[g3.group(1).upper()] + g3.group(2) + THREE_TO_ONE[g3.group(3).upper()]
        g = re.fullmatch(r"([A-Z]?)(\d+)([A-Z])", m.upper())
        if not g:
            raise ValueError(f"I couldn't read the mutation '{m}'. Use a form like L43D (original, position, new).")
        wt, pos, new = g.group(1), int(g.group(2)), g.group(3)
        if not wt:
            if not seq:
                raise ValueError(f"Which residue is at position {pos}? Please write it like L{pos}{new}.")
            if not 1 <= pos <= len(seq):
                raise ValueError(f"Position {pos} is outside the sequence (1-{len(seq)}).")
            wt = seq[pos - 1]
        out.append(f"{wt}{pos}{new}")
    return out


async def run_intent(update: Update, ctx: ContextTypes.DEFAULT_TYPE, it: dict, seqs: dict[str, str]):
    action = str(it.get("action") or "chat").lower().strip()
    reply = str(it.get("reply") or "").strip()
    say = update.message.reply_text

    # work out which protein the user means
    p = it.get("protein")
    p = str(p).strip().strip("[]").upper() if p else ""
    if p in seqs:
        protein = seqs[p]
    elif p in ("CURRENT", "IT", "THIS", "SAME"):
        protein = ctx.user_data.get("subject")
    elif UNIPROT_RE.fullmatch(p):
        protein = p
    elif is_sequence(clean_sequence(p)):
        protein = clean_sequence(p)
    else:
        protein = None
    if protein is None and action in PROTEIN_ACTIONS:
        protein = next(iter(seqs.values()), None) or ctx.user_data.get("subject")
    if action in PROTEIN_ACTIONS and not protein:
        msg = "Which protein should I use? Paste a sequence (or FASTA) or give me a UniProt ID like P69905."
        remember(ctx, "Bot", msg)
        await say(msg)
        return

    remember(ctx, "Bot", (reply or f"(running {action})") + f" [action: {action}]")
    if action == "chat":
        await reply_long(update, reply or "Sorry, I didn't catch that. Could you rephrase?")
        return
    if reply:
        await say(reply)

    try:
        if action in ("mutate", "stability"):
            raw = it.get("mutations")
            if not raw:
                await say("Which mutation(s)? For example: L43D, or \"K11R and E20Q\".")
                return
            seq = None
            if any(not re.match(r"[A-Za-z]", str(m).strip()) for m in (raw if isinstance(raw, list) else [raw])):
                seq = protein if is_sequence(protein) else (await sequence_from_arg(protein))[0]
            muts = normalize_mutations(raw, seq)
            ctx.user_data["subject"] = protein
            ctx.args = [protein, ",".join(muts)]
            await (mutate if action == "mutate" else stability)(update, ctx)
        elif action in ("fold", "afdb"):
            ctx.user_data["subject"] = protein
            if UNIPROT_RE.fullmatch(protein):
                ctx.args = [protein]
                await afdb(update, ctx)
            elif action == "afdb":
                await say("AlphaFold DB lookups need a UniProt ID (like P69905). Want me to fold the sequence instead?")
            else:
                ctx.args = [protein]
                await fold(update, ctx)
        elif action == "analyze":
            ctx.user_data["subject"] = protein
            ctx.args = [protein]
            await analyze(update, ctx)
        elif action in ("design", "benchmark"):
            seq = protein if is_sequence(protein) else (await sequence_from_arg(protein))[0]
            ctx.user_data["subject"] = protein
            if action == "design":
                goal = str(it.get("goal") or "improve predicted stability and confidence")
                ctx.args = [seq] + goal.split()
                await design(update, ctx)
            else:
                ctx.args = [seq]
                await benchmark(update, ctx)
        elif action == "pdb":
            pdb_id = str(it.get("pdb_id") or "").strip().upper() or ctx.user_data.get("last_pdb")
            if not pdb_id or not re.fullmatch(r"[0-9][A-Z0-9]{3}", pdb_id):
                await say("Which PDB entry? Give me a 4-character ID like 1UBQ.")
                return
            ctx.user_data["last_pdb"] = pdb_id
            args = [pdb_id]
            if it.get("chain") and len(str(it["chain"]).strip()) == 1:
                args.append(str(it["chain"]).strip().upper())
            if p in seqs:  # the user pasted their own sequence to compare with the entry
                args.append(seqs[p])
            ctx.args = args
            await pdb_cmd(update, ctx)
        elif action == "method":
            key = str(it.get("method") or "").strip()
            if key in METHODS:
                ok, why = method_available(key)
                if ok:
                    ctx.user_data["method"] = key
                    await say(f"Done, I'll use {METHODS[key]['label']} from now on.")
                else:
                    await say(f"{METHODS[key]['label']} isn't available: {why}")
            else:
                ctx.args = []
                await method_cmd(update, ctx)
        elif action == "status":
            await status(update, ctx)
        elif action == "llmtest":
            await llmtest(update, ctx)
        elif action == "whoami":
            await whoami(update, ctx)
        elif action == "reset":
            await reset(update, ctx)
        else:  # "help" or anything unknown
            await start(update, ctx)
    except ValueError as e:
        await say(f"❌ {e}")
    except Exception as e:  # noqa: BLE001
        log.exception("conversational action %s failed", action)
        await say(f"Something went wrong while running {action}: {str(e)[:300]}")


@restricted
async def chat(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    """Entry point for every plain (non-slash) text message."""
    text = update.message.text or ""
    masked, seqs = extract_sequences(text)
    bare = FASTA_TAG_RE.sub("", masked).strip()
    remember(ctx, "User", masked)
    if bare == "[SEQ1]":  # just a sequence/FASTA: fold it straight away, like before
        it = {"action": "fold", "protein": "SEQ1"}
    elif UNIPROT_RE.fullmatch(bare.upper()):  # just a UniProt ID
        it = {"action": "fold", "protein": bare.upper()}
    else:
        await ctx.bot.send_chat_action(update.effective_chat.id, ChatAction.TYPING)
        it = None
        if not await llm_ready():
            try:
                it = await llm_intent(update, ctx, masked, seqs)
            except Exception:  # noqa: BLE001
                log.exception("intent parsing failed, using keyword rules")
        if it is None:
            it = rule_intent(masked, seqs, ctx)
    log.info("intent: %s", {k: v for k, v in it.items() if k != "reply"})
    await run_intent(update, ctx, it, seqs)


@restricted
async def reset(update: Update, ctx: ContextTypes.DEFAULT_TYPE):
    for k in ("subject", "history", "last_pdb"):
        ctx.user_data.pop(k, None)
    _recent_output.pop(update.effective_user.id, None)
    await update.message.reply_text("Fresh start! I've forgotten the current protein and our conversation.")


async def post_init(app):
    for name, q in queues.items():
        q.start(METHODS[name]["workers"])
    if not ALLOWED:
        log.warning("ALLOWED_USER_IDS is empty: everyone will be refused. Use /whoami to find your ID.")


def main():
    t_req = HTTPXRequest(
        connection_pool_size=8,
        connect_timeout=30.0,
        read_timeout=60.0,
        write_timeout=60.0,
        media_write_timeout=120.0,
        pool_timeout=5.0,
    )
    app = ApplicationBuilder().token(TOKEN).request(t_req).post_init(post_init).build()
    app.add_handler(CommandHandler("start", start))
    app.add_handler(CommandHandler("help", start))
    app.add_handler(CommandHandler("whoami", whoami))
    app.add_handler(CommandHandler("status", status))
    app.add_handler(CommandHandler("llmtest", llmtest))
    app.add_handler(CommandHandler("method", method_cmd))
    app.add_handler(CallbackQueryHandler(method_cb, pattern=r"^method:"))
    app.add_handler(CommandHandler("fold", fold))
    app.add_handler(CommandHandler("afdb", afdb))
    app.add_handler(CommandHandler("mutate", mutate))
    app.add_handler(CommandHandler("design", design))
    app.add_handler(CommandHandler("benchmark", benchmark))
    app.add_handler(CommandHandler("stability", stability))
    app.add_handler(CommandHandler("pdb", pdb_cmd))
    app.add_handler(CommandHandler("analyze", analyze))
    app.add_handler(CommandHandler("reset", reset))
    app.add_handler(MessageHandler(filters.TEXT & ~filters.COMMAND, chat))
    app.run_polling()


if __name__ == "__main__":
    main()