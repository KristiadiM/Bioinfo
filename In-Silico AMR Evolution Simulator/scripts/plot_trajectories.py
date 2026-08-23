#!/usr/bin/env python3
"""
plot_trajectories.py

Plots the sequencing-derived resistance allele frequency trajectory for
every simulated lineage, faceted by condition (drug vs no_drug), with the
cross-replicate mean overlaid and the fixation threshold marked.
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trajectories", nargs="+", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--outfile", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    fix_thresh = cfg["calling"]["fixation_threshold"]

    dfs = [pd.read_csv(f) for f in args.trajectories]
    data = pd.concat(dfs, ignore_index=True)

    conditions = sorted(data["condition"].unique())
    fig, axes = plt.subplots(1, len(conditions), figsize=(6 * len(conditions), 4.5), sharey=True)
    if len(conditions) == 1:
        axes = [axes]

    colors = {"drug": "#B3202C", "no_drug": "#2E6F95"}

    for ax, cond in zip(axes, conditions):
        sub = data[data["condition"] == cond]
        for rep, g in sub.groupby("replicate"):
            g = g.sort_values("generation")
            ax.plot(g["generation"], g["af_obs"], color=colors.get(cond, "gray"), alpha=0.35, lw=1.2)

        mean_traj = sub.groupby("generation")["af_obs"].mean()
        ax.plot(mean_traj.index, mean_traj.values, color="black", lw=2.2, label="mean across replicates")

        ax.axhline(fix_thresh, color="gray", ls="--", lw=1, label=f"fixation threshold ({fix_thresh})")
        ax.set_title(cond.replace("_", " "))
        ax.set_xlabel("Generation (passage)")
        ax.set_ylim(-0.02, 1.02)
        ax.legend(fontsize=8, loc="lower right")

    axes[0].set_ylabel("Observed resistance allele frequency")
    fig.suptitle("In-silico serial-passage evolution: resistance allele frequency across generations")
    fig.tight_layout()

    outpath = Path(args.outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=160)
    print(f"[plot_trajectories] wrote {outpath}", file=sys.stderr)


if __name__ == "__main__":
    main()
