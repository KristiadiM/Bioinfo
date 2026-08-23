#!/usr/bin/env python3
"""
generate_reads.py

Simulates targeted/amplicon pooled sequencing of the resistance locus at
each generation of a serial-passage lineage. Takes the *true* population
allele frequency (from simulate_population.py) and draws finite-depth,
error-prone read counts from it -- exactly the noisy observation a real
variant-calling pipeline would work with, rather than the ground truth
itself.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sim", required=True, help="generations.csv from simulate_population.py")
    ap.add_argument("--config", required=True)
    ap.add_argument("--outfile", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    seq_cfg = cfg["sequencing"]
    depth = seq_cfg["depth"]
    err = seq_cfg["error_rate"]

    sim = pd.read_csv(args.sim)
    base_seed = cfg.get("random_seed", 0)
    condition = sim["condition"].iloc[0]
    replicate = int(sim["replicate"].iloc[0])
    rng = np.random.default_rng(base_seed + 500_000 + hash((condition, replicate, "reads")) % 1_000_000)

    rows = []
    for _, r in sim.iterrows():
        true_af = r["true_freq"]
        if pd.isna(true_af):
            alt_count, ref_count = 0, 0
        else:
            # probability a read reports "alt" = true carrier frequency,
            # adjusted for symmetric sequencing error in both directions
            p_alt_read = true_af * (1 - err) + (1 - true_af) * err
            alt_count = rng.binomial(depth, p_alt_read)
            ref_count = depth - alt_count

        rows.append(
            {
                "condition": r["condition"],
                "replicate": r["replicate"],
                "generation": r["generation"],
                "true_freq": true_af,
                "depth": depth,
                "ref_count": ref_count,
                "alt_count": alt_count,
            }
        )

    out = pd.DataFrame(rows)
    outpath = Path(args.outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(outpath, index=False)
    print(f"[generate_reads] wrote {len(out)} generation samples -> {outpath}", file=sys.stderr)


if __name__ == "__main__":
    main()
