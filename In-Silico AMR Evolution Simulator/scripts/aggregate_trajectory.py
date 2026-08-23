#!/usr/bin/env python3
"""
aggregate_trajectory.py

Takes the per-generation variant calls for one lineage and derives the
allele frequency trajectory plus the fixation time: the first generation at
which the *called* (sequencing-derived) allele frequency reaches the
fixation threshold defined in the config. If the lineage never reaches
fixation within the observed generations, fixation_generation is left null.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calls", required=True, help="calls.tsv from call_variants.py")
    ap.add_argument("--config", required=True)
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--fixation-outfile", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    fix_thresh = cfg["calling"]["fixation_threshold"]

    calls = pd.read_csv(args.calls).sort_values("generation")
    outpath = Path(args.outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    calls.to_csv(outpath, index=False)  # trajectory == annotated calls, kept as its own artifact

    fixed = calls[calls["af_obs"] >= fix_thresh]
    fixation_gen = int(fixed["generation"].iloc[0]) if not fixed.empty else None

    condition = calls["condition"].iloc[0]
    replicate = int(calls["replicate"].iloc[0])
    final_af = float(calls["af_obs"].iloc[-1])

    summary = pd.DataFrame(
        [
            {
                "condition": condition,
                "replicate": replicate,
                "fixation_generation": fixation_gen,
                "fixed": fixation_gen is not None,
                "final_af_obs": final_af,
            }
        ]
    )
    fixation_path = Path(args.fixation_outfile)
    fixation_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(fixation_path, index=False)

    print(
        f"[aggregate_trajectory] {condition} rep{replicate}: "
        f"fixation_generation={fixation_gen} final_af={final_af:.3f} -> {outpath}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
