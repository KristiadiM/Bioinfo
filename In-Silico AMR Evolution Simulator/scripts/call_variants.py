#!/usr/bin/env python3
"""
call_variants.py

A deliberately simple low-frequency variant caller, applied per generation
to the simulated read counts. This stands in for tools like LoFreq/VarScan
used to call resistance-conferring SNPs across longitudinal bacterial WGS
samples: for each generation it estimates the alt (resistance) allele
frequency, a Wilson confidence interval, and a binomial-test-based
significance call against the background sequencing error rate.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml
from scipy import stats


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    half = (z * ((p * (1 - p) / n + z**2 / (4 * n**2)) ** 0.5)) / denom
    return max(0.0, center - half), min(1.0, center + half)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reads", required=True, help="reads.tsv from generate_reads.py")
    ap.add_argument("--config", required=True)
    ap.add_argument("--outfile", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    call_cfg = cfg["calling"]
    err = cfg["sequencing"]["error_rate"]
    min_alt = call_cfg["min_alt_reads"]
    af_thresh = call_cfg["af_call_threshold"]

    reads = pd.read_csv(args.reads)

    rows = []
    for _, r in reads.iterrows():
        depth = int(r["depth"])
        alt = int(r["alt_count"])
        af_obs = alt / depth if depth > 0 else 0.0
        ci_low, ci_high = wilson_ci(alt, depth)

        # one-sided binomial test: is alt_count higher than expected from
        # sequencing error alone?
        pval = stats.binomtest(alt, depth, err, alternative="greater").pvalue if depth > 0 else 1.0

        variant_called = bool(
            (alt >= min_alt) and (af_obs >= af_thresh) and (pval < 0.01)
        )

        rows.append(
            {
                "condition": r["condition"],
                "replicate": r["replicate"],
                "generation": r["generation"],
                "true_freq": r["true_freq"],
                "depth": depth,
                "alt_count": alt,
                "af_obs": af_obs,
                "af_ci_low": ci_low,
                "af_ci_high": ci_high,
                "pval_vs_error": pval,
                "variant_called": variant_called,
            }
        )

    out = pd.DataFrame(rows)
    outpath = Path(args.outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(outpath, index=False)
    print(f"[call_variants] wrote {len(out)} calls -> {outpath}", file=sys.stderr)


if __name__ == "__main__":
    main()
