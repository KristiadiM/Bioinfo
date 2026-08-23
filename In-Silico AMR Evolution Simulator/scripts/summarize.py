#!/usr/bin/env python3
"""
summarize.py

Combines the per-lineage fixation summaries into condition-level statistics:
fixation probability (fraction of replicates that reached the fixation
threshold within the observed passages) and, among lineages that fixed, the
mean/median/std fixation time (in generations).
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fixation-files", nargs="+", required=True)
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    dfs = [pd.read_csv(f) for f in args.fixation_files]
    data = pd.concat(dfs, ignore_index=True)

    rows = []
    for cond, g in data.groupby("condition"):
        n = len(g)
        n_fixed = int(g["fixed"].sum())
        fixed_times = g.loc[g["fixed"], "fixation_generation"]
        rows.append(
            {
                "condition": cond,
                "n_replicates": n,
                "n_fixed": n_fixed,
                "fixation_probability": n_fixed / n,
                "mean_fixation_generation": fixed_times.mean() if n_fixed else None,
                "median_fixation_generation": fixed_times.median() if n_fixed else None,
                "std_fixation_generation": fixed_times.std() if n_fixed else None,
                "mean_final_af": g["final_af_obs"].mean(),
            }
        )

    summary = pd.DataFrame(rows).sort_values("condition")
    outpath = Path(args.outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outpath, index=False)

    lines = ["# In-silico evolution experiment: summary report", ""]
    lines.append("| Condition | Replicates | Fixed | Fixation probability | Mean fixation gen | Median fixation gen |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in summary.iterrows():
        mean_g = f"{r['mean_fixation_generation']:.1f}" if pd.notna(r["mean_fixation_generation"]) else "n/a"
        med_g = f"{r['median_fixation_generation']:.1f}" if pd.notna(r["median_fixation_generation"]) else "n/a"
        lines.append(
            f"| {r['condition']} | {int(r['n_replicates'])} | {int(r['n_fixed'])} | "
            f"{r['fixation_probability']:.2f} | {mean_g} | {med_g} |"
        )
    lines.append("")
    lines.append(
        "Fixation is called on the *sequencing-derived* allele frequency "
        "(post variant-calling), not the underlying true simulated frequency, "
        "so it reflects what a real genomic surveillance pipeline would observe."
    )

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines))

    print(f"[summarize] wrote {outpath} and {report_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
