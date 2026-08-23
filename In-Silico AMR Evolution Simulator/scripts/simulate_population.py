#!/usr/bin/env python3
"""
simulate_population.py

Simulates a serial-passage in-vitro evolution experiment for a single
bacterial lineage under (optional) antibiotic selection at one resistance
locus. Each passage ("generation") consists of:

  1. Drug exposure / kill step
       - resistant cells survive antibiotic exposure
       - susceptible cells mostly die (kill_prob_susceptible), EXCEPT for a
         small phenotypically tolerant/persister subpopulation that survives
         at a higher rate independent of genotype (the "non-genomic" axis
         of resistance evolution: transient phenotypic tolerance, not an
         inherited mutation).
  2. Mutation step
       - susceptible survivors can acquire the resistance allele de novo
         (mutation_rate)
       - resistant survivors can revert (reversion_rate)
  3. Regrowth to carrying capacity K
       - modeled as binomial resampling from the post-kill/mutation
         bottleneck, which is the source of genetic drift (effective
         population size ~= bottleneck size, standard for serial-passage /
         Wright-Fisher-with-bottleneck models)
  4. Serial dilution back down to the founder size (N0) for the next passage

The true (population-level) resistance allele frequency is recorded at the
post-regrowth timepoint of every generation -- this is the state that would
be sampled if you pulled DNA and sequenced the culture at the end of each
passage, before diluting into fresh medium/host for the next passage.
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


def simulate_lineage(cond_cfg, pop_cfg, tol_cfg, seed):
    rng = np.random.default_rng(seed)

    N0 = pop_cfg["N0"]
    K = pop_cfg["K"]
    generations = pop_cfg["generations"]
    mu = pop_cfg["mutation_rate"]
    nu = pop_cfg["reversion_rate"]

    kill_prob_sus = cond_cfg["kill_prob_susceptible"]
    tol_frac = tol_cfg["tolerant_fraction"]
    tol_surv = tol_cfg["tolerant_survival_prob"]

    # founder population starts fully susceptible
    n_res = 0
    n_sus = N0

    records = []
    for gen in range(1, generations + 1):
        # ---- 1. drug exposure / kill step ----
        # resistant cells are unaffected by the antibiotic
        survivors_res = n_res

        if n_sus > 0:
            n_tol = rng.binomial(n_sus, tol_frac)
            n_normal_sus = n_sus - n_tol
            survivors_tol = rng.binomial(n_tol, tol_surv)
            survivors_normal = rng.binomial(n_normal_sus, 1.0 - kill_prob_sus)
            survivors_sus = survivors_tol + survivors_normal
        else:
            survivors_sus = 0

        # ---- 2. mutation step ----
        new_res_from_mut = rng.binomial(survivors_sus, mu) if survivors_sus > 0 else 0
        new_sus_from_rev = rng.binomial(survivors_res, nu) if survivors_res > 0 else 0

        res_bottleneck = survivors_res - new_sus_from_rev + new_res_from_mut
        sus_bottleneck = survivors_sus - new_res_from_mut + new_sus_from_rev
        res_bottleneck = max(res_bottleneck, 0)
        sus_bottleneck = max(sus_bottleneck, 0)
        bottleneck_total = res_bottleneck + sus_bottleneck

        # ---- 3. regrowth to carrying capacity (drift step) ----
        if bottleneck_total == 0:
            # lineage went extinct
            n_res_grown, n_sus_grown = 0, 0
        else:
            p_res = res_bottleneck / bottleneck_total
            n_res_grown = rng.binomial(K, p_res)
            n_sus_grown = K - n_res_grown

        N_total = n_res_grown + n_sus_grown
        true_freq = (n_res_grown / N_total) if N_total > 0 else np.nan

        records.append(
            {
                "generation": gen,
                "bottleneck_size": bottleneck_total,
                "n_res": n_res_grown,
                "n_sus": n_sus_grown,
                "N_total": N_total,
                "true_freq": true_freq,
            }
        )

        # ---- 4. serial dilution: founder for next passage ----
        if N_total == 0:
            n_res, n_sus = 0, 0
        else:
            p_res = n_res_grown / N_total
            n_res = rng.binomial(N0, p_res)
            n_sus = N0 - n_res

    return pd.DataFrame.from_records(records)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--condition", required=True, choices=["drug", "no_drug"])
    ap.add_argument("--replicate", required=True, type=int)
    ap.add_argument("--outfile", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    cond_cfg = cfg["conditions"][args.condition]
    pop_cfg = cfg["population"]
    tol_cfg = cfg["non_genomic"]

    base_seed = cfg.get("random_seed", 0)
    # distinct, reproducible seed per (condition, replicate)
    seed = base_seed + hash((args.condition, args.replicate)) % 1_000_000

    df = simulate_lineage(cond_cfg, pop_cfg, tol_cfg, seed)
    df.insert(0, "replicate", args.replicate)
    df.insert(0, "condition", args.condition)

    outpath = Path(args.outfile)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(outpath, index=False)
    print(f"[simulate_population] wrote {len(df)} generations -> {outpath}", file=sys.stderr)


if __name__ == "__main__":
    main()
