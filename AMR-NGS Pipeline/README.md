# amr-ngs-pipeline

A small, complete, **runnable** Snakemake pipeline that takes bacterial
whole-genome sequencing reads from raw FASTQ through to antimicrobial
resistance (AMR) gene and variant reporting.

Built as a portfolio piece for a PhD application in bioinformatics
(Leibniz ScienceCampus **AMR-PLAS**, Research Center Borstel), which
asks for exactly this combination: basic NGS analysis, Snakemake/Nextflow,
Linux, Python, and microbial genomics in an AMR context. Rather than
describe those skills in a cover letter, this repo demonstrates them.

```
raw reads ─▶ FastQC ─▶ fastp (trim) ─▶ FastQC ─┬─▶ SPAdes assembly ─▶ QUAST
                                                 │                     │
                                                 │                     ▼
                                                 │              BLAST vs AMR
                                                 │              marker genes
                                                 │              (acquired resistance)
                                                 │
                                                 └─▶ BWA-MEM ─▶ bcftools
                                                     (map to reference) variant calling
                                                     (chromosomal target-mutation
                                                      resistance)
                                                                       │
                                              all QC/assembly logs ────┴──▶ MultiQC report
```

## Why two detection strategies?

This is the one design choice worth reading before the rest: real AMR
genomics almost always needs **both** an assembly-based and a
mapping-based route, because the two resistance mechanisms look
completely different at the sequence level:

- **Horizontally acquired resistance genes** (e.g. plasmid-borne
  beta-lactamases, efflux pump genes) are often simply *absent* from
  whatever reference genome you map to. Mapping-based variant calling
  can't detect the presence of DNA that isn't in the reference at all —
  you have to assemble the reads and search the assembly against a
  curated gene database (ABRicate / AMRFinderPlus / RGI against
  CARD/ResFinder in production; a small local BLAST database here).
- **Chromosomal target-mutation resistance** (e.g. gyrA/rpoB-type SNPs
  behind fluoroquinolone/rifampicin resistance) is a point mutation *on*
  the shared backbone. It's invisible to gene-presence screening (the
  gene is present either way) and is exactly what mapping + variant
  calling is built to find.

The synthetic test data in this repo (see below) is deliberately
constructed so each method catches something the other one misses —
you can see this directly in `example_results/`.

## Repo layout

```
Snakefile                  the pipeline: every rule, in dependency order
config/
  config.yaml              paths, tool parameters, thresholds
  samples.tsv              sample name -> R1/R2 FASTQ paths
envs/                      one conda env per pipeline stage (bioconda)
scripts/
  build_synthetic_genomes.py   builds the offline demo dataset (see below)
  simulate_reads.py            paired-end read simulator w/ sequencing errors
  filter_amr_hits.py           BLAST hits -> presence/absence report
resources/
  reference_genome.fasta       "reference strain" reads are mapped to
  amr_marker_db.fasta          local stand-in for CARD/ResFinder
  sim_genomes/                 truth genomes used only to simulate reads
example_results/            small curated outputs from a real run, so a
                             reviewer can see it works without running it
docs/
  rulegraph.png                rule-level DAG (this README's diagram, rendered)
  dag_full.png                 full per-sample job DAG
```

## Running it

```bash
conda install -n base -c conda-forge mamba   # once
snakemake --use-conda -c4                    # full run, 4 cores
snakemake -n                                 # dry run / inspect the DAG
snakemake --dag | dot -Tpng > dag.png        # regenerate the diagram
```

No conda? Every tool here is also installable via `apt` on Ubuntu
(`fastqc bwa samtools bcftools spades ncbi-blast+`) or `pip`
(`fastp` is apt-installable too; `quast`, `multiqc`, `snakemake` are on
PyPI) — that's how this repo was actually developed and tested.

### Using real data

Swap the demo for a real dataset in two edits:

1. **`config/samples.tsv`** — point `r1`/`r2` at real `fastq.gz` files
   (e.g. downloaded with `prefetch`+`fasterq-dump` from an SRA/ENA run
   accession for an ESKAPE or *M. tuberculosis* isolate).
2. **`config/config.yaml`** — remove the sample from `demo_samples`
   (or just leave it; the simulate rule only fires when its output
   files don't already exist, so real files simply take priority),
   point `reference` at a real reference genome (e.g. an NCBI RefSeq
   assembly for the same species), and point `amr_db` at a real curated
   database, e.g.:
   ```bash
   # ResFinder DB
   git clone https://bitbucket.org/genomicepidemiology/resfinder_db
   cat resfinder_db/*.fsa > resources/amr_marker_db.fasta
   ```
   For production use I'd actually swap the custom BLAST step for
   **ABRicate** or **AMRFinderPlus** run against CARD — they add curated
   resistance-mechanism annotation (not just sequence identity) that a
   raw BLAST search doesn't give you. The custom script here exists so
   the pipeline runs with nothing but BLAST+, useful for this demo but I would use a maintained tool for real work.

## The synthetic demo dataset

There's no internet access to SRA/ENA in the sandbox this was built in,
so `scripts/build_synthetic_genomes.py` constructs a small (~12-14 kb)
synthetic genome pair instead of using a real bacterial genome — this
is clearly labeled everywhere it's used and is **not** real biological
data. It exists purely so the pipeline can be smoke-tested end-to-end,
offline, in seconds, which is also just good practice for CI testing of
a real pipeline.

Two "isolates" are simulated:

| | `isolate01_resistant` | `isolate02_susceptible` |
|---|---|---|
| Shared chromosome backbone | ✓ | ✓ |
| `blaSYN` (beta-lactamase-like) | **present** | absent |
| `mefSYN` (efflux-pump-like) | **present** | absent |
| `gyrSYN` target gene | mutant allele (1 SNP) | wild-type allele |

Reads are simulated at 40x coverage with a 1% per-base error rate
(`scripts/simulate_reads.py`) from each isolate's true genome, then fed
into the pipeline exactly as if they'd come off a sequencer.

### What the pipeline actually recovered (`example_results/`)

AMR gene screening (`combined_amr_report.tsv`) — assembly + BLAST:

| sample | blaSYN | mefSYN | vanSYN (decoy, never present) |
|---|---|---|---|
| isolate01_resistant | PRESENT (100% id / 100% cov) | PRESENT (100% id / 100% cov) | absent |
| isolate02_susceptible | absent | absent | absent |

Variant calling (`isolate0{1,2}_*.filtered.vcf`) — mapping to reference:

- `isolate01_resistant`: 2 adjacent SNPs called inside `gyrSYN`, at the
  exact codon position the simulator mutated — QUAL 225, DP 37×,
  homozygous alt.
- `isolate02_susceptible`: **zero** variants, as expected for a strain
  identical to the reference.

That's the two-methods point made concretely: the acquired genes are
invisible to variant calling (they're not in the reference), and the
target mutation is invisible to gene-presence screening (the gene's
there in both samples) — each isolate's resistance profile only comes
together once both branches of the pipeline are combined, which is
what `combine_amr_reports` + the per-sample VCFs do.

## Design notes / what I'd change for a real production version

- **Reproducible environments**: each rule has a dedicated conda env
  in `envs/` pinned to bioconda package versions, following current
  Snakemake best practice, so results — not just code — are
  reproducible with `--use-conda`.
- **Sample scaling**: adding a sample is a one-line edit to
  `samples.tsv`; Snakemake's wildcard-based rules mean every stage
  scales to N samples without touching pipeline logic — the same
  reason the AMR-PLAS group's HPC resources are useful here (this
  would trivially move to a Slurm/HPC executor via a
  `--profile` config for a cohort of real isolates).
- **AMR database**: swap the local BLAST search for ABRicate/
  AMRFinderPlus + CARD/ResFinder/NCBI's curated database, which encode
  resistance mechanism and confidence, not just % identity.
- **Variant calling**: for real mycobacterial work I'd add
  lineage-aware calling and a curated resistance-SNP catalogue (e.g.
  the WHO TB mutation catalogue) rather than generic `bcftools call`.
- **Testing**: the synthetic dataset here is effectively a minimal CI
  fixture — I'd keep it in the repo as a `pytest`-driven
  `snakemake --use-conda -n` / small end-to-end smoke test even after
  wiring in real data, so regressions in the pipeline logic are caught
  without needing a real dataset on every CI run.

## Why Snakemake over Nextflow here

Both are reasonable choices and the position lists either. I picked
Snakemake because its rules are plain Python-plus-shell, which sits
closer to the "Python/R + Linux command line" profile the position
describes than Nextflow's Groovy-based DSL — and because Snakemake's
wildcard system made the "same pipeline, two detection strategies,
scales to N samples" structure above very direct to express. I'd be
equally happy picking up Nextflow/nf-core on the job, since the
underlying ideas (DAG of rules/processes, per-step environments,
resumable execution) transfer directly.
