# AMR NGS Pipeline

A small, complete, **runnable** Snakemake pipeline that takes bacterial whole-genome sequencing reads from raw FASTQ through to antimicrobial resistance (AMR) gene and variant reporting.
It demonstrates basic NGS analysis, reproducible environments, and microbial genomics in an AMR context using two distinct detection strategies.

## Features

* Process raw reads through assembly and mapping
* Quality control and read trimming using FastQC and fastp
* Assembly-based detection for **horizontally acquired resistance genes** (SPAdes + BLAST)
* Mapping-based variant calling for **chromosomal target-mutation resistance** (BWA-MEM + bcftools)
* Scalable to N samples using Snakemake wildcards
* Fully reproducible environments via Bioconda
* Includes a synthetic offline dataset for immediate testing

## Requirements

```bash
conda install -n base -c conda-forge mamba

```

Alternatively, tools can be installed via `apt` or `pip`: `fastqc`, `fastp`, `bwa`, `samtools`, `bcftools`, `spades`, `ncbi-blast+`, `quast`, `multiqc`, `snakemake`.

## How to Run

1. Clone or download the repository.
2. Open terminal in the project directory.
3. Perform a dry run to inspect the DAG:
```bash
snakemake -n

```

4. Run the full pipeline (using 4 cores):
```bash
snakemake --use-conda -c4

```

5. (Optional) Generate a visual DAG diagram:
```bash
snakemake --dag | dot -Tpng > dag.png

```

## This pipeline will:

1. Perform quality control and trim raw reads
2. Assemble reads and evaluate the assembly quality
3. Identify horizontally acquired resistance genes via a local database
4. Map reads to a reference genome and call variants to find point mutations
5. Combine outputs into a final MultiQC report and a comprehensive AMR report
   
### Detection Strategies

Real AMR genomics requires **both** methods. After running, the pipeline highlights:
* **Horizontally acquired genes:** Detected via assembly-based screening (captures genes absent from the reference genome).
* **Chromosomal mutations:** Detected via mapping-based variant calling (captures point mutations on the shared backbone).


### Output Files

The pipeline generates results in the `example_results/` directory, including:
* `combined_amr_report.tsv` for gene presence/absence.
* VCF files (`*.filtered.vcf`) for identified point mutations.
* Aggregated QC logs in a MultiQC report.

## Configuration Notes

You can customize:

* **Samples**: Add real dataset FASTQ paths in `config/samples.tsv`.
* **Reference Genome**: Update `reference` path in `config/config.yaml` to point to a real reference (e.g., NCBI RefSeq).
* **AMR Database**: Swap the local BLAST DB in `config/config.yaml` with a curated database like **CARD** or **ResFinder**.
* **Execution Profile**: Easily adapt for HPC environments using Slurm by configuring Snakemake profiles.

## Contributing

Contributions are welcome! Feel free to open issues or pull requests for:

* Integrating **ABRicate** or **AMRFinderPlus** for curated resistance annotation
* Adding lineage-aware variant calling for specific pathogens (e.g., *M. tuberculosis*)
* Enhancing pipeline testing using `pytest`
* Nextflow/nf-core adaptation

📬 Contact
For questions or suggestions, reach out at:
📧 kristiadimikael@gmail.com
