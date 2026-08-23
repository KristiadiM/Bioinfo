# In-Silico AMR Evolution Simulator

A small simulation and Snakemake analysis pipeline modeling how a bacterial population acquires and fixes an antibiotic-resistance mutation under serial drug exposure.
Built as a project pitch for the **AMR-PLAS ScienceCampus**, it demonstrates both **genomic** and **non-genomic mechanisms of antimicrobial resistance evolution** in a single longitudinal genomics workflow.

**1. Simulates evolution, not just genetics.** Each of two arms — an
antibiotic-exposed lineage and an untreated control — is propagated through
serial passages (`scripts/simulate_population.py`). Every passage:

1. **Drug exposure / kill step** — resistant cells survive; susceptible cells mostly die, *except* for a small phenotypically tolerant/persister subpopulation that survives at an elevated rate independent of genotype. This is the **non-genomic** axis: transient phenotypic tolerance that buffers a lineage from extinction without changing allele frequency.
2. **Mutation step** — susceptible cells can acquire the resistance allele de novo; resistant cells can revert. This is the **genomic** axis.
3. **Regrowth to carrying capacity** via binomial resampling from the post-kill bottleneck — the source of genetic drift (effective population size ≈ bottleneck size, standard for serial-passage / Wright–Fisher models with bottlenecks).
4. **Serial dilution** back to founder size for the next passage, exactly as in a real serial-passage evolution assay.

**2. Runs it through a pipeline shaped like a real AMR genomics workflow.**
Rather than reading the "true" simulated allele frequency directly, the Snakemake pipeline treats each generation like a sequenced isolate/timepoint.

## Features

* Models **both axes** of resistance evolution (genomic mutations and phenotypic tolerance)
* Simulates true serial-passage population dynamics (bottlenecks, genetic drift, selection)
* Mimics a real AMR surveillance workflow using simulated pooled/amplicon sequencing reads
* Implements a low-frequency variant caller (Wilson CI + binomial test against error rates)
* Calculates allele frequency trajectories and **fixation time** automatically
* Highly parameterized and reproducible using Snakemake

## Requirements

```bash
pip install -r requirements.txt
# Alternatively:
conda env create -f environment.yml

```

## How to Run

1. Clone or download the repository.
2. Open terminal in the project directory.
3. Install the dependencies using pip or conda.
4. Run the full end-to-end pipeline:
```bash
snakemake --cores 4

```

## This pipeline will:

1. Simulate bacterial populations undergoing serial drug exposure alongside untreated controls.
2. Generate sequencing reads based on finite-depth, error-prone sampling of true frequencies.
3. Call variants against the background error rate to mimic longitudinal WGS calling.
4. Aggregate point-in-time calls into allele frequency trajectories.
5. Plot fixation probability and time between the drug and no-drug arms.
### The Biological Model

* **Genomic Axis:** Includes de novo mutation, reversion, genetic drift, and selection. Susceptible cells acquire a resistance allele over time.
* **Non-genomic Axis:** A transient, phenotypically tolerant subpopulation survives drug exposure at an elevated rate independent of genotype, buffering the lineage against immediate extinction.

### Output Files

Running the pipeline populates the `results/` directory, including:
* `results/plots/fixation_trajectories.png` showing the mutation sweeping and fixing.
* `results/report/summary.md` detailing cross-lineage fixation statistics.

## Configuration Notes

You can easily explore the parameter space by editing `config/config.yaml` without touching the code:

* **Simulation Dynamics**: Adjust carrying capacity, standard effective population size, and bottleneck severity.
* **Evolutionary Forces**: Modify the selection coefficient, mutation/reversion rates, and the fraction of the population exhibiting non-genomic tolerance.
* **Sequencing Metrics**: Alter sequencing depth, background error rates, and variant-calling thresholds.

## Contributing

Natural extensions are welcome! Feel free to open issues or pull requests for:

* Implementing multiple resistance loci with epistasis or a complex fitness landscape.
* Adding population structure or horizontal transfer of resistance (e.g., plasmid-borne).
* Replacing simulated reads with real longitudinal WGS from serial clinical or in-vitro isolates.
* Swapping the native python script for established callers like LoFreq or BCFtools inside the Snakemake DAG.

📬 Contact
For questions or suggestions, reach out at:
📧 kristiadimikael@gmail.com
