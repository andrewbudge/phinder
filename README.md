# phinder

[![CI](https://github.com/andrewbudge/phinder/actions/workflows/ci.yml/badge.svg)](https://github.com/andrewbudge/phinder/actions/workflows/ci.yml)

A Nextflow pipeline for phage discovery from metagenomic assemblies.

Takes a combined contig FASTA and runs viral identification, quality assessment,
annotation, and classification end-to-end.

```
contigs.fasta
    │
    ├─ geNomad          viral classification + provirus detection
    ├─ R filter         DTR topology always kept; Provirus filtered by score
    ├─ CheckV           completeness + quality estimation
    ├─ R filter         quality-based filtering with DTR bypass
    ├─ Pharokka         phage genome annotation (MMseqs2 + PyHMMER)
    │
    └─ PhaBOX*          taxonomy, lifestyle, host, vOTU clustering, phylogenetic tree
                        (*optional — requires separate installation)
```

---

## Requirements

- [Nextflow](https://nextflow.io/) >= 23.04
- One of:
  - [conda](https://docs.conda.io/) / [mamba](https://mamba.readthedocs.io/), or
  - [Docker](https://www.docker.com/) (laptop/workstation), or
  - [Singularity / Apptainer](https://apptainer.org/) (HPC)

---

## Quick start

phinder needs two things: the reference **databases** (downloaded once) and a way
to **provision the tools** (a `-profile` — containers or conda).

**1. Download the databases**

```bash
curl -O https://raw.githubusercontent.com/andrewbudge/phinder/main/setup.sh
bash setup.sh --skip-envs        # databases only → ~/.phinder_dbs
```

Drop `--skip-envs` to *also* build pinned conda environments (`genomad_phinder`,
`checkv_phinder`, `pharokka_phinder`) — only needed for the conda profile below.
Add `--with-phabox2` to set up the optional PhaBOX step.

**2. Run**

With **Docker** — pinned images, nothing to build (recommended):

```bash
nextflow run andrewbudge/phinder \
    --input contigs.fasta \
    --genomad_db  ~/.phinder_dbs/genomad_db \
    --checkv_db   ~/.phinder_dbs/checkv_db \
    --pharokka_db ~/.phinder_dbs/pharokka_db \
    -profile docker
```

With **conda** — uses the envs built by `setup.sh` (run it without `--skip-envs`):

```bash
nextflow run andrewbudge/phinder \
    --input contigs.fasta \
    --genomad_db  ~/.phinder_dbs/genomad_db \
    --checkv_db   ~/.phinder_dbs/checkv_db \
    --pharokka_db ~/.phinder_dbs/pharokka_db \
    --genomad_env  $(conda info --base)/envs/genomad_phinder \
    --checkv_env   $(conda info --base)/envs/checkv_phinder \
    --pharokka_env $(conda info --base)/envs/pharokka_phinder \
    -profile conda
```

Use `-resume` on reruns to skip completed steps. For HPC / Singularity / Apptainer,
see [Execution profiles](#execution-profiles).

---

## Verify your install

Before running on real data, confirm Nextflow and your tool profile are wired up
correctly — **no databases required**. This runs the whole pipeline on a tiny
bundled dataset in *stub* mode, where every step emits placeholder outputs in
seconds:

```bash
nextflow run andrewbudge/phinder -profile test -stub-run
```

A `[SUCCESS]` line with all 8 processes completed means your setup is good. To
also check that your engine pulls images / builds envs, add it to the profile:

```bash
nextflow run andrewbudge/phinder -profile test,docker -stub-run   # or test,conda
```

This is the same check phinder's [CI](https://github.com/andrewbudge/phinder/actions/workflows/ci.yml)
runs on every change.

---

## Input

A single combined contig FASTA — uncompressed or gzipped. If working with
multiple samples, merge and rename headers before passing to phinder to avoid
ID collisions:

```bash
# Example: prefix each sample's contigs with its ID
zcat sample1_contigs.fasta.gz | sed 's/>/>sample1_/' >> combined.fasta
zcat sample2_contigs.fasta.gz | sed 's/>/>sample2_/' >> combined.fasta
```

Coverage filtering, circular contig extraction, and other pre-processing are
left to the user — phinder is assembler-agnostic and makes no assumptions about
header format. If your assembler embeds `coverage=N` in headers (e.g. metaMDBG),
phinder will parse it and flag low-coverage candidates.

---

## Output

```
results/
├── genomad/
│   ├── filtered_genomad.tsv       filtered geNomad hits
│   └── output/                    full geNomad output
├── checkv/
│   ├── potential_phage.tsv        final candidate table
│   └── output/                    full CheckV output
├── candidates/
│   └── candidate_phages.fna       candidate phage sequences
├── pharokka/
│   └── output/                    per-contig annotation files
└── phabox/                        (if --phabox2_env provided)
    ├── end_to_end/                 taxonomy + lifestyle + host predictions
    ├── votu/                       AAI-based vOTU clusters
    └── tree/                       phylogenetic tree (terl + portal markers)
```

---

## All parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Combined contig FASTA (.fa or .fa.gz) |
| `--outdir` | `results` | Output directory |
| `--threads` | `8` | Threads per process |
| `--genomad_db` | required | Path to geNomad database |
| `--checkv_db` | required | Path to CheckV database |
| `--pharokka_db` | required | Path to Pharokka database |
| `--phabox_db` | required if `--phabox2_env` set | Path to PhaBOX database |
| `--genomad_env` | builds from `envs/genomad.yml` | Path to existing geNomad conda env |
| `--checkv_env` | builds from `envs/checkv.yml` | Path to existing CheckV conda env |
| `--pharokka_env` | builds from `envs/pharokka.yml` | Path to existing Pharokka conda env |
| `--phabox2_env` | unset (PhaBOX skipped) | Path to existing phabox2 conda env |
| `--min_provirus_score` | `0.9` | geNomad Provirus minimum virus_score |
| `--checkv_quality_keep` | `High-quality,Complete` | Comma-separated CheckV quality tiers to keep |
| `--min_coverage` | `3` | Coverage threshold for low-coverage flagging |
| `--genomad_splits` | `20` | geNomad MMseqs2 database splits (lower = faster, higher RAM) |
| `--pharokka_gene_predictor` | `prodigal-gv` | Gene predictor for Pharokka |
| `--phabox_skip_phamer` | `true` | Skip Phamer (sequences already confirmed viral) |
| `--phabox_votu_mode` | `AAI` | vOTU clustering mode (AAI or ANI) |
| `--phabox_tree_markers` | `terl,portal` | Marker genes for phylogenetic tree |

---

## Filtering logic

**geNomad filter:**
- DTR topology → always kept (hallmark of complete linear phage genome)
- Provirus topology → kept if `virus_score >= --min_provirus_score`
- All other topologies → dropped

**CheckV filter:**
- Kept if `checkv_quality` in `--checkv_quality_keep`
- DTR contigs bypass CheckV quality thresholds (CheckV under-scores DTRs
  due to absent host flanking regions)
- Low coverage is flagged (`low_coverage = TRUE`) but not dropped

---

## Execution profiles

Pick how tools are provisioned with `-profile`. Profiles are composable
(comma-separated):

| Profile | Tools provided by | Use when |
|---------|-------------------|----------|
| `conda` | conda envs built from `envs/*.yml` | local conda/mamba install |
| `mamba` | same, resolved with mamba | faster conda solves |
| `docker` | pinned biocontainer images | laptop / workstation |
| `singularity` | same images, via Singularity | HPC without root |
| `apptainer` | same images, via Apptainer | HPC without root |
| `slurm` | (executor only) | submit processes as SLURM jobs |

Containers pull pinned, frozen tool images — no conda solve, identical on every
machine. The reference **databases are still downloaded separately** (via
`setup.sh`) regardless of profile, and passed with the `--*_db` flags; Nextflow
mounts them into the container automatically.

```bash
# Laptop, with Docker
nextflow run andrewbudge/phinder --input contigs.fasta \
    --genomad_db ... --checkv_db ... --pharokka_db ... \
    -profile docker

# HPC, Singularity images submitted as SLURM jobs
nextflow run andrewbudge/phinder ... -profile singularity,slurm
```

> **PhaBOX** (optional) currently has no container image and runs via conda only.
> Combine `-profile docker` with `--phabox2_env <env>` if you need it, or omit
> PhaBOX under the container profiles.

---

## Running on HPC

Add `-profile slurm` to submit processes as SLURM jobs (compose with a tool
profile, e.g. `conda` or `singularity`):

```bash
nextflow run andrewbudge/phinder ... -profile singularity,slurm
```

---

## Setup options

```
bash setup.sh [--db-dir DIR] [--skip-envs] [--skip-dbs] [--with-phabox2]

  --db-dir DIR      Database directory (default: ~/.phinder_dbs)
  --skip-envs       Skip conda environment creation
  --skip-dbs        Skip database downloads
  --with-phabox2    Also install the phabox2 conda environment
```

---

## Citation

If you use phinder in your work, please cite the underlying tools:

- **geNomad**: Camargo, A.P., Roux, S., Schulz, F. et al. Identification of mobile genetic elements with geNomad. Nat Biotechnol 42, 1303–1312 (2024). https://doi.org/10.1038/s41587-023-01953-y
- **CheckV**: Nayfach, S., Camargo, A.P., Schulz, F. et al. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nat Biotechnol 39, 578–585 (2021). https://doi.org/10.1038/s41587-020-00774-7
- **Pharokka**: George Bouras, Roshan Nepal, Ghais Houtak, Alkis James Psaltis, Peter-John Wormald, Sarah Vreugde, Pharokka: a fast scalable bacteriophage annotation tool, Bioinformatics, Volume 39, Issue 1, January 2023, btac776, https://doi.org/10.1093/bioinformatics/btac776
- **PhaBOX**: Shang, J., Peng, C., Guan, J., Cai, D., Wang, D., & Sun, Y. (2026). PhaBOX2: an enhanced web server for discovering and analyzing viral contigs in metagenomic data. Nucleic Acids Research, gkag382, https://doi.org/10.1093/nar/gkag382
