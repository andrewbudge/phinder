# phinder

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
- [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/)

---

## Quick start

**1. Download and run setup**

```bash
curl -O https://raw.githubusercontent.com/andrewcbudge/phinder/main/setup.sh
bash setup.sh
```

This creates pinned conda environments (`genomad_phinder`, `checkv_phinder`,
`pharokka_phinder`) and downloads all required databases to `~/.phinder_dbs`.
The exact run command is printed at the end.

To also set up PhaBOX (optional):
```bash
bash setup.sh --with-phabox2
```

**2. Run**

```bash
nextflow run andrewcbudge/phinder \
    --input contigs.fasta \
    --genomad_db  ~/.phinder_dbs/genomad_db \
    --checkv_db   ~/.phinder_dbs/checkv_db \
    --pharokka_db ~/.phinder_dbs/pharokka_db \
    --genomad_env  $(conda info --base)/envs/genomad_phinder \
    --checkv_env   $(conda info --base)/envs/checkv_phinder \
    --pharokka_env $(conda info --base)/envs/pharokka_phinder \
    -profile conda
```

Use `-resume` on reruns to skip completed steps.

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
| `--genomad_env` | builds from bioconda | Path to existing geNomad conda env |
| `--checkv_env` | builds from bioconda | Path to existing CheckV conda env |
| `--pharokka_env` | builds from bioconda | Path to existing Pharokka conda env |
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

## Running on HPC

Add `-profile slurm` to submit processes as SLURM jobs:

```bash
nextflow run andrewcbudge/phinder ... -profile conda,slurm
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

- **geNomad**: Camargo et al. (2023) *Nature Biotechnology*
- **CheckV**: Nayfach et al. (2021) *Nature Biotechnology*
- **Pharokka**: Bouras et al. (2023) *Bioinformatics*
- **PhaBOX**: Shang et al. (2023) *Briefings in Bioinformatics*
