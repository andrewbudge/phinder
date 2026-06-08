# Bundled test data

`test_contigs.fasta` is a tiny (~64 KB) input for verifying a phinder install.
The four contigs are **real metaMDBG assembly contigs from the biocrust dataset
phinder was developed on** (SRA run accessions, underlying reads public), chosen
to exercise distinct branches of the filter logic:

| Contig | Length | Coverage | geNomad call | Expected fate |
|--------|--------|----------|--------------|---------------|
| `SRR35808553_ctg1561` | 15,576 | 138.3 | DTR (virus_score 0.86) | **Kept** — DTR bypasses the geNomad score gate *and* the CheckV quality gate |
| `SRR35808779_ctg381` | 20,746 | 90.5 | DTR | **Kept** |
| `SRR35808900_ctg9697` | 17,385 | 2.89 | DTR | **Kept, flagged** `low_coverage=TRUE` (coverage < `--min_coverage` 3) |
| `SRR35808553_ctg1953` | 10,125 | 64.8 | plasmid | **Dropped** — non-viral, never enters the viral path |

So a correct real run keeps the three DTR phages (one flagged low-coverage) and
drops the plasmid.

## Verifying an install

Fast wiring check (no databases needed — every process is stubbed):

```bash
nextflow run . -profile test -stub-run
```

This validates that the whole workflow is wired correctly end to end. It does
**not** run the real tools, so it does not assert the fates above.

## What this set does not yet cover

Every *non-DTR* positive in the source data is a **provirus embedded in a
multi-Mb host contig**, which would blow the size budget — so the provirus
excision path and the "passes geNomad, fails CheckV" branch are not represented
here. Those belong in a heavier, on-demand fixture (real tools + host contigs +
databases), not this lightweight bundle.
