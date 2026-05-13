#!/usr/bin/env nextflow
// ============================================================================
// phinder — phage finder for metagenomic assemblies
//
// Input:  a single combined contig FASTA (one or many samples, pre-merged)
// Output: viral candidates classified, annotated (Pharokka), and clustered /
//         placed phylogenetically (PhaBOX).
// ============================================================================

nextflow.enable.dsl = 2

// ============================================================================
// Processes
// ============================================================================

process GENOMAD {
    conda params.genomad_env ?: 'bioconda::genomad'
    publishDir "${params.outdir}/genomad", mode: 'copy'

    input:
    path contigs

    output:
    path 'output/**/*_virus_summary.tsv', emit: summary
    path 'output/**/*_virus.fna',         emit: fasta
    path 'output'

    script:
    """
    genomad end-to-end --cleanup \\
        -t ${task.cpus} --splits ${params.genomad_splits} \\
        ${contigs} output ${params.genomad_db}
    """
}

process FILTER_GENOMAD {
    conda 'conda-forge::r-base conda-forge::r-tidyverse'
    publishDir "${params.outdir}/genomad", mode: 'copy'

    input:
    path summary_tsv

    output:
    path 'filtered_genomad.tsv'

    script:
    """
    filter_genomad.R \\
        ${summary_tsv} \\
        filtered_genomad.tsv \\
        ${params.min_provirus_score}
    """
}

process SUBSET_GENOMAD_FASTA {
    conda 'bioconda::seqkit'

    input:
    path filtered_tsv
    path virus_fna

    output:
    path 'hq_viral_hits.fna'

    script:
    """
    tail -n +2 ${filtered_tsv} | cut -f1 > ids.txt
    seqkit grep -n -f ids.txt ${virus_fna} > hq_viral_hits.fna
    """
}

process CHECKV {
    conda params.checkv_env ?: 'bioconda::checkv'
    publishDir "${params.outdir}/checkv", mode: 'copy'

    input:
    path hq_fna

    output:
    path 'output/quality_summary.tsv', emit: summary
    path 'output/viruses.fna',         emit: viruses
    path 'output/proviruses.fna',      emit: proviruses
    path 'output'

    script:
    """
    checkv end_to_end ${hq_fna} output \\
        -d ${params.checkv_db} \\
        --remove_tmp \\
        -t ${task.cpus}
    """
}

process FILTER_CHECKV {
    conda 'conda-forge::r-base conda-forge::r-tidyverse'
    publishDir "${params.outdir}/checkv", mode: 'copy'

    input:
    path filtered_genomad_tsv
    path checkv_summary
    path input_contigs

    output:
    path 'potential_phage.tsv'

    script:
    """
    filter_checkv.R \\
        ${filtered_genomad_tsv} \\
        ${checkv_summary} \\
        ${input_contigs} \\
        potential_phage.tsv \\
        ${params.min_coverage} \\
        '${params.checkv_quality_keep}'
    """
}

process CLEAN_PROVIRUS_HEADERS {
    // CheckV renames trimmed proviruses with a `_1 start-end/total` suffix,
    // which breaks seqkit ID matching. Strip the suffix back to a stable form.
    input:
    path proviruses

    output:
    path 'proviruses_clean.fna'

    script:
    """
    sed 's/|provirus_\\([0-9]*\\)_\\([0-9]*\\)_[0-9]* .*/|provirus_\\1_\\2/' \\
        ${proviruses} > proviruses_clean.fna
    """
}

process BUILD_CANDIDATES {
    conda 'bioconda::seqkit'
    publishDir "${params.outdir}/candidates", mode: 'copy'

    input:
    path potential_tsv
    path viruses
    path proviruses_clean
    path genomad_fasta

    output:
    path 'candidate_phages.fna'

    script:
    """
    cat ${viruses} ${proviruses_clean} > all_checkv.fna

    tail -n +2 ${potential_tsv} | cut -f1 > wanted_ids.txt
    seqkit grep -n -f wanted_ids.txt all_checkv.fna > candidate_phages.fna

    # Fall back to the geNomad FASTA for any IDs CheckV dropped.
    grep ">" candidate_phages.fna | sed 's/>//g' | awk '{print \$1}' | sort > have_ids.txt
    sort wanted_ids.txt > wanted_sorted.txt
    comm -23 wanted_sorted.txt have_ids.txt > missing_ids.txt

    if [ -s missing_ids.txt ]; then
        echo "[WARN] missing IDs from CheckV, falling back to geNomad:"
        cat missing_ids.txt
        seqkit grep -f missing_ids.txt ${genomad_fasta} >> candidate_phages.fna
    fi

    # Normalize headers for downstream tools (| breaks some parsers).
    sed -i 's/|/_/g' candidate_phages.fna
    """
}

process PHAROKKA {
    conda params.pharokka_env ?: "${projectDir}/envs/pharokka.yml"
    publishDir "${params.outdir}/pharokka", mode: 'copy'

    input:
    path candidates

    output:
    path 'output'

    script:
    """
    pharokka.py \\
        -i ${candidates} \\
        -o output \\
        -d ${params.pharokka_db} \\
        -t ${task.cpus} \\
        -m -s \\
        --dnaapler \\
        --meta_hmm \\
        -g ${params.pharokka_gene_predictor} \\
        -f
    """
}

process PHABOX_END_TO_END {
    conda params.phabox2_env
    publishDir "${params.outdir}/phabox/end_to_end", mode: 'copy'

    input:
    path candidates

    output:
    path 'output'

    script:
    def skip_arg = params.phabox_skip_phamer ? '--skip Y' : ''
    """
    phabox2 --task end_to_end \\
        --dbdir ${params.phabox_db} \\
        --outpth output \\
        --contigs ${candidates} \\
        --threads ${task.cpus} \\
        ${skip_arg}
    """
}

process PHABOX_VOTU {
    conda params.phabox2_env
    publishDir "${params.outdir}/phabox/votu", mode: 'copy'

    input:
    path candidates

    output:
    path 'output'

    script:
    """
    phabox2 --task votu \\
        --dbdir ${params.phabox_db} \\
        --outpth output \\
        --contigs ${candidates} \\
        --threads ${task.cpus} \\
        --mode ${params.phabox_votu_mode}
    """
}

process PHABOX_TREE {
    conda params.phabox2_env
    publishDir "${params.outdir}/phabox/tree", mode: 'copy'

    input:
    path candidates

    output:
    path 'output'

    script:
    def markers = params.phabox_tree_markers
        .tokenize(',')
        .collect { "--marker ${it.trim()}" }
        .join(' ')
    """
    phabox2 --task tree \\
        --dbdir ${params.phabox_db} \\
        --outpth output \\
        --contigs ${candidates} \\
        --threads ${task.cpus} \\
        ${markers} \\
        --tree Y --msa Y
    """
}

// ============================================================================
// Workflow
// ============================================================================

workflow {
    // --- Sanity checks ------------------------------------------------------
    if (!params.input)       error "Missing --input (combined contig FASTA)"
    if (!params.genomad_db)  error "Missing --genomad_db"
    if (!params.checkv_db)   error "Missing --checkv_db"
    if (!params.pharokka_db) error "Missing --pharokka_db"
    if (params.phabox2_env && !params.phabox_db) error "PhaBOX enabled (--phabox2_env set) but --phabox_db not provided"

    contigs_ch = Channel.fromPath(params.input, checkIfExists: true)

    // --- geNomad classification + R filter ---------------------------------
    GENOMAD(contigs_ch)
    FILTER_GENOMAD(GENOMAD.out.summary)
    SUBSET_GENOMAD_FASTA(FILTER_GENOMAD.out, GENOMAD.out.fasta)

    // --- CheckV completeness + R filter ------------------------------------
    CHECKV(SUBSET_GENOMAD_FASTA.out)
    FILTER_CHECKV(FILTER_GENOMAD.out, CHECKV.out.summary, contigs_ch)

    // --- Build the candidate FASTA -----------------------------------------
    CLEAN_PROVIRUS_HEADERS(CHECKV.out.proviruses)
    BUILD_CANDIDATES(
        FILTER_CHECKV.out,
        CHECKV.out.viruses,
        CLEAN_PROVIRUS_HEADERS.out,
        GENOMAD.out.fasta
    )

    // --- Annotation --------------------------------------------------------
    PHAROKKA(BUILD_CANDIDATES.out)

    // --- Classification (optional — requires --phabox2_env) ----------------
    if (params.phabox2_env) {
        PHABOX_END_TO_END(BUILD_CANDIDATES.out)
        PHABOX_VOTU(BUILD_CANDIDATES.out)
        PHABOX_TREE(BUILD_CANDIDATES.out)
    }
}
