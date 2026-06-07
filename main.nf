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
    conda params.genomad_env ?: "${projectDir}/envs/genomad.yml"
    container 'quay.io/biocontainers/genomad:1.11.2--pyhdfd78af_0'
    publishDir "${params.outdir}/genomad", mode: 'copy'

    input:
    path contigs

    output:
    path 'output/**/*_virus_summary.tsv', emit: summary
    path 'output/**/*_virus.fna',         emit: fasta
    path 'output'
    path 'versions.yml',                  emit: versions

    script:
    """
    genomad end-to-end --cleanup \\
        -t ${task.cpus} --splits ${params.genomad_splits} \\
        ${contigs} output ${params.genomad_db}

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    genomad: \$(genomad --version | sed 's/geNomad, version //')
END_VERSIONS
    """

    stub:
    """
    mkdir -p output/summary
    touch output/summary/test_virus_summary.tsv output/summary/test_virus.fna

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    genomad: 1.11.2
END_VERSIONS
    """
}

process FILTER_GENOMAD {
    conda "${projectDir}/envs/rfilter.yml"
    container 'rocker/tidyverse:4.3.3'
    publishDir "${params.outdir}/genomad", mode: 'copy'

    input:
    path summary_tsv

    output:
    path 'filtered_genomad.tsv', emit: tsv
    path 'versions.yml',         emit: versions

    script:
    """
    filter_genomad.R \\
        ${summary_tsv} \\
        filtered_genomad.tsv \\
        ${params.min_provirus_score}

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    r-base: \$(Rscript -e 'cat(paste0(R.version\$major,".",R.version\$minor))')
    r-tidyverse: \$(Rscript -e 'cat(as.character(packageVersion("tidyverse")))')
END_VERSIONS
    """

    stub:
    """
    touch filtered_genomad.tsv

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    r-base: 4.3.3
    r-tidyverse: 2.0.0
END_VERSIONS
    """
}

process SUBSET_GENOMAD_FASTA {
    conda "${projectDir}/envs/seqkit.yml"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_1'

    input:
    path filtered_tsv
    path virus_fna

    output:
    path 'hq_viral_hits.fna', emit: fasta
    path 'versions.yml',      emit: versions

    script:
    """
    tail -n +2 ${filtered_tsv} | cut -f1 > ids.txt
    seqkit grep -n -f ids.txt ${virus_fna} > hq_viral_hits.fna

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    seqkit: \$(seqkit version | sed 's/seqkit v//')
END_VERSIONS
    """

    stub:
    """
    touch hq_viral_hits.fna

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    seqkit: 2.8.2
END_VERSIONS
    """
}

process CHECKV {
    conda params.checkv_env ?: "${projectDir}/envs/checkv.yml"
    container 'quay.io/biocontainers/checkv:1.0.3--pyhdfd78af_0'
    publishDir "${params.outdir}/checkv", mode: 'copy'

    input:
    path hq_fna

    output:
    path 'output/quality_summary.tsv', emit: summary
    path 'output/viruses.fna',         emit: viruses
    path 'output/proviruses.fna',      emit: proviruses
    path 'output'
    path 'versions.yml',               emit: versions

    script:
    """
    checkv end_to_end ${hq_fna} output \\
        -d ${params.checkv_db} \\
        --remove_tmp \\
        -t ${task.cpus}

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    checkv: \$(python -c 'import checkv; print(checkv.__version__)')
END_VERSIONS
    """

    stub:
    """
    mkdir -p output
    touch output/quality_summary.tsv output/viruses.fna output/proviruses.fna

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    checkv: 1.0.3
END_VERSIONS
    """
}

process FILTER_CHECKV {
    conda "${projectDir}/envs/rfilter.yml"
    container 'rocker/tidyverse:4.3.3'
    publishDir "${params.outdir}/checkv", mode: 'copy'

    input:
    path filtered_genomad_tsv
    path checkv_summary
    path input_contigs

    output:
    path 'potential_phage.tsv', emit: tsv
    path 'versions.yml',        emit: versions

    script:
    """
    filter_checkv.R \\
        ${filtered_genomad_tsv} \\
        ${checkv_summary} \\
        ${input_contigs} \\
        potential_phage.tsv \\
        ${params.min_coverage} \\
        '${params.checkv_quality_keep}'

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    r-base: \$(Rscript -e 'cat(paste0(R.version\$major,".",R.version\$minor))')
    r-tidyverse: \$(Rscript -e 'cat(as.character(packageVersion("tidyverse")))')
END_VERSIONS
    """

    stub:
    """
    touch potential_phage.tsv

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    r-base: 4.3.3
    r-tidyverse: 2.0.0
END_VERSIONS
    """
}

process CLEAN_PROVIRUS_HEADERS {
    // CheckV renames trimmed proviruses with a `_1 start-end/total` suffix,
    // which breaks seqkit ID matching. Strip the suffix back to a stable form.
    conda "${projectDir}/envs/seqkit.yml"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_1'

    input:
    path proviruses

    output:
    path 'proviruses_clean.fna', emit: fasta
    path 'versions.yml',         emit: versions

    script:
    """
    sed 's/|provirus_\\([0-9]*\\)_\\([0-9]*\\)_[0-9]* .*/|provirus_\\1_\\2/' \\
        ${proviruses} > proviruses_clean.fna

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    seqkit: \$(seqkit version | sed 's/seqkit v//')
END_VERSIONS
    """

    stub:
    """
    touch proviruses_clean.fna

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    seqkit: 2.8.2
END_VERSIONS
    """
}

process BUILD_CANDIDATES {
    conda "${projectDir}/envs/seqkit.yml"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_1'
    publishDir "${params.outdir}/candidates", mode: 'copy'

    input:
    path potential_tsv
    path viruses
    path proviruses_clean
    path genomad_fasta

    output:
    path 'candidate_phages.fna', emit: fasta
    path 'versions.yml',         emit: versions

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

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    seqkit: \$(seqkit version | sed 's/seqkit v//')
END_VERSIONS
    """

    stub:
    """
    touch candidate_phages.fna

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    seqkit: 2.8.2
END_VERSIONS
    """
}

process PHAROKKA {
    conda params.pharokka_env ?: "${projectDir}/envs/pharokka.yml"
    container 'quay.io/biocontainers/pharokka:1.8.2--pyhdfd78af_0'
    publishDir "${params.outdir}/pharokka", mode: 'copy'

    input:
    path candidates

    output:
    path 'output'
    path 'versions.yml', emit: versions

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

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    pharokka: \$(pharokka.py --version 2>&1 | sed 's/.*pharokka //I; s/^v//')
END_VERSIONS
    """

    stub:
    """
    mkdir -p output

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    pharokka: 1.8.2
END_VERSIONS
    """
}

// NOTE: PhaBOX (optional, gated on --phabox2_env) has no `container` directive
// yet — there is no biocontainer for the pinned 2.2 release. Under -profile
// docker/singularity these steps would need an image; build one once a 2.2+
// biocontainer is published (see envs/phabox2.yml). The conda path is unaffected.
process PHABOX_END_TO_END {
    conda params.phabox2_env ?: "${projectDir}/envs/phabox2.yml"
    publishDir "${params.outdir}/phabox/end_to_end", mode: 'copy'

    input:
    path candidates

    output:
    path 'output'
    path 'versions.yml', emit: versions

    script:
    def skip_arg = params.phabox_skip_phamer ? '--skip Y' : ''
    """
    phabox2 --task end_to_end \\
        --dbdir ${params.phabox_db} \\
        --outpth output \\
        --contigs ${candidates} \\
        --threads ${task.cpus} \\
        ${skip_arg}

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    phabox2: \$(python -c "import importlib.metadata as m; print(m.version('phabox2'))")
END_VERSIONS
    """

    stub:
    """
    mkdir -p output

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    phabox2: 2.2
END_VERSIONS
    """
}

process PHABOX_VOTU {
    conda params.phabox2_env ?: "${projectDir}/envs/phabox2.yml"
    publishDir "${params.outdir}/phabox/votu", mode: 'copy'

    input:
    path candidates

    output:
    path 'output'
    path 'versions.yml', emit: versions

    script:
    """
    phabox2 --task votu \\
        --dbdir ${params.phabox_db} \\
        --outpth output \\
        --contigs ${candidates} \\
        --threads ${task.cpus} \\
        --mode ${params.phabox_votu_mode}

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    phabox2: \$(python -c "import importlib.metadata as m; print(m.version('phabox2'))")
END_VERSIONS
    """

    stub:
    """
    mkdir -p output

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    phabox2: 2.2
END_VERSIONS
    """
}

process PHABOX_TREE {
    conda params.phabox2_env ?: "${projectDir}/envs/phabox2.yml"
    publishDir "${params.outdir}/phabox/tree", mode: 'copy'

    input:
    path candidates

    output:
    path 'output'
    path 'versions.yml', emit: versions

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

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    phabox2: \$(python -c "import importlib.metadata as m; print(m.version('phabox2'))")
END_VERSIONS
    """

    stub:
    """
    mkdir -p output

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    phabox2: 2.2
END_VERSIONS
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

    // Tool versions accumulate here, one small versions.yml per process.
    ch_versions = Channel.empty()

    // --- geNomad classification + R filter ---------------------------------
    GENOMAD(contigs_ch)
    FILTER_GENOMAD(GENOMAD.out.summary)
    SUBSET_GENOMAD_FASTA(FILTER_GENOMAD.out.tsv, GENOMAD.out.fasta)
    ch_versions = ch_versions.mix(GENOMAD.out.versions, FILTER_GENOMAD.out.versions, SUBSET_GENOMAD_FASTA.out.versions)

    // --- CheckV completeness + R filter ------------------------------------
    CHECKV(SUBSET_GENOMAD_FASTA.out.fasta)
    FILTER_CHECKV(FILTER_GENOMAD.out.tsv, CHECKV.out.summary, contigs_ch)
    ch_versions = ch_versions.mix(CHECKV.out.versions, FILTER_CHECKV.out.versions)

    // --- Build the candidate FASTA -----------------------------------------
    CLEAN_PROVIRUS_HEADERS(CHECKV.out.proviruses)
    BUILD_CANDIDATES(
        FILTER_CHECKV.out.tsv,
        CHECKV.out.viruses,
        CLEAN_PROVIRUS_HEADERS.out.fasta,
        GENOMAD.out.fasta
    )
    ch_versions = ch_versions.mix(CLEAN_PROVIRUS_HEADERS.out.versions, BUILD_CANDIDATES.out.versions)

    // --- Annotation --------------------------------------------------------
    PHAROKKA(BUILD_CANDIDATES.out.fasta)
    ch_versions = ch_versions.mix(PHAROKKA.out.versions)

    // --- Classification (optional — requires --phabox2_env) ----------------
    if (params.phabox2_env) {
        PHABOX_END_TO_END(BUILD_CANDIDATES.out.fasta)
        PHABOX_VOTU(BUILD_CANDIDATES.out.fasta)
        PHABOX_TREE(BUILD_CANDIDATES.out.fasta)
        ch_versions = ch_versions.mix(PHABOX_END_TO_END.out.versions, PHABOX_VOTU.out.versions, PHABOX_TREE.out.versions)
    }

    // --- Provenance: one combined versions.yml for the run -----------------
    ch_versions
        .collectFile(name: 'versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)
}
