process POST_PROCESS {
    tag "${id}"
    publishDir { output_dir }, mode: 'link', pattern: "*_contigs_min*.fa"
    conda "bioconda::seqkit"

    input:
    tuple val(id), path(assembly), path(megahit_log), val(output_dir)

    output:
    tuple val(id), path("${id}_contigs_min500.fa"), path("${id}_contigs_min2000.fa"), val(output_dir), emit: ok, optional: true

    script:
    """
    set -e
    set -o pipefail

    seqkit seq -m 500 ${assembly} | \\
        seqkit replace -p .+ -r "${id}_{nr}" > ${id}_contigs_min500.fa

    seqkit seq -m 2000 ${id}_contigs_min500.fa > ${id}_contigs_min2000.fa
    """
}
