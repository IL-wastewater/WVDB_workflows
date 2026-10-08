process GENOMAD {
    tag "${id}"
    conda "bioconda::genomad"
    publishDir { output_dir }, mode: 'link', pattern: 'genomad_output'

    input:
    tuple val(id), path(assembly_500), val(output_dir)

    output:
    tuple val(id), path("genomad_output"), val(output_dir), emit: ok

    script:
    """
    set -e
    set -o pipefail

    genomad end-to-end \
        --cleanup \
        ${assembly_500} \
        genomad_output \
        ${params.genomad_db} \
        -t ${task.cpus}
    """
}
