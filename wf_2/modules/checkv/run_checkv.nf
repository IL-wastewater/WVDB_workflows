process CHECKV {
    tag "${id}"
    conda "bioconda::prodigal-gv bioconda::hmmer bioconda::checkv"
    publishDir { output_dir }, mode: 'copy', pattern: 'checkv_out'

    input:
    tuple val(id), path(assembly_500), val(output_dir)

    output:
    tuple val(id), path("checkv_out"), val(output_dir), emit: ok

    script:
    """
    set -e
    set -o pipefail

    checkv end_to_end \
        ${assembly_500} \
        checkv_out \
        -d ${params.checkv_db} \
        -t ${task.cpus}
    """
}
