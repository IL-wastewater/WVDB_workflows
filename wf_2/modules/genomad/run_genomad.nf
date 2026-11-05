process GENOMAD {
    tag "${id}"
    conda "${params.genomad_env}"
    publishDir { output_dir }, mode: 'move', pattern: { "${id}*" }

    beforeScript """
        export PATH=${params.genomad_env}/bin:\$PATH
    """

    input:
    tuple val(id), path(assembly_500), val(output_dir)

    output:
    tuple val(id), path("${id}*"), val(output_dir), emit: ok, optional: true
    tuple val(id), val("GENOMAD"), path("${id}.GENOMAD.failed.txt"), val(output_dir), emit: failed, optional: true

    script:
    """
    mkdir -p ${output_dir}
    set +e
    set -o pipefail

    write_failed() {
        reason="\$1"
        exit_status="\$2"

        rm -rf genomad_output ${id}_*

        {
            echo "sample: ${id}"
            echo "process: GENOMAD"
            echo "reason: \$reason"
            echo "exit_status: \$exit_status"
            echo "workdir: \$PWD"
            echo "timestamp: \$(date -Is)"
            echo
            echo "stderr_tail:"
            if [ -s ${id}.stderr.log ]; then
                tail -n 80 ${id}.stderr.log
            else
                echo "[no stderr captured]"
            fi
        } > ${id}.GENOMAD.failed.txt

        cp ${id}.GENOMAD.failed.txt ${output_dir}/${id}.GENOMAD.failed.txt
    }

    ${params.genomad_env}/bin/genomad end-to-end \
        --cleanup \
        ${assembly_500} \
        genomad_output \
        ${params.genomad_db} \
        -t ${task.cpus} \
        > ${id}.stdout.log 2> ${id}.stderr.log

    status=\$?

    if [ \$status -ne 0 ]; then
        write_failed "tool exited nonzero" "\$status"
        exit 0
    fi

    if [ -d genomad_output ]; then
        mv genomad_output/* .
    fi

    if ! ls ${id}* >/dev/null 2>&1; then
        write_failed "no results: geNomad outputs missing" "\$status"
        exit 0
    fi
    """
}
