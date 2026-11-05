process CHECKV {
    tag "${id}"
    conda "${params.checkv_env}"
    publishDir { output_dir }, mode: 'move', pattern: 'checkv_out'

    input:
    tuple val(id), path(assembly_500), val(output_dir)

    output:
    tuple val(id), path("checkv_out"), val(output_dir), emit: ok, optional: true
    tuple val(id), val("CHECKV"), path("${id}.CHECKV.failed.txt"), val(output_dir), emit: failed, optional: true

    script:
    """
    mkdir -p ${output_dir}
    set +e
    set -o pipefail

    write_failed() {
        reason="\$1"
        exit_status="\$2"

        rm -rf checkv_out

        {
            echo "sample: ${id}"
            echo "process: CHECKV"
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
        } > ${id}.CHECKV.failed.txt

        cp ${id}.CHECKV.failed.txt ${output_dir}/${id}.CHECKV.failed.txt
    }

    ${params.checkv_env}/bin/checkv end_to_end \
        ${assembly_500} \
        checkv_out \
        -d ${params.checkv_db} \
        -t ${task.cpus} \
        > ${id}.stdout.log 2> ${id}.stderr.log

    status=\$?

    if [ \$status -ne 0 ]; then
        write_failed "tool exited nonzero" "\$status"
        exit 0
    fi

    if [ ! -d checkv_out ]; then
        write_failed "no results: checkv_out directory missing" "\$status"
        exit 0
    fi
    """
}
