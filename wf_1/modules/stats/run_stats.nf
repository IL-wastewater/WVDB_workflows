process STATS {
    tag "${id}"
    publishDir "${output_dir}", mode: 'move'

    input:
    tuple val(id), path(no_rrna_final), path(no_rrna_human), path(rrna), val(output_dir)

    output:
    tuple val(id), path("${id}_combined_stats.txt"), val(output_dir), emit: ok, optional: true
    tuple val(id), val("STATS"), path("${id}.STATS.failed.txt"), val(output_dir), emit: failed, optional: true

    script:
    """
    mkdir -p ${output_dir}
    set +e

    write_failed() {
        reason="\$1"
        exit_status="\$2"

        rm -f ${id}_combined_stats.txt

        {
            echo "sample: ${id}"
            echo "process: STATS"
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
        } > ${id}.STATS.failed.txt

        cp ${id}.STATS.failed.txt ${output_dir}/${id}.STATS.failed.txt
    }

    seqkit stats -j ${task.cpus} ${rrna} ${no_rrna_human} ${no_rrna_final} \\
        > ${id}_combined_stats.txt 2> ${id}.stderr.log

    status=\$?

    if [ \$status -ne 0 ]; then
        write_failed "tool exited nonzero" "\$status"
        exit 0
    fi

    if [ ! -s ${id}_combined_stats.txt ]; then
        write_failed "no results: combined stats file missing or empty" "\$status"
        exit 0
    fi
    """
}
