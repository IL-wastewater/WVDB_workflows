process BBDUK {
    tag "${id}"
    publishDir "${output_dir}", mode: 'link', pattern: "{${id}_rrna.fastq.gz,${id}_stats.txt}"

    input:
    tuple val(id), path(clean_ilv), val(output_dir)

    output:
    tuple val(id), path("${id}_no_rrna.fastq"), path("${id}_rrna.fastq.gz"), path("${id}_stats.txt"), val(output_dir), emit: ok, optional: true
    tuple val(id), val("BBDUK"), path("${id}.BBDUK.failed.txt"), val(output_dir), emit: failed, optional: true

    script:
    """
    mkdir -p ${output_dir}
    set +e
    set -o pipefail

    write_failed() {
        reason="\$1"
        exit_status="\$2"

        rm -f ${id}_no_rrna.fastq ${id}_rrna.fastq.gz ${id}_stats.txt

        {
            echo "sample: ${id}"
            echo "process: BBDUK"
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
        } > ${id}.BBDUK.failed.txt

        cp ${id}.BBDUK.failed.txt ${output_dir}/${id}.BBDUK.failed.txt
    }

    mkfifo rrna_pipe

    {
        pigz -p 64 < rrna_pipe > ${id}_rrna.fastq.gz &
        PIGZ_PID=\$!

        bbduk.sh t=${task.cpus} ref=${params.bbduk_db} \\
                 in=${clean_ilv} \\
                 out=${id}_no_rrna.fastq \\
                 outm=rrna_pipe \\
                 stats=${id}_stats.txt
        BBDUK_STATUS=\$?

        wait \$PIGZ_PID
        PIGZ_STATUS=\$?

        rm -f rrna_pipe

        if [ \$BBDUK_STATUS -ne 0 ]; then
            exit \$BBDUK_STATUS
        fi
        if [ \$PIGZ_STATUS -ne 0 ]; then
            exit \$PIGZ_STATUS
        fi
    } > ${id}.stdout.log 2> ${id}.stderr.log

    status=\$?
    rm -f rrna_pipe

    if [ \$status -ne 0 ]; then
        write_failed "tool exited nonzero" "\$status"
        exit 0
    fi

    if [ ! -s ${id}_no_rrna.fastq ]; then
        write_failed "no results: non-rRNA FASTQ missing or empty" "\$status"
        exit 0
    fi

    read_count=\$(awk 'END {print int(NR/4)}' ${id}_no_rrna.fastq)

    if [ "\${read_count:-0}" -eq 0 ]; then
        write_failed "no results: non-rRNA FASTQ contains zero reads" "\$status"
        exit 0
    fi
    """
}
