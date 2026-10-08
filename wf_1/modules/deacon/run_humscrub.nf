process DEACON {
    tag "${id}"
    publishDir "${output_dir}", mode: 'link', pattern: "*_no_rrna_{final,human}_reads.fastq.gz"
    scratch true

    input:
    tuple val(id), path(no_rrna_fastq), val(output_dir)

    output:
    tuple val(id), path("${id}_no_rrna_final_reads.fastq.gz"), path("${id}_no_rrna_human_reads.fastq.gz"), val(output_dir), emit: ok, optional: true
    tuple val(id), val("DEACON"), path("${id}.DEACON.failed.txt"), val(output_dir), emit: failed, optional: true

    script:
    """
    set -euo pipefail

    mkdir -p ${output_dir}
    mkdir -p tmp

    write_failed() {
        reason="\$1"
        exit_status="\$2"

        rm -f ${id}_no_rrna_final_reads.fastq.gz ${id}_no_rrna_human_reads.fastq.gz

        {
            echo "sample: ${id}"
            echo "process: DEACON"
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
        } > ${id}.DEACON.failed.txt

        cp ${id}.DEACON.failed.txt ${output_dir}/${id}.DEACON.failed.txt
    }

    LOCAL_DB=/dev/shm/sra_db_${id}_\${SLURM_JOB_ID:-\$\$}

    echo "[INFO] Copying database to RAM..." >&2
    cp -r ${params.sra_hum_scrub_db} \$LOCAL_DB

    trap "echo '[INFO] Cleaning up RAM database...'; rm -rf \$LOCAL_DB; exit" EXIT SIGTERM SIGINT

    echo "[INFO] Starting no_rrna processing..." >&2

    mkfifo tmp/no_rrna_final.pipe tmp/no_rrna_human.pipe

    pigz -p 32 < tmp/no_rrna_final.pipe > ${id}_no_rrna_final_reads.fastq.gz &
    PIGZ1=\$!

    pigz -p 32 < tmp/no_rrna_human.pipe > ${id}_no_rrna_human_reads.fastq.gz &
    PIGZ2=\$!

    set +e

    scrub.sh -d \$LOCAL_DB -x -p ${params.scrub_threads} -s \\
             -i ${no_rrna_fastq} \\
             -o tmp/no_rrna_final.pipe \\
             -u tmp/no_rrna_human.pipe \
             > ${id}.stdout.log 2> ${id}.stderr.log

    SCRUB_STATUS=\$?

    wait \$PIGZ1
    PIGZ1_STATUS=\$?

    wait \$PIGZ2
    PIGZ2_STATUS=\$?

    set -e

    rm -f tmp/no_rrna_final.pipe tmp/no_rrna_human.pipe
    rm -rf tmp

    if [ \$SCRUB_STATUS -ne 0 ]; then
        write_failed "scrub.sh exited nonzero" "\$SCRUB_STATUS"
        exit 0
    fi

    if [ \$PIGZ1_STATUS -ne 0 ]; then
        write_failed "pigz final output exited nonzero" "\$PIGZ1_STATUS"
        exit 0
    fi

    if [ \$PIGZ2_STATUS -ne 0 ]; then
        write_failed "pigz human output exited nonzero" "\$PIGZ2_STATUS"
        exit 0
    fi

    if [ ! -s ${id}_no_rrna_final_reads.fastq.gz ]; then
        write_failed "no results: final non-human FASTQ missing or empty" "0"
        exit 0
    fi

    set +e
    read_count=\$(zcat ${id}_no_rrna_final_reads.fastq.gz 2>/dev/null | awk 'END {print int(NR/4)}')
    zcat_status=\${PIPESTATUS[0]}
    set -e

    if [ \$zcat_status -ne 0 ]; then
        write_failed "no results: final non-human FASTQ is not readable gzip" "\$zcat_status"
        exit 0
    fi

    if [ "\${read_count:-0}" -eq 0 ]; then
        write_failed "no results: final non-human FASTQ contains zero reads" "0"
        exit 0
    fi

    echo "[INFO] Processing complete for ${id}" >&2
    """
}
