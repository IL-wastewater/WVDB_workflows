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

    echo "[INFO] Starting DEACON depletion for ${id}" >&2

    set +e

    deacon filter -d \\
        -t ${params.deacon_threads} \\
        ${params.deacon_db} \\
        ${no_rrna_fastq} \\
        -o ${id}_no_rrna_final_reads.fastq.gz \\
        > ${id}.deacon_final.stdout.log 2> ${id}.deacon_final.stderr.log

    FINAL_STATUS=\$?

    deacon filter \\
        -t ${params.deacon_threads} \\
        ${params.deacon_db} \\
        ${no_rrna_fastq} \\
        -o ${id}_no_rrna_human_reads.fastq.gz \\
        > ${id}.deacon_human.stdout.log 2> ${id}.deacon_human.stderr.log

    HUMAN_STATUS=\$?

    cat ${id}.deacon_final.stderr.log ${id}.deacon_human.stderr.log > ${id}.stderr.log

    set -e

    if [ \$FINAL_STATUS -ne 0 ]; then
        write_failed "deacon final depletion command exited nonzero" "\$FINAL_STATUS"
        exit 0
    fi

    if [ \$HUMAN_STATUS -ne 0 ]; then
        write_failed "deacon human extraction command exited nonzero" "\$HUMAN_STATUS"
        exit 0
    fi

    if [ ! -s ${id}_no_rrna_final_reads.fastq.gz ]; then
        write_failed "no results: final non-human FASTQ missing or empty" "0"
        exit 0
    fi

    if [ ! -s ${id}_no_rrna_human_reads.fastq.gz ]; then
        write_failed "human FASTQ missing or empty" "0"
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

    echo "[INFO] DEACON processing complete for ${id}" >&2
    """
}
