process FASTP {
    tag "${id}"
    publishDir "${output_dir}", mode: 'link', pattern: "*.{json,html}"

    input:
    tuple val(id), val(read_type), path(reads), val(output_dir)

    output:
    tuple val(id), path("${id}_clean_ilv.fastq.gz"), path("${id}_clean.json"), path("${id}_clean.html"), val(output_dir), emit: ok, optional: true
    tuple val(id), val("FASTP"), path("${id}.FASTP.failed.txt"), val(output_dir), emit: failed, optional: true

    script:
    def fastp_common = "fastp -w ${task.cpus} -l 50 --dont_eval_duplication --trim_poly_g --poly_g_min_len=5 --poly_x_min_len=10 --low_complexity_filter 30 -y"

    def cmd = read_type == 'paired'
        ? """
          ${fastp_common} \\
              -i ${reads[0]} -I ${reads[1]} \\
              -j ${id}_clean.json -h ${id}_clean.html \\
              --stdout | pigz -p ${task.cpus} > ${id}_clean_ilv.fastq.gz
          """
        : """
          ${fastp_common} \\
              --interleaved_in \\
              -i ${reads} \\
              -j ${id}_clean.json -h ${id}_clean.html \\
              --stdout | pigz -p ${task.cpus} > ${id}_clean_ilv.fastq.gz
          """

    """
    mkdir -p ${output_dir}
    set +e
    set -o pipefail

    write_failed() {
        reason="\$1"
        exit_status="\$2"

        rm -f ${id}_clean_ilv.fastq.gz ${id}_clean.json ${id}_clean.html

        {
            echo "sample: ${id}"
            echo "process: FASTP"
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
        } > ${id}.FASTP.failed.txt

        cp ${id}.FASTP.failed.txt ${output_dir}/${id}.FASTP.failed.txt
    }

    {
        ${cmd}
    } > ${id}.stdout.log 2> ${id}.stderr.log

    status=\$?

    if [ \$status -ne 0 ]; then
        write_failed "tool exited nonzero" "\$status"
        exit 0
    fi

    if [ ! -s ${id}_clean_ilv.fastq.gz ]; then
        write_failed "no results: cleaned FASTQ missing or empty" "\$status"
        exit 0
    fi

    read_count=\$(zcat ${id}_clean_ilv.fastq.gz 2>/dev/null | awk 'END {print int(NR/4)}')
    zcat_status=\${PIPESTATUS[0]}

    if [ \$zcat_status -ne 0 ]; then
        write_failed "no results: cleaned FASTQ is not readable gzip" "\$zcat_status"
        exit 0
    fi

    if [ "\${read_count:-0}" -eq 0 ]; then
        write_failed "no results: cleaned FASTQ contains zero reads" "\$status"
        exit 0
    fi
    """
}
