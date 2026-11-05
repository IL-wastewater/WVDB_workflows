process POST_PROCESS {
    tag "${id}"
    publishDir { output_dir }, mode: 'link', pattern: { "{${id}_contigs_min500.fa,${id}_contigs_min2000.fa}" }

    input:
    tuple val(id), path(assembly), path(megahit_log), val(output_dir)

    output:
    tuple val(id), path("${id}_contigs_min500.fa"), path("${id}_contigs_min2000.fa"), val(output_dir), emit: ok, optional: true
    tuple val(id), val("POST_PROCESS"), path("${id}.POST_PROCESS.failed.txt"), val(output_dir), emit: failed, optional: true

    script:
    """
    mkdir -p ${output_dir}
    set +e
    set -o pipefail

    rm -f ${id}_contigs_min500.fa ${id}_contigs_min2000.fa ${id}.stdout.log ${id}.stderr.log ${id}.POST_PROCESS.failed.txt

    write_failed() {
        reason="\$1"
        exit_status="\$2"

        rm -f ${id}_contigs_min500.fa ${id}_contigs_min2000.fa

        {
            echo "sample: ${id}"
            echo "process: POST_PROCESS"
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
        } > ${id}.POST_PROCESS.failed.txt

        cp ${id}.POST_PROCESS.failed.txt ${output_dir}/${id}.POST_PROCESS.failed.txt
    }

    {
        ${params.seqkit_env}/bin/seqkit seq -m 500 ${assembly} | \\
            ${params.seqkit_env}/bin/seqkit replace -p .+ -r "${id}_{nr}" > ${id}_contigs_min500.fa

        ${params.seqkit_env}/bin/seqkit seq -m 2000 ${id}_contigs_min500.fa > ${id}_contigs_min2000.fa
    } > ${id}.stdout.log 2> ${id}.stderr.log

    status=\$?

    if [ \$status -ne 0 ]; then
        write_failed "tool exited nonzero" "\$status"
        exit 0
    fi

    if [ ! -s ${id}_contigs_min500.fa ]; then
        write_failed "no results: min500 contigs file missing or empty" "\$status"
        exit 0
    fi

    if [ ! -e ${id}_contigs_min2000.fa ]; then
        write_failed "no results: min2000 contigs file was not created" "\$status"
        exit 0
    fi
    """
}
