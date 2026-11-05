process MEGAHIT {
    tag "${id}"
    publishDir { output_dir }, mode: 'link', pattern: { "{${id}.contigs.fa,${id}.log}" }

    input:
    tuple val(id), path(reads), val(output_dir)

    output:
    tuple val(id), path("${id}.contigs.fa"), path("${id}.log"), val(output_dir), emit: ok, optional: true
    tuple val(id), val("MEGAHIT"), path("${id}.MEGAHIT.failed.txt"), val(output_dir), emit: failed, optional: true

    script:
    """
    mkdir -p ${output_dir}
    set +e
    set -o pipefail

    rm -rf megahit_out
    rm -f ${id}.contigs.fa ${id}.log ${id}.stdout.log ${id}.stderr.log ${id}.MEGAHIT.failed.txt

    write_failed() {
        reason="\$1"
        exit_status="\$2"

        rm -f ${id}.contigs.fa ${id}.log
        rm -rf megahit_out

        {
            echo "sample: ${id}"
            echo "process: MEGAHIT"
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
        } > ${id}.MEGAHIT.failed.txt

        cp ${id}.MEGAHIT.failed.txt ${output_dir}/${id}.MEGAHIT.failed.txt
    }

    ${params.megahit_env}/megahit \
        --12 ${reads} \
        --num-cpu-threads ${task.cpus} \
        --out-dir megahit_out \
        --out-prefix ${id} \
        > ${id}.stdout.log 2> ${id}.stderr.log

    status=\$?

    if [ \$status -ne 0 ]; then
        write_failed "tool exited nonzero" "\$status"
        exit 0
    fi

    if [ -s megahit_out/${id}.contigs.fa ]; then
        mv megahit_out/${id}.contigs.fa .
    fi

    if [ -s megahit_out/${id}.log ]; then
        mv megahit_out/${id}.log .
    elif [ -s megahit_out/log ]; then
        cp megahit_out/log ${id}.log
    elif [ -s ${id}.stderr.log ]; then
        cp ${id}.stderr.log ${id}.log
    else
        touch ${id}.log
    fi

    if [ ! -s ${id}.contigs.fa ]; then
        write_failed "no results: contigs file missing or empty" "\$status"
        exit 0
    fi
    """
}
