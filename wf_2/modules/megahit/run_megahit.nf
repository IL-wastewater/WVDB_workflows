process MEGAHIT {
    tag "${id}"
    publishDir [path: { output_dir }, mode: 'copy', pattern: "*.contigs.fa"]
    conda "bioconda::megahit"

    time 2.d

    input:
    tuple val(id), path(reads), val(output_dir)

    output:
    tuple val(id), path("${id}.contigs.fa"), path("${id}.log"), val(output_dir), emit: ok

    script:
    """
    set -e
    set -o pipefail

    megahit \
        --12 ${reads} \
        --num-cpu-threads ${task.cpus} \
        --out-dir megahit_out \
        --out-prefix ${id}

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
    """
}
