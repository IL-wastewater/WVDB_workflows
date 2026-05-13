reads_ch = channel.fromFilePairs(params.reads) 

process fastp {
    tag "${id}"
    publishDir "${params.output_dir}/${id}", mode : 'link', pattern: "*.{json,html}"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_clean_ilv.fastq.gz"), path("${id}_clean.json"), path("${id}_clean.html")

    script:
    """
    mkdir -p ${params.output_dir}/${id}

    fastp -w ${task.cpus} -l 50 --dont_eval_duplication --trim_poly_g --poly_g_min_len=5 \
          --poly_x_min_len=10 --low_complexity_filter 30 -y \
          -i ${reads[0]} -I ${reads[1]} \
          -j ${id}_clean.json -h ${id}_clean.html \
          --stdout | pigz -p ${task.cpus} > ${id}_clean_ilv.fastq.gz
    """
}

process bbduk {
    tag "${id}"
    publishDir "${params.output_dir}/${id}", mode : 'link', pattern: "{${id}_rrna.fastq.gz,${id}_stats.txt}"

    input:
    tuple val(id), path(clean_ilv)

    output:
    tuple val(id), path("${id}_no_rrna.fastq"), path("${id}_rrna.fastq.gz"), path("${id}_stats.txt")

    script:
    """
    # Create named pipe only for rrna (which we're keeping)
    mkfifo rrna_pipe
    
    # Start compression for rrna in background
    pigz -p params.scrub_threads < rrna_pipe > ${id}_rrna.fastq.gz &
    PIGZ_PID=\$!
    
    # Run bbduk with uncompressed no_rrna output
    bbduk.sh t=${task.cpus} ref=$params.bbduk_db \
             in=${clean_ilv} \
             out=${id}_no_rrna.fastq \
             outm=rrna_pipe \
             stats=${id}_stats.txt
    
    # Wait for rrna compression to finish
    wait \$PIGZ_PID
    
    # Cleanup pipe
    rm -f rrna_pipe
    """
}

process sra_hum_scrub {
    tag "${id}"
    publishDir "${params.output_dir}/${id}", mode : 'link', pattern: "*_no_rrna_{final,human}_reads.fastq.gz"
    scratch true

    input:
    tuple val(id), file(no_rrna_fastq)  // Now uncompressed

    output:
    tuple val(id), file("${id}_no_rrna_final_reads.fastq.gz"), file("${id}_no_rrna_human_reads.fastq.gz")

    script:
    """
    set -euo pipefail
    
    # Setup temp directory in scratch space
    mkdir -p tmp
    
    # Copy database to RAM with unique name per job
    LOCAL_DB=/dev/shm/sra_db_${id}_\${SLURM_JOB_ID:-\$\$}
    echo "[INFO] Copying database to RAM..." >&2
    cp -r $params.sra_hum_scrub_db \$LOCAL_DB
    
    # CRITICAL: Trap ensures cleanup even if job fails/cancelled
    trap "echo '[INFO] Cleaning up RAM database...'; rm -rf \$LOCAL_DB; exit" EXIT SIGTERM SIGINT
    
    echo "[INFO] Starting no_rrna processing..." >&2
    
    mkfifo tmp/no_rrna_final.pipe tmp/no_rrna_human.pipe
    
    pigz -p 32 < tmp/no_rrna_final.pipe > ${id}_no_rrna_final_reads.fastq.gz &
    PIGZ1=\$!
    pigz -p 32 < tmp/no_rrna_human.pipe > ${id}_no_rrna_human_reads.fastq.gz &
    PIGZ2=\$!
    
    # No decompression needed! Direct input
    scrub.sh -d \$LOCAL_DB -x -p ${params.scrub_threads} -s \
             -i ${no_rrna_fastq} \
             -o tmp/no_rrna_final.pipe \
             -u tmp/no_rrna_human.pipe
    
    wait \$PIGZ1 \$PIGZ2
    rm -f tmp/no_rrna_final.pipe tmp/no_rrna_human.pipe
    rm -rf tmp
    
    echo "[INFO] Processing complete for ${id}" >&2
    """
}

process stats {
    tag "${id}"
    publishDir "${params.output_dir}/${id}", mode : 'move'

    input:
    tuple val(id), file(no_rrna_final), file(no_rrna_human), file(rrna)

    output:
    tuple val(id), file("${id}_combined_stats.txt")

    script:
    """
    # Calculate stats for the 3 key output files
    seqkit stats -j ${task.cpus} ${rrna} ${no_rrna_human} ${no_rrna_final} \
        > ${id}_combined_stats.txt
    """
}

workflow {
    fastp_output = fastp(reads_ch)
    bbduk_output = bbduk(fastp_output.map { id, clean_ilv, json, html -> tuple(id, clean_ilv) })
    
    // Only process no_rrna reads through human scrubber
    sra_hum_scrub_output = sra_hum_scrub(bbduk_output.map { id, no_rrna, rrna, stats_file -> tuple(id, no_rrna) })
    
    // Combine outputs for stats: scrubbed files + original rrna from bbduk
    stats_input = sra_hum_scrub_output
        .map { id, no_rrna_final, no_rrna_human -> tuple(id, no_rrna_final, no_rrna_human) }
        .join(bbduk_output.map { id, no_rrna, rrna, stats_file -> tuple(id, rrna) })
    
    stats_output = stats(stats_input)
}

