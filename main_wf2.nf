#!/usr/bin/env nextflow

params.reads    = params.reads    ?: "${params.read_dir}/*/*_no_rrna_final_reads.fastq.gz"

reads_ch = Channel
  .fromPath(params.reads)
  .map { file ->
    def sample_id = file.name.replaceFirst(/_no_rrna_final_reads\.fastq\.gz$/, '')
    tuple(sample_id, file)
  }

process megahit_assemble {
    tag "$sample_id"
    time params.time
    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'link'
    conda params.megahit_env

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), file("${sample_id}.contigs.fa"), file("${sample_id}.log")

    script:
    """
    /path/to/envs/megahit/megahit --12 ${reads}  --num-cpu-threads ${params.threads} --out-prefix ${sample_id}
    mv megahit_out/* .
    """
}

process post_process {
    tag "$sample_id"
    time params.time
    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'link'
    conda params.seqkit_env    

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id),
          file("${sample_id}_contigs_min500.fa"),
          file("${sample_id}_contigs_min2000.fa"),
          path("has_min500.txt")

    script:
    """
    seqkit seq -m 500 ${assembly} | seqkit replace -p .+ -r "${sample_id}_{nr}" > ${sample_id}_contigs_min500.fa
    seqkit seq -m 2000 ${sample_id}_contigs_min500.fa > ${sample_id}_contigs_min2000.fa

    if [ -s ${sample_id}_contigs_min500.fa ]; then
        echo true > has_min500.txt
    else
        echo false > has_min500.txt
    fi
    """
}

process checkv {
    tag "$sample_id"
    time params.time
    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'move'
    conda params.checkv_env

    input:
    tuple val(sample_id), path(assembly_500)

    output:
    tuple val(sample_id), path("checkv_out")

    script:
    """
    ${params.checkv_env}/bin/checkv end_to_end ${assembly_500} checkv_out -d ${params.checkv_db} -t ${params.threads}
    """
}

process genomad {
    tag "$sample_id"
    time params.time
    publishDir { "${params.out_dir}/${sample_id}" }, mode: 'move'

    beforeScript = """
        export PATH=${params.genomad_env}/bin:\$PATH
    """

    input:
    tuple val(sample_id), path(assembly_500)

    output:
    tuple val(sample_id), path("${sample_id}*")

    script:
    """
    ${params.genomad_env}/bin/genomad end-to-end --cleanup ${assembly_500} genomad_output ${params.genomad_db} -t ${params.threads}
    mv genomad_output/* .
    """
}

workflow {
    megahit_output = megahit_assemble(reads_ch)
    post_process_output = post_process(megahit_output.map { sample_id, assembly, log -> tuple(sample_id, assembly) })
    post_process_min500_nonempty = post_process_output.filter { sample_id, assembly_500, assembly_2000, has_min500_file -> has_min500_file.text.trim() == 'true' }
    checkv_output = checkv(post_process_min500_nonempty.map { sample_id, assembly_500, assembly_2000, has_min500_file -> tuple(sample_id, assembly_500) })
    genomad_output = genomad(post_process_min500_nonempty.map { sample_id, assembly_500, assembly_2000, has_min500_file -> tuple(sample_id, assembly_500) })
}
