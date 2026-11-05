nextflow.enable.dsl = 2

include { MEGAHIT } from './modules/megahit/run_megahit.nf'
include { POST_PROCESS } from './modules/post_process/run_post_process.nf'
include { CHECKV } from './modules/checkv/run_checkv.nf'
include { GENOMAD } from './modules/genomad/run_genomad.nf'
include { CHMOD_OUTPUT } from './modules/permissions/run_chmod.nf'
include { WRITE_FAILURE_LIST } from './modules/failures/write_failure_list.nf'

workflow {
    def metadata_path = params.m ?: params.metadata
    def metadata_file = file(metadata_path, checkIfExists: true)
    def metadata_rows = metadata_file
        .readLines()
        .findAll { it?.trim() }

    if (metadata_rows.size() < 2) {
        throw new IllegalArgumentException("Metadata file must contain a header and at least one data row: ${metadata_path}")
    }

    def metadata_header = metadata_rows[0].split('\t', -1)
    def experiment_idx = metadata_header.findIndexOf { it == 'Experiment' }

    if (experiment_idx < 0) {
        throw new IllegalArgumentException('Metadata file is missing required column: Experiment')
    }

    def unique_experiments = metadata_rows
        .drop(1)
        .collect { line -> line.split('\t', -1)[experiment_idx]?.trim() }
        .findAll { it && it != 'nan' && it != 'NaN' }
        .unique()

    if (unique_experiments.size() != 1) {
        throw new IllegalArgumentException("Metadata file must contain exactly one non-empty Experiment per run. Found: ${unique_experiments}")
    }

    samples_ch = channel
        .fromPath(metadata_path, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .flatMap { row ->
            def sample_name = row.'specimen collector sample ID'?.toString()?.trim()
            def batch_name = row.Experiment?.toString()?.trim()

            if (!sample_name || sample_name == 'nan' || sample_name == 'NaN') {
                return []
            }

            if (!batch_name || batch_name == 'nan' || batch_name == 'NaN') {
                batch_name = sample_name
            }

            return [tuple(batch_name, sample_name)]
        }
        .unique()

    reads_ch = samples_ch.flatMap { batch_name, sample_name ->
        def sample_dir = "${params.read_root}/${batch_name}/${params.read_subdir}/${sample_name}"
        def output_dir = "${params.read_root}/${batch_name}/${params.output_subdir}/${sample_name}"
        def reads = file("${sample_dir}/${sample_name}${params.final_reads_suffix}")

        if (!reads.exists()) {
            return []
        }

        return [tuple(sample_name, reads, output_dir)]
    }

    output_dirs_ch = samples_ch
        .map { batch_name, sample_name ->
            "${params.read_root}/${batch_name}/${params.output_subdir}"
        }
        .unique()
    MEGAHIT(reads_ch)

    POST_PROCESS(MEGAHIT.out.ok)

    viral_input_ch = POST_PROCESS.out.ok.map { id, assembly_500, assembly_2000, output_dir ->
        tuple(id, assembly_500, output_dir)
    }

    CHECKV(viral_input_ch)

    GENOMAD(viral_input_ch)

    failure_ch = MEGAHIT.out.failed
        .mix(POST_PROCESS.out.failed)
        .mix(CHECKV.out.failed)
        .mix(GENOMAD.out.failed)

    output_dirs_list_ch = output_dirs_ch.collect()

    completion_ch = CHECKV.out.ok
        .mix(CHECKV.out.failed)
        .mix(GENOMAD.out.ok)
        .mix(GENOMAD.out.failed)
        .collect()

    WRITE_FAILURE_LIST(failure_ch.collect(), output_dirs_ch)

    CHMOD_OUTPUT(output_dirs_list_ch, completion_ch)
}
