nextflow.enable.dsl = 2

include { FASTP } from './modules/fastp/run_fastp'
include { BBDUK } from './modules/bbduk/run_bbduk'
include { DEACON } from './modules/deacon/run_deacon'
include { STATS } from './modules/stats/run_stats'
include { CHMOD_OUTPUT } from './modules/permissions/run_chmod'
include { WRITE_FAILURE_LIST } from './modules/failures/write_failure_list'

/*
 * Metadata input:
 *   LabKey-style tab-delimited metadata file.
 *
 * Required columns:
 *   Experiment                    -> batch directory name, e.g. 31157
 *   specimen collector sample ID  -> sample prefix, e.g. NegativeControl_20250217_S23
 *
 * If specimen collector sample ID is blank, all samples in that Experiment batch
 * are discovered based on the selected input mode.
 *
 * Input mode:
 *   default: interleaved
 *   paired:  --paired_end true
 *
 * Metadata aliases:
 *   --metadata metadat.csv
 *   --m metadat.csv
 */

def metadata_path = params.m ?: params.metadata

def requireFile(path, label, id) {
    def f = file(path)

    if (!f.exists()) {
        throw new IllegalArgumentException("Missing ${label} file for '${id}': ${path}")
    }

    return f
}

workflow {

    samples_ch = Channel
        .fromPath(metadata_path, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .flatMap { row ->

            def batch_name  = row.Experiment?.toString()?.trim()
            def sample_name = row.'specimen collector sample ID'?.toString()?.trim()
            def read_type   = params.paired_end ? 'paired' : 'interleaved'

            if (!batch_name) {
                throw new IllegalArgumentException("Metadata row is missing Experiment, which is used as the batch directory name")
            }

            def batch_dir = file("${params.read_root}/${batch_name}/${params.read_subdir}")

            if (!batch_dir.exists()) {
                throw new IllegalArgumentException("Batch directory does not exist: ${batch_dir}")
            }

            /*
             * If sample_name is present, run one sample.
             * If sample_name is blank, discover all samples in this batch.
             */
            if (sample_name) {
                return [tuple(batch_name, sample_name, read_type)]
            }

            if (read_type == 'paired') {
                return batch_dir
                    .listFiles()
                    .findAll { it.name.endsWith(params.paired_r1_suffix) }
                    .collect { r1 ->
                        def inferred_sample = r1.name.replaceFirst(/_R1\.clumped\.fastq\.gz$/, '')
                        tuple(batch_name, inferred_sample, read_type)
                    }
            }

            if (read_type == 'interleaved') {
                return batch_dir
                    .listFiles()
                    .findAll { it.name.endsWith(params.interleaved_suffix) }
                    .findAll { !(it.name.endsWith(params.paired_r1_suffix) || it.name.endsWith(params.paired_r2_suffix)) }
                    .collect { fq ->
                        def inferred_sample = fq.name.replaceFirst(/\.clumped\.fastq\.gz$/, '')
                        tuple(batch_name, inferred_sample, read_type)
                    }
            }
        }

    reads_ch = samples_ch.map { batch_name, sample_name, read_type ->

        def sample_dir = "${params.read_root}/${batch_name}/${params.read_subdir}"
        def output_dir = "${params.read_root}/${batch_name}/${params.output_subdir}/${sample_name}"

        if (read_type == 'paired') {
            def r1 = requireFile("${sample_dir}/${sample_name}${params.paired_r1_suffix}", 'R1', sample_name)
            def r2 = requireFile("${sample_dir}/${sample_name}${params.paired_r2_suffix}", 'R2', sample_name)

            return tuple(sample_name, read_type, [r1, r2], output_dir)
        }

        if (read_type == 'interleaved') {
            def reads = requireFile("${sample_dir}/${sample_name}${params.interleaved_suffix}", 'interleaved', sample_name)

            return tuple(sample_name, read_type, reads, output_dir)
        }
    }

    fastp_out = FASTP(reads_ch)

    bbduk_input = fastp_out.ok.map { id, clean_ilv, json, html, output_dir ->
        tuple(id, clean_ilv, output_dir)
    }

    bbduk_out = BBDUK(bbduk_input)

    deacon_input = bbduk_out.ok.map { id, no_rrna, rrna, stats_file, output_dir ->
        tuple(id, no_rrna, output_dir)
    }

    deacon_out = DEACON(deacon_input)

    stats_input = deacon_out.ok
        .map { id, no_rrna_final, no_rrna_human, output_dir ->
            tuple(id, no_rrna_final, no_rrna_human, output_dir)
        }
        .join(
            bbduk_out.ok.map { id, no_rrna, rrna, stats_file, output_dir ->
                tuple(id, rrna)
            }
        )
        .map { id, no_rrna_final, no_rrna_human, output_dir, rrna ->
            tuple(id, no_rrna_final, no_rrna_human, rrna, output_dir)
        }

    stats_out = STATS(stats_input)

    all_failures_ch = Channel
    .empty()
    .mix(fastp_out.failed)
    .mix(bbduk_out.failed)
    .mix(deacon_out.failed)
    .mix(stats_out.failed)

    failure_files_ch = all_failures_ch
        .map { id, process_name, failed_file, output_dir -> failed_file }
        .collect()
    
    batch_output_dir_ch = samples_ch
        .map { batch_name, sample_name, read_type ->
            "${params.read_root}/${batch_name}/${params.output_subdir}"
        }
        .unique()
        .first()
    
    WRITE_FAILURE_LIST(failure_files_ch, batch_output_dir_ch)

    sample_output_dirs_ch = stats_out.ok
        .map { id, combined_stats, output_dir -> output_dir }
        .unique()
        .collect()

    stats_done_signal_ch = stats_out.ok.collect()

    CHMOD_OUTPUT(sample_output_dirs_ch, stats_done_signal_ch)

}
