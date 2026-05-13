#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
import pysam

'''
get high-quality viral genomes from single assemblies, including trimmed proviruses
1) checkV summary: filter contigs on contig length >2kb and completeness >90 -> virus and provirus lists
2) geNomad summary: filter contigs on virus_score >0.7 -> virus and provirus lists
    * note: this filter should not usually remove anything more because checkV >90 is more stringent
3) get list intersections and pull contigs:
    * virus by checkV and virus/provirus by genomad: pull from checkv_out/virus.fa
    * provirus by checkV, already trimmed: pull from checkv_out/provirus.fa
'''

def filter_genomad(genomad_summary):
    genomad_df = pd.read_csv(genomad_summary, sep='\t', low_memory=False)

    # Guard for missing required columns
    required_cols = {'virus_score', 'topology', 'seq_name'}
    if not required_cols.issubset(genomad_df.columns):
        missing = required_cols - set(genomad_df.columns)
        raise ValueError(f"genoMad summary missing required columns: {', '.join(missing)}")

    genomad_df = genomad_df[genomad_df.virus_score > 0.7]

    # contig names of proviruses have additional text, fix this
    genomad_virus_set = set(
        genomad_df[genomad_df.topology != 'Provirus'].seq_name
    )

    genomad_provirus_df = genomad_df[genomad_df.topology == 'Provirus'].copy()
    # Extract contig_id from provirus seq_name; keep only the part before '|'
    genomad_provirus_df['contig_id'] = (
        genomad_provirus_df.seq_name
        .astype(str)
        .str.extract(r'^(.*)\|', expand=False)
    )

    genomad_provirus_set = set(
        genomad_provirus_df['contig_id'].dropna()
    )
    genomad_virus_set = genomad_virus_set.union(genomad_provirus_set)

    return genomad_virus_set


def filter_checkv(checkv_summary):
    # If file exists but is empty, treat as no viruses found
    if os.path.getsize(checkv_summary) == 0:
        print(
            f"Warning: CheckV quality summary '{checkv_summary}' is empty, "
            "treating as no high-quality viruses detected.",
            file=sys.stderr
        )
        return set(), set()

    checkv_df = pd.read_csv(checkv_summary, sep='\t')

    # Guard for missing required columns
    required_cols = {'contig_length', 'completeness', 'contig_id', 'provirus'}
    if not required_cols.issubset(checkv_df.columns):
        missing = required_cols - set(checkv_df.columns)
        raise ValueError(f"CheckV summary missing required columns: {', '.join(missing)}")

    checkv_df = checkv_df[
        (checkv_df.contig_length > 2000) & (checkv_df.completeness > 90)
    ]

    if checkv_df.empty:
        print(
            "Warning: After filtering on contig_length > 2000 and completeness > 90, "
            "no CheckV contigs remain.",
            file=sys.stderr
        )
        return set(), set()

    checkv_virus_set = set(checkv_df[checkv_df.provirus == 'No'].contig_id)
    checkv_provirus_set = set(checkv_df[checkv_df.provirus == 'Yes'].contig_id)

    return checkv_virus_set, checkv_provirus_set


def pull_viruses(checkv_virus_fasta, filtered_virus_fasta, filtered_virus_set):
    if not filtered_virus_set:
        print(
            "Info: No high-quality non-provirus viruses to write.",
            file=sys.stderr
        )
        # Create an empty file so downstream steps can rely on its existence
        open(filtered_virus_fasta, "w").close()
        return

    if not os.path.isfile(checkv_virus_fasta):
        print(
            f"Warning: CheckV virus FASTA '{checkv_virus_fasta}' not found. "
            "Cannot write filtered virus sequences.",
            file=sys.stderr
        )
        # Still create the output file as empty
        open(filtered_virus_fasta, "w").close()
        return

    fasta = pysam.FastaFile(checkv_virus_fasta)
    try:
        with open(filtered_virus_fasta, "w") as out:
            for ref in filtered_virus_set:
                if ref not in fasta.references:
                    print(
                        f"Warning: Virus contig '{ref}' not found in '{checkv_virus_fasta}'.",
                        file=sys.stderr
                    )
                    continue
                seq = fasta.fetch(ref)
                out.write(f">{ref}\n{seq}\n")
    finally:
        fasta.close()


def pull_proviruses(checkv_provirus_fasta, filtered_provirus_fasta, filtered_provirus_set):
    """
    Pull filtered_provirus_set from checkv_provirus_fasta.

    Example:
        quality_summary.tsv contig_id:  LosAngeles_CA_20231203_322
        provirus.fa header:            >LosAngeles_CA_20231203_322_1 1-552/962

    pysam will see the reference name as "LosAngeles_CA_20231203_322_1".
    We match by removing the last underscore segment:
        "LosAngeles_CA_20231203_322_1" -> "LosAngeles_CA_20231203_322"

    If that prefix is in filtered_provirus_set, we write the full sequence
    using the full reference name.
    """
    if not filtered_provirus_set:
        print(
            "Info: No high-quality proviruses to write.",
            file=sys.stderr
        )
        open(filtered_provirus_fasta, "w").close()
        return

    if not os.path.isfile(checkv_provirus_fasta):
        print(
            f"Warning: CheckV provirus FASTA '{checkv_provirus_fasta}' not found. "
            "Cannot write filtered provirus sequences.",
            file=sys.stderr
        )
        open(filtered_provirus_fasta, "w").close()
        return

    fasta = pysam.FastaFile(checkv_provirus_fasta)
    try:
        with open(filtered_provirus_fasta, "w") as out:
            for ref in fasta.references:
                # ref is something like "LosAngeles_CA_20231203_322_1"
                parts = ref.split('_')
                if len(parts) <= 1:
                    prefix = ref
                else:
                    prefix = '_'.join(parts[:-1])

                if prefix in filtered_provirus_set:
                    seq = fasta.fetch(ref)
                    out.write(f">{ref}\n{seq}\n")
    finally:
        fasta.close()


def main():
    parser = argparse.ArgumentParser(
        description="Filter high-quality viral genomes from CheckV and geNomad outputs."
    )
    parser.add_argument(
        '-g', '--genomad_virus_summary', required=True, type=str,
        help='geNomad virus_summary.tsv'
    )
    parser.add_argument(
        '-c', '--checkv_dir', required=True, type=str,
        help='CheckV output directory containing quality_summary.tsv, virus.fa, and provirus.fa'
    )
    parser.add_argument(
        '-o', '--outdir', required=True, type=str,
        help='Output directory for filtered FASTA files'
    )

    args = parser.parse_args()
    genomad_summary = args.genomad_virus_summary
    checkv_dir = args.checkv_dir
    outdir = args.outdir

    # Construct paths inside the CheckV directory
    checkv_summary = os.path.join(checkv_dir, 'quality_summary.tsv')
    checkv_virus_fasta = os.path.join(checkv_dir, 'viruses.fna')
    checkv_provirus_fasta = os.path.join(checkv_dir, 'proviruses.fna')

    filtered_virus_fasta = os.path.join(outdir, 'filtered_virus.fasta')
    filtered_provirus_fasta = os.path.join(outdir, 'filtered_provirus.fasta')

    # Check output directory
    if not os.path.isdir(outdir):
        print(f"Error: Output directory '{outdir}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Check that genoMad summary exists
    if not os.path.isfile(genomad_summary):
        print(
            f"Error: geNomad summary file '{genomad_summary}' does not exist.",
            file=sys.stderr
        )
        sys.exit(1)

    # Check that CheckV directory exists
    if not os.path.isdir(checkv_dir):
        print(
            f"Error: CheckV directory '{checkv_dir}' does not exist.",
            file=sys.stderr
        )
        sys.exit(1)

    # quality_summary.tsv is required, even if empty
    if not os.path.isfile(checkv_summary):
        print(
            f"Error: CheckV quality summary '{checkv_summary}' does not exist.",
            file=sys.stderr
        )
        sys.exit(1)

    # Now run the filtering logic
    try:
        genomad_virus_set = filter_genomad(genomad_summary)
    except Exception as e:
        print(f"Error while filtering geNomad summary: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        checkv_virus_set, checkv_provirus_set = filter_checkv(checkv_summary)
    except Exception as e:
        print(f"Error while filtering CheckV summary: {e}", file=sys.stderr)
        sys.exit(1)

    # High quality non-proviruses
    filtered_virus_set = genomad_virus_set & checkv_virus_set
    # High quality proviruses
    filtered_provirus_set = genomad_virus_set & checkv_provirus_set

    # Write FASTA outputs; these functions already handle missing FASTA files
    pull_viruses(checkv_virus_fasta, filtered_virus_fasta, filtered_virus_set)
    pull_proviruses(checkv_provirus_fasta, filtered_provirus_fasta, filtered_provirus_set)


if __name__ == '__main__':
    main()