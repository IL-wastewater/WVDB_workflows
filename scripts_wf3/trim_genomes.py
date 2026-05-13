#!/usr/bin/env python3

import argparse
import concurrent.futures
from glob import glob
import os
import pandas as pd
import re
import subprocess
import shutil
import sys

def initialize(query_fasta, outdir):
    '''
    Check for dependencies and required files. Make output subdirectories.
    '''
    REQUIRED_PROGS = ["seqkit", "nucmer", "show-coords"]

    missing = [prog for prog in REQUIRED_PROGS if shutil.which(prog) is None]

    if missing:
        print(
        "Error: required program(s) not found in PATH: {}\n"
        "Please install them or add them to your PATH.\n".format(", ".join(missing)),
        file=sys.stderr
        )
        sys.exit(1)

    if not os.path.isfile(query_fasta):
        print(f"Error: fasta file not found: {query_fasta}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(f'{outdir}/pairs', exist_ok=True)
    os.makedirs(f'{outdir}/aln', exist_ok=True)


def get_alns(blastani_file, qcov=85, tcov=85):
    '''
    In BLAST-based trimming mode:
    Load BLAST ANI results and derive a qname/tname pairs DataFrame.
    Filters rows where qcov >= qcov or tcov >= tcov, then for each qname
    keeps the row with maximum qcov.
    ALTERNATIVELY - run pull_aln on just one row.
    '''
    ani_df = pd.read_csv(blastani_file, sep='\t')
    # allow either query cov or target cov to be > specified value 
    # this is because the query could be a chimeric misassembly, resulting in a much longer sequence with <85% covered by the target
    ani_df = ani_df[(ani_df.qcov >= qcov) | (ani_df.tcov >= tcov)]
    # from the filtered hits to each query, choose only the single target with max query coverage to use for trimming
    if ani_df.empty:
        print("Warning: no BLAST ANI hits passed coverage filters, no pairs will be generated.", file=sys.stderr)
        return pd.DataFrame(columns=["qname", "tname"])
    
    else:
        idx = ani_df.groupby('qname')['qcov'].idxmax()
        pairs_df = ani_df.loc[idx][['qname', 'tname']]
        return pairs_df

def build_seqkitfaidx(pairs_df, query_fasta, target_fasta, fastadir, threads):
    '''
    Run a single seqkit faidx command for both the query_fasta and target_fasta, to generate the faidx file.
    This command will pull a single entry from the multi-fasta file target_fasta
    This should speed up subsequent seqkit faidx runs called in pull_aln().
    Also, if this fails, we won't try to run pull_aln().
    '''
    if pairs_df.empty:
        print("No pairs to align", file=sys.stderr)
        return None

    query = pairs_df['qname'].values[0]
    target = pairs_df['tname'].values[0]
    target_newname = re.split(r'\||\ |,', target)[0]

    # 1) Test pulling one query
    cmd1 = [
        "seqkit", "faidx", query_fasta, query,
        "--threads", str(threads)
    ]
    q_out = os.path.join(fastadir, f"{query}.fasta")

    # 2) Test pulling one target
    cmd2 = [
        "seqkit", "faidx", target_fasta, target,
        "--threads", str(threads)
    ]
    t_out = os.path.join(fastadir, f"{target_newname}.fasta")

    try:
        with open(q_out, "w") as fq:
            r1 = subprocess.run(cmd1, stdout=fq, stderr=subprocess.PIPE, text=True)
        with open(t_out, "w") as ft:
            r2 = subprocess.run(cmd2, stdout=ft, stderr=subprocess.PIPE, text=True)
    except OSError as e:
        print(f"Error creating test fasta files: {e}", file=sys.stderr)
        return None

    if r1.returncode != 0:
        print(f"Error running seqkit faidx on query: {r1.stderr}", file=sys.stderr)
        return r1
    if r2.returncode != 0:
        print(f"Error running seqkit faidx on target: {r2.stderr}", file=sys.stderr)
        return r2

    return r2  # or r1, both have returncode

def pull_aln(row_id, query, target, query_fasta, target_fasta, fastadir, alndir):
    '''
    This function takes a single pair of sequences that should be aligned for use in trimming the query sequence.
    Seqkit is used to pull the target and query sequences from their respective multi-fasta files into a single fasta file. 
    Then nucmer is used to align them and get alignment coordinates.
    Note that show-coords commands are hard-coded to require 1000 nt overlapping alignment at >=90% ID.
    
    :param row: row from pairs_df
    '''
    # get a target name that can be used as a fasta file name by removing bad characters
    target_newname = re.split(r'\||\ |,', target)[0]

    qsuffix = f"{query}_{row_id}"
    tsuffix = f"{target_newname}_{row_id}"

    cmd = (
        f'seqkit faidx "{query_fasta}" "{query}" > "{fastadir}/{qsuffix}.query.fasta"; '
        f'seqkit faidx "{target_fasta}" "{target}" > "{fastadir}/{tsuffix}.target.fasta"; '
        f'nucmer -p "{alndir}/query_{qsuffix}" "{fastadir}/{qsuffix}.query.fasta" "{fastadir}/{tsuffix}.target.fasta"; '
        f'show-coords -r -c -l -L 1000 -I 90 -T "{alndir}/query_{qsuffix}.delta" > "{alndir}/query_{qsuffix}.coords"'
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error for {query}: {result.stderr}", file=sys.stderr)

    return result.returncode


def parse_nucmer(alndir):
    '''
    Read individual nucmer coords files found in the alndir
    Generate bed file that can be used to correctly trim the query sequences downstream
    '''

    files = glob(f'{alndir}/*.coords')

    if not files:
        print(f"No .coords files found in {alndir}, cannot create BED.", file=sys.stderr)
        return pd.DataFrame(columns=['ref_name', 'start1', 'end1'])
    else:
        dfs = []
        headers = ['start1','end1', 'start2', 'end2', 'len1', 'len2', 'pid',
                   'len_ref', 'len_query', 'cov_ref', 'cov_query', 'ref_name', 'query_name']
        for file in files:
            df = pd.read_csv(file, sep='\t', names=headers, skiprows=4)
            if not df.empty:          # filter out empty coords files
                dfs.append(df)

        if not dfs:
            # All coords files were empty
            print(f"Warning: all coords files in {alndir} are empty, no alignments found.", file=sys.stderr)
            return pd.DataFrame(columns=['ref_name', 'start1', 'end1'])

        aln_df = pd.concat(dfs, ignore_index=True)
        # Ensure start is always less than stop
        aln_df['start1'], aln_df['end1'] = (
            aln_df[['start1', 'end1']].min(axis=1),
            aln_df[['start1', 'end1']].max(axis=1)
        )

        # account for multiple alignments per cluster pair by taking the longest one
        idx = aln_df.groupby('ref_name')['cov_ref'].idxmax()
        aln_max_df = aln_df.loc[idx]

        # make a bed file that can be used by seqkit subseq
        bed_df = aln_max_df[['ref_name', 'start1', 'end1']].copy()

        return bed_df


def run_trimming(query_fasta, outdir, threads):
    """
    Run seqkit subseq with the BED file to trim sequences, then
    clean up contig names by removing coordinate suffixes from fasta headers.
    """
    bed_file = os.path.join(outdir, "trimming.bed")
    temp_fasta = os.path.join(outdir, "trimmed_temp.fasta")
    final_fasta = os.path.join(outdir, "trimmed.fasta")

    if not os.path.isfile(bed_file):
        print(f"Error: BED file '{bed_file}' not found, cannot trim.", file=sys.stderr)
        sys.exit(1)

    cmd_subseq = [
        "seqkit", "subseq",
        "--bed", bed_file,
        query_fasta,
        "-j", str(threads),
    ]
    try:
        with open(temp_fasta, "w") as fout:
            result = subprocess.run(
                cmd_subseq,
                stdout=fout,
                stderr=subprocess.PIPE,
                text=True,
            )
    except OSError as e:
        print(f"Error: unable to write temporary trimmed FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    if result.returncode != 0:
        print(f"Error running seqkit subseq: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Fix contig names
    pattern = re.compile(r'_\d+-\d+:\.')

    try:
        with open(temp_fasta, "r") as fin, open(final_fasta, "w") as fout:
            for line in fin:
                # Apply substitution only to header lines if desired:
                if line.startswith(">"):
                    line = pattern.sub("", line)
                fout.write(line)
    except OSError as e:
        print(f"Error rewriting trimmed FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    # Remove temp file
    try:
        os.remove(temp_fasta)
    except OSError:
        # Not fatal; just warn
        print(f"Warning: could not remove temporary file '{temp_fasta}'", file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run nucmer alignments from either cluster pairs or BLAST-ANI results."
    )

    # Common required arguments
    parser.add_argument(
        "-f", "--fasta", required=True,
        help="FASTA used as input to clustering or BLASTn"
    )
    parser.add_argument(
        "-o", "--outdir", required=True,
        help="Output directory"
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of threads to use"
    )

    # Mode specific arguments
    group = parser.add_mutually_exclusive_group(required=True)

    # Branch A: cluster mode
    group.add_argument(
        "-c", "--cluster_pairs",
        help="Cluster pairs file output from parse_clusters.py (enables cluster mode). Either -a or -c is required."
    )

    # Branch B: BLAST mode
    group.add_argument(
        "-a", "--blast_ani",
        help="BLAST ANI output file (enables BLAST mode, requires --blast_db_fasta. Either -a or -c is required.)"
    )

    parser.add_argument(
        "-d", "--blast_db_fasta",
        help="FASTA file corresponding to the BLAST database (required in BLAST mode)"
    )

    args = parser.parse_args()

    # Enforce branch B requirement: if blast_ani is given, blast_db_fasta must be given
    if args.blast_ani is not None and args.blast_db_fasta is None:
        parser.error("Argument --blast_db_fasta is required when --blast_ani is used")

    # Optionally, enforce that blast_db_fasta is NOT given when in cluster mode
    if args.cluster_pairs is not None and args.blast_db_fasta is not None:
        parser.error("--blast_db_fasta should only be used together with --blast_ani")

    return args

def main():

    args = parse_args()
    threads = args.threads
    query_fasta = args.fasta
    outdir = args.outdir

    initialize(query_fasta, outdir) # check dependencies and files; mkdirs
    # 1. Get pairs of sequences to align for downstream trimming
    if args.cluster_pairs is not None:
        # Branch A: cluster-based trimming mode
        target_fasta = query_fasta
        cluster_pairs = args.cluster_pairs
        if not os.path.isfile(cluster_pairs):
            print(f"Error: File '{cluster_pairs}' does not exist.", file=sys.stderr)
            sys.exit(1)
        # get cluster pairs, provided in file generated by parse_clusters.py
        pairs_df = pd.read_csv(cluster_pairs, sep='\t', names=['qname', 'tname'])

    else:
        # Branch B: BLAST-based trimming mode
        blast_ani = args.blast_ani
        target_fasta = args.blast_db_fasta
        files = [blast_ani, target_fasta]
        missing = [f for f in files if not os.path.isfile(f)]
        if missing:
            for f in missing:
                print(f"Error: File '{f}' does not exist.", file=sys.stderr)
            sys.exit(1)
        # get query-target pairs from the BLAST ani results generated by blastani.py
        pairs_df = get_alns(blast_ani)

    fastadir = os.path.join(outdir, "pairs")
    alndir = os.path.join(outdir, "aln")

    # 2. Test pulling of a single sequence from each multifasta and have seqkit generate faidx
    result = build_seqkitfaidx(pairs_df, query_fasta, target_fasta, fastadir, threads)
    if result is None:
        sys.exit(0)
    elif result.returncode != 0:
        print("Error: failed to build seqkit faidx file(s)", file=sys.stderr)
        sys.exit(1)

    # 3. For each pair of sequences in pairs_df, pull single fasta from multi-fasta, then align with nucmer
    ## before pulling fasta files, check that query names (qname) are unique to avoid conflicts from multiple workers trying to create the same file
    if pairs_df["qname"].duplicated().any():
        dups = pairs_df[pairs_df["qname"].duplicated()]["qname"].unique()
        print(f"Error: duplicate qname values in pairs_df: {', '.join(map(str, dups))}",
            file=sys.stderr)
        sys.exit(1)
    with concurrent.futures.ProcessPoolExecutor(max_workers=min(threads, len(pairs_df))) as executor:
        futures = []
        for idx, r in pairs_df.iterrows():
            qname = r["qname"]
            tname = r["tname"]
            futures.append(
                executor.submit(pull_aln, idx, qname, tname, query_fasta, target_fasta, fastadir, alndir)
            )
        # wait for all to finish
        for future in concurrent.futures.as_completed(futures):
            try:
                rc = future.result()
                if rc != 0:
                    print(f"Alignment task failed with code {rc}", file=sys.stderr)
            except Exception as e:
                print(f"Alignment task raised an exception: {e}", file=sys.stderr)

    # 4. Parse the nucmer coords output into a bed file for rapid downstream trimming with seqkit subseq
    bed_df = parse_nucmer(alndir)
    if bed_df.empty:
        # error message
        print(f"Error: No alignments created. BED file will not be written.", file=sys.stderr)
        sys.exit(1)
    else:
        bed_path = os.path.join(outdir, "trimming.bed")
        bed_df.to_csv(bed_path, sep='\t', index=False, header=False)

    # 5. Run trimming with seqkit subseq and clean contig names
    run_trimming(query_fasta, outdir, threads)

if __name__ == "__main__":
    main()
