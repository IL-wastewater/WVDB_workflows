#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def read_checkv_ids(quality_summary_tsv: str, threshold: float):
    """
    Returns:
      complete_ids: set[str] where completeness > threshold and not NA
      incomplete_ids: set[str] where completeness is NA/missing or <= threshold
    """
    df = pd.read_csv(quality_summary_tsv, sep="\t", dtype=str)

    required_cols = {"contig_id", "completeness"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing columns in CheckV TSV: {sorted(missing)}. "
            f"Found columns: {list(df.columns)}"
        )

    # Normalize completeness to numeric, NA stays NaN
    comp = pd.to_numeric(df["completeness"], errors="coerce")

    contig_id = df["contig_id"].astype(str)

    complete_mask = comp.notna() & (comp > threshold)
    incomplete_mask = comp.isna() | (comp <= threshold)

    complete_ids = set(contig_id[complete_mask].tolist())
    incomplete_ids = set(contig_id[incomplete_mask].tolist())

    return complete_ids, incomplete_ids


def write_fasta_subset(in_fasta: str, id_set: set, out_fasta: str):
    """
    Streams through in_fasta and writes records whose record.id is in id_set.
    Returns number written.
    """
    written = 0
    out_path = Path(out_fasta)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w") as out_handle:
        for rec in SeqIO.parse(in_fasta, "fasta"):
            # Biopython rec.id is header up to first whitespace, matches typical tools
            if rec.id in id_set:
                SeqIO.write(rec, out_handle, "fasta")
                written += 1

    return written


def main():
    p = argparse.ArgumentParser(
        description=(
            "Select complete trimmed genomes using CheckV quality_summary.tsv, "
            "optionally write incomplete genomes from the untrimmed FASTA."
        )
    )
    p.add_argument("-q", "--quality-summary", required=True,
                   help="CheckV quality_summary.tsv")
    p.add_argument("--trimmed-fasta", required=True,
                   help="FASTA used as input to CheckV (trimmed)")
    p.add_argument("--untrimmed-fasta", required=True,
                   help="Original untrimmed FASTA")
    p.add_argument("--threshold", type=float, default=90.0,
                   help="Completeness threshold, complete if > threshold (default 90)")
    p.add_argument("-o", "--complete-trimmed-fasta", required=True,
                   help="Output FASTA of trimmed genomes with completeness > threshold")

    p.add_argument("--incomplete-untrimmed-fasta", default=None,
                   help="Optional output FASTA of untrimmed genomes that are incomplete (<= threshold or NA)")
    p.add_argument("--incomplete-ids", default=None,
                   help="Optional output text file with incomplete contig IDs (one per line)")

    args = p.parse_args()

    complete_ids, incomplete_ids = read_checkv_ids(args.quality_summary, args.threshold)

    n_complete = write_fasta_subset(args.trimmed_fasta, complete_ids, args.complete_trimmed_fasta)

    n_incomplete_fa = None
    if args.incomplete_untrimmed_fasta:
        n_incomplete_fa = write_fasta_subset(
            args.untrimmed_fasta, incomplete_ids, args.incomplete_untrimmed_fasta
        )

    if args.incomplete_ids:
        out_path = Path(args.incomplete_ids)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w") as handle:
            for cid in sorted(incomplete_ids):
                handle.write(cid + "\n")

    # Basic sanity reporting to stderr
    print(f"[checkv_select] threshold: {args.threshold}", file=sys.stderr)
    print(f"[checkv_select] complete IDs: {len(complete_ids)}", file=sys.stderr)
    print(f"[checkv_select] incomplete IDs: {len(incomplete_ids)}", file=sys.stderr)
    print(f"[checkv_select] wrote complete trimmed FASTA records: {n_complete}", file=sys.stderr)
    if n_incomplete_fa is not None:
        print(f"[checkv_select] wrote incomplete untrimmed FASTA records: {n_incomplete_fa}", file=sys.stderr)

    # Hard fail if we wrote zero complete records (often indicates ID mismatch)
    if len(complete_ids) > 0 and n_complete == 0:
        raise RuntimeError(
            "Found complete IDs in quality_summary.tsv, but wrote 0 records from trimmed_fasta. "
            "Possible ID mismatch (whitespace in FASTA headers, or different ID normalization)."
        )


if __name__ == "__main__":
    main()