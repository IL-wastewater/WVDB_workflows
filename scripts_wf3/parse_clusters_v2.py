#!/usr/bin/env python3
import argparse
from collections import defaultdict

"""
script to get contig pairs from each cluster

Modes:
  rank12: longest and second-longest contigs (default)
  rank23: second-longest and third-longest contigs, only if cluster has >= 3 contigs

Enhancement:
  In rank23 mode, optionally restrict clusters to those whose representative (rank1, rep)
  appears in a provided ID list file. This avoids upstream grep-based filtering.
"""

def read_id_set(path):
    """Read newline-delimited IDs into a set (ignores blank lines)."""
    ids = set()
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            ids.add(s)
    return ids

def get_lengths(length_file):
    """Read tab-delimited table of contigs and lengths."""
    contig_lengths = {}
    with open(length_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            contig, length = line.split('\t')
            if contig == 'name':  # drop header row if it exists
                continue
            contig_lengths[contig] = int(length)

    return contig_lengths

def get_clusters(cluster_file, contig_lengths, restrict_reps=None):
    """
    Read clusters into a dict:
      clusters[rep][0] = [list of target_name]
      clusters[rep][1] = [list of target_len]
    Note that the cluster rep is the key and is also one of the values.

    If restrict_reps is provided, only keep clusters whose rep is in the set.
    """
    clusters = defaultdict(lambda: [[], []])
    with open(cluster_file, 'r') as f:
        for line in f:
            if line.startswith('object'):
                continue
            line = line.strip()
            if not line:
                continue

            target, rep = line.split('\t')

            if restrict_reps is not None and rep not in restrict_reps:
                continue

            if target not in contig_lengths:
                raise KeyError(
                    f"Contig {target} found in cluster file but missing from length file: {cluster_file} vs {contig_lengths}"
                )

            target_len = contig_lengths[target]
            clusters[rep][0].append(target)
            clusters[rep][1].append(target_len)

    return clusters

def get_sorted_by_length(names, lengths):
    """Return a list of (name, length) sorted by length descending."""
    return sorted(zip(names, lengths), key=lambda x: x[1], reverse=True)

def make_contig_list(clusters, mode):
    """
    Build singletons list and pairs list based on mode.

    mode:
      "rank12" -> pairs are [longest, second-longest]
      "rank23" -> pairs are [second-longest, third-longest], only if cluster has >= 3 contigs
    """
    if mode == 'rank12':
        ranks = [0, 1]
        min_size_for_pair = 2
    elif mode == 'rank23':
        ranks = [1, 2]
        min_size_for_pair = 3
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    singletons = []
    pairs = []

    for rep, (names, lengths) in clusters.items():
        n = len(names)

        # singletons are always tracked, regardless of mode
        if n == 1:
            singletons.append(rep)

        # check if we have enough contigs for this mode
        if n < min_size_for_pair:
            continue

        sorted_contigs = get_sorted_by_length(names, lengths)
        # sorted_contigs is a list of (name, length), sorted descending

        # make sure requested ranks exist
        if max(ranks) >= len(sorted_contigs):
            continue

        selected_names = [sorted_contigs[r][0] for r in ranks]
        pairs.append(selected_names)

    return singletons, pairs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--cluster_file', type=str, required=True,
        help='cluster tsv file from vclust (object<TAB>rep)'
    )
    parser.add_argument(
        '-l', '--length_file', type=str, required=True,
        help='tsv with contig name and length (e.g., seqkit fx2tab -nl output)'
    )
    parser.add_argument(
        '-o', '--out_pairs', type=str, required=True,
        help='output file for contig pairs within a cluster'
    )
    parser.add_argument(
        '-s', '--out_singletons', type=str, required=True,
        help='output file for singleton contigs'
    )
    parser.add_argument(
        '-m', '--pair_mode', type=str,
        choices=['rank12', 'rank23'],
        default='rank12',
        help='pairing mode: "rank12" (longest and second-longest, default) '
             'or "rank23" (second-longest and third-longest, clusters with >=3 contigs)'
    )
    parser.add_argument(
        '--restrict-reps', type=str, default=None,
        help='Optional file with contig IDs, only clusters whose representative is in this list are kept. '
             'Intended for -m rank23 (dropout reps from branch A).'
    )

    args = parser.parse_args()

    # Guard against misuse
    if args.pair_mode != 'rank23' and args.restrict_reps:
        raise SystemExit('--restrict-reps is intended for -m rank23; omit it for rank12')

    restrict_reps = None
    if args.pair_mode == 'rank23' and args.restrict_reps:
        restrict_reps = read_id_set(args.restrict_reps)

    contig_lengths = get_lengths(args.length_file)
    clusters = get_clusters(args.cluster_file, contig_lengths, restrict_reps=restrict_reps)

    singletons, pairs = make_contig_list(clusters, args.pair_mode)

    print(f'total clusters: {len(clusters)}')
    print(f'singletons: {len(singletons)}')
    print(f'pairs ({args.pair_mode}): {len(pairs)}')
    if args.pair_mode == 'rank23' and restrict_reps is not None:
        print(f'restrict reps: {len(restrict_reps)} (clusters kept if rep in list)')

    with open(args.out_singletons, 'w') as f:
        if singletons:
            f.write('\n'.join(singletons) + '\n')

    with open(args.out_pairs, 'w') as f:
        for pair in pairs:
            f.write('\t'.join(pair) + '\n')

if __name__ == '__main__':
    main()