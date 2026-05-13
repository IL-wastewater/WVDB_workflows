#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -o stdout_wf3_%j
#SBATCH -e stderr_wf3_%j
#SBATCH --job-name wf3_cluster
##SBATCH -c 100 # required if the system is not 1-job-per-node

# usage:
# sbatch your_script.sbatch /path/to/params_wf3.sh

set -euo pipefail
IFS=$'\n\t'

# load params file
PARAMS_FILE="$1"

if [[ -z "$PARAMS_FILE" ]]; then
    echo "ERROR: Please provide a params file"
    echo "Usage: sbatch your_script.sbatch /path/to/params_wf3.sh"
    exit 1
fi

if [[ ! -f "$PARAMS_FILE" ]]; then
    echo "ERROR: Params file not found: $PARAMS_FILE"
    exit 1
fi

source "$PARAMS_FILE"

# check dependencies
# ensure the following are in your path (versions included for reference): 
#   seqkit v2.9.0
#   vclust v1.3.1
#   blastn v2.16.0+
#   checkv v1.0.3
#   nucmer v4.0.1
#   /path/to/scripts_wf3 (this repo)

# check for dependencies and fail fast
command -v python >/dev/null
command -v vclust >/dev/null
command -v seqkit >/dev/null
command -v nucmer >/dev/null
command -v blastn >/dev/null
command -v checkv >/dev/null

# record versions to stdout
python --version
vclust -v
seqkit -h | head -n 2
blastn -version
checkv -h | head -n 1

# python helper scripts must be located within the scripts_dir: 
#   blastani_nayfach.py, parse_clusters_v2.py, trim_genomes.py, completeness_filter.py
blastani=blastani_nayfach.py
parse_clusters=parse_clusters_v2.py
trim_genomes=trim_genomes.py
completeness_filter=completeness_filter.py

# create dirs
filterdir="$workdir/1_filtered"
clusterdir="$workdir/2_clustered"
branchAdir="$workdir/3_branchA"
branchBdir="$workdir/4_branchB"
branchCdir="$workdir/5_branchC"
reclusterdir="$workdir/6_reclustered"

mkdir -p "$filterdir"
mkdir -p "$clusterdir"
mkdir -p "$branchAdir"
mkdir -p "$branchBdir"
mkdir -p "$branchCdir"
mkdir -p "$reclusterdir"

# outputs
fasta="$filterdir"/filtered_all.fasta

### 0. Run get_HQ_viruses_v1.py on assemblies (backfill for all existing assemblies, do this before running)

### 1. Collect all input genomes ###
shopt -s nullglob # safety in case no matches to the glob
virus_files=( "$fastqdir"/*/assembly/*/filtered_virus.fasta )
provirus_files=( "$fastqdir"/*/assembly/*/filtered_provirus.fasta )
((${#virus_files[@]})) || { echo "ERROR: no filtered_virus.fasta found" >&2; exit 1; }
((${#provirus_files[@]})) || { echo "ERROR: no filtered_provirus.fasta found" >&2; exit 1; }
cat "${virus_files[@]}" > "$filterdir/filtered_virus_all.fasta"
cat "${provirus_files[@]}" > "$filterdir/filtered_provirus_all.fasta"

cat "$filterdir/filtered_virus_all.fasta" "$filterdir/filtered_provirus_all.fasta" > "$fasta"

# get lengths
seqkit fx2tab -nl "$filterdir/filtered_all.fasta" > "$filterdir/filtered_all.length.txt"

# check
[[ -s "$fasta" ]] || { echo "ERROR: merged FASTA is empty: $fasta" >&2; exit 1; }
[[ -s "$filterdir/filtered_all.length.txt" ]] || { echo "ERROR: length file empty" >&2; exit 1; }

### 2. Cluster input genomes ###
vclust prefilter -i "$fasta" -o "$clusterdir/vclust_prefilter.txt" --min-ident "$ani" --threads "$threads"
vclust align --filter "$clusterdir/vclust_prefilter.txt" -i "$fasta" \
--out-ani "$ani" --out-qcov "$qcov" -o "$clusterdir/vclust_ani.tsv" \
--out-aln "$clusterdir/vclust_ani.aln.tsv" --threads "$threads"
vclust cluster -i "$clusterdir/vclust_ani.tsv" -o "$clusterdir/vclust_clusters.tsv" \
--ids "$clusterdir/vclust_ani.ids.tsv" --algorithm leiden \
--metric ani --ani "$ani" --qcov "$qcov" --out-repr

### 3. For each cluster, trim longest contig against 2nd longest contig (branch A) ###
### 3.1 Parse clusters into longest contig pairs for trimming ###
"$parse_clusters" \
-c "$clusterdir/vclust_clusters.tsv" -l "$filterdir/filtered_all.length.txt" \
-o "$clusterdir/cluster_pairs.tsv" -s "$clusterdir/cluster_singletons.txt" \
-m rank12

### 3.2. Trim centroid genomes against cluster mates ###
"$trim_genomes" -c "$clusterdir/cluster_pairs.tsv" -f "$fasta" -o "$branchAdir" -t "$threads"

### 3.3. Check quality of trimmed genomes ###
checkv end_to_end -d "checkvdb" -t "$threads" --remove_tmp "$branchAdir/trimmed.fasta" "$branchAdir/checkv_out"

### 3.4. collect genomes based on completeness ###
"$completeness_filter" \
  -q "$branchAdir/checkv_out/quality_summary.tsv" \
  --trimmed-fasta "$branchAdir/trimmed.fasta" \
  --untrimmed-fasta "$fasta" \
  --threshold "$completeness" \
  -o "$branchAdir/cluster_reps_complete.fasta" \
  --incomplete-ids "$branchAdir/cluster_reps_incomplete_or_na.txt"

### 4. Perform second round of trimming 2nd longest against 3rd longest contig (branch B) ###

### 4.1 Parse clusters into contig pairs for trimming  
# use mode "rank23" to get second and third longest contigs from each cluster that dropped out of branch A
"$parse_clusters" \
  -c "$clusterdir/vclust_clusters.tsv" \
  -l "$filterdir/filtered_all.length.txt" \
  -o "$branchBdir/cluster_pairs23.tsv" \
  -s "$branchBdir/cluster_singletons23.tsv" \
  -m rank23 \
  --restrict-reps "$branchAdir/cluster_reps_incomplete_or_na.txt"

### 4.2. Trim centroid genomes against cluster mates ###
"$trim_genomes" -c "$branchBdir/cluster_pairs23.tsv" -f "$fasta" -o "$branchBdir" -t "$threads"

### 4.3. Check quality of trimmed genomes ###
checkv end_to_end -d "checkvdb" -t "$threads" --remove_tmp "$branchBdir/trimmed.fasta" "$branchBdir/checkv_out"

### 4.4. collect genomes based on completeness ###
"$completeness_filter" \
  -q "$branchBdir/checkv_out/quality_summary.tsv" \
  --trimmed-fasta "$branchBdir/trimmed.fasta" \
  --untrimmed-fasta "$fasta" \
  --threshold "$completeness" \
  -o "$branchBdir/cluster_reps23_complete.fasta" \
  --incomplete-ids "$branchBdir/cluster_reps_notcomplete_or_na_after23.txt" \
  --incomplete-untrimmed-fasta "$branchBdir/cluster_reps_notcomplete_after23.fasta"

### 5. Trimming of singletons and genomes incomplete after cluster-mode trimming (branch C) ###
### 5.1 collect contigs that still need trimming

seqkit grep -f "$clusterdir/cluster_singletons.txt" "$fasta" -o "$branchCdir/singletons.fasta"

cat "$branchCdir/singletons.fasta" "$branchBdir/cluster_reps_notcomplete_after23.fasta" > "$branchCdir/branchC_genomes.fasta"

# 5.2 BLASTN against RefSeq and curated EsViritu ref db
out_tsv="$branchCdir/branchC_genomes-vs-esviritu_plus_refseq_virus.tsv"
ani_tsv="$branchCdir/branchC_genomes-vs-esviritu_plus_refseq_virus.ani.tsv"

blastn -query "$branchCdir/branchC_genomes.fasta" -db "$refseq_ev_blastdb" \
-max_target_seqs 10 -perc_identity 90 \
-num_threads "$threads" -out "$out_tsv" -outfmt "6 std qlen slen"
"$blastani" -i "$out_tsv" -o "$ani_tsv"

mkdir -p "$branchCdir/refseq_esviritu"

"$trim_genomes" -a $ani_tsv -f "$fasta" -d "$refseq_ev_blastdb"  -o "$branchCdir/refseq_esviritu" -t "$threads"

# 5.3 BLASTN against IMG-VR v5 (metaVR)
out_tsv="$branchCdir/branchC_genomes-vs-metaVR.tsv"
ani_tsv="$branchCdir/branchC_genomes-vs-metaVR.ani.tsv"

blastn -query "$branchCdir/branchC_genomes.fasta" -db "$IMGVR_blastdb" \
-max_target_seqs 10 -perc_identity 90 \
-num_threads "$threads" -out "$out_tsv" -outfmt "6 std qlen slen"
"$blastani" -i "$out_tsv" -o "$ani_tsv"

mkdir -p "$branchCdir/metavr"

"$trim_genomes" -a "$ani_tsv" -f "$fasta" -d "$IMGVR_blastdb"  -o "$branchCdir/metavr" -t "$threads"

### 5.4 collect: keep all trimmed against refseq+esviritu dbs plus any additional against metavr ###
awk '{print $1}' "$branchCdir/refseq_esviritu/trimming.bed" | LC_ALL=C sort > "$branchCdir/refseq_esviritu/trimmed_hits_sorted.txt"
awk '{print $1}' "$branchCdir/metavr/trimming.bed" | LC_ALL=C sort > "$branchCdir/metavr/trimmed_hits_sorted.txt"
LC_ALL=C comm -13 "$branchCdir/refseq_esviritu/trimmed_hits_sorted.txt" "$branchCdir/metavr/trimmed_hits_sorted.txt" > "$branchCdir/metavr/trimmed_hits_metavr-only.txt"
seqkit grep -f "$branchCdir/metavr/trimmed_hits_metavr-only.txt" "$branchCdir/metavr/trimmed.fasta" > "$branchCdir/metavr/trimmed_metavr-only.fasta"
cat "$branchCdir/refseq_esviritu/trimmed.fasta" "$branchCdir/metavr/trimmed_metavr-only.fasta" > "$branchCdir/trimmed.fasta"

### 5.5. Check quality of trimmed genomes ###
checkv end_to_end -d "checkvdb" -t "$threads" --remove_tmp "$branchCdir/trimmed.fasta" "$branchCdir/checkv_out"

### 5.6. collect genomes based on completeness after blast-mode trimming ###
"$completeness_filter" \
  -q "$branchCdir/checkv_out/quality_summary.tsv" \
  --trimmed-fasta "$branchCdir/trimmed.fasta" \
  --untrimmed-fasta "$fasta" \
  --threshold "$completeness" \
  -o "$branchCdir/blastmode_complete.fasta" \
  --incomplete-ids "$branchCdir/blastmode_incomplete_or_na.txt"

### 6. Recluster ###

# note: manual checkpoint, not coding for now
# did extra trimming of two contigs (one that was a complete duplicate and one that was huge with provirus only detected by geNomad)

cat "$branchAdir/cluster_reps_complete.fasta" "$branchBdir/cluster_reps23_complete.fasta" "$branchCdir/blastmode_complete.fasta" > "$reclusterdir/recluster_input.fasta"

vclust prefilter -i "$reclusterdir/recluster_input.fasta" -o "$reclusterdir/vclust_prefilter.txt" --min-ident "$ani" --threads "$threads"
vclust align --filter "$reclusterdir/vclust_prefilter.txt" -i "$reclusterdir/recluster_input.fasta" \
--out-ani "$ani" --out-qcov "$qcov" -o "$reclusterdir/vclust_ani.tsv" \
--out-aln "$reclusterdir/vclust_ani.aln.tsv" --threads "$threads"
vclust cluster -i "$reclusterdir/vclust_ani.tsv" -o "$reclusterdir/vclust_clusters.tsv" \
--ids "$reclusterdir/vclust_ani.ids.tsv" --algorithm leiden \
--metric ani --ani "$ani" --qcov "$qcov" --out-repr

# get centroid fasta
awk -F'\t' 'NR>1 {print $2}' "$reclusterdir/vclust_clusters.tsv" | LC_ALL=C sort -u > "$reclusterdir/vclust_centroids.txt"
seqkit grep -f "$reclusterdir/vclust_centroids.txt" "$reclusterdir/recluster_input.fasta" > "$reclusterdir/vclust_centroids.fasta"
