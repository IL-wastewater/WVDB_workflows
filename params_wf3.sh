# params and thresholds

threads=100

completeness=90
ani=0.95
qcov=0.85

# paths
## specify top-level data dir used by wf2 (fastq/)
fastqdir="/path/to/fastq"
## specify workdir for intermediate and final outputs
workdir="/path/to/workdir"

# add scripts to path
export PATH=$PATH:/path/to/repo/scripts_wf3

# db paths
## checkv database
checkvdb="/path/to/checkv_db/checkv-db-v1.5"

## Concatenate NCBI RefSeq virus and EsViritu reference database, then make blastn db located at the path below
## Alternatively, provide a high-priority reference database of your choice for use in reference-based trimming
refseq_ev_blastdb="/path/to/esviritu_plus_refseq_virus.fna"

## download MetaVR (IMG-VR5) and make a blastn db located at the path below:
IMGVR_blastdb="/path/to/IMG-VR_2025-12-02/IMGVR5_UViG.fna"
