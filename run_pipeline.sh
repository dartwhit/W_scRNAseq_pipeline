#!/bin/bash
#SBATCH --job-name=scRNAseq_submit_GSE264508
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --output=logs/GSE264508_submit_%j.out
#SBATCH --error=logs/GSE264508_submit_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=f005c3n@dartmouth.edu

set -euo pipefail

# Defaults (can be overridden via CLI args)
DATASET="GSE264508"
MAX_JOBS=6

usage() {
	echo "Usage: $0 [-d DATASET] [-j MAX_JOBS]"
	echo "  -d DATASET   GEO dataset ID (default: ${DATASET})"
	echo "  -j MAX_JOBS  Max concurrent Snakemake jobs (default: ${MAX_JOBS})"
	echo "  -h           Show help"
}

while getopts ":d:j:h" opt; do
	case ${opt} in
		d)
			DATASET="${OPTARG}"
			;;
		j)
			MAX_JOBS="${OPTARG}"
			;;
		h)
			usage
			exit 0
			;;
		:) 
			echo "Error: Option -${OPTARG} requires an argument."
			usage
			exit 1
			;;
		\?)
			echo "Error: Invalid option -${OPTARG}"
			usage
			exit 1
			;;
	esac
done

if ! [[ "${MAX_JOBS}" =~ ^[0-9]+$ ]] || [ "${MAX_JOBS}" -lt 1 ]; then
	echo "Error: MAX_JOBS must be a positive integer (got: ${MAX_JOBS})"
	exit 1
fi

echo "Running dataset: ${DATASET}"
echo "Max parallel jobs: ${MAX_JOBS}"

mkdir -p logs logs/slurm

source ~/.bashrc
conda activate scrna_pipeline

DATASET=${DATASET} snakemake \
	--jobs ${MAX_JOBS} \
	--keep-going \
	--rerun-incomplete \
	--latency-wait 60 \
	--default-resources mem_mb=8000 runtime=240 \
	--cluster "sbatch --partition=standard --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.runtime} --output=logs/slurm/%x-%j.out --error=logs/slurm/%x-%j.err"
