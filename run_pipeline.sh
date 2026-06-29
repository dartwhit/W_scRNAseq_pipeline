#!/bin/bash
#SBATCH --job-name=scRNAseq_GSE264508
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/GSE264508_%j.out
#SBATCH --error=logs/GSE264508_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=f005c3n@dartmouth.edu

mkdir -p logs

source ~/.bashrc
conda activate snakemake

DATASET=GSE264508 snakemake --cores 16
