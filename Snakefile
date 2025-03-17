
import os
import yaml

# Get dataset name from command-line argument
dataset = os.getenv("DATASET", None)
if dataset is None:
    raise ValueError("Please specify a dataset name using: DATASET=GSEXXXX snakemake --cores 8")

# Load config.yaml
with open("config.yaml", "r") as f:
    config_data = yaml.safe_load(f)

# Ensure the dataset exists in config.yaml
if dataset not in config_data["datasets"]:
    raise ValueError(f"Dataset '{dataset}' not found in config.yaml.")

# Extract sample information
samples_dict = config_data["datasets"][dataset]["samples"]  # Should be a dictionary {sample_name: SRA_ID}
sra_ids = list(samples_dict.values())  # Extract SRA IDs directly

rule all:
    input:
        expand("datasets/{dataset}/cellranger_output/{sra_id}/outs/filtered_feature_bc_matrix",
               dataset=dataset,
               sra_id=sra_ids),
        # Ensure final annotated Seurat object exists
        expand("datasets/{dataset}/seurat/{dataset}_annotated.rds",
               dataset=dataset),
        # Ensure final plots PDF exists
        expand("datasets/{dataset}/seurat/{dataset}_plots.pdf",
               dataset=dataset)

rule download_and_rename_fastq:
    output:
        "datasets/{dataset}/fastq/{sra_id}_S1_L001_R1_001.fastq.gz",
        "datasets/{dataset}/fastq/{sra_id}_S1_L001_R2_001.fastq.gz"
    params:
        fastq_dir="datasets/{dataset}/fastq/"
    shell:
        """
        # Ensure the FASTQ output directory exists
        mkdir -p {params.fastq_dir}

        # Download FASTQ files using fastq-dl
        fastq-dl -a {wildcards.sra_id} --cpus 8 --sra-lite --outdir {params.fastq_dir} || (echo "fastq-dl failed" && exit 1)

        # Rename downloaded files to match Cell Ranger's expected format
        mv {params.fastq_dir}/{wildcards.sra_id}_1.fastq.gz {params.fastq_dir}/{wildcards.sra_id}_S1_L001_R1_001.fastq.gz
        mv {params.fastq_dir}/{wildcards.sra_id}_2.fastq.gz {params.fastq_dir}/{wildcards.sra_id}_S1_L001_R2_001.fastq.gz
        """

rule cellranger:
    input:
        "datasets/{dataset}/fastq/{sra_id}_S1_L001_R1_001.fastq.gz",
        "datasets/{dataset}/fastq/{sra_id}_S1_L001_R2_001.fastq.gz"
    output:
        barcode="datasets/{dataset}/cellranger_output/{sra_id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        features="datasets/{dataset}/cellranger_output/{sra_id}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        matrix="datasets/{dataset}/cellranger_output/{sra_id}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
    params:
        reference="references/refdata-gex-GRCh38-2024-A",
        raw_output="{sra_id}"  # This is where Cell Ranger originally saves output
    shell:
        """
        module load cellranger/9.0.0

        # Run Cell Ranger
        cellranger count --id={wildcards.sra_id} \
                         --transcriptome={params.reference} \
                         --fastqs=datasets/{wildcards.dataset}/fastq/ \
                         --create-bam false \
                         --sample={wildcards.sra_id}

        # Ensure the final output directory exists
        mkdir -p datasets/{wildcards.dataset}/cellranger_output/{wildcards.sra_id}

        # Move results to the correct directory
        mv {params.raw_output}/outs datasets/{wildcards.dataset}/cellranger_output/{wildcards.sra_id}/

        # Cleanup: Remove the old directory if necessary
        rm -rf {params.raw_output}
        """

rule process_seurat:
    input:
        expand("datasets/{dataset}/cellranger_output/{sra_id}/outs/filtered_feature_bc_matrix",
               dataset=dataset,
               sra_id=sra_ids)
    output:
        seurat_rds = "datasets/{dataset}/seurat/{dataset}_annotated.rds",
        plots_pdf = "datasets/{dataset}/seurat/{dataset}_plots.pdf"
    params:
        seurat_script="scripts/process_seurat.R"
    shell:
        """
        Rscript {params.seurat_script} {wildcards.dataset}
        """

