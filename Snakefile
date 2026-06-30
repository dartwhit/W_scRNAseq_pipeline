
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

# Extract sample information: {GSM: [SRR1, SRR2, ...]}
samples_dict = config_data["datasets"][dataset]["samples"]

# Normalize to lists (supports both string and list values in config)
for gsm in samples_dict:
    if isinstance(samples_dict[gsm], str):
        samples_dict[gsm] = [samples_dict[gsm]]

# Build flat lookup: (gsm, lane_str) -> sra_id
# Lane is 1-indexed and zero-padded to 3 digits (e.g. "001", "002")
gsm_lane_to_sra = {
    (gsm, f"{i:03d}"): sra_id
    for gsm, srr_list in samples_dict.items()
    for i, sra_id in enumerate(srr_list, start=1)
}

gsm_ids = list(samples_dict.keys())


rule all:
    input:
        expand("datasets/{dataset}/cellranger_output/{gsm}/outs/filtered_feature_bc_matrix",
               dataset=dataset,
               gsm=gsm_ids),
        expand("datasets/{dataset}/seurat/{dataset}_annotated.rds",
               dataset=dataset),
        expand("datasets/{dataset}/seurat/{dataset}_plots.pdf",
               dataset=dataset)


rule download_and_rename_fastq:
    threads: 4
    resources:
        mem_mb=16000,
        runtime=180
    output:
        r1="datasets/{dataset}/fastq/{gsm}_S1_L{lane}_R1_001.fastq.gz",
        r2="datasets/{dataset}/fastq/{gsm}_S1_L{lane}_R2_001.fastq.gz"
    params:
        fastq_dir="datasets/{dataset}/fastq/",
        sra_id=lambda wildcards: gsm_lane_to_sra[(wildcards.gsm, wildcards.lane)]
    shell:
        """
        # Ensure the FASTQ output directory exists
        mkdir -p {params.fastq_dir}

        # Download FASTQ files using fastq-dl
        fastq-dl -a {params.sra_id} --cpus {threads} --outdir {params.fastq_dir} || (echo "fastq-dl failed" && exit 1)

        # Rename downloaded files to match Cell Ranger's expected format
        mv {params.fastq_dir}/{params.sra_id}_1.fastq.gz {output.r1}
        mv {params.fastq_dir}/{params.sra_id}_2.fastq.gz {output.r2}
        """


rule cellranger:
    threads: 16
    resources:
        mem_mb=128000,
        runtime=1440
    input:
        lambda wildcards: expand(
            "datasets/{dataset}/fastq/{gsm}_S1_L{lane}_R{read}_001.fastq.gz",
            dataset=wildcards.dataset,
            gsm=wildcards.gsm,
            lane=[f"{i:03d}" for i in range(1, len(samples_dict[wildcards.gsm]) + 1)],
            read=["1", "2"]
        )
    output:
        directory("datasets/{dataset}/cellranger_output/{gsm}/outs/filtered_feature_bc_matrix")
    params:
        reference="references/refdata-gex-GRCh38-2024-A",
    shell:
        """
        module load cellranger/9.0.0

        # Run Cell Ranger (fastqs dir contains all lanes; --sample filters by GSM prefix)
        cellranger count --id={wildcards.gsm} \
                         --transcriptome={params.reference} \
                         --fastqs=datasets/{wildcards.dataset}/fastq/ \
                         --localcores={threads} \
                         --localmem=120 \
                         --create-bam false \
                         --sample={wildcards.gsm}

        # Move results to the expected output location
        mkdir -p datasets/{wildcards.dataset}/cellranger_output/{wildcards.gsm}
        mv {wildcards.gsm}/outs datasets/{wildcards.dataset}/cellranger_output/{wildcards.gsm}/

        # Cleanup Cell Ranger working directory
        rm -rf {wildcards.gsm}
        """


rule process_seurat:
    threads: 8
    resources:
        mem_mb=64000,
        runtime=720
    input:
        expand("datasets/{dataset}/cellranger_output/{gsm}/outs/filtered_feature_bc_matrix",
               dataset=dataset,
               gsm=gsm_ids)
    output:
        seurat_rds="datasets/{dataset}/seurat/{dataset}_annotated.rds",
        plots_pdf="datasets/{dataset}/seurat/{dataset}_plots.pdf"
    params:
        seurat_script="scripts/process_seurat.R"
    shell:
        """
        Rscript {params.seurat_script} {wildcards.dataset}
        """
