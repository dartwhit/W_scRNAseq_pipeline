# W_scRNAseq_pipeline

A comprehensive Snakemake pipeline for processing single-cell RNA sequencing (scRNA-seq) data. This pipeline automates the entire workflow from raw FASTQ files to annotated Seurat objects with cell type identification.

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline Workflow](#pipeline-workflow)
- [Output](#output)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)

## 🔬 Overview

This pipeline processes scRNA-seq data from public repositories (e.g., GEO datasets) through the following stages:
1. Download raw FASTQ files from SRA (Sequence Read Archive)
2. Process with Cell Ranger for alignment and quantification
3. Quality control and doublet detection with scDblFinder
4. Cell type annotation using SingleR and CAMML
5. Generate comprehensive visualization and analysis reports

## ✨ Features

- **Automated Data Retrieval**: Downloads FASTQ files directly from SRA using fastq-dl
- **Cell Ranger Integration**: Automatic alignment and feature counting
- **Quality Control**: Filters cells based on feature counts and mitochondrial content
- **Doublet Detection**: Uses scDblFinder for identifying and removing doublets
- **Cell Type Annotation**: Dual annotation approach with SingleR and CAMML
- **Comprehensive Visualization**: Generates UMAP plots, feature plots, and heatmaps
- **Scalable**: Handles multiple samples across different datasets
- **Reproducible**: Snakemake ensures reproducible and parallelizable execution

## 📦 Requirements

### Software Dependencies

- **Python 3.7+** with packages:
  - `snakemake`
  - `pyyaml`
  - `pandas`
  - `pysradb`
  
- **R 4.0+** with packages:
  - `Seurat`
  - `SingleR`
  - `celldex`
  - `CAMML`
  - `scDblFinder`
  - `SingleCellExperiment`
  - `pheatmap`
  - `ggplot2`
  - `yaml`
  - `dplyr`

- **External Tools**:
  - `Cell Ranger` (v9.0.0 or compatible)
  - `fastq-dl` (for SRA data download)

### Reference Genome

- Human reference genome (GRCh38) for Cell Ranger
  - Default path: `references/refdata-gex-GRCh38-2024-A`
  - Download from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads)

## 🚀 Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/dartwhit/W_scRNAseq_pipeline.git
   cd W_scRNAseq_pipeline
   ```

2. **Create and activate a Conda environment**:
   ```bash
   conda env create -f environment.yml
   conda activate scrna_pipeline
   ```

3. **Install CAMML (not available in Conda channels)**:
   ```R
   devtools::install_github("Chenmengpin/CAMML")
   ```

4. **Install Cell Ranger**:
   - Download from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads)
   - Follow installation instructions for your system
   - Ensure Cell Ranger is available in your PATH or load as a module

5. **Download Reference Genome**:
   ```bash
   cd references
   wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
   tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
   ```

## 🏃 Quick Start

```bash
# 1. Configure your dataset in config.yaml (or generate it)
python make_config.py GSE249279

# 2. Run the pipeline
DATASET=GSE249279 snakemake --cores 8

# 3. Check results
ls -l datasets/GSE249279/seurat/
```

## ⚙️ Configuration

### Using config.yaml

The `config.yaml` file defines datasets and their associated samples. Each dataset contains a mapping of sample names to SRA accession IDs.

**Example config.yaml**:
```yaml
datasets:
  GSE249279:
    samples:
      GSM7932779: SRR27061746
      GSM7932780: SRR27061745
      # ... more samples
```

### Generating Configuration Automatically

Use the `make_config.py` script to automatically generate a configuration file from a GEO dataset:

```bash
python make_config.py GSE249279
```

This script:
- Fetches metadata from GEO using pysradb
- Extracts sample-to-SRA mappings
- Creates a properly formatted `config.yaml` file

**Script usage**:
```bash
python make_config.py <GEO_ACCESSION>
```

## 💻 Usage

Before running commands below, activate the Conda environment:

```bash
conda activate scrna_pipeline
```

### Basic Usage

Run the pipeline for a specific dataset:

```bash
DATASET=GSE249279 snakemake --cores 8
```

### Advanced Options

**Dry run** (see what will be executed):
```bash
DATASET=GSE249279 snakemake -n
```

**Run specific rules**:
```bash
DATASET=GSE249279 snakemake --cores 8 datasets/GSE249279/seurat/GSE249279_annotated.rds
```

**Use more cores**:
```bash
DATASET=GSE249279 snakemake --cores 16
```

**Resume after failure**:
```bash
DATASET=GSE249279 snakemake --cores 8 --rerun-incomplete
```

### Pipeline Rules

The pipeline consists of three main rules:

1. **download_and_rename_fastq**: Downloads FASTQ files from SRA
2. **cellranger**: Runs Cell Ranger count for alignment and quantification
3. **process_seurat**: Performs QC, doublet detection, clustering, and annotation

## 🔄 Pipeline Workflow

```
┌─────────────────────────────────────────────┐
│  Input: SRA Accession IDs (config.yaml)    │
└──────────────────┬──────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────┐
│  1. Download FASTQ files (fastq-dl)        │
│     - Downloads paired-end reads            │
│     - Renames to Cell Ranger format        │
└──────────────────┬──────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────┐
│  2. Cell Ranger Processing                 │
│     - Alignment to reference genome         │
│     - Gene quantification                   │
│     - Output: filtered feature matrices    │
└──────────────────┬──────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────┐
│  3. Seurat Processing                      │
│     - Load and merge samples               │
│     - QC filtering (features, MT%)         │
│     - Doublet detection (scDblFinder)      │
│     - Normalization and scaling            │
│     - PCA and UMAP                         │
│     - Cell type annotation:                │
│       * SingleR (HumanPrimaryCellAtlas)    │
│       * CAMML (skin immune cells)          │
│     - Generate visualization plots         │
└──────────────────┬──────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────┐
│  Output: Annotated Seurat object + Plots   │
└─────────────────────────────────────────────┘
```

### Quality Control Filters

- **Minimum features**: > 200 genes per cell
- **Maximum features**: < 2500 genes per cell
- **Mitochondrial content**: < 25%
- **Minimum cells per feature**: 3 cells

### Cell Type Annotation

The pipeline uses two complementary annotation methods:

1. **SingleR**: Uses HumanPrimaryCellAtlasData reference
   - Annotates major cell types (B cells, T cells, macrophages, etc.)
   
2. **CAMML**: Uses skin immune cell gene sets
   - Specialized annotation for skin-related cell types

## 📊 Output

The pipeline generates the following outputs for each dataset:

### Directory Structure

```
datasets/
└── GSE249279/
    ├── fastq/
    │   ├── SRR27061746_S1_L001_R1_001.fastq.gz
    │   ├── SRR27061746_S1_L001_R2_001.fastq.gz
    │   └── ... (more samples)
    ├── cellranger_output/
    │   └── SRR27061746/
    │       └── outs/
    │           └── filtered_feature_bc_matrix/
    │               ├── barcodes.tsv.gz
    │               ├── features.tsv.gz
    │               └── matrix.mtx.gz
    └── seurat/
        ├── GSE249279_annotated.rds          # Final Seurat object
        ├── GSE249279_plots.pdf              # Comprehensive plots
        └── GSE249279_doublet_report.csv     # Doublet detection report
```

### Output Files

1. **`{dataset}_annotated.rds`**: 
   - Annotated Seurat object with all metadata
   - Contains normalized data, UMAP coordinates, clusters, and cell type labels
   
2. **`{dataset}_plots.pdf`**:
   - UMAP plots (by clusters, samples, annotations)
   - Feature plots for cell type markers
   - Violin plots for marker expression
   - Heatmap comparing annotation methods
   
3. **`{dataset}_doublet_report.csv`**:
   - Summary of doublet detection results
   - Counts of singlets vs doublets

### Metadata in Seurat Object

The final Seurat object contains:
- `orig.ident`: Original sample ID
- `seurat_clusters`: Cluster assignments
- `SingleR_label`: SingleR cell type annotation
- `Base_CAMML`: CAMML cell type annotation
- `scDblFinder.class`: Doublet classification
- `percent.mt`: Mitochondrial percentage

## 🔧 Troubleshooting

### Common Issues

**Issue**: Module not found errors for Cell Ranger
```bash
# Solution: Load Cell Ranger module or ensure it's in PATH
module load cellranger/9.0.0
# OR
export PATH=/path/to/cellranger:$PATH
```

**Issue**: fastq-dl fails to download
```bash
# Solution: Check SRA accession ID and internet connection
# Verify the accession exists:
fastq-dl -a SRR27061746 --verbose
```

**Issue**: Out of memory during Seurat processing
```bash
# Solution: Increase memory or adjust filtering parameters
# In process_seurat.R, modify the filtering criteria
# or process fewer samples at once
```

**Issue**: Reference genome not found
```bash
# Solution: Update the reference path in Snakefile
# Line 62: params: reference="references/refdata-gex-GRCh38-2024-A"
```

**Issue**: R packages not installed
```R
# Solution: Install missing packages
BiocManager::install(c("package_name"))
```

### Getting Help

If you encounter issues:
1. Check the Snakemake log files
2. Review Cell Ranger output in `{dataset}/cellranger_output/`
3. Examine R script errors in the console output
4. Open an issue on GitHub with:
   - Error message
   - Dataset being processed
   - System information

### Check Previous Run Status (Per Sample)

Use the reusable status script to see which stages completed and which failed/missed outputs.

Generate a report for a dataset:

```bash
python scripts/report_run_status.py --dataset GSE264508
```

This writes a TSV report to:

```text
datasets/GSE264508/status/GSE264508_run_status.tsv
```

You can also override paths:

```bash
python scripts/report_run_status.py \
   --dataset GSE264508 \
   --config config.yaml \
   --base-dir . \
   --snakemake-log-dir .snakemake/log \
   --output datasets/GSE264508/status/custom_status.tsv
```

Report columns:
- `download_fastq_status`: `ok`, `partial`, or `missing` based on expected lane FASTQs from `config.yaml`
- `cellranger_status`: `ok` when `outs/filtered_feature_bc_matrix` exists for a sample
- `process_seurat_status`: dataset-level stage status (`ok`, `partial`, `missing`) based on final Seurat outputs
- `failure_hints`: best-effort hints parsed from recent Snakemake logs in `.snakemake/log`

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is available for academic and research purposes.

## 📧 Contact

For questions or support, please open an issue on the [GitHub repository](https://github.com/dartwhit/W_scRNAseq_pipeline).

## 🙏 Acknowledgments

This pipeline integrates several excellent tools and packages:
- [Snakemake](https://snakemake.readthedocs.io/) - Workflow management
- [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger) - 10x Genomics scRNA-seq processing
- [Seurat](https://satijalab.org/seurat/) - Single-cell analysis
- [SingleR](https://bioconductor.org/packages/SingleR/) - Cell type annotation
- [scDblFinder](https://bioconductor.org/packages/scDblFinder/) - Doublet detection
- [CAMML](https://github.com/Chenmengpin/CAMML) - Cell type annotation
- [fastq-dl](https://github.com/rpetit3/fastq-dl) - SRA data download
