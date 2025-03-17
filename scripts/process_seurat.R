# Load required libraries
library(Seurat)
library(yaml)
# library(future)
library(SingleR)
library(celldex)
library(CAMML)
library(scDblFinder)
library(SingleCellExperiment)
library(pheatmap)
# library(BiocParallel)

library(ggplot2)

# # **Make BiocParallel use future**
# bp_param <- MulticoreParam(workers = 8)
# register(bp_param)

# options(future.globals.maxSize = 10 * 1024^3)  # Increase to 10 GB


# Set up parallel backend (use all but one core)
#plan("multisession", workers = 8)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
config_file <- "config.yaml"  # Path to the config file
# Load config.yaml
config <- yaml.load_file(config_file)
# Get the list of samples from config.yaml
sample_ids <- config$datasets[[dataset]]$samples
if (is.null(sample_ids)) {
  stop("Error: No samples found in config.yaml for dataset ", dataset)
}

output_file <- paste0("datasets/", dataset, "/seurat/", dataset, "_annotated.rds")
output_dir <- dirname(output_file)
dir.create(output_dir, showWarnings = FALSE)
report_file <- file.path(output_dir, paste0(dataset, "_doublet_report.csv"))

plots_pdf <- paste0("datasets/", dataset, "/seurat/", dataset, "_plots.pdf")

# Path to Cell Ranger outputs (for specified samples only)
data_dir <- paste0("datasets/", dataset, "/cellranger_output/")
sample_paths <- file.path(data_dir, sample_ids, "outs/filtered_feature_bc_matrix")

# Check if all samples exist
missing_samples <- sample_paths[!file.exists(sample_paths)]
if (length(missing_samples) > 0) {
  stop("Error: The following sample directories do not exist:\n", paste(missing_samples, collapse = "\n"))
}


seurat_list <- lapply(sample_ids, function(sample_name) {
    message("Processing: ", sample_name)
	sample_path <- file.path(data_dir, sample_name, "outs/filtered_feature_bc_matrix")
	
	seurat_obj <- Read10X(data.dir = sample_path)
    seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = dataset, min.cells = 3, min.features = 200)
	
	# Assign sample identity
    seurat_obj$orig.ident <- sample_name
	# QC filtering
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
	
	  return(seurat_obj)
})


# Merge all samples into a single Seurat object
if (length(seurat_list) > 1) {
  seurat_combined <- merge(seurat_list[[1]], y = seurat_list[-1])
} else {
  seurat_combined <- seurat_list[[1]]
}

seurat_combined <- JoinLayers(seurat_combined)

# **Retain Sample Identity**
message("Checking sample identity retention:")
table(seurat_combined$orig.ident)

# Normalize & Identify Variable Features
seurat_combined <- NormalizeData(seurat_combined)
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)

# Scale Data & PCA
seurat_combined <- ScaleData(seurat_combined)
seurat_combined <- RunPCA(seurat_combined)

# **Convert Seurat object to SingleCellExperiment for scDblFinder**
sce_obj <- as.SingleCellExperiment(seurat_combined)
# Ensure `orig.ident` is a factor (important for grouping samples)
sce_obj$orig.ident <- as.factor(sce_obj$orig.ident)

# **Check for empty samples before running scDblFinder**
num_cells_per_sample <- table(sce_obj$orig.ident)
empty_samples <- names(num_cells_per_sample[num_cells_per_sample == 0])

if (length(empty_samples) > 0) {
    message("Warning: The following samples have zero cells and will be skipped: ", paste(empty_samples, collapse=", "))
    sce_obj <- sce_obj[, !sce_obj$orig.ident %in% empty_samples]  # Remove empty samples
}

# If no valid cells remain, stop execution
if (ncol(sce_obj) == 0) {
    stop("Error: No valid cells remaining after filtering. Check input data.")
}

# **Doublet Detection and Removal using scDblFinder**
message("Running scDblFinder to detect doublets...")
# Reduce the number of PCA dimensions if dataset is large
num_cells <- ncol(sce_obj)
print(num_cells)
dims_to_use <- if (num_cells > 50000) 20 else 50  # Adjust dimensions for large datasets
artificialDoublets <- ifelse(num_cells > 50000, FALSE, TRUE)

sce_obj <- scDblFinder(sce_obj, samples = "orig.ident", dims = dims_to_use)


# sce_obj <- scDblFinder(sce_obj,
                        # samples = "orig.ident",
                        # dims = dims_to_use,
                        # artificialDoublets = artificialDoublets,  # Disable artificial doublets for large datasets
						# svd.solver = "auto",
                        # BPPARAM = bp_param)

# **Convert back to Seurat**
seurat_combined <- as.Seurat(sce_obj)
seurat_combined$doublet_class <- seurat_combined@meta.data$scDblFinder.class
# Count doublets and singlets
doublet_counts <- table(seurat_combined$scDblFinder.class)
doublet_summary <- as.data.frame(doublet_counts)
colnames(doublet_summary) <- c("Cell_Type", "Count")

# Save doublet detection report
write.csv(doublet_summary, report_file, row.names = FALSE)
message("Doublet report saved: ", report_file)

# Remove doublets, keeping only singlets
seurat_combined <- subset(seurat_combined, subset = scDblFinder.class == "singlet")

# Clustering & UMAP
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:10,reduction="PCA")
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)
seurat_combined <- RunUMAP(seurat_combined, dims = 1:10,reduction="PCA")

# Load reference for SingleR
message("Loading HumanPrimaryCellAtlasData reference...")
ref <- celldex::HumanPrimaryCellAtlasData()
coldata_df <- as.data.frame(colData(ref))

# Keep only relevant cell types
ref_keep_idx <- coldata_df$label.main %in% c("B_cell", "DC", "Endothelial_cells",
                                             "Epithelial_cells", "Fibroblasts",
                                             "Keratinocytes", "Macrophage",
                                             "Monocyte", "Neutrophils", "NK_cell",
                                             "Smooth_muscle_cells", "T_cells")
ref <- ref[, ref_keep_idx]

# Run SingleR
message("Running SingleR annotation...")
log_counts <- GetAssayData(seurat_combined, layer = "data")
predictions <- SingleR(test = log_counts, ref = ref, labels = ref$label.main)
seurat_combined$SingleR_label <- predictions$labels

# Run CAMML Annotation
message("Running CAMML annotation...")
complete_gene_set <- GetGeneSets("skin.immune.cells")
seurat_combined <- CAMML(seurat_combined, complete_gene_set)
results <- GetCAMMLLabels(seurat_combined, labels = "top1")
seurat_combined$Base_CAMML <- unlist(results)

# **Define Marker Genes for Feature Plots**
markers_list <- list(
  "Keratinocytes" = c("KRT2", "KRT14", "KRT17", "LOR"),
  "Fibroblasts" = c("COL1A1", "FN1", "SFRP2", "DCN"),
  "Endothelial" = c("PECAM1", "VWF", "ACKR1", "CDH5"),
  "Smooth Muscle" = c("MYH11", "DES", "CALD1", "TAGLN"),
  "T Cells" = c("CD3E", "CXCR4", "TOX", "PTPRC"),
  "Myeloid Cells" = c("MRC1", "LYZ", "CXCL8", "CCR1"),
  "B Cells" = c("IGHG1", "JCHAIN", "BANK1", "IGHD"),
  "Melanocytes" = c("PMEL", "TYRP1", "DCT", "MLANA"),
  "Eccrine Cells" = c("MUC1", "DCD", "SCGB2A2", "PIP"),
  "Lymphatic Endothelium" = c("CCL21", "TFF3", "MMRN1", "ECSCR"),
  "Mast Cells" = c("TPSAB1", "KIT", "HPGD", "SLCO2B1"),
  "Neuronal Cells" = c("MPZ", "PMP22", "NRXN1", "S100B")
)

# Generate All Plots in One PDF
pdf(plots_pdf, width = 8, height = 6)

# UMAP - Seurat clusters
p1 <- DimPlot(seurat_combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
      ggtitle("UMAP by Seurat Clusters")
print(p1)

# UMAP - Sample origin
p2 <- DimPlot(seurat_combined, reduction = "umap", group.by = "orig.ident", label = FALSE) +
      ggtitle("UMAP by Sample (orig.ident)")
print(p2)

# UMAP - SingleR annotation
p3 <- DimPlot(seurat_combined, reduction = "umap", group.by = "SingleR_label") +
      ggtitle("UMAP by SingleR Annotation")
print(p3)

# UMAP - CAMML annotation
p4 <- DimPlot(seurat_combined, reduction = "umap", group.by = "Base_CAMML") +
      ggtitle("UMAP by CAMML Annotation")
print(p4)

# Heatmap Comparing SingleR and CAMML
message("Generating heatmap of annotation agreement...")
cont_calls <- table(seurat_combined$Base_CAMML, seurat_combined$SingleR_label)
pheatmap::pheatmap(log10(cont_calls + 10), color = colorRampPalette(c('white','blue'))(10),
                   main = "Comparison of SingleR and CAMML Annotations")


# Feature and Violin Plots for Cell Markers
for (cell_type in names(markers_list)) {
    markers <- markers_list[[cell_type]]
    message("Generating FeaturePlot and VlnPlot for: ", cell_type)

    feature_plot <- FeaturePlot(seurat_combined, features = markers, reduction = "umap") +
                    ggtitle(paste("FeaturePlot -", cell_type))
    print(feature_plot)

    vln_plot <- VlnPlot(seurat_combined, features = markers, ncol = 2) +
                ggtitle(paste("Violin Plot -", cell_type))
    print(vln_plot)
}

dev.off()
message("All plots saved in: ", plots_pdf)

# Save Annotated Seurat Object
saveRDS(seurat_combined, file = output_file)
message("Final Seurat object saved: ", output_file)
