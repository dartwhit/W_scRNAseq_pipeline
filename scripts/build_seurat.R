library(Seurat)
library(scDblFinder)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
output_file <- args[2]

#output_file <- "../datasets/GSE209635/seurat/GSE209635_seurat.rds"
#dataset <- "GSE209635"

# Define output directory
output_dir <- dirname(output_file)
dir.create(output_dir, showWarnings = FALSE)
report_file <- file.path(output_dir, paste0(dataset, "_doublet_report.csv"))
umap_pdf <- file.path(output_dir, paste0(dataset, "_umap_plots.pdf"))

# Path to all filtered feature matrices
data_dir <- paste0("datasets/", dataset, "/cellranger_output/")
samples <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)

# Load each sample as a Seurat object and assign `orig.ident`
seurat_list <- lapply(samples, function(sample) {
  sample_name <- basename(sample)  # Extracts sample name (e.g., SRR20614738)
  message("Processing: ", sample)

  seurat_obj <- Read10X(data.dir = file.path(sample, "outs/filtered_feature_bc_matrix"))
  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = dataset, min.cells = 3, min.features = 200)

  # Assign sample identity before merging
  seurat_obj$orig.ident <- sample_name

  # Quality control filtering
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

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

# **Doublet Detection and Removal using scDblFinder**
message("Running scDblFinder to detect doublets...")

sce_obj <- scDblFinder(sce_obj)

# **Convert back to Seurat**
seurat_combined <- as.Seurat(sce_obj)

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


# Save as RDS file
saveRDS(seurat_combined, file = output_file)

message("Seurat object saved: ", output_file)

# **Save UMAP plots to PDF**
pdf(umap_pdf, width = 8, height = 6)

# UMAP colored by clusters
p1 <- DimPlot(seurat_combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
      ggtitle("UMAP by Seurat Clusters")
print(p1)

# UMAP colored by samples
p2 <- DimPlot(seurat_combined, reduction = "umap", group.by = "orig.ident", label = FALSE) +
      ggtitle("UMAP by Sample (orig.ident)")
print(p2)

dev.off()
message("UMAP plots saved: ", umap_pdf)

