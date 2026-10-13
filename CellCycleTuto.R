# scRNAseq analysis - Cell cycle scoring tutorial
#
# Dataset: brain_E18_Filtre.rds (must be in the same folder as this script).
# Cells from a combined cortex, hippocampus and sub ventricular zone of an E18
# mouse:
# https://www.10xgenomics.com/datasets/10-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0 # nolint
#
# This script starts from the filtered dataset produced by the preprocessing
# (quality control, doublet removal, normalization). It computes the PCA, the
# clusters and the UMAP, then scores the cell cycle of each cell.
#
# Run step by step: each '# %%' line starts a new cell (Ctrl+Enter / Cmd+Enter
# runs the current line, select a cell to run it). Run everything: Rscript
# analyse.R

# %% Set up ----------------------------------------------------------------
#
# fix the random seed for reproducibility of the results
the_seed <- 42
set.seed(the_seed)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

sample_name <- "brain_E18"

# filtered dataset (output of the preprocessing): Seurat object
sobj_filtrd <- readRDS(paste0(sample_name, "_Filtre.rds"))

# %% PCA, clusters and UMAP ------------------------------------------------
#
# The filtered dataset contains the normalized data, but no dimension reduction
# and no clusters yet. We compute them on the log-normalized `RNA` assay, with
# the parameters of the clustering course (2000 variable genes, 30 PCs,
# resolution 0.2).
DefaultAssay(sobj_filtrd) <- "RNA"

# number of PC kept for the analysis
n_pc <- 30

sobj_filtrd <- FindVariableFeatures(sobj_filtrd,
  selection.method = "vst",
  nfeatures = 2000, verbose = FALSE
)
sobj_filtrd <- ScaleData(sobj_filtrd, verbose = FALSE)
sobj_filtrd <- RunPCA(sobj_filtrd, verbose = FALSE)
sobj_filtrd <- FindNeighbors(sobj_filtrd, dims = 1:n_pc, verbose = FALSE)
sobj_filtrd <- FindClusters(sobj_filtrd, resolution = 0.2, verbose = FALSE)
sobj_filtrd <- RunUMAP(sobj_filtrd,
  dims = 1:n_pc, n.neighbors = 20,
  min.dist = 0.8, verbose = FALSE
)

# number of cells per cluster
table(Idents(sobj_filtrd))

DimPlot(sobj_filtrd, reduction = "umap", label = TRUE)

# %% Cell cycle scoring ----------------------------------------------------
#
# The cell cycle can have a strong effect on the transcriptome and the cluster
# formation. We check that our clustering was not too influenced by it using
# the Seurat method `CellCycleScoring` (another option is the Scran package
# Cyclone).
#
# `CellCycleScoring` gives each cell a score for the S phase and a score for the
# G2/M phase, based on the expression of known marker genes (Tirosh et al.,
# 2016). Each score is computed like `AddModuleScore`: the mean expression of
# the marker genes minus the mean expression of random control genes with
# similar expression levels. The cell is then assigned to a phase:
#
# - **G2M** if its G2M score is the highest and positive
# - **S** if its S score is the highest and positive
# - **G1** if both scores are negative (the cell does not express S or G2/M
#   genes: G1 or non-cycling G0 cells)

# %% Cell cycle genes ------------------------------------------------------
#
# Seurat provides the list of marker genes in `cc.genes.updated.2019`. These are
# **human** gene symbols (e.g. `MCM5`), while our data is from **mouse** (e.g.
# `Mcm5`). For these genes, the mouse symbol is simply the human symbol with
# only the first letter in upper case.
#
# human gene symbols -> mouse gene symbols (MCM5 -> Mcm5)
to_mouse <- function(genes) {
  paste0(substr(genes, 1, 1), tolower(substring(genes, 2)))
}

s_genes <- to_mouse(cc.genes.updated.2019$s.genes)
g2m_genes <- to_mouse(cc.genes.updated.2019$g2m.genes)

# check that the genes are found in the count matrix
length(intersect(s_genes, rownames(sobj_filtrd))) # out of 43
length(intersect(g2m_genes, rownames(sobj_filtrd))) # out of 54

# For other genes, the simple case conversion may fail. A more robust option is
# to convert the genes with orthology databases, e.g.
# `gprofiler2::gorth(cc.genes.updated.2019$s.genes, source_organism =
# "hsapiens", target_organism = "mmusculus")$ortholog_name` (requires an
# internet connection).

# %% Score the cells -------------------------------------------------------
#
# The scores are computed on the log-normalized data of the `RNA` assay. They
# are added to the metadata in 3 new columns: `S.Score`, `G2M.Score` and
# `Phase`.
sobj_filtrd <- CellCycleScoring(sobj_filtrd,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = FALSE
)

head(sobj_filtrd@meta.data[, c("S.Score", "G2M.Score", "Phase")])

# number of cells in each phase
table(sobj_filtrd$Phase)

# %% Visualize the phases --------------------------------------------------
DimPlot(sobj_filtrd, reduction = "umap", group.by = "Phase")

FeaturePlot(sobj_filtrd,
  reduction = "umap",
  features = c("S.Score", "G2M.Score"), order = TRUE
)

# proportion of cells in each phase for each cluster
round(prop.table(table(Idents(sobj_filtrd), sobj_filtrd$Phase), margin = 1), 2)

phase_by_cluster <- data.frame(
  cluster = Idents(sobj_filtrd),
  Phase = sobj_filtrd$Phase
)
ggplot(phase_by_cluster, aes(x = cluster, fill = Phase)) +
  geom_bar(position = "fill") +
  labs(x = "Cluster", y = "Proportion of cells") +
  theme_minimal()

# The scores should match the expression of classical proliferation markers:
FeaturePlot(sobj_filtrd, reduction = "umap", features = c("Mki67", "Top2a"))

# **Question**: Are the cycling cells (S and G2M) spread across all the
# clusters, or grouped in some clusters? Do they match the cells expressing
# `Mki67` and `Top2a`?

# %% Does the cell cycle drive the variability ? ---------------------------
#
# To check how strongly the cell cycle structures the data, we run a PCA using
# only the cell cycle genes. If the cells separate by phase on this PCA, the
# cell cycle is a strong source of variation in the data.
#
# By default, `ScaleData` only scaled the highly variable genes, and most cell
# cycle genes are not among them. We first scale the variable genes **and** the
# cell cycle genes (the values of the variable genes do not change).
sobj_filtrd <- ScaleData(sobj_filtrd,
  features = union(
    VariableFeatures(sobj_filtrd),
    c(s_genes, g2m_genes)
  )
)

sobj_filtrd <- RunPCA(sobj_filtrd,
  features = c(s_genes, g2m_genes),
  reduction.name = "pca_cc", reduction.key = "PCcc_",
  verbose = FALSE
)

DimPlot(sobj_filtrd, reduction = "pca_cc", group.by = "Phase")

# %% Should we remove the cell cycle effect ? ------------------------------
#
# If the cell cycle drives the clustering, you can remove its effect by
# regressing out the scores while scaling the data (`vars.to.regress`). There
# are two options:
#
# - regress `S.Score` and `G2M.Score`: removes **all** the cell cycle signal.
#   Cycling and non-cycling cells are then mixed.
# - regress the difference `S.Score - G2M.Score`: removes only the differences
#   **between** the cycling phases (S vs G2M), but keeps the separation between
#   cycling and non-cycling cells.
#
# At E18, the brain is in active neurogenesis: dividing progenitors (radial
# glia, intermediate progenitors) are a biological population we want to keep.
# Removing all the cell cycle signal would merge them with the other cells. The
# second option is therefore the most appropriate here.
#
# To keep the clusters unchanged for the rest of the analysis, we test the
# regression on a copy of the object (this step takes a few minutes):
sobj_cc <- sobj_filtrd
sobj_cc$CC.Difference <- sobj_cc$S.Score - sobj_cc$G2M.Score

sobj_cc <- ScaleData(sobj_cc,
  vars.to.regress = "CC.Difference",
  features = VariableFeatures(sobj_cc)
)
sobj_cc <- RunPCA(sobj_cc,
  features = VariableFeatures(sobj_cc),
  verbose = FALSE
)
sobj_cc <- RunUMAP(sobj_cc,
  dims = 1:n_pc, n.neighbors = 20, min.dist = 0.8,
  verbose = FALSE
)

# compare before / after regression
DimPlot(sobj_filtrd, reduction = "umap", group.by = "Phase") +
  ggtitle("Before regression")
DimPlot(sobj_cc, reduction = "umap", group.by = "Phase") +
  ggtitle("After regression of S - G2M")

# **Question**: After regression, are the S and G2M cells still separated on the
# UMAP ? Are the cycling cells still distinct from the non-cycling ones ? Would
# you keep the regression for the rest of the analysis ?
#
# Note: if you normalize the data with `SCTransform`, the regression is done
# directly in this function: `SCTransform(sobj, vars.to.regress =
# "CC.Difference")`.

# %% Save ------------------------------------------------------------------
#
# The scores and phases (`S.Score`, `G2M.Score`, `Phase`) are kept in the
# metadata of `sobj_filtrd` for the rest of the analysis.
save(sobj_filtrd, file = paste0(sample_name, "_CELLCYCLE.rda"))
