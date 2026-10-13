# Introduction to scRNAseq analysis - Clustering
# Author: Lorette Noiret
#
# Dataset: brain_E18_FILTERED.rda (must be in the same folder as this script)
# Run step by step: each '# %%' line starts a new cell (Ctrl+Enter / Cmd+Enter
# runs the current line, select a cell to run it). Run everything: Rscript
# analyse.R

# %% Data ------------------------------------------------------------------
#
# Cells from a combined cortex, hippocampus and sub ventricular zone of an E18
# mouse.
#
# https://www.10xgenomics.com/datasets/10-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0 # nolint

# %% Identify cell populations - Introduction ------------------------------
#
# In this section, we will analyse further the brain data and try to identify
# sub-populations of cells. The objective is to create clusters of cells, i.e to
# make groups of cells that share similar expression profile.
#
# The mains steps are as follow :
# - Select a subset of genes to perfom the downstream analyses: highly variable
#   genes
# - Scale the normalized data
# - Perform a dimension reduction (PCA)
# - Cluster the cells
# - Annotate the clusters
#
# Important : We will apply this pipeline to log-normalized data. If we used
# `SCTransform`, the steps selection of highly variables genes and scaling are
# not necessary. The `SCTransform` command replaces NormalizeData(),
# ScaleData(), and FindVariableFeatures().

# %% Set up ----------------------------------------------------------------
#
# fix the random seed for reproducibility of the results
the_seed <- 42
set.seed(the_seed)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(gprofiler2)
  library(scales)
})

sample_name <- "brain_E18"
species <- "mouse"

load(paste0(sample_name, "_FILTERED.rda"))

# %% Set the default assay -------------------------------------------------
DefaultAssay(sobj_filtrd) <- "RNA"

# %% Select highly variables (hvg) -----------------------------------------
#
# The downstream analyses (dimension reduction and then clustering) will be
# performed on the subset of genes. The aims is to exclude genes with minor
# biological relevance and decrease the computational load.
#
# Find Variable Genes (usually between 500 to 3000):
sobj_filtrd <- FindVariableFeatures(sobj_filtrd,
  selection.method = "vst",
  nfeatures = 2000, verbose = FALSE
)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sobj_filtrd), 10)

# you can also access the complete list using VariableFeatures()
list_hvg <- VariableFeatures(sobj_filtrd)

# %% Dimension Reduction with PCA ------------------------------------------
#
# Before performing PCA, we scale the data.

# %% Scaling ---------------------------------------------------------------
#
# Expression of each gene is scaled (centered around its mean with a unit
# standard deviation). This step gives equal weight to all genes while
# performing the in the dimension reduction, and insures that the most expressed
# genes do not dominate.
sobj_filtrd <- ScaleData(sobj_filtrd)

# By default, scaling is only performed on the variable genes.
#
# Since Seurat version 3, you can perform the normalization, selection of highly
# variable genes and scaling in a single step with the function `sctransform`.
# This approach use regularized negative binomial regression to normalize UMI
# count data.

# %% PCA -------------------------------------------------------------------
#
# Even though, the gene selection step reduces the data considerably, there is a
# further need to reduce dimension to keep the most important variability in the
# data.There are several methods for dimension reduction, the most commonly used
# is PCA, but you can also check ICA or more single-cell specific methods that
# tries to cope with over-dispersion of the data such as ZIFA and pCMF.
sobj_filtrd <- RunPCA(sobj_filtrd,
  features = VariableFeatures(object = sobj_filtrd)
)

# Exercice: Seurat object manipulation. Open your Seurat object again and notice
# which fields have changed. Perform the following tasks without using IA help !
# - Extract the normalize data (counting matrix)
# - Extract the scaled data (counting matrix)
# - Extract the highly variable gene
# - Extract the first 2 PC components

# %% Explore the PCA results -----------------------------------------------
#
# It is worth analyzing the genes that contribute the most to each axes. It may
# help you to identify genes whose expression decrease/increase together. The
# `DimHeatmap` function in Seurat helps you to visualize the genes that are
# driving the components and allows to get some insight about the heterogeneity
# of the data.
#
# Examine and visualize PCA results :
print(sobj_filtrd[["pca"]], dims = 1:6, nfeatures = 5)
VizDimLoadings(sobj_filtrd, dims = 3, reduction = "pca")

# Plot expression of gene driving principal component
DimHeatmap(sobj_filtrd, dims = 1:12)

# %% Choose the number of axes ---------------------------------------------
#
# After PCA, you need to decide how many components you want to keep for the
# clustering. Using too few PCs might merge distinct populations, while too many
# PCs could introduce noise and spurious clusters. In practice, keep between
# 10-50 (up to 100) principal components.
#
# An elbow plot show you the inertia associated with each axe (how much of the
# data variation is captured by each PC axe). It helps determine the "elbow
# point," where the percentage of explained variance by each successive PC
# begins to level off. PCs before this point contribute significantly to the
# variability in the data, while PCs after this point likely represent noise.
DimPlot(sobj_filtrd, reduction = "pca")

ElbowPlot(sobj_filtrd, ndims = 40)

# There is "elbow" is at 2, 7,14 components. The inertia stabilizes around 30
# components. So what do we choose ?
#
# I choose to go further with n_pc = 14, but I could also explore the results
# with other values (10 ,20 or ben 30 !) to see the impact on the downstream
# results.

# %% Cluster cells ---------------------------------------------------------
#
# There are two main steps to cluster the cells:
#
# - Make a **Shared Nearest Neighbor graph**:
#   - identify the k-nearest neighbours of each cell (by default Seurat takes
#     k=20 neighbours).
#   - The distance between cells is calculated using the coordinates obtained
#     with the PCA.
# - Make clusters (community of cells) using **Leiden or Louvain algorithm**.
#   - The number of clusters depends on the **resolution parameter** chosen: the
#     higher the value of the parameter is, the more groups you will get.
#
# number of PC kept for the analysis
n_pc <- 30

# s-nn graph
sobj_filtrd <- FindNeighbors(sobj_filtrd, dims = 1:n_pc)

# make the clusters (here we try several resolutions)
sobj_filtrd <- FindClusters(sobj_filtrd, resolution = c(0.1, 0.2, 0.5, 1, 2))

# the clusters are stored in `metadata`:
head(sobj_filtrd@meta.data)

# The various clusterings are stored in the last columns of the layer
# `metadata`.
#
# check the number of cells per cluster
table(sobj_filtrd[["RNA_snn_res.0.1"]])
table(sobj_filtrd[["RNA_snn_res.0.2"]])
table(sobj_filtrd[["RNA_snn_res.0.5"]])

# Usually the clustering results are presented using UMAP (dimension reduction)
# which visually separates the group more nicely than PCA.
#
# create umap
sobj_filtrd <- RunUMAP(sobj_filtrd, dims = 1:15)

# Visualize clusters on UMAP
#
# specify the clustering you want to see
Idents(sobj_filtrd) <- "RNA_snn_res.0.2"

DimPlot(sobj_filtrd, reduction = "umap")

# **Important** : can you interpret the distance between cluster on a UMAP ?
#
# change the UMAP appearance by playing with the parameters
sobj_filtrd <- RunUMAP(sobj_filtrd,
  dims = 1:n_pc, n.neighbors = 20,
  min.dist = 0.8
)

DimPlot(sobj_filtrd, reduction = "umap")

# %% Group exercice --------------------------------------------------------
#
# Try to play with the number of PCs, and the parameters affecting the
# clustering. How could you decide which set of parameters is the best ? Comment
# on the PC results. Could you extract any useful information ?

# %% Check for biais -------------------------------------------------------
#
# Check if you can detect a biais in the clustering.
#
# Metrics
FeaturePlot(
  object = sobj_filtrd, reduction = "umap", pt.size = .2,
  features = "percent_mt", cols = c("gold", "blue"), order = TRUE
)

FeaturePlot(
  object = sobj_filtrd, reduction = "umap", pt.size = .2,
  features = "nCount_RNA", cols = c("gold", "blue"), order = TRUE
)

FeaturePlot(
  object = sobj_filtrd, reduction = "umap", pt.size = .2,
  features = "nFeature_RNA", cols = c("gold", "blue"),
  order = TRUE
)

# The cell cycle can have a strong effect on the transcriptome and the cluster
# formation. We could also check that our clustering was not too influence by it
# using Seurat method `CellCycleScoring` or the Scran package Cyclone.

# %% Cell cycle scoring ----------------------------------------------------
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
  g2m.features = g2m_genes
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

ggplot(
  sobj_filtrd@meta.data,
  aes(x = .data[["RNA_snn_res.0.2"]], fill = Phase)
) +
  geom_bar(position = "fill") +
  labs(x = "Cluster", y = "Proportion of cells") +
  theme_minimal()

# **Question**: Are the cycling cells (S and G2M) spread across all the
# clusters, or grouped in some clusters? Compare with the expression of the
# `Cycling` markers (`Mki67`, `Top2a`) in the annotation part below.

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
# To keep the clusters of the next sections unchanged, we test the regression on
# a copy of the object (this step takes a few minutes):
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

# %% Annotate your clusters ------------------------------------------------
#
# Hurray, you succeeded in clustering your cells!
#
# Now your biological work starts. You need to annotate each cluster based on
# marker gene expression, and check that your annotation makes biological sense.
#
# At E18, the brain tissue is in the middle of active neurogenesis. This has an
# important consequence for how you should interpret your results: instead of a
# set of fully differentiated, discrete cell types, you are mostly looking at a
# **differentiation continuum** — radial glia (stem-like progenitors)
# progressively give rise to intermediate progenitors, then to immature,
# migrating neurons, which only later commit to an excitatory or inhibitory
# identity. Don't expect these stages to form clearly separated islands on the
# UMAP: a gradient-like structure is the expected biological result here, not a
# sign that something went wrong. In contrast, non-neuronal populations
# (vascular cells, microglia, blood contaminants) are unrelated lineages and
# should appear as distinct, well-isolated clusters.
#
# We provide a curated list of marker genes covering both the neurogenic lineage
# and the non-neuronal populations you can expect to find in this tissue:
marker_genes <- list(
  Radial_glia      = c("Pax6", "Sox2", "Hes5", "Fabp7", "Nes"),
  Interm_prog      = c("Eomes", "Neurog2", "Insm1"),
  Cycling          = c("Mki67", "Top2a", "Pcna"),
  Immature_neuron  = c("Dcx", "Stmn2", "Neurod2", "Neurod6"),
  Excitatory       = c("Satb2", "Bcl11b", "Slc17a7"),
  Inhibitory       = c("Gad1", "Gad2", "Dlx1", "Dlx2"),
  Cajal_Retzius    = c("Reln", "Trp73"),
  Astro_precursor  = c("Gfap", "Aldh1l1", "Apoe"),
  OPC              = c("Olig1", "Olig2", "Pdgfra"),
  Microglia        = c("Aif1", "Cx3cr1", "C1qa"),
  Vascular         = c("Cldn5", "Pdgfrb", "Col1a1"),
  Erythroid        = c("Hba-a1", "Hba-a2", "Car2")
)

## keep only genes also found in the count matrix
marker_genes <- lapply(
  marker_genes,
  function(vec) intersect(vec, rownames(sobj_filtrd))
)

VlnPlot(sobj_filtrd, features = marker_genes$Radial_glia)
FeaturePlot(sobj_filtrd, features = marker_genes$Radial_glia)

# While VlnPlot/FeaturePlot let you inspect one marker at a time, DotPlot
# summarizes all markers across all clusters at once — useful once you have
# several candidate cell types to compare.
#
# **How to read the DotPlot**: the **size** of each dot shows the percentage of
# cells in the cluster expressing the gene; the **color** shows the average
# scaled expression level. A good marker for a given cluster is one that is both
# large and dark within that cluster, and small/pale elsewhere. Genes marking
# rare populations (vascular cells, microglia) will always have small dots
# overall — that reflects the small size of the population, not a weak marker.
DotPlot(sobj_filtrd, features = marker_genes) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
    strip.text = element_text(angle = 45, hjust = 0.5)
  )

# **Question**: For each cluster, which marker set(s) show the clearest signal?
# Are there clusters where no marker set stands out — and if so, could that
# reflect a genuine transitional/undecided cell state rather than a technical
# issue?
#
# **Question**: Which clusters would you expect to form well-separated islands
# on the UMAP, and which would you expect to form a continuum? Check this
# against the UMAP you generated earlier — does it match your expectation?
#
# Individual markers can be noisy or lowly expressed. `AddModuleScore` combines
# several markers of the same cell type into a single per-cell score, which is
# more robust than relying on one gene alone — particularly useful for
# progenitor populations where no single gene is perfectly specific.
#
# The `AddModuleScore()` function calculates a score for each cell based on the
# expression of a list of genes. This score, automatically added to the
# metadata, can then be plotted on the UMAP.
sobj_filtrd <- AddModuleScore(sobj_filtrd,
  features = marker_genes, pool = NULL,
  nbin = 5, seed = 1, k = FALSE, name = "program"
)

## score for each list of gene is named program1, program2...
## rename the columns in metadata to get more interpretable names
vec_rename <- paste0("program", seq_along(marker_genes))
names(vec_rename) <- names(marker_genes)
sobj_filtrd@meta.data <- sobj_filtrd@meta.data |>
  dplyr::rename(all_of(vec_rename))
head(sobj_filtrd@meta.data)

for (f in names(marker_genes)) {
  print(VlnPlot(sobj_filtrd, features = f))
  print(FeaturePlot(sobj_filtrd, features = f))
}

# Exercise: - Which resolution would you choose ?
# - Are the marker specific ?
#
# Once you've interpreted the DotPlot and module scores, modify the code below
# and rename the cluster identities based on your conclusions:
#
# Exercise: complete the code below, then set `run_exercise <- TRUE`
run_exercise <- FALSE
if (run_exercise) {
  head(x = sobj_filtrd@active.ident)
  levels(sobj_filtrd@active.ident)

  # modify the code below and rename ALL the clusters based on your
  # DotPlot interpretation
  sobj_filtrd <- RenameIdents(
    object = sobj_filtrd,
    "0" = ,
    "1" = ,
    "2" = ,
    "3" = ,
    "4" = ,
    "5" = ,
    "6" = ,
    "7" = ,
    "8" = ,
    "9" = ,
    "10" = ,
    "11" =
    )

  head(x = sobj_filtrd@active.ident)
  DimPlot(sobj_filtrd, reduction = "umap", label = TRUE)
}

# For some tissues, the different cell types have already been largely described
# and databases exist with lists of referenced marker genes or datasets. This is
# particularly true for human and mouse.
#
# Instead of manual annotation, these datasets (e.g. Panglao, celldex) can be
# used, along with published tools (e.g. CelliD, SingleR) to call the cell types
# in your data. This process is called automatic annotation.
#
# Your annotation above was based on a handful of marker genes chosen from prior
# knowledge. This is a good starting point, but it's inherently limited: you
# only looked at genes you already expected to be relevant.
#
# Differential expression analysis reverses the logic — instead of asking "does
# this cluster express the genes I expect?", it asks "what genes actually
# distinguish this cluster from all the others, regardless of whether you
# already know them?" This is a useful way to **validate** your manual
# annotation (do the top DE genes for a cluster match the markers you used to
# name it?) and to **discover** genes you hadn't considered, including markers
# for populations you may have left ambiguous.

# %% Differential expression -----------------------------------------------
#
# Differential expression will help you to validate your clusters. By default
# you can use a Wilcoxon test, but if you have several replicates, packages such
# as DESeq2 are more appropriate.
#
# 1 group versus the other min.pct = 0.2 : statistical test is only performed on
# genes expressed in at least 20% of cells in either of the two populations
result_de <- FindAllMarkers(sobj_filtrd, min.pct = 0.2)

head(result_de, n = 5)
tail(result_de)

# Visualize the results. To properly check whether differential expression
# supports your manual annotation, we need the top markers for **each** cluster
# — not just the first rows of the table, which would mostly reflect a single
# cluster.
top_de <- result_de |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 5)

DoHeatmap(subset(sobj_filtrd),
  slot = "data", features = unique(top_de$gene),
  size = 3
)

# **Question**: Compare the top DE genes for each cluster to the marker genes
# you used earlier for annotation. Do they overlap? For clusters where they
# don't overlap well, what could explain the difference — a weak manual marker
# choice, or a genuinely new biological insight from the unbiased DE analysis?

# %% Save ------------------------------------------------------------------
#
# Save seurat object
save(sobj_filtrd, file = paste0(sample_name, "_FILTERED_CLUSTERS.rda"))

# %% Session info ----------------------------------------------------------
sessionInfo()
