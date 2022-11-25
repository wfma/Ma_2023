#### Library ####
# load other libraries i like to have
suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(Seurat) # CRAN
  library(patchwork) # CRAN
  library(readr) # CRAN
  library(SingleR) # BIOCONDUCTOR
  library(tidyverse) # CRAN
  library(monocle3) # SPECIFIC INSTALLATION ON WEBSITE
  library(SeuratData) # satijalab/seurat-data
  library(magrittr)# CRAN
  library(ggrepel)# CRAN
  # library(dyno) # devtools::install_github("dynverse/dyno")
  library(SeuratDisk) # remotes::install_github("mojaveazure/seurat-disk")
  library(celldex) # BiocManager::install("celldex")
  library(data.table) # CRAN
  library(matrixStats)# CRAN
  library(Matrix)# CRAN
  # library(bayNorm) # for transposition of sparase matrix
  library(future)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(AUCell) # BiocManager::install("AUCell")
  # library(STRINGdb) # BiocManager::install("STRINGdb")
  # library(DoubletFinder) #remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  library(scCustomize) # devtools::install_github(repo = "samuel-marsh/scCustomize")
  
  set.seed(1234)
})

# color definitions
original_color_list <-
  {c("rosybrown2",
     "cadetblue1",
     "lemonchiffon3",
     "darkseagreen",
     "skyblue3",
     "thistle3",
     "cadetblue3",
     "darkseagreen1",
     "palevioletred3",
     "palevioletred1",
     "darkseagreen2",
     "rosybrown3",
     "thistle2",
     "lightsteelblue3",
     "salmon1",
     "palevioletred4",
     "lemonchiffon4",
     "cadetblue2"
  )}

color_function <- colorRampPalette(original_color_list)
# color_function <- colorRampPalette(metcolors)

manual_color_list_extended <- color_function(22) # change this if clusters >40

manual_color_list <-
  {c("rosybrown2",
     "cadetblue1",
     "lemonchiffon3",
     "darkseagreen",
     "skyblue3",
     "thistle3",
     "cadetblue3",
     "darkseagreen1",
     "palevioletred3",
     "palevioletred1",
     "darkseagreen2",
     "rosybrown3",
     "thistle2",
     "lightsteelblue3",
     "salmon1",
     "palevioletred4",
     "lemonchiffon4",
     "cadetblue2"
  )}


#### READ DATA ####
Ma_2023 <- readRDS(file = "Ma_2023_RNA_Obj_Pt2_Annotated_Rivanna.rds")

#### Run Monocle on the Dataset ####
expressiondata <- Ma_2023@assays[["RNA"]]@data
cellmd <- Ma_2023@meta.data
genemd <- data.frame(gene_short_name = row.names(expressiondata), 
                     row.names = row.names(expressiondata))

# Construct monocle cds
Ma_2023.cds <- new_cell_data_set(expression_data = expressiondata,
                                     cell_metadata = cellmd,
                                     gene_metadata = genemd)

Ma_2023.cds <- preprocess_cds(Ma_2023.cds, num_dim = 30) # we used 30 in earlier seurat scripts

# 
# run clustering again (didnt transfer from seurat)
Ma_2023.cds <- reduce_dimension(Ma_2023.cds, reduction_method = "UMAP")
Ma_2023.cds <- cluster_cells(Ma_2023.cds, reduction_method = "UMAP")


#### STEP4A: TRANSFER SEURAT EMBEDDINGS ###
# Note that these may be calculated on the Integrated object, not the counts
#   and thus will involve fewer genes
temp.cds <- ProjectDim(Ma_2023, reduction = "pca") # this will be removed
reducedDim(Ma_2023.cds, type = "PCA") <- temp.cds@reductions$pca@cell.embeddings
Ma_2023.cds@reduce_dim_aux$prop_var_expl <- temp.cds@reductions$pca@stdev
plot_pc_variance_explained(Ma_2023.cds)

# Transfer Seurat UMAP embeddings
Ma_2023.cds@int_colData@listData$reducedDims$UMAP <- temp.cds@reductions$umap@cell.embeddings

## transfer seurat labels to moncle3 object
colData(Ma_2023.cds)$assigned_cell_type <- 
  Ma_2023$Author_Provided_Clusters_LVL1 # call this by opening the object

most.common.cell <- names(sort(table(Ma_2023@meta.data[["Author_Provided_Clusters_LVL1"]]),decreasing = T)[1])

#### MONOCLE3 CONT. ---
# now learn the PATH (trajectory)
Ma_2023.cds <- learn_graph(Ma_2023.cds)

Ma_2023.cds <- order_cells(Ma_2023.cds, reduction_method = "UMAP")

plot_cells(Ma_2023.cds, color_cells_by = "Author_Provided_Clusters_LVL1",reduction_method = "UMAP", group_label_size = 3, alpha = 0.8,
           label_groups_by_cluster=FALSE)
plot_cells(Ma_2023.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

#### Monocle3 just the fibroblasts ####
cds_sub <- choose_cells(Ma_2023.cds) # will need manual selection 
cds_sub <- cluster_cells(cds_sub) 
cds_sub <- learn_graph(cds_sub)
cds_sub <- order_cells(cds_sub, reduction_method = "UMAP")

plot_cells(cds_sub, color_cells_by = "Author_Provided_Clusters_LVL1",reduction_method = "UMAP", group_label_size = 3, alpha = 0.8,
           label_groups_by_cluster=FALSE)

pdf(file = "OUTPUT-RNA-Trajectory-FB_Pseudotime_Subset.pdf", width = 8, height = 5)
plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()
#### Monocle3 just the endothelial ####
cds_sub <- choose_cells(Ma_2023.cds) # will need manual selection 
cds_sub <- cluster_cells(cds_sub) 
cds_sub <- learn_graph(cds_sub)
cds_sub <- order_cells(cds_sub, reduction_method = "UMAP")

plot_cells(cds_sub, color_cells_by = "Author_Provided_Clusters_LVL1",reduction_method = "UMAP", group_label_size = 3, alpha = 0.8,
           label_groups_by_cluster=FALSE)

pdf(file = "OUTPUT-RNA-Trajectory-Endothelial_Pseudotime_Subset.pdf", width = 6, height = 5)
plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()
#### Monocle3 just the macrophages ####
cds_sub <- choose_cells(Ma_2023.cds) # will need manual selection 
cds_sub <- cluster_cells(cds_sub) 
cds_sub <- learn_graph(cds_sub)
cds_sub <- order_cells(cds_sub, reduction_method = "UMAP")

plot_cells(cds_sub, color_cells_by = "Author_Provided_Clusters_LVL1",reduction_method = "UMAP", group_label_size = 3, alpha = 0.8,
           label_groups_by_cluster=FALSE)

pdf(file = "OUTPUT-RNA-Trajectory-Macrophages_Pseudotime_Subset.pdf", width = 6, height = 5)
plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()
#### Monocle3 just the CM ####
cds_sub <- choose_cells(Ma_2023.cds) # will need manual selection 
cds_sub <- cluster_cells(cds_sub) 
cds_sub <- learn_graph(cds_sub)
cds_sub <- order_cells(cds_sub, reduction_method = "UMAP")

plot_cells(cds_sub, color_cells_by = "Author_Provided_Clusters_LVL1",reduction_method = "UMAP", group_label_size = 3, alpha = 0.8,
           label_groups_by_cluster=FALSE)

pdf(file = "OUTPUT-RNA-Trajectory-CMs_Pseudotime_Subset.pdf", width = 6, height = 5)
plot_cells(cds_sub,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()


#### LIGER: FB ####
Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
subset.cells <- subset(Ma_2023, idents = "Fibroblast")

Ma_2023_liger <- seuratToLiger(subset.cells, combined.seurat = T, meta.var = "Status")
Ma_2023_liger <- normalize(Ma_2023_liger)
Ma_2023_liger <- selectGenes(Ma_2023_liger)
Ma_2023_liger <- scaleNotCenter(Ma_2023_liger)
Ma_2023_liger <- optimizeALS(Ma_2023_liger, k = 20)
Ma_2023_liger <- quantile_norm(Ma_2023_liger)
Ma_2023_liger <- louvainCluster(Ma_2023_liger, resolution = 0.25)
Ma_2023_liger <- runUMAP(Ma_2023_liger, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)

pdf(file = "OUTPUT-RNA-LIGER-Intg-Compare-FBs_Gene_Loadings.pdf", width = 8, height = 6)
plotGeneLoadings(Ma_2023_liger, return.plots = FALSE)
dev.off()

gene_loadings <- plotGeneLoadings(Ma_2023_liger, do.spec.plot = FALSE, return.plots = TRUE)
gene_loadings[[7]]

pdf(file = "OUTPUT-RNA-LIGER-Intg-Compare-FBs.pdf", width = 28, height = 8)
DimPlot(Ma_2023, group.by = c("Status", "Donor", "Author_Provided_Clusters_LVL2"), label = T, ncol = 3)
dev.off()
#### GSEA with irGSEA ####
# see installation here https://github.com/chuiqin/irGSEA
library(UCell)
library(irGSEA)

test <- irGSEA.score(object = Ma_2023, assay = "SCT", 
                             slot = "data", seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens", category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

pdf(file = "OUTPUT-RNA-irGSEA-HeatMap-by-Author_Provided_Clusters_LVL1.pdf", width = 9, height = 8)
result.dge <- irGSEA.integrate(object = test, 
                               group.by = "Author_Provided_Clusters_LVL1",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot
dev.off()

pdf(file = "OUTPUT-RNA-irGSEA-HeatMap-by-Author_Provided_Clusters_LVL2.pdf", width = 9, height = 8)
result.dge <- irGSEA.integrate(object = test, 
                               group.by = "Author_Provided_Clusters_LVL2",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot
dev.off()


pdf(file = "OUTPUT-RNA-irGSEA-HeatMap-by-Status.pdf", width = 13, height = 8)
result.dge <- irGSEA.integrate(object = test, 
                               group.by = "Status",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot
dev.off()

#### Assigning additional level by status ####
test$Author_Provided_Clusters_LVL2_with_Status <- paste(Ma_2023$Author_Provided_Clusters_LVL2, Ma_2023$Status, Sep = "-")

pdf(file = "OUTPUT-RNA-irGSEA-HeatMap-by-Author_Provided_Clusters_LVL2_by_Status.pdf", width = 20, height = 8)
result.dge <- irGSEA.integrate(object = test, 
                               group.by = "Author_Provided_Clusters_LVL2_with_Status",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL) + guides(fill=guide_legend(ncol=2))

irGSEA.heatmap.plot
dev.off()

 


#### Cluster Profiler on our FBs####
# prep the data into a DEG
Ma_2023$Author_Provided_Clusters_LVL1_with_Status <- paste(Ma_2023$Status, Ma_2023$Author_Provided_Clusters_LVL1, sep = "-")
Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1_with_Status"

pbmc.markers <- FindAllMarkers(Ma_2023, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top100 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)

# ClusterProfiler
library("clusterProfiler") # BiocManager::install("clusterProfiler")
library("org.Hs.eg.db") # BiocManager::install("org.Hs.eg.db")
library("AnnotationHub") # BiocManager::install("AnnotationHub")
library(ReactomePA) # BiocManager::install("ReactomePA")


df <- top100pval[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample)

dfsample$`Normal-Pericyte` = bitr(dfsample$`Normal-Pericyte`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Normal-Endothelial` = bitr(dfsample$`Normal-Endothelial`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Normal-T` = bitr(dfsample$`Normal-T`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Normal-Cardiomyocyte` = bitr(dfsample$`Normal-Cardiomyocyte`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Normal-Macrophage` = bitr(dfsample$`Normal-Macrophage`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Normal-Fibroblast` = bitr(dfsample$`Normal-Fibroblast`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Normal-VSMC` = bitr(dfsample$`Normal-VSMC`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Normal-Unknown` = bitr(dfsample$`Normal-Unknown`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Normal-Neuronal` = bitr(dfsample$`Normal-Neuronal`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Normal-Lymphatic` = bitr(dfsample$`Normal-Lymphatic`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Normal-Adipocyte` = bitr(dfsample$`Normal-Adipocyte`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

dfsample$`Chemo-Pericyte` = bitr(dfsample$`Chemo-Pericyte`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Chemo-Endothelial` = bitr(dfsample$`Chemo-Endothelial`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Chemo-T` = bitr(dfsample$`Chemo-T`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Chemo-Cardiomyocyte` = bitr(dfsample$`Chemo-Cardiomyocyte`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Chemo-Macrophage` = bitr(dfsample$`Chemo-Macrophage`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Chemo-Fibroblast` = bitr(dfsample$`Chemo-Fibroblast`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Chemo-VSMC` = bitr(dfsample$`Chemo-VSMC`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Chemo-Unknown` = bitr(dfsample$`Chemo-Unknown`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Chemo-Neuronal` = bitr(dfsample$`Chemo-Neuronal`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Chemo-Lymphatic` = bitr(dfsample$`Chemo-Lymphatic`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$`Chemo-Adipocyte` = bitr(dfsample$`Chemo-Adipocyte`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")



#### Plot for Fibroblasts ####
genelist <- list("CH-Fibroblast" = dfsample$`Chemo-Fibroblast`$ENTREZID, 
                 "N-Fibroblast" = dfsample$`Normal-Fibroblast`$ENTREZID
                 
)

pdf(file = "OUTPUT-RNA-ClusterProfiler_Fibroblasts_by_Status_enrichGO.pdf", width = 6, height = 6)
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
dotplot(GOclusterplot)
dev.off()

pdf(file = "OUTPUT-RNA-ClusterProfiler_Fibroblasts_by_Status_KEGG.pdf", width = 6, height = 6)
KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
dotplot(KEGGclusterplot)
dev.off()

pdf(file = "OUTPUT-RNA-ClusterProfiler_Fibroblasts_by_Status_Pathway.pdf", width = 6, height = 6)
RXlusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway")
dotplot(RXlusterplot)
dev.off()

#### Cluster Profiler on our for CMs ####
genelist <- list("CH-CM" = dfsample$`Chemo-Cardiomyocyte`$ENTREZID, 
                 "N-CM" = dfsample$`Normal-Cardiomyocyte`$ENTREZID
                 
)

pdf(file = "OUTPUT-RNA-ClusterProfiler_Cardiomyocyte_by_Status_enrichGO.pdf", width = 6, height = 6)
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
dotplot(GOclusterplot)
dev.off()

pdf(file = "OUTPUT-RNA-ClusterProfiler_Cardiomyocyte_by_Status_KEGG.pdf", width = 6, height = 6)
KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
dotplot(KEGGclusterplot)
dev.off()

pdf(file = "OUTPUT-RNA-ClusterProfiler_Cardiomyocyte_by_Status_Pathway.pdf", width = 6, height = 6)
RXlusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway")
dotplot(RXlusterplot)
dev.off()
