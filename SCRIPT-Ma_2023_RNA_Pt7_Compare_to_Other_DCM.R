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
  library(CellChat) # devtools::install_github("sqjin/CellChat")
  #library(STRINGdb) # BiocManager::install("STRINGdb")
  # library(DoubletFinder) #remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  
  set.seed(1234)
})

library(CIPR)

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

library(CellChat)
#### READ DATA ####
Ma_2023 <- readRDS(file = "Ma_2023_RNA_Obj_Pt2_Annotated_Rivanna.rds")
koenig <- readRDS(file = "../../../../Desktop/DataProcessing/data/Koenig_2022/source_files/UNPROCESSED.RDS") # this one i special reach each time bc its sooo big. i dont put it in the same directory https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183852


# #### Subset Koenig for FB DCM Cell Types ####
# Idents(koenig) <- "Author_Provided"
# table(koenig$Names)
# 
# koenig.fibroblasts <- subset(koenig, idents = "Fibroblasts")
# # koenig.cardiomyocoytes <- subset(koenig, idents = "Macrophages")
# # koenig.endothelium <- subset(koenig, idents = "Endothelium")
# # koenig.macrophages <- subset(koenig, idents = c("Macrophages", "Monocytes"))
# # koenig.pericytes <- subset(koenig, idents = "Pericytes")
# # koenig.smc <- subset(koenig, idents = "Smooth_Muscle")
# # koenig.t <- subset(koenig, idents = "T/NK_Cells")
# 
# 
# Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
# table(Ma_2023$Author_Provided_Clusters_LVL1)
# Ma_2023.fibroblasts <- subset(Ma_2023, idents = "Fibroblast")
# # Ma_2023.cardiomyocoytes <- subset(Ma_2023, idents = "Macrophage")
# # Ma_2023.endothelium <- subset(Ma_2023, idents = "Endothelial")
# # Ma_2023.macrophages <- subset(Ma_2023, idents = "Macrophage")
# # Ma_2023.pericytes <- subset(Ma_2023, idents = "Pericyte")
# # Ma_2023.smc <- subset(Ma_2023, idents = "VSMC")
# # Ma_2023.t <- subset(Ma_2023, idents = "T")
# 
# 
# #### Integrate fibroblasts ####
# koenig.fibroblasts$Status <- koenig.fibroblasts$Condition
# 
# koenig.fibroblasts$combined.idents <- paste(koenig.fibroblasts$Condition, koenig.fibroblasts$Author_Provided, sep = "-")
# Ma_2023.fibroblasts$combined.idents <- paste(Ma_2023.fibroblasts$Status, Ma_2023.fibroblasts$Author_Provided_Clusters_LVL1, sep = "-")
# 
# 
# 
# ## integrative analysis with LIGER 
# tointegrate <- list(Ma_2023.fibroblasts, koenig.fibroblasts)
# tointegrate <- lapply(X = tointegrate, FUN = SCTransform)
# features <- SelectIntegrationFeatures(object.list = tointegrate, nfeatures = 3000)
# tointegrate <- PrepSCTIntegration(object.list = tointegrate, anchor.features = features)
# anchors <- FindIntegrationAnchors(object.list = tointegrate, anchor.features = features, normalization.method = "SCT")
# integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# integrated <- ScaleData(integrated, verbose = FALSE)
# integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
# integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
# integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
# integrated <- FindClusters(integrated, resolution = 0.5)
# 
# pdf(file = "OUTPUT-RNA-iDCM_Compare_fibroblast_by_Status.pdf", width = 8, height = 3)
# DimPlot(integrated, reduction = "umap", group.by = "combined.idents", cols = manual_color_list, 
#         split.by = "combined.idents") +
#   labs(title = "Fibroblasts by Status") +
#   theme_bw() +
#   theme(legend.position="none") 
# dev.off()
# 
# Idents(integrated) <- "combined.idents"
# markers <- FindAllMarkers(integrated)
# write.csv(markers, file = "OUTPUT-RNA-iDCM_Compare_fibroblast_DEGS.csv")
# 
# library(rliger) # install.packages('rliger')
# ligerobj <- seuratToLiger(list(Chemo = Ma_2023.fibroblasts, DCM = koenig.fibroblasts))
# ligerobj <- normalize(ligerobj)
# ligerobj <- selectGenes(ligerobj)
# ligerobj <- scaleNotCenter(ligerobj)
# ligerobj <- optimizeALS(ligerobj, k = 20)
# ligerobj <- quantile_norm(ligerobj)
# ligerobj <- louvainCluster(ligerobj, resolution = 0.25)
# ligerobj <- runUMAP(ligerobj, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
# all.plots <- plotByDatasetAndCluster(ligerobj, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
# all.plots[[1]] + all.plots[[2]]
# 
# 
# 
# #### FB ClusterProfiler ####
# 
# Idents(integrated) <- "combined.idents"
# 
# pbmc.markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
# write.csv(pbmc.markers, file = "OUTPUT-RNA-Difflist_iDCM_vs_Chemo_FBs.csv")
# top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# write.csv(top5, file = "OUTPUT-RNA-Difflist_iDCM_vs_Chemo_FBs_TOP5.csv")
# 
# top100 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC)
# top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)
# 
# library("clusterProfiler") # BiocManager::install("clusterProfiler")
# library("org.Hs.eg.db") # BiocManager::install("org.Hs.eg.db")
# library("AnnotationHub") # BiocManager::install("AnnotationHub")
# library(ReactomePA) # BiocManager::install("ReactomePA")
# 
# 
# df <- top100pval[,7:6]
# dfsample <- split(df$gene,df$cluster)
# length(dfsample)
# 
# dfsample$`Normal-Fibroblast` = bitr(dfsample$`Normal-Fibroblast`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# dfsample$`Donor-Fibroblasts` = bitr(dfsample$`Donor-Fibroblasts`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# dfsample$`Chemo-Fibroblast` = bitr(dfsample$`Chemo-Fibroblast`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# dfsample$`DCM-Fibroblasts` = bitr(dfsample$`DCM-Fibroblasts`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# 
# 
# genelist <- list("Chemo-FB" = dfsample$`Chemo-Fibroblast`$ENTREZID, 
#                  # "Normal-FB" = dfsample$`Normal-Fibroblast`$ENTREZID,
#                  "DCM-FB" = dfsample$`DCM-Fibroblast`$ENTREZID
#                  # "Donor-FB" = dfsample$`Donor-Fibroblast`$ENTREZID
#                  
#                  
# )
# 
# pdf(file = "OUTPUT-RNA-ClusterProfiler_Fibroblasts_with-iDCM_enrichGO.pdf", width = 6, height = 6)
# GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
# dotplot(GOclusterplot)
# dev.off()
# 
# pdf(file = "OUTPUT-RNA-ClusterProfiler_Fibroblasts_with-iDCM_KEGG.pdf", width = 6, height = 6)
# KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
# dotplot(KEGGclusterplot)
# dev.off()
# 
# pdf(file = "OUTPUT-RNA-ClusterProfiler_Fibroblasts_with-iDCM_Pathway.pdf", width = 6, height = 6)
# RXlusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway")
# dotplot(RXlusterplot)
# dev.off()
# 
# 
# #### FB Oxidative Pathway ####
# pdf(file = "OUTPUT-RNA-Feature_Oxidative_NOX1_by_Author_Provided.pdf", width = 5, height = 5)
# Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
# FeaturePlot(Ma_2023, feature = "NOX1", label = T, repel = T)
# dev.off()
# 
# pdf(file = "OUTPUT-RNA-Feature_Oxidative_NOX2_by_Author_Provided.pdf", width = 5, height = 5)
# Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
# FeaturePlot(Ma_2023, feature = "CYBB", label = T, repel = T)
# dev.off()
# 
# pdf(file = "OUTPUT-RNA-Feature_Oxidative_NOX4_by_Author_Provided.pdf", width = 5, height = 5)
# Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
# FeaturePlot(Ma_2023, feature = "NOX4", label = T, repel = T)
# dev.off()
# 
# pdf(file = "OUTPUT-RNA-Feature_Oxidative_NOX5_by_Author_Provided.pdf", width = 5, height = 5)
# Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
# FeaturePlot(Ma_2023, feature = "NOX5", label = T, repel = T)
# dev.off()
# 
# pdf(file = "OUTPUT-RNA-Feature_Oxidative_CYBA_by_Author_Provided.pdf", width = 5, height = 5)
# Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
# FeaturePlot(Ma_2023, feature = "CYBA", label = T, repel = T)
# dev.off()
# 
# 
# pdf(file = "OUTPUT-RNA-Feature_Oxidative_MPO_by_Author_Provided.pdf", width = 5, height = 5)
# Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
# FeaturePlot(Ma_2023, feature = "MPO", label = T, repel = T)
# dev.off()
# 
# pdf(file = "OUTPUT-RNA-Feature_Oxidative_NOX4_by_Author_Provided_splitby_status.pdf", width = 8, heighot = 5)
# VlnPlot(Ma_2023, features = "NOX4", split.by = "Status", log = F)
# dev.off()
# 
# pdf(file = "OUTPUT-RNA-Feature_Oxidative_CYBA_by_Author_Provided_splitby_status.pdf", width = 8, height = 5)
# VlnPlot(Ma_2023, features = "CYBA", split.by = "Status", log = F)
# dev.off()
# 
# ## compare all fibroblasts oxidative pathways 
# pdf(file = "OUTPUT-RNA-Feature_Oxidative_ALL_by_Author_Provided_FB_ONLY.pdf", width = 10, height = 7)
# Ma_2023.fibroblasts <- subset(Ma_2023, idents = "Fibroblast")
# features <- c("NOX1", "CYBB", "MPO", "NOS1", "NOS2", "NOS3")
# Idents(Ma_2023.fibroblasts) <- "Author_Provided_Clusters_LVL1_with_Status"
# VlnPlot(
#   Ma_2023.fibroblasts,
#   features = features,
#   split.by = "Status") 
# dev.off()
# 
# #### FB GTPase Activity ####
# pdf(file = "OUTPUT-RNA-Feature_GPCR_GLUT2_by_Author_Provided.pdf", width = 10, height = 5)
# Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
# FeaturePlot(Ma_2023, feature = "SLC2A2", label = T, repel = T, split.by = "Status")
# dev.off()
# 
# pdf(file = "OUTPUT-RNA-Feature_GPCR_GLUT1_by_Author_Provided.pdf", width = 10, height = 5)
# Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
# FeaturePlot(Ma_2023, feature = "SLC2A1", label = T, repel = T, split.by = "Status")
# dev.off()
# 
# # co-expression for GPCR
# FeaturePlot(Ma_2023, features = c("MTERF1", "SLC2A1"), blend = TRUE, order = T)
# 
##### Cluster Profiler Cardiomyocytes
library("clusterProfiler") # BiocManager::install("clusterProfiler")
library("org.Hs.eg.db") # BiocManager::install("org.Hs.eg.db")
library("AnnotationHub") # BiocManager::install("AnnotationHub")
library(ReactomePA) # BiocManager::install("ReactomePA")

