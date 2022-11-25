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
  # library(CellChat) # devtools::install_github("sqjin/CellChat")
  # library(STRINGdb) # BiocManager::install("STRINGdb")
  # library(DoubletFinder) #remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  
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


#### Compare FB Clusters ####
# chemo/dx by cluster divided
calcpercluster <- function(cluster = NULL){
  sub <- subset(Ma_2023, subset = Seurat_Clusters == cluster)

  Idents(object = sub) <- "Status"
  difflist <-  tryCatch(FindMarkers(sub, ident.1 = "Normal", ident.2 = "Chemo"),  error=function(e) NULL)
  difflist$GeneName <- rownames(difflist)
  tryCatch(write_csv(difflist, 
            file = paste("OUTPUT-RNA-Difflist_1NORMAL_2CHEMO_Cluster_", 
                         cluster, ".csv", sep = "")), error=function(e) NULL)
}

lapply(1:28, calcpercluster)

# chemo/dx by author-provided divided
calcpercluster <- function(cluster = NULL){
  sub <- subset(Ma_2023, subset = Author_Provided_Clusters_LVL1 == cluster)
  
  Idents(object = sub) <- "Status"
  difflist <-  tryCatch(FindMarkers(sub, ident.1 = "Normal", ident.2 = "Chemo"),  error=function(e) NULL)
  difflist$GeneName <- rownames(difflist)
  tryCatch(write_csv(difflist, 
                     file = paste("OUTPUT-RNA-Difflist_1NORMAL_2CHEMO_Cluster_", 
                                  cluster, "_LVL1.csv", sep = "")), error=function(e) NULL)
}

clusterstocalculate <- c("Fibroblast", "Macrophage", "Pericyte", "Cardiomyocyte", "Endothelial")
lapply(clusterstocalculate, calcpercluster)

## comparision fb 1 vs 2,3,4

Idents(object = Ma_2023) <- "Author_Provided_Clusters_LVL2"
difflist <-  FindMarkers(Ma_2023, ident.1 = "FB_1", ident.2 = c("FB_2", "FB_3", "FB_4"))
difflist$GeneName <- rownames(difflist)
write.csv(difflist, file = "OUTPUT-RNA-Difflist_FB1_vs_FB2-4.csv")
#### Compare Cell Exp. Changes to Dx ####
library(cowplot)
theme_set(theme_cowplot())
Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
CM_subset <- subset(Ma_2023, idents = c("Cardiomyocyte"))
Idents(Ma_2023) <- "Status"
avg_expression <- as.data.frame(log1p(AverageExpression(Ma_2023, verbose = FALSE)$RNA))
avg_expression$gene <- rownames(avg_expression)

pdf(file = "OUTPUT-RNA-Volcano_CMs.pdf", width = 5, height = 5)
p1 <- ggplot(avg_expression, aes(Normal, Chemo)) + geom_point(alpha = 0.6) + ggtitle("Cardiomyocye Response")
genes.to.label = c("ANKRD1", "MYH9", "DLG1", "COX6A2", "ACTA1", "CRYAB", "NPPA")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 2, color = "blue")
p1
dev.off()

## add cell type + dx label ###
Ma_2023$Author_Provided_Clusters_DX <- paste(Ma_2023$Author_Provided_Clusters, Ma_2023$Status, sep = "_")
Idents(Ma_2023) <- "Author_Provided_Clusters_DX"

CM_normal <- paste(c("CM_1", "CM_2", "CM_4", "CM_5"), "_Normal", sep = "") # theres no CM3_normal
CM_chemo <- paste(c("CM_2", "CM_3", "CM_4", "CM_5"), "_Chemo", sep = "") # theres no CM1_chemo

CM_response_markers <- FindMarkers(Ma_2023, ident.1 = CM_normal, ident.2 = CM_chemo, verbose = T)
head(CM_response_markers, n = 15)
write.csv(CM_response_markers, file = "OUTPUT-RNA_ONLY/OUTPUT-RNA-Difflist-CMs_response_1Normal2Chemo.csv")


#### Save Seurat Obj ####
saveRDS(Ma_2023, file = "Ma_2023_RNA_Obj_Pt4_Condition_Compared_Rivanna.rds")

