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

#### CowPlot for CM Response ####

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
CM <- subset(Ma_2023, idents = "Cardiomyocyte")

Idents(CM) <- "Status"
avg.CM <- as.data.frame(log1p(AverageExpression(CM, verbose = FALSE)$RNA))
avg.CM$gene <- rownames(avg.CM)

pdf(file = "OUTPUT-RNA-Cowplots-CM-Response.pdf", height = 5, width = 5)
genes.to.label = c("NPPA", "ACTA1", "ANKRD1", "XIST", "CHRM2", "SLC27A6", "TNNC1", "MYL2", "COX6A2", "MT-ND3")
p1 <- ggplot(avg.CM, aes(Normal, Chemo)) + geom_point(alpha = 0.6, color = "black") + ggtitle("Cardimyocyte Changes") + theme_bw()
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, color = "Red" ) 
p1
dev.off()




