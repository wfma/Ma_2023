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






#### Data Loading ####
# Batch1
counts.n1 <- Read10X_h5( "Data/CellRanger_BASE_COUNTS_N1/outs/filtered_feature_bc_matrix.h5")
counts.n3 <- Read10X_h5( "Data/CellRanger_BASE_COUNTS_N3/outs/filtered_feature_bc_matrix.h5")
counts.ch1 <- Read10X_h5( "Data/CellRanger_BASE_COUNTS_CH1/outs/filtered_feature_bc_matrix.h5")
counts.ch3 <- Read10X_h5( "Data/CellRanger_BASE_COUNTS_CH3/outs/filtered_feature_bc_matrix.h5")

#Batch 2
counts.ch1B <- Read10X_h5( "Data/CellRanger_BASE_COUNTS_CH1_B/outs/filtered_feature_bc_matrix.h5")
counts.ch3B <- Read10X_h5( "Data/CellRanger_BASE_COUNTS_CH3_B/outs/filtered_feature_bc_matrix.h5")
counts.ch4 <- Read10X_h5( "Data/CellRanger_BASE_COUNTS_CH4/outs/filtered_feature_bc_matrix.h5")
counts.ch5 <- Read10X_h5( "Data/CellRanger_BASE_COUNTS_CH5/outs/filtered_feature_bc_matrix.h5")

# create a Seurat object containing the RNA adata
N1 <- CreateSeuratObject(
  counts = counts.n1,
  assay = "RNA"
)

N3 <- CreateSeuratObject(
  counts = counts.n3,
  assay = "RNA"
)

CH1 <- CreateSeuratObject(
  counts = counts.ch1,
  assay = "RNA",
)

CH3 <- CreateSeuratObject(
  counts = counts.ch3,
  assay = "RNA"
)

CH1B <- CreateSeuratObject(
  counts = counts.ch1B,
  assay = "RNA"
)

CH3B <- CreateSeuratObject(
  counts = counts.ch3B,
  assay = "RNA"
)

CH4 <- CreateSeuratObject(
  counts = counts.ch4,
  assay = "RNA"
)

CH5 <- CreateSeuratObject(
  counts = counts.ch5,
  assay = "RNA"
)

#### Add Metadata ####
N1$Donor <- "N1"
N3$Donor <- "N3"
CH1$Donor <- "CH1"
CH3$Donor <- "CH3"
CH1B$Donor <- "CH1"
CH3B$Donor <- "CH3"
CH4$Donor <- "CH4"
CH5$Donor <- "CH5"

N1$SampleID <- "N1"
N3$SampleID <- "N3"
CH1$SampleID <- "CH1"
CH3$SampleID <- "CH3"
CH1B$SampleID <- "CH1B"
CH3B$SampleID <- "CH3B"
CH4$SampleID <- "CH4"
CH5$SampleID <- "CH5"

N1$Status <- "Normal"
N3$Status <- "Normal"
CH1$Status <- "Chemo"
CH3$Status <- "Chemo"
CH1B$Status <- "Chemo"
CH3B$Status <- "Chemo"
CH4$Status <- "Chemo"
CH5$Status <- "Chemo"

N1$UVA_ID <- "UVA017"
N3$UVA_ID <- "UVA014"
CH1$UVA_ID <- "UVA109"
CH3$UVA_ID <- "UVA174"
CH1B$UVA_ID <- "UVA109"
CH3B$UVA_ID <- "UVA174"
CH4$UVAID <- "UVA120"
CH5$UVA_ID <- "UVA198"

N1$Age <- "55"
N3$Age <- "49"
CH1$Age <- "53"
CH3$Age <- "55"
CH1B$Age <- "53"
CH3B$Age <- "55"
CH4$Age <- "59"
CH5$Age <- "41"

N1$Sex <- "F"
N3$Sex <- "M"
CH1$Sex <- "M"
CH3$Sex <- "F"
CH1B$Sex <- "M"
CH3B$Sex <- "F"
CH4$Sex <- "F"
CH5$Sex <- "F"

N1$Diag_Code <- "None"
N3$Diag_Code <- "None"
CH1$Diag_Code <- "I50.84"
CH3$Diag_Code <- "I42.0"
CH1B$Diag_Code <- "I50.84"
CH3B$Diag_Code <- "I42.0"
CH4$Diag_Code <- "I42.8"
CH5$Diag_Code <- "I50.2"

N1$Batch <- "First"
N3$Batch <- "First"
CH1$Batch <- "First"
CH3$Batch <- "First"
CH1B$Batch <- "Second"
CH3B$Batch <- "Second"
CH4$Batch <- "Second"
CH5$Batch <- "Second"

run1 <- merge(N1, y = c(N3, CH1, CH3), add.cell.ids = c("N1", "N3", "CH1", "CH3"), project = "Ma_2023_ChemoDCM")
run2 <- merge(CH1B, y = c(CH3B, CH4, CH5), add.cell.ids = c("CH1B", "CH3B", "CH4", "CH5"), project = "Ma_2023_ChemoDCM")

run1 <- NormalizeData(run1)
run1 <- FindVariableFeatures(run1, selection.method = "vst", nfeatures = 2000)

run2 <- NormalizeData(run2)
run2 <- FindVariableFeatures(run2, selection.method = "vst", nfeatures = 2000)

features <- SelectIntegrationFeatures(object.list = list(run1, run2))
anchors <- FindIntegrationAnchors(object.list = list(run1, run2), anchor.features = features)

Ma_2023 <- IntegrateData(anchorset = anchors)

Ma_2023
DefaultAssay(Ma_2023) <-"integrated"
table((Ma_2023$Donor))

#### Save Prelim Obj ####
saveRDS(Ma_2023, file = "Ma_2023_RNA_Obj_UNPROCESSED.rds")

#### Read Prelim Obj ####
Ma_2023 <- readRDS(file = "Ma_2023_RNA_Obj_UNPROCESSED.rds")

#### Preprocessing ####
DefaultAssay(Ma_2023) <-"integrated"

# Run standard cleanup (remove low feature/too many feature/too many mt) (not always needed)
Ma_2023[["percent.mt"]] <- PercentageFeatureSet(Ma_2023, pattern = "^MT-", assay = "RNA")

# this just tells us the distribution of counts
VlnPlot(Ma_2023, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# here we will take off the top 1 percentile nfeatures to rm outlier cells
Ma_2023 <- subset(Ma_2023, subset = nFeature_RNA > 200 & nFeature_RNA < quantile(Ma_2023$nFeature_RNA, .99) & percent.mt < 5)

# Run the standard workflow for visualization and clustering
Ma_2023 <- ScaleData(Ma_2023, verbose = FALSE)
Ma_2023 <- RunPCA(Ma_2023, npcs = 30, verbose = FALSE)
Ma_2023 <- RunUMAP(Ma_2023, reduction = "pca", dims = 1:20)
Ma_2023 <- FindNeighbors(Ma_2023, reduction = "pca", dims = 1:20)
Ma_2023 <- FindClusters(Ma_2023, resolution = 0.5)

Ma_2023

table(Ma_2023$Donor)
table(Ma_2023$Status)

write.csv(table(Ma_2023$Donor), file = "OUTPUT-RNA_ONLY/OUTPUT-RNA-CellPop-by-Donors.csv")
write.csv(table(Ma_2023$Status), file = "OUTPUT-RNA_ONLY/OUTPUT-RNA-CellPop-by-Status.csv.csv")

#### Save filtered Obj ####
saveRDS(Ma_2023, file = "Ma_2023_RNA_Obj_Pt1_Filtered.rds")

