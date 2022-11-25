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

#### Druggability Chemo CM   ####

dgidb <- read_delim("DGIDB_Common_HF_Drug_Targets.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

overexp <- read.csv(file = "OUTPUT-RNA-Difflist_Cardiomyocyte_vs_all.csv")
head(overexp)

overexp$gene <- overexp$X

normal_Cardiomyocyte <- dplyr::filter(overexp, avg_log2FC > 0.5)
chemo_Cardiomyocyte <- dplyr::filter(overexp, avg_log2FC < 0.5)

target_normal_Cardiomyocyte <- inner_join(normal_Cardiomyocyte, dgidb, by = "gene")
target_chemo_Cardiomyocyte <- inner_join(chemo_Cardiomyocyte, dgidb, by = "gene")

write.csv(target_normal_Cardiomyocyte, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Normal-Cardiomyocyte.csv")
write.csv(target_chemo_Cardiomyocyte, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Chemo-Cardiomyocyte.csv")

#### Druggability Chemo FB  ####

dgidb <- read_delim("DGIDB_Common_HF_Drug_Targets.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

overexp <- read.csv(file = "OUTPUT-RNA-Difflist_Fibroblast_vs_all.csv")
head(overexp)

overexp$gene <- overexp$X

normal_Fibroblast <- dplyr::filter(overexp, avg_log2FC > 0.5)
chemo_Fibroblast <- dplyr::filter(overexp, avg_log2FC < 0.5)

target_normal_Fibroblast <- inner_join(normal_Fibroblast, dgidb, by = "gene")
target_chemo_Fibroblast <- inner_join(chemo_Fibroblast, dgidb, by = "gene")

write.csv(target_normal_Fibroblast, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Normal-Fibroblast.csv")
write.csv(target_chemo_Fibroblast, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Chemo-Fibroblast.csv")

#### Druggability Chemo Endothelial   ####

dgidb <- read_delim("DGIDB_Common_HF_Drug_Targets.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

overexp <- read.csv(file = "OUTPUT-RNA-Difflist_Endothelial_vs_all.csv")
head(overexp)

overexp$gene <- overexp$X

normal_Endothelial <- dplyr::filter(overexp, avg_log2FC > 0.5)
chemo_Endothelial <- dplyr::filter(overexp, avg_log2FC < 0.5)

target_normal_Endothelial <- inner_join(normal_Endothelial, dgidb, by = "gene")
target_chemo_Endothelial <- inner_join(chemo_Endothelial, dgidb, by = "gene")

write.csv(target_normal_Endothelial, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Normal-Endothelial.csv")
write.csv(target_chemo_Endothelial, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Chemo-Endothelial.csv")

#### Druggability Chemo Macrophage  ####

dgidb <- read_delim("DGIDB_Common_HF_Drug_Targets.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

overexp <- read.csv(file = "OUTPUT-RNA-Difflist_Macrophage_vs_all.csv")
head(overexp)

overexp$gene <- overexp$X

normal_Macrophage <- dplyr::filter(overexp, avg_log2FC > 0.5)
chemo_Macrophage <- dplyr::filter(overexp, avg_log2FC < 0.5)

target_normal_Macrophage <- inner_join(normal_Macrophage, dgidb, by = "gene")
target_chemo_Macrophage <- inner_join(chemo_Macrophage, dgidb, by = "gene")

write.csv(target_normal_Macrophage, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Normal-Macrophage.csv")
write.csv(target_chemo_Macrophage, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Chemo-Macrophage.csv")


#### Druggability Chemo VSMC ####

dgidb <- read_delim("DGIDB_Common_HF_Drug_Targets.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

overexp <- read.csv(file = "OUTPUT-RNA-Difflist_VSMC_vs_all.csv")
head(overexp)

overexp$gene <- overexp$X

normal_VSMC <- dplyr::filter(overexp, avg_log2FC > 0.5)
chemo_VSMC <- dplyr::filter(overexp, avg_log2FC < 0.5)

target_normal_VSMC <- inner_join(normal_VSMC, dgidb, by = "gene")
target_chemo_VSMC <- inner_join(chemo_VSMC, dgidb, by = "gene")

write.csv(target_normal_VSMC, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Normal-VSMC.csv")
write.csv(target_chemo_VSMC, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Chemo-VSMC.csv")



#### Druggability DCM CM   ####

dgidb <- read_delim("DGIDB_Common_HF_Drug_Targets.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

overexp <- read.csv(file = "OUTPUT-RNA_Koenig_Cardiomyocytes_by_Condition.csv")
head(overexp)

overexp$gene <- overexp$X

Donor_Cardiomyocyte <- dplyr::filter(overexp, avg_log2FC > 0.5)
DCM_Cardiomyocyte <- dplyr::filter(overexp, avg_log2FC < 0.5)

target_Donor_Cardiomyocyte <- inner_join(Donor_Cardiomyocyte, dgidb, by = "gene")
target_DCM_Cardiomyocyte <- inner_join(DCM_Cardiomyocyte, dgidb, by = "gene")

write.csv(target_Donor_Cardiomyocyte, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Donor-Cardiomyocyte.csv")
write.csv(target_DCM_Cardiomyocyte, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-DCM-Cardiomyocyte.csv")

#### Druggability DCM FB  ####

dgidb <- read_delim("DGIDB_Common_HF_Drug_Targets.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

overexp <- read.csv(file = "OUTPUT-RNA_Koenig_Fibroblast_by_Condition.csv")
head(overexp)

overexp$gene <- overexp$X

Donor_Fibroblast <- dplyr::filter(overexp, avg_log2FC > 0.5)
DCM_Fibroblast <- dplyr::filter(overexp, avg_log2FC < 0.5)

target_Donor_Fibroblast <- inner_join(Donor_Fibroblast, dgidb, by = "gene")
target_DCM_Fibroblast <- inner_join(DCM_Fibroblast, dgidb, by = "gene")

write.csv(target_Donor_Fibroblast, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Donor-Fibroblast.csv")
write.csv(target_DCM_Fibroblast, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-DCM-Fibroblast.csv")

#### Druggability DCM Endothelial   ####

dgidb <- read_delim("DGIDB_Common_HF_Drug_Targets.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

overexp <- read.csv(file = "OUTPUT-RNA_Koenig_Endothelium_by_Condition.csv")
head(overexp)

overexp$gene <- overexp$X

Donor_Endothelial <- dplyr::filter(overexp, avg_log2FC > 0.5)
DCM_Endothelial <- dplyr::filter(overexp, avg_log2FC < 0.5)

target_Donor_Endothelial <- inner_join(Donor_Endothelial, dgidb, by = "gene")
target_DCM_Endothelial <- inner_join(DCM_Endothelial, dgidb, by = "gene")

write.csv(target_Donor_Endothelial, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Donor-Endothelial.csv")
write.csv(target_DCM_Endothelial, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-DCM-Endothelial.csv")

#### Druggability DCM Macrophage  ####

dgidb <- read_delim("DGIDB_Common_HF_Drug_Targets.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

overexp <- read.csv(file = "OUTPUT-RNA_Koenig_Macrophage_by_Condition.csv")
head(overexp)

overexp$gene <- overexp$X

Donor_Macrophage <- dplyr::filter(overexp, avg_log2FC > 0.5)
DCM_Macrophage <- dplyr::filter(overexp, avg_log2FC < 0.5)

target_Donor_Macrophage <- inner_join(Donor_Macrophage, dgidb, by = "gene")
target_DCM_Macrophage <- inner_join(DCM_Macrophage, dgidb, by = "gene")

write.csv(target_Donor_Macrophage, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Donor-Macrophage.csv")
write.csv(target_DCM_Macrophage, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-DCM-Macrophage.csv")


#### Druggability DCM VSMC ####

dgidb <- read_delim("DGIDB_Common_HF_Drug_Targets.tsv", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

overexp <- read.csv(file = "OUTPUT-RNA_Koenig_Smooth_Muscle_by_Condition.csv")
head(overexp)

overexp$gene <- overexp$X

Donor_VSMC <- dplyr::filter(overexp, avg_log2FC > 0.5)
DCM_VSMC <- dplyr::filter(overexp, avg_log2FC < 0.5)

target_Donor_VSMC <- inner_join(Donor_VSMC, dgidb, by = "gene")
target_DCM_VSMC <- inner_join(DCM_VSMC, dgidb, by = "gene")

write.csv(target_Donor_VSMC, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-Donor-VSMC.csv")
write.csv(target_DCM_VSMC, file = "OUTPUT-RNA_Targeted_by_HF_Drugs-DCM-VSMC.csv")

