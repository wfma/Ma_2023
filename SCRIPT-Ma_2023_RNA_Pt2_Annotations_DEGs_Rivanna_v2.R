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
  library(sctransform)
  
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

manual_color_list_extended <- color_function(28) # change this if clusters >40

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








#### Load Data ####
Ma_2023 <- readRDS(file = "Ma_2023_RNA_Obj_Pt1_Filtered.rds")

#### Future: Setting Parallel Computing ####
plan("multicore", workers = 8)
options(future.globals.maxSize= 25000 * 1024^2)

#### scTransform the data ####
# store mitochondrial percentage in object meta data
Ma_2023 <- PercentageFeatureSet(Ma_2023, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
# run sctransform
Ma_2023 <- SCTransform(Ma_2023, vars.to.regress = "percent.mt", verbose = FALSE)

Ma_2023 <- RunPCA(Ma_2023, verbose = FALSE)
Ma_2023 <- RunUMAP(Ma_2023, dims = 1:30, verbose = FALSE)

Ma_2023 <- FindNeighbors(Ma_2023, dims = 1:30, verbose = FALSE)
Ma_2023 <- FindClusters(Ma_2023, verbose = FALSE)

# #### SingleR ####
# 
# # point this FROM THE ROOT "./" to DataProcessing's 'data' folder
# root.to.data <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/data/"
# root.to.TS.ref <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/references/Tabula_sapiens_reference/TS_Vasculature.h5seurat"
# root.to.TS.mouse.ref <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/references/Tabula_muris_reference/updated.TS.muris.RDS"
# root.to.mastermetadata <- "~/Documents/My Drive/PlaqView_Master/DataProcessing/Summary-Master_Metadata.csv"
# hpca.se <- HumanPrimaryCellAtlasData()
# 
# # singleR requires that it be in a 'singlecellexperiment' format
# # they are workout agnostic
# 
# for_singleR_input <- GetAssayData(Ma_2023)
# pred.Ma_2023 <- SingleR(test = for_singleR_input,
#                             ref = hpca.se,
#                             label = hpca.se$label.main) # reference cell types
# pred.Ma_2023
# # summarize distribution
# table(pred.Ma_2023$labels)
# 
# # to show annotation confidence map
# plotScoreHeatmap(pred.Ma_2023)
# 
# # to show # that are pruned due to low score
# summary(is.na(pred.Ma_2023$pruned.labels))
# 
# ### to place the singleR predictions into Seurat as a sep unit ###
# # seurat.obj[["SingleR.labels"]] <- singler.results$labels
# Ma_2023[["SingleR.labels"]] <- pred.Ma_2023$labels # this nest under metadata
# 
# # Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
# Ma_2023$SingleR.pruned.calls <- pred.Ma_2023$pruned.labels
# Ma_2023$SingleR.calls <- pred.Ma_2023$labels
# 
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], Smooth_muscle_cells = "SMC")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], Endothelial_cells = "EC")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], NK_cell = "NK")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], Chondrocytes = "CH")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], Fibroblasts = "FB")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], Monocyte = "Mono")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], B_cell = "B_Cells")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], Macrophage = "Mø")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], Tissue_stem_cells = "Stem")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], T_cells = "T_Cells")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], 'Pre-B_cell_CD34-' = "PreB_CD34-")
# Ma_2023@meta.data[["SingleR.calls"]] <- recode(Ma_2023@meta.data[["SingleR.calls"]], 'Pro-B_cell_CD34+' = "ProB_CD34+")
# 
# 
# table(Ma_2023@meta.data[["SingleR.calls"]])
# 
# 
# saveRDS(Ma_2023, file = "Ma_2023_RNA_Obj_Pt2_Annotated_Rivanna_singlROnly.rds")
# #### Label Transfer with Livtinokova et al (submit as local job) ####
# TSref <- readRDS(file = "../../PlaqView/PlaqView/data/Litvinukova_2020/Litvinukova_2020.rds")
# 
# #### preprocess references  
# TSref <- NormalizeData(TSref, verbose = T)
# TSref <- FindVariableFeatures(TSref, selection.method = "vst", verbose = T)
# 
# DefaultAssay(Ma_2023) <- 'RNA'
# DefaultAssay(TSref) <- 'RNA'
# 
# 
# anchors <- FindTransferAnchors(reference = TSref, query = Ma_2023, 
#                                dims = 1:30)
# 
# predictions <- TransferData(anchorset = anchors, refdata = TSref$cell_type, 
#                             dims = 1:30)
# 
# Ma_2023 <- AddMetaData(Ma_2023, metadata = predictions)
# 
# #### rename transferred column metadata 
# Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]] <- Ma_2023@meta.data[["predicted.id"]]
# 
# # capitalize the lettering
# Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]] <-str_to_title(Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]], locale = "en")
# 
# # set to active idents
# Idents(Ma_2023) <- Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]]
# 
# Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]], 
#                                                          'Smooth Muscle Cell' = "SMCs")
# Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]], 
#                                                          'Pancreatic Acinar Cell' = "Panc Acinar Cell")
# Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]], 
#                                                          'Fibroblast' = "FB")
# Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]], 
#                                                          'Endothelial Cell' = "EC")
# Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]], 
#                                                          'Macrophage' = "Mø")
# Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]], 
#                                                          'Natural Killer Cell' = "NK")
# Idents(Ma_2023) <- Ma_2023@meta.data[["Seurat_with_Litvinokova_Ref"]]
# 
# saveRDS(Ma_2023, file = "Ma_2023_RNA_Obj_Pt2_Annotated_Rivanna_singleR_Litv.rds")
#### Label Transfer with Chaffin et al (submit as local job) ####
gc()
TSref <- readRDS(file = "../../PlaqView/PlaqView/data/Chaffin_2022/Chaffin_2022.rds")
gc()
#### preprocess references  
TSref <- NormalizeData(TSref, verbose = T)
TSref <- FindVariableFeatures(TSref, selection.method = "vst", verbose = T)

DefaultAssay(Ma_2023) <- 'RNA'
DefaultAssay(TSref) <- 'RNA'


anchors <- FindTransferAnchors(reference = TSref, query = Ma_2023, 
                               dims = 1:30)

predictions <- TransferData(anchorset = anchors, refdata = TSref$Author_Provided, 
                            dims = 1:30)

Ma_2023 <- AddMetaData(Ma_2023, metadata = predictions)

#### rename transferred column metadata 
Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]] <- Ma_2023@meta.data[["predicted.id"]]

# capitalize the lettering
Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]] <-str_to_title(Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]], locale = "en")

# set to active idents
Idents(Ma_2023) <- Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]]

Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]], 
                                                            'Smooth Muscle Cell' = "SMCs")
Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]], 
                                                            'Pancreatic Acinar Cell' = "Panc Acinar Cell")
Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]], 
                                                            'Fibroblast' = "FB")
Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]], 
                                                            'Endothelial Cell' = "EC")
Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]], 
                                                            'Macrophage' = "Mø")
Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]] <- recode(Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]], 
                                                            'Natural Killer Cell' = "NK")
Idents(Ma_2023) <- Ma_2023@meta.data[["Seurat_with_Chaffin_Ref"]]



saveRDS(Ma_2023, file = "Ma_2023_RNA_Obj_Pt2_Annotated_Rivanna_singleR_Litv_Chaffin.rds")

#### Cluster-based Markers Plots ####
## VSMC
pdf(file = "OUTPUT-RNA-Feature_VSMC_ACTA2_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ACTA2", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_VSMC_MYH11_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "MYH11", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_VSMC_PDGFRB_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "PDGFRB", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_VSMC_SEMA3D_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "SEMA3D", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_WT1_ACTA2_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ACTA2", label = T)
dev.off()

## LYMPHOCYTES
pdf(file = "OUTPUT-RNA-Feature_T_CD3D_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CD3D", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_T_PECAM1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CD4", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_T_CD19_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CD19", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_T_CD20_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "MS4A1", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_NK_CD56_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "NCAM1", label = T)
dev.off()

## ENDOTHELIAL
pdf(file = "OUTPUT-RNA-Feature_Endo_PECAM1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "PECAM1", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Endo_VWF_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "VWF", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Endo_GNG11_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "GNG11", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Endo_ADGRF5_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ADGRF5", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Endo_NRP2_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "NRP2", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Endo_CDH5_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CDH5", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Endo_VCAM1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "VCAM1", label = T)
dev.off()

## CARDIAC STEM CELL
pdf(file = "OUTPUT-RNA-Feature_Stem_TBX18_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "TBX18", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Stem_TBX4_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "TBX4", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Stem_WNT5A_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "WNT5A", label = T)
dev.off()

## MACROPHAGES
pdf(file = "OUTPUT-RNA-Feature_MO_CD16_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CD16", label = T)
dev.off()
pdf(file = "OUTPUT-RNA-Feature_MO_CD14_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CD14", label = T)
dev.off()
pdf(file = "OUTPUT-RNA-Feature_MO_CD68_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CD68", label = T)
dev.off()
pdf(file = "OUTPUT-RNA-Feature_MO_CCR5_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CCR5", label = T)
dev.off()

## SMCS 
# TCF21 is not that specific
pdf(file = "OUTPUT-RNA-Feature_SMC_TCF21_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "TCF21", label = T, order = T)
dev.off()


pdf(file = "OUTPUT-RNA-Feature_SMC_FAS_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "FAS", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_SMC_NOX4_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "NOX4", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_SMC_CNN1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CNN1", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_SMC_Desmin_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "DES", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_SMC_MYH11_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "MYH11", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_SMC_PDGFRB_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "PDGFRB", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_SMC_SEMA3D_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "SEMA3D", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_SMC_TBX18_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "TBX18", label = T, order = T)
dev.off()

## cardiomyocyte markers
pdf(file = "OUTPUT-RNA-Feature_CM_ANKRD1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ANKRD1", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_CM_NPPB_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "NPPB", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_CM_TNNT2_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "TNNT2", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_CM_CD31_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CD31", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_CM_MCAM_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "MCAM", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_CM_GNG2_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "GNG2", label = T)
dev.off()

## PERICYTE "SPECIFIC" lol not that specific honestly
pdf(file = "OUTPUT-RNA-Feature_Pericyte_ANGPT_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ANGPT1", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Pericyte_CSPG4_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CSPG4", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Pericyte_ACTG1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ACTG1", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Pericyte_STEAP4_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "STEAP4", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Pericyte_ANGPT2_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ANGPT2", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Pericyte_MYO1B_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "MYO1B", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Pericyte_VEGFA_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "VEGFA", label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_Pericyte_ECM1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ECM1", label = T)
dev.off()

## markers for activated fibroblasts
pdf(file = "OUTPUT-RNA-Feature_FBa_CDH9_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CDH9", label = T, order = T)
dev.off()
pdf(file = "OUTPUT-RNA-Feature_FBa_CD248_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CD248", label = T, order = T)
dev.off()
pdf(file = "OUTPUT-RNA-Feature_FBa_COL1A1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "COL1A1", label = T, order = T)
dev.off()
pdf(file = "OUTPUT-RNA-Feature_FBa_COL1A2_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "COL1A2", label = T, order = T)
dev.off()
pdf(file = "OUTPUT-RNA-Feature_FBa_CCN2_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CCN2", label = T, order = T)
dev.off()
pdf(file = "OUTPUT-RNA-Feature_FBa_FAP_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "FAP", label = T, order = T)
dev.off()
pdf(file = "OUTPUT-RNA-Feature_FBa_S100A4_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "S100A4", label = T, order = T)
dev.off()


pdf(file = "OUTPUT-RNA-Feature_FBa_ACTA2_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ACTA2", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FBa_FAP_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "FAP", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FBa_ENG_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ENG", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FBa_TAGLN_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "TAGLN", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FBa_TNC_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "TNC", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FBa_SPP1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "SPP1", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FBa_POSTN_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "POSTN", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FBa_NOX4_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "NOX4", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FB_DES_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "DES", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FB_VIM_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "VIM", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FB_ACTA2_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "ACTA2", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FB_CDH1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "CDH1", label = T, order = T)
dev.off()

pdf(file = "OUTPUT-RNA-Feature_FB_FN1_by_Seurat_Cluster.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
FeaturePlot(Ma_2023, feature = "FN1", label = T, order = T)
dev.off()
#### Author-Labels LEVEL1 by Clusters ####
Ma_2023$Seurat_Clusters <- Ma_2023$seurat_clusters
Idents(Ma_2023) <- "Seurat_Clusters"

topmarkers <- FindAllMarkers(Ma_2023, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topmarkers %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC) -> top10

pdf(file = "OUTPUT-RNA-Heatmap-Seurat_Clusters-topmarkers.pdf", width = 16, height = 30)
DefaultAssay(Ma_2023) <- "SCT" 
DoHeatmap(Ma_2023, features = top10$gene) 
dev.off()

write_csv(topmarkers, file = "OUTPUT-RNA-Difflist_TopMarkers.csv")

# BASED ON THE DIFF LIST AND ENRICHR CHECKING, I AM TESTING THESE
Ma_2023$Author_Provided_Clusters <- Ma_2023$Seurat_Clusters

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "14" = "Endothelial")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "2" = "Endothelial")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "11" = "Endothelial")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "13" = "Endothelial")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "5" = "Endothelial")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "7" = "Endothelial")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "3" = "Pericyte")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "16" = "Pericyte")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "8" = "Fibroblast")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "18" = "Fibroblast")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "24" = "Fibroblast")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "1" = "Fibroblast")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "4" = "Macrophage")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "21" = "Macrophage")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "20" = "Macrophage")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "12" = "Macrophage")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "0" = "Cardiomyocyte")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "17" = "Cardiomyocyte")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "9" = "Cardiomyocyte")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "6" = "Cardiomyocyte")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "15" = "Cardiomyocyte")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "22" = "Cardiomyocyte")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "10" = "T")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "19" = "VSMC")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "23" = "Unknown")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "25" = "Neuronal")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "26" = "Lymphatic")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "27" = "Adipocyte")

# recode it into level1
Ma_2023@meta.data[["Author_Provided_Clusters_LVL1"]] <- Ma_2023@meta.data[["Author_Provided_Clusters"]] 

Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
Ma_2023$Author_Provided_Clusters_LVL1 <- factor(x = Idents(Ma_2023), levels = sort(levels(Ma_2023)))

pdf(file = "OUTPUT-RNA-Dimplot-Author_Provided_by_Cluster_LVL1.pdf", width = 7, height = 5)
DimPlot(object = Ma_2023, cols = manual_color_list_extended, label = T, group.by = "Author_Provided_Clusters_LVL1", 
        repel = T) +  
  guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
dev.off()


#### Author-Labels LEVEL2 by Clusters ####
Ma_2023$Seurat_Clusters <- Ma_2023$seurat_clusters
Idents(Ma_2023) <- "Seurat_Clusters"

Ma_2023$Author_Provided_Clusters <- Ma_2023$Seurat_Clusters

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "14" = "Endothelial_1")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "2" = "Endothelial_2")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "11" = "Endothelial_3")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "13" = "Endothelial_4")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "5" = "Endothelial_5")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "7" = "Endothelial_6")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "3" = "Pericyte_1")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "16" = "Pericyte_2")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "8" = "FB_1")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "18" = "FB_2")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "24" = "FB_3")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "1" = "FB_4")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "4" = "Mø_1")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "21" = "Mø_2")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "20" = "Mø_3")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "12" = "Mø_4")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "0" = "CM_1")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "17" = "CM_2")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "9" = "CM_3")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "6" = "CM_4")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "15" = "CM_5")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "22" = "CM_6")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "10" = "T")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "19" = "VSMC")

Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "23" = "Unknown")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "25" = "Neuronal")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "26" = "Lymphatic")
Ma_2023@meta.data[["Author_Provided_Clusters"]] <- recode(Ma_2023@meta.data[["Author_Provided_Clusters"]], "27" = "Adipocyte")

# recode it into level1
Ma_2023@meta.data[["Author_Provided_Clusters_LVL2"]] <- Ma_2023@meta.data[["Author_Provided_Clusters"]] 

Idents(Ma_2023) <- "Author_Provided_Clusters_LVL2"
Ma_2023$Author_Provided_Clusters_LVL2 <- factor(x = Idents(Ma_2023), levels = sort(levels(Ma_2023)))

pdf(file = "OUTPUT-RNA-Dimplot-Author_Provided_by_Cluster_LVL2.pdf", width = 9, height = 5)
DimPlot(object = Ma_2023, cols = manual_color_list_extended, 
        label = T, group.by = "Author_Provided_Clusters_LVL2", repel = T) +  
  guides(color = guide_legend(override.aes = list(size=4), ncol=2) )
dev.off()



#### What the heck is going on in cluster 23??? ####
CellPop_Seurat_Cluster_by_Donor <- table(Ma_2023$Donor, Ma_2023$Seurat_Clusters)


#### Output Annotated Object ####

saveRDS(Ma_2023, file = "Ma_2023_RNA_Obj_Pt2_Annotated_Rivanna.rds")

#### Cluster-based markers ####
Idents(object = Ma_2023) <- "seurat_clusters"
difflist <- Seurat::FindAllMarkers(Ma_2023)
write_csv(difflist, file = "OUTPUT-RNA-Difflist_by_Seurat_Clusters.csv")

Idents(object = Ma_2023) <- "Status"
difflist <- Seurat::FindAllMarkers(Ma_2023)
write_csv(difflist, file = "OUTPUT-RNA-Difflist_by_Status.csv")

Idents(object = Ma_2023) <- "Seurat_with_Chaffin_Ref"
difflist <- Seurat::FindAllMarkers(Ma_2023)
write_csv(difflist, file = "OUTPUT-RNA-Difflist_by_Seurat_Chaffin.csv")

Idents(object = Ma_2023) <- "Author_Provided_Cells"
difflist <- Seurat::FindAllMarkers(Ma_2023)
write_csv(difflist, file = "OUTPUT-RNA-Difflist_by_Author_Provided_Cells.csv")

Idents(object = Ma_2023) <- "Author_Provided_Clusters_LVL1"
difflist <- Seurat::FindAllMarkers(Ma_2023)
write_csv(difflist, file = "OUTPUT-RNA-Difflist_by_Author_Provided_Clusters_LVL1.csv")


Idents(object = Ma_2023) <- "Author_Provided_Clusters_LVL2"
difflist <- Seurat::FindAllMarkers(Ma_2023)
write_csv(difflist, file = "OUTPUT-RNA-Difflist_by_Author_Provided_Clusters_LVL2.csv")


#### Basic Plots ####
pdf(file = "OUTPUT-RNA-Dimplot-Seurat_Clusters-labels.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Seurat_Clusters"
DimPlot(object = Ma_2023, cols = manual_color_list_extended, label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Dimplot-Status-labels.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Status"
DimPlot(object = Ma_2023, cols = manual_color_list_extended, label = T)
dev.off()

pdf(file = "OUTPUT-RNA-Dimplot-Seurat_Clusters-labels-SPLIT-by-status.pdf", width = 10, height = 5)
Idents(Ma_2023) <- "seurat_clusters"
DimPlot(object = Ma_2023, cols = manual_color_list_extended, label = T, split.by = "Status")
dev.off()

pdf(file = "OUTPUT-RNA-Dimplot-Author_Provided-labels-SPLIT-by-status_LVL1.pdf", width = 10, height = 5)
Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
DimPlot(object = Ma_2023, cols = manual_color_list_extended, label = T, split.by = "Status") +
  ggtitle("Author_Provided_Clusters_LVL1_by_Status") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


pdf(file = "OUTPUT-RNA-Dimplot-Author_Provided-labels-SPLIT-by-status_LVL2.pdf", width = 10, height = 5)
Idents(Ma_2023) <- "Author_Provided_Clusters_LVL2"
DimPlot(object = Ma_2023, cols = manual_color_list_extended, label = T, split.by = "Status") +
  ggtitle("Author_Provided_Clusters_LVL2_by_Status")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


pdf(file = "OUTPUT-RNA-Dimplot-Donors-labels.pdf", width = 7, height = 5)
Idents(Ma_2023) <- "Donor"
DimPlot(object = Ma_2023, cols = manual_color_list, label = F) +
  ggtitle("Cellular_Contribution_by_Donors")+
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

pdf(file = "OUTPUT-RNA-Dimplot-Chaffin-labels.pdf", width = 9, height = 5)
DimPlot(object = Ma_2023, cols = manual_color_list, label = T, group.by = "Seurat_with_Chaffin_Ref", repel = T)
dev.off()

pdf(file = "OUTPUT-RNA-Dimplot-SingleR_calls-labels.pdf", width = 9, height = 5)
DimPlot(object = Ma_2023, cols = manual_color_list_extended, label = T, group.by = "SingleR.calls", repel = T)
dev.off()

pdf(file = "OUTPUT-RNA-Dimplot-Batch-labels.pdf", width = 9, height = 5)
DimPlot(object = Ma_2023, cols = manual_color_list_extended, label = T, group.by = "Author_Provided_Clusters_LVL1", repel = T,
        split.by = "Batch")
dev.off()

#### Additional Genes to Plot after Finalization of Annotation ####
pdf(file = "OUTPUT-RNA-Feature_SMC_TCF21_by_Author_Provided.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
FeaturePlot(Ma_2023, feature = "TCF21", label = T, order = T) 
dev.off()

pdf(file = "OUTPUT-RNA-Feature_SMC_NOX4_by_Author_Provided.pdf", width = 5, height = 5)
Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
FeaturePlot(Ma_2023, feature = "NOX4", label = T, order = T) 
dev.off()
