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

#### MAKE CELLCHAT OBJ ####
# makes the actual obj
cellchat <- createCellChat(object = Ma_2023, # this is my seurat obj
                           group.by = "Author_Provided_Clusters_LVL1")

#### PREPROCESS CELLCHAT OBJ ####
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

#### COMPUTE CELL-CELL NETWORK ####
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# compute pathways involved
cellchat <- computeCommunProbPathway(cellchat)

# compute aggregate networks 
cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# compute functional 
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional") # may need to run reticulate::py_install(packages = 'umap-learn')

#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space


####-#-#-# PRELIM ANALYSIS FOR ENTIRE DATASET  ####
### print out signal matrices
df.net <- subsetCommunication(cellchat) #ligand/receptor
df.path <- subsetCommunication(cellchat, slot.name = "netP") #pathway level

write.csv(df.net, file = "OUTPUT-RNA-CellChat_all-enriched-comm-by-ligand-by-Author_Provided_Clusters_LVL1.csv")
write.csv(df.path, file = "OUTPUT-RNA-CellChat_all-enriched-comm-by-pathway-by-Author_Provided_Clusters_LVL1.csv")

### dominant senders and receivers
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)

pdf(file = "strength.pdf", width = 5, height = 5)
gg1
dev.off()

### circos plot- aggregate networks
pdf(file = "OUTPUT-RNA-CellChat_Circos_Author_Provided_Clusters_LVL1.pdf", width = 4, height = 4) # this says to make a pdf
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off() # this say stop making pdf


pdf(file = "OUTPUT-RNA-CellChat_circos-strength of interactions.pdf", width = 4, height = 4) # this says to make a pdf
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off() # this say stop making pdf

### plot all outgoings 
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
## here you need to choose which one may be informative


### signal of certain genes
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)

### signal of certain pathways
pathways.show <- c("CXCL") # pathway

netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# create a new metadata column to store EGFP status
# cant used merged object... easier to use original
full <- readRDS(file = "merge.rds") # this merged.rds file is provided by caitlin, not to be confused with my own merging

Idents(full) <- full$manual_clusters
# take out the ones we dont wanna work with 

full <- subset(full, idents = c("Macrophage", "Ki67+ Mac", "Tumor", "S100a9+ Mac", "DC", "SMC"))

# the next two lines will label egfp + and neg and put them into 'active.ident'
Idents(full, WhichCells(object = full, expression = eGFP > 0., slot = 'data')) <- 'e+'
Idents(full, WhichCells(object = full, expression = eGFP <= 0, slot = 'data')) <- 'e-'

### create a new identity group that tells you the cell type and GFP 
full$eGFP_detail <- str_c(full$manual_clusters, full@active.ident, sep = " ")

# next two lines reformats macrphage for consistency
full@meta.data[["eGFP_detail"]] <- recode(full@meta.data[["eGFP_detail"]], 
                                          'Macrophage e-' = "Mac e-")
full@meta.data[["eGFP_detail"]] <- recode(full@meta.data[["eGFP_detail"]], 
                                          'Macrophage e+' = "Mac e+")



unique(full$eGFP_detail)
table(full$eGFP_detail)

DimPlot(full)
DimPlot(full, group.by = "eGFP_detail", cols = manual_color_list, label = T, repel = T) # if you get error like viewpoint has zero dim make ur window bigger


#### Compare Chemo vs. Normal #### 
Idents(Ma_2023) <- "Status"
chemo <- subset(x = Ma_2023, idents = "Chemo")
normal <- subset(x = Ma_2023, idents = "Normal")

chemo.cellchat <- createCellChat(object = chemo, # this is my seurat obj
                           group.by = "Author_Provided_Clusters_LVL1")
normal.cellchat <- createCellChat(object = normal, # this is my seurat obj
                           group.by = "Author_Provided_Clusters_LVL1")

# process normal
normal.cellchat@DB <- CellChatDB
normal.cellchat <- subsetData(normal.cellchat) # This step is necessary even if using the whole database
normal.cellchat <- identifyOverExpressedGenes(normal.cellchat)
normal.cellchat <- identifyOverExpressedInteractions(normal.cellchat)
normal.cellchat <- projectData(normal.cellchat, PPI.human)
normal.cellchat <- computeCommunProb(normal.cellchat)
normal.cellchat <- filterCommunication(normal.cellchat, min.cells = 10)
normal.cellchat <- computeCommunProbPathway(normal.cellchat)
normal.cellchat <- aggregateNet(normal.cellchat)
normal.cellchat <- netAnalysis_computeCentrality(normal.cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
normal.cellchat <- computeNetSimilarity(normal.cellchat, type = "functional")
normal.cellchat <- netEmbedding(normal.cellchat, type = "functional") # may need to run reticulate::py_install(packages = 'umap-learn')
normal.cellchat <- netClustering(normal.cellchat, type = "functional")

# process chemo
chemo.cellchat@DB <- CellChatDB
chemo.cellchat <- subsetData(chemo.cellchat) # This step is necessary even if using the whole database
chemo.cellchat <- identifyOverExpressedGenes(chemo.cellchat)
chemo.cellchat <- identifyOverExpressedInteractions(chemo.cellchat)
chemo.cellchat <- projectData(chemo.cellchat, PPI.human)
chemo.cellchat <- computeCommunProb(chemo.cellchat)
chemo.cellchat <- filterCommunication(chemo.cellchat, min.cells = 10)
chemo.cellchat <- computeCommunProbPathway(chemo.cellchat)
chemo.cellchat <- aggregateNet(chemo.cellchat)
chemo.cellchat <- netAnalysis_computeCentrality(chemo.cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
chemo.cellchat <- computeNetSimilarity(chemo.cellchat, type = "functional")
chemo.cellchat <- netEmbedding(chemo.cellchat, type = "functional") # may need to run reticulate::py_install(packages = 'umap-learn')
chemo.cellchat <- netClustering(chemo.cellchat, type = "functional")

#### combine the processed sub cellchat groups...
object.list <- list(Normal = normal.cellchat, Chemo = chemo.cellchat) 

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

pdf(file = "OUTPUT-RNA-CellChat_Comparative_Total_Int.pdf", width = 5, height = 5)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

pdf(file = "OUTPUT-RNA-CellChat_Comparative_LVL1_Circos_Weighted_by_CellCounts.pdf", width = 15, height = 8)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

pdf(file = "OUTPUT-RNA-CellChat_Comparative_LVL1_2D_Weighted_by_CellCounts.pdf", width = 15, height = 8)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
dev.off()