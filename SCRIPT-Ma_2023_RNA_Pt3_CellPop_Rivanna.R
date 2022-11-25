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


#### READ DATA ####
Ma_2023 <- readRDS(file = "Ma_2023_RNA_Obj_Pt2_Annotated_Rivanna.rds")

#### Cell Count Tables ####
# counts by donor
CellPop_Donors <- table(Ma_2023$Donor)
write.csv(CellPop_Donors, file = "OUTPUT-RNA-CellPop-by-Donors.csv")

# counts by disease state
CellPop_Status <- table(Ma_2023$Status)
write.csv(CellPop_Status, file = "OUTPUT-RNA-CellPop-by-Status.csv")

# counts by author per cluster
CellPop_Author_Provided_Cluster <- table(Ma_2023$Author_Provided_Clusters_LVL1)
write.csv(CellPop_Author_Provided_Cluster, file = "OUTPUT-RNA-CellPop-by-Author_Provided_Clusters_LVL1.csv")

# counts by author per cluster2
CellPop_Author_Provided_Cluster <- table(Ma_2023$Author_Provided_Clusters_LVL2)
write.csv(CellPop_Author_Provided_Cluster, file = "OUTPUT-RNA-CellPop-by-Author_Provided_Clusters_LVL2.csv")


# counts per cluster by donors
CellPop_Seurat_Cluster_by_Donor <- table(Ma_2023$Donor, Ma_2023$Seurat_Clusters)
write.csv(CellPop_Seurat_Cluster_by_Donor, file = "OUTPUT-RNA-CellPop-by-Seurat_Clusters-by-Donor.csv")

CellPop_Seurat_Cluster_PROP <- prop.table(CellPop_Seurat_Cluster_by_Donor, 2)
write.csv(CellPop_Seurat_Cluster_PROP, file = "OUTPUT-RNA-CellPop-by-Seurat_Clusters-by-Donor-Proportions.csv")

# counts per cluster by status
CellPop_Seurat_Cluster_by_Status <- table(Ma_2023$Status, Ma_2023$Seurat_Clusters)
write.csv(CellPop_Seurat_Cluster_by_Status, file = "OUTPUT-RNA-CellPop-by-Seurat_Clusters-by-Status.csv")

CellPop_Seurat_Cluster_by_Status_PROP <- prop.table(CellPop_Seurat_Cluster_by_Status, 2)
write.csv(CellPop_Seurat_Cluster_by_Status_PROP, file = "OUTPUT-RNA-CellPop-by-Seurat_Clusters-by-Status-Proportions.csv")

# counts per cluster by dx
CellPop_Author_Provided_Cluster_by_Status <- table(Ma_2023$Status, Ma_2023$Author_Provided_Clusters_LVL1)
write.csv(CellPop_Author_Provided_Cluster_by_Status, file = "OUTPUT-RNA-CellPop-by-Author_Provided_Clusters-by-Status.csv")

CellPop_Author_Provided_Cluster_by_Status_Prop <-  prop.table(CellPop_Author_Provided_Cluster_by_Status, 2)
write.csv(CellPop_Author_Provided_Cluster_by_Status_Prop, file = "OUTPUT-RNA-CellPop-by-Author_Provided_Clusters-by-Status-Proportions.csv.csv")

CellPop_Author_Provided_Cluster_by_Donor <- table(Ma_2023$Donor, Ma_2023$Author_Provided_Clusters_LVL1)
write.csv(CellPop_Author_Provided_Cluster_by_Status, file = "OUTPUT-RNA-CellPop-by-Author_Provided_Clusters-by-Donor.csv")


#### Population Comparison Charts ####


pdf(file = "OUTPUT-RNA-PopChart-Author_Provided_by_Cluster_BY_STATUS_LVL1.pdf", width = 5, height = 6)
CellPop_Author_Provided_Cluster_by_Status <- table(Ma_2023$Status, Ma_2023$Author_Provided_Clusters_LVL1)
popchart <- as.data.frame(CellPop_Author_Provided_Cluster_by_Status)
popchart[order(popchart$Var2),]
popchart$Var2 <- as_factor((popchart$Var2))
ggplot(popchart, aes(x = Var2, fill = Var1,
                     y = ifelse(test = Var1 == "Chemo",
                                yes = -Freq, no = Freq))) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  guides(fill=guide_legend(title="Status")) +
  scale_colour_manual(values = c("#f58142", "#407dde"),
                      aesthetics = c("colour", "fill")) +
  xlab("Cell Counts") +
  theme_bw()

dev.off()

pdf(file = "OUTPUT-RNA-PopChart-Author_Provided_by_Cluster_BY_DONOR_LVL1.pdf", width = 5, height = 6)
CellPop_Author_Provided_Cluster_by_Status <- table(Ma_2023$Donor, Ma_2023$Author_Provided_Clusters_LVL1)
popchart <- as.data.frame(CellPop_Author_Provided_Cluster_by_Status)
popchart[order(popchart$Var2),]
popchart$Var2 <- as_factor((popchart$Var2))
ggplot(popchart, aes(x = Var2, fill = Var1,
                     y = ifelse(test = Var1 == "Chemo",
                                yes = -Freq, no = Freq))) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  guides(fill=guide_legend(title="Status")) +
  scale_colour_manual(values = manual_color_list,
                      aesthetics = c("colour", "fill")) +
  xlab("Cell Counts") +
  theme_bw()

dev.off()



pdf(file = "OUTPUT-RNA-PopChart_PROP-Author_Provided_by_Cluster_BY_STATUS_LVL1.pdf", width = 5, height = 6)
CellPop_Author_Provided_Cluster_by_Status_Prop <-  prop.table(CellPop_Author_Provided_Cluster_by_Status, 2)
popchart <- as.data.frame(CellPop_Author_Provided_Cluster_by_Status_Prop)
popchart[order(popchart$Var2),]
popchart$Var2 <- as_factor(popchart$Var2)
popchart %>%
  ggplot(aes(x = Var2, fill = Var1,
           y = ifelse(test = Var1 == "Chemo",
                      yes = -Freq, no = Freq))) + # changes chemo to the negative axis
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Cellular Contribution of Each Cluster", x = "Cell Type", y = "Proportion of Population") +
  scale_colour_manual(values = manual_color_list,
                      aesthetics = c("colour", "fill")) +
  guides(fill=guide_legend(title="Status")) +
  theme_bw()
dev.off()

pdf(file = "OUTPUT-RNA-PopChart_PROP-Author_Provided_by_Cluster_BY_STATUS_LVL2.pdf", width = 5, height = 6)
CellPop_Author_Provided_Cluster_by_Status <- table(Ma_2023$Status, Ma_2023$Author_Provided_Clusters_LVL2)
CellPop_Author_Provided_Cluster_by_Status_Prop <-  prop.table(CellPop_Author_Provided_Cluster_by_Status, 2)
popchart <- as.data.frame(CellPop_Author_Provided_Cluster_by_Status_Prop)
popchart[order(popchart$Var2),]
popchart$Var2 <- as_factor(popchart$Var2)
popchart %>%
  ggplot(aes(x = Var2, fill = Var1,
             y = ifelse(test = Var1 == "Chemo",
                        yes = -Freq, no = Freq))) + # changes chemo to the negative axis
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Cellular Contribution of Each Cluster", x = "Cell Type", y = "Proportion of Population") +
  scale_colour_manual(values = c("#f58142", "#407dde"),
                      aesthetics = c("colour", "fill")) +
  guides(fill=guide_legend(title="Status")) +
  theme_bw()
dev.off()

pdf(file = "OUTPUT-RNA-PopChart_PROP-Cluster-By-Status.pdf", width = 5, height = 6)
popchart <- as.data.frame(CellPop_Seurat_Cluster_by_Status_PROP)
popchart$Var2 <- as_factor(popchart$Var2)
popchart %>%
  mutate(Var2 = fct_reorder(Var2, Freq)) %>% 
  ggplot(aes(x = Var2, fill = Var1,
             y = ifelse(test = Var1 == "Chemo",
                        yes = -Freq, no = Freq))) + # changes chemo to the negative axis
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Cellular Contribution of Each Cluster", x = "Cell Cluster", y = "Proportion of Population") +
  scale_colour_manual(values = c("#f58142", "#407dde"),
                      aesthetics = c("colour", "fill")) +
  guides(fill=guide_legend(title="Status")) +
  theme_bw()
dev.off()

pdf(file = "OUTPUT-RNA-PopChart-Total-Status.pdf", width = 4, height = 6)
popchart <- as.data.frame(CellPop_Status)
ggplot(popchart) + 
  geom_bar(aes(x = Var1, y = Freq, fill = Var1), 
    stat='identity', position = 'dodge') +
  geom_text(
    aes(x = Var1, y = Freq, label = Freq),
    vjust = -0.5, size = 4,
    position = position_dodge(width = 1),
    inherit.aes = TRUE
  ) +
  labs(title = "Cellular Contribution from Each Group", x = "Disease State", y = "Number of Cells") +
  theme_bw() +
  scale_colour_manual(values = c("#f58142", "#407dde"),
                      aesthetics = c("colour", "fill")) +
  theme(legend.position="none") 
dev.off()

pdf(file = "OUTPUT-RNA-PopChart-Donor-by-Status.pdf", width = 4, height = 6)
popchart <- as.data.frame(CellPop_Donors)
ggplot(popchart) + 
  geom_bar(aes(x = Var1, y = Freq, fill = Var1), 
           stat='identity', position = 'dodge') +
  geom_text(
    aes(x = Var1, y = Freq, label = Freq),
    vjust = -0.5, size = 4,
    position = position_dodge(width = 1),
    inherit.aes = TRUE
  ) +
  labs(title = "Cellular Contribution from Each Donor", x = "Donor", y = "Number of Cells") +
  theme_bw() +
  scale_colour_manual(values = manual_color_list,
                      aesthetics = c("colour", "fill")) +
  theme(legend.position="none") 
dev.off()

pdf(file = "OUTPUT-RNA-PopChart-Author_Provided_LVL1_Counts.pdf", width = 5, height = 6)
CellPop_Author_Provided_Cluster <- table(Ma_2023$Author_Provided_Clusters_LVL1)
popchart <- as.data.frame(CellPop_Author_Provided_Cluster)
ggplot(popchart) + 
  geom_bar(aes(x = Var1, y = Freq, fill = Var1), 
           stat='identity', position = 'dodge') +
  geom_text(
    aes(x = Var1, y = Freq, label = Freq),
    vjust = -0.5, size = 4,
    position = position_dodge(width = 1),
    inherit.aes = TRUE
  ) +
  labs(title = "Cell Counts for Each Cell Type", x = "Cell Type", y = "Number of Cells") +
  theme_bw() +
  scale_colour_manual(values = manual_color_list_extended,
                      aesthetics = c("colour", "fill")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="none") 
dev.off()

pdf(file = "OUTPUT-RNA-PopChart-Author_Provided_LVL2_Counts.pdf", width = 10, height = 6)
CellPop_Author_Provided_Cluster <- table(Ma_2023$Author_Provided_Clusters_LVL2)
popchart <- as.data.frame(CellPop_Author_Provided_Cluster)
ggplot(popchart) + 
  geom_bar(aes(x = Var1, y = Freq, fill = Var1), 
           stat='identity', position = 'dodge') +
  geom_text(
    aes(x = Var1, y = Freq, label = Freq),
    vjust = -0.5, size = 4,
    position = position_dodge(width = 1),
    inherit.aes = TRUE
  ) +
  labs(title = "Cell Counts for Each Cell Type", x = "Cell Type", y = "Number of Cells") +
  theme_bw() +
  scale_colour_manual(values = manual_color_list_extended,
                      aesthetics = c("colour", "fill")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="none") 
dev.off()

# cell pop by donor by type
pdf(file = "OUTPUT-RNA-PopChart-Author_Provided-by-Donor.pdf", width = 4, height = 6)
popchart <- as.data.frame((CellPop_Author_Provided_Cluster_by_Donor))
ggplot(popchart, aes(fill=Var1, y=Freq, x=Var2)) + 
  scale_colour_manual(values = manual_color_list,
                      aesthetics = c("colour", "fill")) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Cell Type (Lvl1)") +
  ylab("Number of Cells") +
  guides(fill=guide_legend(title="Donor")) 
  
dev.off()

#### Finding Conserved Cell Markers Across Status ####
Idents(Ma_2023) <- "Author_Provided_Clusters_LVL2"
write.csv(FindConservedMarkers(Ma_2023, ident.1 = "CM_1", grouping.var = "Status", verbose = FALSE),
          file = "OUTPUT-RNA-Difflist_CM1_Conserved_Markers_ignore_status.csv")
write.csv(FindConservedMarkers(Ma_2023, ident.1 = "CM_2", grouping.var = "Status", verbose = FALSE),
          file = "OUTPUT-RNA-Difflist_CM2_Conserved_Markers_ignore_status.csv")
write.csv(FindConservedMarkers(Ma_2023, ident.1 = "CM_3", grouping.var = "Status", verbose = FALSE),
          file = "OUTPUT-RNA-Difflist_CM3_Conserved_Markers_ignore_status.csv")
write.csv(FindConservedMarkers(Ma_2023, ident.1 = "CM_4", grouping.var = "Status", verbose = FALSE),
          file = "OUTPUT-RNA-Difflist_CM4_Conserved_Markers_ignore_status.csv")
write.csv(FindConservedMarkers(Ma_2023, ident.1 = "CM_5", grouping.var = "Status", verbose = FALSE),
          file = "OUTPUT-RNA-Difflist_CM5_Conserved_Markers_ignore_status.csv")
write.csv(FindConservedMarkers(Ma_2023, ident.1 = "CM_6", grouping.var = "Status", verbose = FALSE),
          file = "OUTPUT-RNA-Difflist_CM6_Conserved_Markers_ignore_status.csv")

#### Compare CM Expressions ####
write.csv(FindMarkers(Ma_2023, ident.1 = "CM_1", ident.2 = c("CM_2", "CM_3", "CM_4", "CM_5"), min.pct = 0.5),
          file = "OUTPUT-RNA-Difflist_CM1-vs-other-CMs.csv")
write.csv(FindMarkers(Ma_2023, ident.1 = "CM_2", ident.2 = c("CM_1", "CM_3", "CM_4", "CM_5"), min.pct = 0.5),
          file = "OUTPUT-RNA-Difflist_CM2-vs-other-CMs.csv")
write.csv(FindMarkers(Ma_2023, ident.1 = "CM_3", ident.2 = c("CM_2", "CM_1", "CM_4", "CM_5"), min.pct = 0.5),
          file = "OUTPUT-RNA-Difflist_CM3-vs-other-CMs.csv")
write.csv(FindMarkers(Ma_2023, ident.1 = "CM_4", ident.2 = c("CM_2", "CM_3", "CM_1", "CM_5"), min.pct = 0.5),
          file = "OUTPUT-RNA-Difflist_CM4-vs-other-CMs.csv")
write.csv(FindMarkers(Ma_2023, ident.1 = "CM_5", ident.2 = c("CM_2", "CM_3", "CM_4", "CM_1"), min.pct = 0.5),
          file = "OUTPUT-RNA-Difflist_CM5-vs-other-CMs.csv")
write.csv(FindMarkers(Ma_2023, ident.1 = "CM_6", ident.2 = c("CM_2", "CM_3", "CM_4", "CM_1"), min.pct = 0.5),
          file = "OUTPUT-RNA-Difflist_CM6-vs-other-CMs.csv")

Idents(Ma_2023) <- "Author_Provided_Clusters_LVL1"
write.csv(FindMarkers(Ma_2023, ident.1 = "CM", grouping.var = "Status", min.pct = 0.5),
          file = "OUTPUT-RNA-Difflist_Conserved_Genes_in_CM_by_Status.csv")


#### Compare FB Expression




#### Compare MO Expression ####
write.csv(FindMarkers(Ma_2023, ident.1 = "Macrophage", grouping.var = "Status", min.pct = 0.5),
          file = "OUTPUT-RNA-Difflist_Conserved_Genes_in_MO_by_Status.csv")