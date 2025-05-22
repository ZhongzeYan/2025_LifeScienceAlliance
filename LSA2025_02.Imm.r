library(Seurat)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggplot2)

setwd("result_LSA")
scRNA.combined_datas5_addcelltype_filted <- readRDS( 
    "result_LSA/scRNA.combined_datas5_addcelltype_filted.RDS")

dir.create("result_LSA/02.cell_identify/02.immune_cells")

DefaultAssay(scRNA.combined_datas5_addcelltype_filted) <- "RNA"
Imm_cells <- subset(scRNA.combined_datas5_addcelltype_filted,
                    CellType=="Immune_cells" | 
                    CellType == "Astrocytes" | 
                    CellType == "Astrocytes2*" | 
                    CellType == "Microglia_Macrophages")
Imm.list <- SplitObject(Imm_cells, split.by="orig.ident")
Imm.list <- lapply(X = Imm.list, FUN = function(x){
     x <- NormalizeData(x)
     x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Imm.list)
Imm.list.anchors <- FindIntegrationAnchors(object.list = Imm.list, 
                                           anchor.features = features)
Imm.combined_datas5 <- IntegrateData(anchorset = Imm.list.anchors)
Imm.combined_datas5 <- ScaleData(Imm.combined_datas5)
Imm.combined_datas5 <- RunPCA(Imm.combined_datas5, 
                              features = VariableFeatures(Imm.combined_datas5), 
                              npcs = 50)
Imm.combined_datas5 <- RunUMAP(Imm.combined_datas5, 
                               reduction = "pca", dims = 1:30)
Imm.combined_datas5 <- FindNeighbors(Imm.combined_datas5, 
                                     reduction = "pca", dims = 1:50)
Imm.combined_datas5 <- FindClusters(Imm.combined_datas5, resolution = 0.9)

#pdf("Choose_Neighbor_dims.pdf")
#DimHeatmap(Imm.combined_datas5, dims = 1:50, genes = 20, cells = 500, balanced = TRUE)
#dev.off()

setwd("result_LSA/02.cell_identify/02.immune_cells/")

pdf("00.umap_clusters_res.0.5_Imm.pdf", width = 16, height = 15)
DimPlot(Imm.combined_datas5, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster=FALSE)
dev.off()

pdf("00.umap_clusters_res.0.5_Imm_raster.pdf", width = 16, height = 15)
DimPlot(Imm.combined_datas5, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster=TRUE)
dev.off()


pdf("01.umap_Imm_combined_clusters-split_raster.pdf", width = 40, height = 30)
DimPlot(Imm.combined_datas5, reduction = "umap", split.by = "seurat_clusters", ncol=4, label = TRUE, raster=TRUE)
dev.off()

pdf("01.umap_Imm_combined_origin.pdf", width = 10)
DimPlot(Imm.combined_datas5, reduction = "umap", group.by = "orig.ident", raster=TRUE)
dev.off()

pdf("01.umap_Imm_combined_origin_split.pdf", width = 28, height = 15)
DimPlot(Imm.combined_datas5, reduction = "umap", split.by = "orig.ident", ncol=3, raster=TRUE)
dev.off()

pdf("01.umap_Imm_combined_ref.celltype.pdf", width = 10)
DimPlot(Imm.combined_datas5, reduction = "umap", group.by = "ref_celltype", label=TRUE, raster=TRUE)
dev.off()

pdf("01.umap_Imm_ref.Celltype_2019_AltOlig.pdf", width = 13.5,height=7.5)
DimPlot(subset(Imm.combined_datas5,subset=orig.ident=="2019_AltOligHtg"), reduction = "umap", group.by = "ref_celltype", label=TRUE, raster=TRUE)
dev.off()

pdf("01.umap_Imm_ref.Celltype_split_2019_AltOlig.pdf", width = 21.5, height=7.5)
DimPlot(subset(Imm.combined_datas5,subset=orig.ident=="2019_AltOligHtg"), reduction = "umap", split.by = "ref_celltype", label=TRUE, raster=TRUE, ncol=6)
dev.off()

pdf("01.umap_Imm_ref.Celltype_2018_Lake.pdf", width = 13.5,height=7.5)
DimPlot(subset(Imm.combined_datas5,subset=orig.ident=="2018_Lake"), reduction = "umap", group.by = "ref_celltype",
        label=TRUE, label.size = 10, raster=TRUE)
dev.off()

##  Check if ImOLGs in immune cells
dir.create("02.subtypes")
oligodendrocyte_genes <- c("PLP1","CARNS1", "MBP", "CTNNA3", "PDGFRA", "TMEM144", "OLIG2","ANLN","CNP","APC", "QDPR", "SLAIN1", "MAG", "MOBP")
DefaultAssay(Imm.combined_datas5) <- "RNA"
pdf("02.subtypes/01.oligo_like.selectgenes_VlnPlot.pdf", height = 14, width = 24)
VlnPlot(Imm.combined_datas5, assay="RNA", features = oligodendrocyte_genes, pt.size=0, ncol=2)
dev.off()

pdf("02.subtypes/01.oligo_like.selectgenes_FeaturePlot_lite.pdf", height = 52.5, width = 20)
FeaturePlot(Imm.combined_datas5, features = oligodendrocyte_genes, reduction = "umap", pt.size = 0.2,label=T, ncol=2, raster=TRUE)
dev.off()

ImOLG_genes <-c("APOE", "CD74", "ITPR2")
pdf("02.subtypes/02.ImOLG_marker.selectgenes_VlnPlot.pdf", height = 4, width = 24)
VlnPlot(Imm.combined_datas5, features = ImOLG_genes, pt.size=0, ncol=2)
dev.off()

pdf("02.subtypes/01.ImOLG_like.selectgenes_FeaturePlot_lite.pdf", height = 15, width = 20)
FeaturePlot(Imm.combined_datas5, features = ImOLG_genes,
            reduction = "umap", pt.size = 0.2,label=T, ncol=2, raster=TRUE)
dev.off()

ImOLGs_candidate <- subset(Imm.combined_datas5, seurat_clusters=="8"| 
                           seurat_clusters=="12"| seurat_clusters=="23")

saveRDS(ImOLGs_candidate,"ImOLGs_candidate.Rds")


## Immune cells after filter out oligo
dir.create("../02.immune_cells_filtOLG")
setwd("../02.immune_cells_filtOLG")

Imm.combined_filtOLG <- subset(Imm.combined_datas5,seurat_clusters!="8" & 
                               seurat_clusters!="12" & seurat_clusters!="23")

DefaultAssay(Imm.combined_filtOLG) <- "RNA"
Imm.list <- SplitObject(Imm.combined_filtOLG, split.by="orig.ident")

Imm.list <- lapply(X = Imm.list, FUN = function(x){
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Imm.list)

Imm.list.anchors <- FindIntegrationAnchors(object.list = Imm.list, 
                                           anchor.features = features)

Imm.combined_filtOLG <- IntegrateData(anchorset = Imm.list.anchors)

### Perform integratave analysis
Imm.combined_filtOLG <- ScaleData(Imm.combined_filtOLG)

Imm.combined_filtOLG <- RunPCA(Imm.combined_filtOLG,
                               features = VariableFeatures(Imm.combined_filtOLG), 
                               npcs = 50)
Imm.combined_filtOLG <- RunUMAP(Imm.combined_filtOLG, 
                                reduction = "pca", dims = 1:30)
Imm.combined_filtOLG <- FindNeighbors(Imm.combined_filtOLG, 
                                      reduction = "pca", dims = 1:30)
Imm.combined_filtOLG <- FindClusters(Imm.combined_filtOLG, resolution = 0.5)
Imm.combined_filtOLG <- RunTSNE(Imm.combined_filtOLG, 
                                reduction = "pca", dims = 1:30)

pdf("Imm.filtOLG_Choose_Neighbor_dims.pdf",height=25,width=15)
DimHeatmap(Imm.combined_filtOLG, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

obj <- FindClusters(Imm.combined_filtOLG, resolution = seq(0.8,2.2,by=0.1))
pdf("Imm.filtOLG_Choose_cluster_resolution.pdf",width=25,height=15)
clustree(obj)
dev.off()

for(i in c(30,32,35,40,42,45,46,49,50)){
    for(j in c(0.8,0.9,1,1.2,1.5,2,2.1,2.2)){
        Imm.combined_filtOLG.new <- FindNeighbors(object = Imm.combined_filtOLG, 
                                                  dims = 1:i, verbose = T)
        Imm.combined_filtOLG.new <- FindClusters(object = Imm.combined_filtOLG.new, 
                                                 resolution = j, verbose = T)
        Imm.combined_filtOLG.new <- RunUMAP(object = Imm.combined_filtOLG.new, 
                                            dims = 1:i, verbose = T)
        Imm.combined_filtOLG.new <- RunTSNE(object = Imm.combined_filtOLG.new, 
                                            dims = 1:i, verbose = T, 
                                            check_duplicates = FALSE)
        gd1 <- DimPlot(object = Imm.combined_filtOLG.new, reduction = "umap",label = T) 
        gd2 <- DimPlot(object = Imm.combined_filtOLG.new, reduction = "tsne",label = T) 
        DefaultAssay(Imm.combined_filtOLG.new) <- "RNA"
          #homeostatic microglia
        gd3 <- FeaturePlot(Imm.combined_filtOLG.new, features = "P2RY12", 
                           reduction = "umap", label=T, raster=TRUE)
        gd4 <- FeaturePlot(Imm.combined_filtOLG.new, features = "P2RY12", 
                           reduction = "tsne", label=T, raster=TRUE)
        gd5 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CX3CR1", 
                           reduction = "umap", label=T, raster=TRUE)
        gd6 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CX3CR1", 
                           reduction = "tsne", label=T, raster=TRUE)
          #inflamed microglia
        gd7 <- FeaturePlot(Imm.combined_filtOLG.new, features = "TREM2",
                           reduction = "umap", label=T, raster=TRUE)
        gd8 <- FeaturePlot(Imm.combined_filtOLG.new, features = "TREM2", 
                           reduction = "tsne", label=T, raster=TRUE)
        gd9 <- FeaturePlot(Imm.combined_filtOLG.new, features = "APOE",
                           reduction = "umap", label=T, raster=TRUE)
        gd10 <- FeaturePlot(Imm.combined_filtOLG.new, features = "APOE",
                            reduction = "tsne", label=T, raster=TRUE)
        gd11 <- FeaturePlot(Imm.combined_filtOLG.new, features = "FTL", 
                            reduction = "umap", label=T, raster=TRUE)
        gd12 <- FeaturePlot(Imm.combined_filtOLG.new, features = "FTL",
                            reduction = "tsne", label=T, raster=TRUE)
          #mature dendritic cells
        gd13 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD83",
                            reduction = "umap", label=T, raster=TRUE)
        gd14 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD83",
                            reduction = "tsne", label=T, raster=TRUE)
        gd15 <- FeaturePlot(Imm.combined_filtOLG.new, features = "NFKB1",
                            reduction = "umap", label=T, raster=TRUE)
        gd16 <- FeaturePlot(Imm.combined_filtOLG.new, features = "NFKB1",
                            reduction = "tsne", label=T, raster=TRUE)
        #T cells
        gd21 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD2", 
                            reduction = "umap", label=T, raster=TRUE)
        gd22 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD2", 
                            reduction = "tsne", label=T, raster=TRUE)
        gd23 <- FeaturePlot(Imm.combined_filtOLG.new, features = "SKAP1", 
                            reduction = "umap", label=T, raster=TRUE)
        gd24 <- FeaturePlot(Imm.combined_filtOLG.new, features = "SKAP1", 
                            reduction = "tsne", label=T, raster=TRUE)
          #monocyte_genes
        gd25 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD14",
                            reduction = "umap", label=T, raster=TRUE)
        gd26 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD14",
                            reduction = "tsne", label=T, raster=TRUE)
          #Macrophage
        gd27 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD68",
                            reduction = "umap", label=T, raster=TRUE)
        gd28 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD68",
                            reduction = "tsne", label=T, raster=TRUE)
          #astrocytes_genes "GFAP", "AQP4"
        gd31 <- FeaturePlot(Imm.combined_filtOLG.new, features = "GFAP",
                            reduction = "umap", label=T, raster=TRUE)
        gd32 <- FeaturePlot(Imm.combined_filtOLG.new, features = "GFAP",
                            reduction = "tsne", label=T, raster=TRUE)

        gd33 <- FeaturePlot(Imm.combined_filtOLG.new, features = "AQP4",
                            reduction = "umap", label=T, raster=TRUE)
        gd34 <- FeaturePlot(Imm.combined_filtOLG.new, features = "AQP4",
                            reduction = "tsne", label=T, raster=TRUE)
          #nonreactive_astrocyte
        gd35 <- FeaturePlot(Imm.combined_filtOLG.new, features = "SLC1A2",
                            reduction = "umap", label=T, raster=TRUE)
        gd36 <- FeaturePlot(Imm.combined_filtOLG.new, features = "SLC1A2",
                            reduction = "tsne", label=T, raster=TRUE)
  
        gd37 <- FeaturePlot(Imm.combined_filtOLG.new, features = "SLC4A4",
                            reduction = "umap", label=T, raster=TRUE)
        gd38 <- FeaturePlot(Imm.combined_filtOLG.new, features = "SLC4A4",
                            reduction = "tsne", label=T, raster=TRUE)
          #reactive_astrocyte
        gd39 <- FeaturePlot(Imm.combined_filtOLG.new, features = "HSP90AA1",
                            reduction = "umap", label=T, raster=TRUE)
        gd40 <- FeaturePlot(Imm.combined_filtOLG.new, features = "HSP90AA1",
                            reduction = "tsne", label=T, raster=TRUE)
        gd41 <- FeaturePlot(Imm.combined_filtOLG.new, features = "SERPINH1",
                            reduction = "umap", label=T, raster=TRUE)
        gd42 <- FeaturePlot(Imm.combined_filtOLG.new, features = "SERPINH1",
                            reduction = "tsne", label=T, raster=TRUE)
          #senescent_astrocyte
        gd43 <- FeaturePlot(Imm.combined_filtOLG.new, features = "DNAH9",
                            reduction = "umap", label=T, raster=TRUE)
        gd44 <- FeaturePlot(Imm.combined_filtOLG.new, features = "DNAH9",
                            reduction = "tsne", label=T, raster=TRUE)
        gd45 <- FeaturePlot(Imm.combined_filtOLG.new, features = "SPAG17",
                            reduction = "umap", label=T, raster=TRUE)
        gd46 <- FeaturePlot(Imm.combined_filtOLG.new, features = "SPAG17",
                            reduction = "tsne", label=T, raster=TRUE)
          # AIMS
        gd47 <- FeaturePlot(Imm.combined_filtOLG.new, features = "VIM",
                            reduction = "umap", label=T, raster=TRUE)
        gd48 <- FeaturePlot(Imm.combined_filtOLG.new, features = "VIM",
                            reduction = "tsne", label=T, raster=TRUE)
        gd49 <- FeaturePlot(Imm.combined_filtOLG.new, features = "S100B",
                            reduction = "umap", label=T, raster=TRUE)
        gd50 <- FeaturePlot(Imm.combined_filtOLG.new, features = "S100B",
                            reduction = "tsne", label=T, raster=TRUE)
          #perinodal astrocytes
        gd51 <- FeaturePlot(Imm.combined_filtOLG.new, features = "IL1RAPL1",
                            reduction = "umap", label=T, raster=TRUE)
        gd52 <- FeaturePlot(Imm.combined_filtOLG.new, features = "IL1RAPL1",
                            reduction = "tsne", label=T, raster=TRUE)
        gd53 <- FeaturePlot(Imm.combined_filtOLG.new, features = "KIRREL3",
                            reduction = "umap", label=T, raster=TRUE)
        gd54 <- FeaturePlot(Imm.combined_filtOLG.new, features = "KIRREL3",
                            reduction = "tsne", label=T, raster=TRUE)
          #macrophage
        gd55 <- FeaturePlot(Imm.combined_filtOLG.new, features = "MRC1",
                            reduction = "umap", label=T, raster=TRUE)
        gd56 <- FeaturePlot(Imm.combined_filtOLG.new, features = "MRC1",
                            reduction = "tsne", label=T, raster=TRUE)
        gd57 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD163",
                            reduction = "umap", label=T, raster=TRUE)
        gd58 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD163",
                            reduction = "tsne", label=T, raster=TRUE)
        gd59 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD163L1",
                            reduction = "umap", label=T, raster=TRUE)
        gd60 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD163L1",
                            reduction = "tsne", label=T, raster=TRUE)
        gd61 <- FeaturePlot(Imm.combined_filtOLG.new, features = "MSR1",
                            reduction = "umap", label=T, raster=TRUE)
        gd62 <- FeaturePlot(Imm.combined_filtOLG.new, features = "MSR1",
                            reduction = "tsne", label=T, raster=TRUE)
          #endothelial_genes 
        gd63 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CLDN5", 
                            reduction = "umap", label=T, raster=TRUE)
        gd64 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CLDN5", 
                            reduction = "tsne", label=T, raster=TRUE)     
          #plasmablasts  CD38", "CD79A", "IGHG", "IGHA", "IGHM")
        gd69 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD38",
                            reduction = "umap", label=T, raster=TRUE)
        gd70 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD38",
                            reduction = "tsne", label=T, raster=TRUE)
        gd71 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD79A",
                            reduction = "umap", label=T, raster=TRUE)
        gd72 <- FeaturePlot(Imm.combined_filtOLG.new, features = "CD79A",
                            reduction = "tsne", label=T, raster=TRUE)
        gd73 <- FeaturePlot(Imm.combined_filtOLG.new, features = "IGHM",
                            reduction = "umap", label=T, raster=TRUE)
        gd74 <- FeaturePlot(Imm.combined_filtOLG.new, features = "IGHM",
                            reduction = "tsne", label=T, raster=TRUE)
        gd75 <- FeaturePlot(Imm.combined_filtOLG.new, features = "IGHG1",
                            reduction = "umap", label=T, raster=TRUE)
        gd76 <- FeaturePlot(Imm.combined_filtOLG.new, features = "IGHG1",
                            reduction = "tsne", label=T, raster=TRUE)

    #Since TSPO has been the most frequently implemented PET imaging target of neuroinflammation in vivo,19 including of chronic active MS lesions
        gd77 <- FeaturePlot(Imm.combined_filtOLG.new, features = "TSPO",
                            reduction = "umap", label=T, raster=TRUE)
        gd78 <- FeaturePlot(Imm.combined_filtOLG.new, features = "TSPO",
                            reduction = "tsne", label=T, raster=TRUE)
        CombinePlots(plots = list(gd1,gd2,gd3,gd4,gd5,gd6,gd7,gd8,gd9,gd10,
                                  gd11,gd12,gd13,gd14,gd15,gd16,
                                  gd21,gd22,gd23,gd24,gd25,gd26,gd27,gd28,
                                  gd31,gd32,gd33,gd34,gd35,gd36,gd37,gd38,gd39,gd40,
                                  gd41,gd42,gd43,gd44,gd45,gd46,gd47,gd48,gd49,gd50,
                                  gd51,gd52,gd53,gd54,gd55,gd56,gd57,gd58,gd59,gd60,
                                  gd61,gd62,gd63,gd64,gd69,gd70,
                                  gd71,gd72,gd73,gd74,gd75,gd76,gd77,gd78), ncol = 2)
        ggsave(filename = paste(
            'All_sample_Micoglia_UMap_tSNE_by_cluster_and_sample_DimPlot_pc' 
            ,i, 'res', j,'.pdf',sep=''), 
               width=20, height=340, limitsize=FALSE)
        }
    }

macrophage_genes <- c("MRC1","MRC1","CD163","CD163L1","MSR1")

Imm.combined_filtOLG <- FindNeighbors(Imm.combined_filtOLG, reduction = "pca", 
                                      dims = 1:30)
Imm.combined_filtOLG <- FindClusters(Imm.combined_filtOLG, resolution = 0.9)
Imm.combined_filtOLG <- RunUMAP(Imm.combined_filtOLG, reduction = "pca", 
                                dims = 1:30)
Imm.combined_filtOLG <- RunTSNE(Imm.combined_filtOLG, reduction = "pca", 
                                dims = 1:30)

#A1 astrocytes 
A1_astro_genes <- c("C3","CD68")
pdf("03.A1_astrocytes.selectgenes_VlnPlot.pdf",  height = 2, width = 24)
VlnPlot(Imm.combined_filtOLG, features = A1_astro_genes, pt.size=0, ncol=2)
dev.off()

pdf("03.A1_astrocytes.selectgenes_FeaturePlot.pdf", height = 7.5, width = 20)
FeaturePlot(Imm.combined_filtOLG, features = A1_astro_genes, 
            reduction = "umap", label=T, ncol=2, raster=TRUE)
dev.off()

## Some failed to pass the cluster re-identification of the Marker gene
### Astrocytes
Astro_11.markers <- FindMarkers(Imm.combined_filtOLG, 
                                ident.1 = 11, 
                                ident.2 = c(0,1,3,4,6,10,13,15,17,18,23,25,26), 
                                min.pct = 0.25, only.pos = TRUE)

Astro_11.all.markers = Astro_11.markers %>% dplyr::select(, everything()
                                                          ) %>% subset(p_val<0.05)

Astro_11.top10 = Astro_11.all.markers %>% group_by(cluster) %>% top_n(
    n = 10, wt = avg_log2FC)

write.csv(Astro_11.markers, "Astro11_diff_genes_wilcox.csv", row.names = T)
write.csv(Imm_filtOLG.top10, "01-2.top10_diff_genes_wilcox.csv", row.names = T)

Astro_11.vs.Nonreactive.markers <- FindMarkers(Imm.combined_filtOLG, ident.1 = 11, ident.2 = c(1,3,6,10,15,18,25), min.pct = 0.25, only.pos = TRUE)

write.csv(Astro_11.vs.Nonreactive.markers, "Astro11.vs.Nonreactive_diff_genes_wilcox.csv", row.names = T)

#
cluster23.markers <- FindMarkers(Imm.combined_filtOLG, ident.1 = 23, 
                                 min.pct = 0.25, only.pos = TRUE)
cluster23.all.markers = cluster23.markers %>% dplyr::select(, everything()) %>% subset(p_val<0.05)
write.csv(cluster23.all.markers, "cluster23.all_diff_genes_wilcox.csv", row.names = T)

cluster26.markers <- FindMarkers(Imm.combined_filtOLG, 
                                 ident.1 = 26, min.pct = 0.25, only.pos = TRUE)
cluster26.all.markers = cluster26.markers %>% dplyr::select(, everything()
                                                           ) %>% subset(p_val<0.05)

write.csv(cluster26.all.markers, "cluster26.all_diff_genes_wilcox.csv", row.names = T)

cluster20.markers <- FindMarkers(Imm.combined_filtOLG, ident.1 = 20, min.pct = 0.25, only.pos = TRUE)
cluster20.all.markers = cluster20.markers %>% dplyr::select(, everything()) %>% subset(p_val<0.05)
write.csv(cluster20.all.markers, "cluster20.all_diff_genes_wilcox.csv", row.names = T)

cluster28.markers <- FindMarkers(Imm.combined_filtOLG, ident.1 = 28, min.pct = 0.25, only.pos = TRUE)
cluster28.all.markers = cluster28.markers %>% dplyr::select(, everything()) %>% subset(p_val<0.05)
write.csv(cluster28.all.markers, "cluster28.all_diff_genes_wilcox.csv", row.names = T)

cluster29.markers <- FindMarkers(Imm.combined_filtOLG, ident.1 = 29, min.pct = 0.25, only.pos = TRUE)
cluster29.all.markers = cluster29.markers %>% dplyr::select(, everything()) %>% subset(p_val<0.05)
write.csv(cluster29.all.markers, "cluster29.all_diff_genes_wilcox.csv", row.names = T)

Astro0.markers <- FindMarkers(Imm.combined_filtOLG, ident.1 = 0, min.pct = 0.25, only.pos = TRUE)
Astro0.all.markers = Astro0.markers %>% dplyr::select(, everything()) %>% subset(p_val<0.05)
write.csv(Astro0.all.markers, "Astro0.all_diff_genes_wilcox.csv", row.names = T)

Astro.markers <- FindMarkers(Imm.combined_filtOLG, ident.1 = , 
                             min.pct = 0.25, only.pos = TRUE)
Astro.all.markers = Astro.markers %>% dplyr::select(,everything()
                                                    ) %>% subset(p_val<0.05)
write.csv(Astro.all.markers, "Astro.all_diff_genes_wilcox.csv", row.names = T)

# perinodal astrocytes
perinodal_astrocyte_genes <- c("IL1RAPL1", "KIRREL3", "UNC5C", "ANK3", 
                               "CNTN2", "MAG", "LRP2", "MYO1D", "MBP", "SHTN1", 
                               "BIN1", "DNM3", "GRM3", "APBB2", "AGTPBP1", 
                               "SPOCK1", "MAP7", "ALCAM", "DLG1")
pdf("01.Imm_astrocytes_perinodal.selectgenes_VlnPlot.pdf", height = 20, width = 24) 
VlnPlot(Imm.combined_filtOLG, features = perinodal_astrocyte_genes, 
        pt.size=0, ncol=2)
dev.off()

pdf("01.Imm_astrocytes_perinodal.selectgenes_FeaturePlot.pdf",
    height = 50, width = 14) # height = 5 x lines
FeaturePlot(Imm.combined_filtOLG, features = perinodal_astrocyte_genes, 
            reduction = "umap", label=T, ncol=2, raster=TRUE)
dev.off()

## Rename immune celltypes 
Immune_celltype <- c("Astrocyte","Astrocyte","Microglia","Astrocyte","Astrocyte","Microglia","Astrocyte","Microglia","Dendritic cells","Microglia","Astrocyte","Astrocyte","Microglia","Astrocyte","T cells","Astrocyte","Astrocyte","Astrocyte","Astrocyte","Macrophage","Filtered","Astrocyte","Microglia","Astrocyte","Macrophage","Astrocyte","Astrocyte","B cells","Filtered","Filtered","Filtered","Filtered")
names(Immune_celltype) <- levels(Imm.combined_filtOLG)
write.csv(Immune_celltype,"02.subtypes/Immune_celltype.csv")

Imm.combined_filtOLG_AddSubtype <- RenameIdents(Imm.combined_filtOLG, Immune_celltype)
CellsMeta <- Imm.combined_filtOLG@meta.data
CellType <- Imm.combined_filtOLG_AddSubtype@active.ident
CellsMeta["Immune_celltype"] <- CellType
CellsMetaTrim <- subset(CellsMeta, select = Immune_celltype)
Imm.combined_filtOLG_AddSubtype <- AddMetaData(Imm.combined_filtOLG_AddSubtype, CellsMetaTrim)

write.csv(CellType, "02.cell_identify/CellType.csv")

Imm.combined_filtOLG_AddSubtype@meta.data$Immune_subtype <- ""

cluster2celltype <- c("0"="Nonreactive astrocyte", "1"="Nonreactive astrocyte",
                      "2"="Homeostatic microglia", "3"="Nonreactive astrocyte",
                      "4"="Reactive astrocyte", "5"="Inflamed microglia",
                      "6"="Nonreactive astrocyte", "7"="Inflamed microglia",
                      "8"="Dendritic cells", "9"="Homeostatic microglia",
                      "10"="Nonreactive astrocyte", "11"="Nonreactive astrocyte",
                      "12"="Inflamed microglia", "13"="Reactive astrocyte",
                      "14"="T cells", "15"="Nonreactive astrocyte",
                      "16"="Senescent astrocytes", "17"="Reactive astrocyte",
                      "18"="Nonreactive astrocyte", "19"="Macrophage",
                      "20"="-", "21"="Senescent astrocytes",
                      "22"="Inflamed microglia", "23"="Nonreactive astrocyte",
                      "24"="Macrophage", "25"="Nonreactive astrocyte",
                      "26"="Nonreactive astrocyte", "27"="B cells",
                      "28"="-", "29"="-",
                      "30"="-", "31"="-")

Imm.combined_filtOLG_AddSubtype[['Immune_subtype']] = unname(
    cluster2celltype[Imm.combined_filtOLG_AddSubtype@meta.data$seurat_clusters])

saveRDS(Imm.combined_filtOLG_AddSubtype, 
        "02.immune_cells_filtOLG/02.subtypes/Imm.combined_filtOLG_AddSubtype.RDS")

## plot作图
setwd("02.immune_cells_filtOLG/")

pdf("01.umap_Imm_Immune.celltype-lite.pdf", width = 8.5,height=7.5)
DimPlot(subset(Imm.combined_filtOLG_AddSubtype,Immune_celltype!="Filtered"),
        reduction = "umap", group.by = "Immune_celltype", 
        label=TRUE, raster=TRUE)
dev.off()

pdf("01.umap_Imm_Immune.celltype.pdf", width = 8.5, height=7.5)
DimPlot(subset(Imm.combined_filtOLG_AddSubtype, 
               Immune_celltype!="Filtered"),
        reduction = "umap", group.by = "Immune_celltype", 
        label=TRUE, raster=FALSE)
dev.off()

pdf("01.umap_Imm_Immune.subtype.pdf", width = 8.5, height=7.5)
DimPlot(subset(Imm.combined_filtOLG_AddSubtype,Immune_celltype!="-"), 
        reduction = "umap", group.by = "Immune_subtype", label=TRUE, raster=FALSE)
dev.off()

pdf("01.umap_Imm_Immune.subtype_nolabel.pdf", width = 8.5, height=7.5)
DimPlot(subset(Imm.combined_filtOLG_AddSubtype,Immune_subtype!="-"), 
        reduction = "umap", group.by = "Immune_subtype", 
        label=FALSE, raster=FALSE)
dev.off()

pdf("01.umap_Imm_Immune.subtype_nolabel_Diagnosis-split.pdf", width = 8.5, height=7.5)
DimPlot(subset(Imm.combined_filtOLG_AddSubtype,Immune_subtype!="-"), 
        reduction = "umap", group.by = "Immune_subtype", 
        split.by="Diagnosis", label=FALSE, raster=FALSE)
dev.off()


pdf("01.umap_Imm_Immune.subtype.region.pdf", width = 8.5, height=7.5)
DimPlot(subset(subset(Imm.combined_filtOLG_AddSubtype,Immune_celltype!="Filtered"), 
               Lregion_sym=="Control"| Lregion_sym=="CA"| Lregion_sym=="CA edge"| Lregion_sym=="CI"| Lregion_sym=="CI edge"| Lregion_sym=="NAWM"), 
        reduction = "umap", group.by = "Lregion_sym", label=TRUE, raster=FALSE)
dev.off()

pdf("01.umap_Imm_Immune.subtype.region-split.pdf", width = 23.5, height=15)
DimPlot(subset(subset(Imm.combined_filtOLG_AddSubtype,Immune_subtype!="-"), 
               Lregion_sym=="Control"| Lregion_sym=="CA"| Lregion_sym=="CA edge"| Lregion_sym=="CI"| Lregion_sym=="CI edge"| Lregion_sym=="NAWM"), 
        reduction = "umap", split.by = "Lregion_sym", label=TRUE, raster=FALSE, ncol=3)
dev.off()

pdf("01.umap_Imm_Immune.ATF3_FeaturePlot.pdf", width = 8.5, height=7.5)
FeaturePlot(subset(subset(Imm.combined_filtOLG_AddSubtype,Immune_subtype!="-"), 
               Lregion_sym=="Control"| Lregion_sym=="CA"| Lregion_sym=="CA edge"| Lregion_sym=="CI"| Lregion_sym=="CI edge"| Lregion_sym=="NAWM"), 
            features = "ATF3", label.by="Immune_subtype",
            reduction = "umap", label=T, raster=TRUE)
dev.off()

pdf("01.Imm_Immune.ATF3_VlnPlot.pdf", height = 4, width = 12)
VlnPlot(subset(subset(Imm.combined_filtOLG_AddSubtype,Immune_subtype!="-"), 
               Lregion_sym=="Control"| Lregion_sym=="CA"| Lregion_sym=="CA edge"| Lregion_sym=="CI"| Lregion_sym=="CI edge"| Lregion_sym=="NAWM"), 
        features = "ATF3", pt.size=0)
dev.off()

pdf("01.umap_Imm_Immune.HES6_FeaturePlot.pdf", width = 8.5, height=7.5)
FeaturePlot(subset(subset(Imm.combined_filtOLG_AddSubtype,Immune_subtype!="-"), 
               Lregion_sym=="Control"| Lregion_sym=="CA"| Lregion_sym=="CA edge"| Lregion_sym=="CI"| Lregion_sym=="CI edge"| Lregion_sym=="NAWM"), 
            features = "HES6", 
            reduction = "umap", label=T, raster=TRUE)
dev.off()

pdf("01.Imm_Immune.HES6_VlnPlot.pdf", height = 4, width = 12)
VlnPlot(subset(subset(Imm.combined_filtOLG_AddSubtype,Immune_subtype!="-"), 
               Lregion_sym=="Control"| Lregion_sym=="CA"| Lregion_sym=="CA edge"| Lregion_sym=="CI"| Lregion_sym=="CI edge"| Lregion_sym=="NAWM"), 
        features = "HES6", pt.size=0)
dev.off()

## Immune scoring
library(tidyverse)
library(Matrix)
library(cowplot)
RK <- Imm.combined_filtOLG_AddSubtype

## immune - cytokine score 
setwd("result_LSA/")
cytokine_gene <- read.csv("reference/03.immune_associate/cytokines_gene.xlsx", header=FALSE)
# 
gene <- as.list(cytokine_gene)

## immune inflammatory genes
inflammatory_gene <- read.csv("reference/03.immune_associate/inflammatory_gene.xlsx", header=FALSE)
# 
gene_inflam <- as.list(inflammatory_gene)

RK_copy <- AddModuleScore(RK, features = gene, 
                           ctrl=100, # 默认是100
                           name = "CD_Features")

RK_inflam <- AddModuleScore(RK, features = gene_inflam, 
                            ctrl=100, 
                            name = "CD_Features")

colnames(RK_copy@meta.data)
colnames(RK_copy@meta.data)[41] <- "cytokine_Score"

colnames(RK_inflam@meta.data)
colnames(RK_inflam@meta.data)[41] <- "Inflammatory_Score"

dir.create("./02.cell_identify/02.immune_cells_filtOLG/03.immune_score/")
setwd("./02.cell_identify/02.immune_cells_filtOLG/03.immune_score/")

pdf("01.Cytokine_score_immune.celltype.pdf")
VlnPlot(RK_copy,features = 'cytokine_Score', pt.size=0, group.by="Immune_celltype")
dev.off()
pdf("01.Inflammatory_score_immune.celltype.pdf")
VlnPlot(RK_inflam,features = 'Inflammatory_Score', pt.size=0, group.by="Immune_celltype")
dev.off()

# box plot
data<- FetchData(RK_inflam,vars = c("Immune_celltype","Inflammatory_Score"))
pdf("01.Inflammatory_score_barplot_immune.celltype.pdf")
ggplot(data,aes(Immune_celltype,Inflammatory_Score)) + 
    geom_boxplot()+theme_bw()+RotatedAxis()
dev.off()

data_inflam_allinfo <- FetchData(RK_inflam, vars = 
                                 c("UMAP_1","UMAP_2","Inflammatory_Score",
                                   "Diagnosis", "seurat_clusters", 
                                   "Immune_celltype","Immune_subtype", "Lregion_sym", 
                                   "Clinical.disease.course","Age","Sex"))
data_inflam_filtered <- subset(data_inflam_allinfo,Immune_celltype!="Filtered"&Immune_subtype!="-")

pdf("Fig.2B.Boxplot.inflam.Score_Immune.celltype_Diagnosis.pdf",width=14, height=5)
    ggplot(data_inflam_filtered, aes(Immune_celltype,Inflammatory_Score,fill=Diagnosis)) + 
    geom_boxplot() + 
    facet_grid(.~Immune_celltype,scales = "free",space="free_x") + 
    theme_bw()+theme(axis.text.x = element_blank())
dev.off()

pdf("Fig.2B.Boxplot.inflam.Score_Immune.celltype_LRegion.pdf",width=14)
ggplot(data_inflam_filtered, aes(Immune_celltype,
                                 Inflammatory_Score,fill=Lregion_sym)) + 
    geom_boxplot() + 
    facet_grid(.~Immune_celltype,scales = "free",space="free_x") + 
    theme_bw()+theme(axis.text.x = element_blank())
dev.off()


