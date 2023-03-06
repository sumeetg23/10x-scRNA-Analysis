library(plyr)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(yaml)
library(AUCell)
library(GSEABase)
library(RColorBrewer)
library(msigdbr)
library(presto)
library(fgsea)
library(scCustomize)
library(org.Mm.eg.db)

parameters <- read_yaml("/home/sgupta/RClass/ESR1-Analysis/paramters.yaml")

markerstoplot <- read.table("/home/sgupta/RClass/ESR1-Analysis/markers-UMAP.txt")

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1","brown","blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4","skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F" # lt orange
)

setwd(parameters$workingdir)
source("ReadData.R")
source("QCPlots.R")
source("Markers.R")
source("Enrichment.R")

allmarkers <- getPanglaoDBmarkers()

allcolors <- getallcolors()

sampleinfo <- read.table(parameters$Input, header = T)

data <- readalldata(sampleinfo)
gc()

#qcplots(data)

for (i in unique(data$Experiment)){

  if(!is.null(parameters$params[[i]])){
    MIN_NCOUNT_RNA = parameters$params[[i]]$nCount_RNA
    MIN_FEATURE_RNA = parameters$params[[i]]$nFeature_RNA
    MIN_PERCENTMT = parameters$params[[i]]$percent.mt
    if(!is.null(parameters$params[[i]]$max_nCount_RNA)){
      MAX_NCOUNT_RNA = parameters$params[[i]]$max_nCount_RNA
    }
    else {
      MAX_NCOUNT_RNA = max(data$nCount_RNA)
    }
  }else{
    MIN_NCOUNT_RNA = parameters$params[["Default"]]$nCount_RNA
    MIN_FEATURE_RNA = parameters$params[["Default"]]$nFeature_RNA
    MIN_PERCENTMT = parameters$params[["Default"]]$percent.mt
    if(!is.null(parameters$params[["Default"]]$max_nCount_RNA)){
      MAX_NCOUNT_RNA = parameters$params[["Default"]]$max_nCount_RNA
    }
    else{
      MAX_NCOUNT_RNA = max(data$nCount_RNA)
    }
  }

  cat("Filtering = ",i," with following paramters\n")
  cat("MIN NCOUNT RNA = ",MIN_NCOUNT_RNA,"\n")
  cat("MIN FEATURE RNA = ",MIN_FEATURE_RNA,"\n")
  cat("MIN PERCENTMT = ",MIN_PERCENTMT,"\n")
  cat("MAX NCOUNT RNA = ",MAX_NCOUNT_RNA,"\n","\n")
  
  sub_dir = paste(parameters$workingdir,i,sep = "/")
  creatdir(sub_dir)
  
  data_filtered <- subset(data, subset = nCount_RNA > MIN_NCOUNT_RNA &
                                nFeature_RNA > MIN_FEATURE_RNA &
                                percent.mt < MIN_PERCENTMT &
                                Experiment == i)
  counts <- GetAssayData(object = data_filtered, slot = "counts")
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  filtered_counts <- counts[keep_genes, ]
  
  combined_filtered <- CreateSeuratObject(filtered_counts, meta.data = data_filtered@meta.data)
  combined_filtered <- NormalizeData(combined_filtered)
  combined_filtered <- FindVariableFeatures(combined_filtered, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(combined_filtered)
  combined_filtered <- ScaleData(combined_filtered, features = all.genes)
  
  combined_filtered <- RunPCA(object = combined_filtered, features = VariableFeatures(object = combined_filtered))
  combined_filtered <- RunUMAP(object = combined_filtered, dims = 1:20)
  combined_filtered <- FindNeighbors(combined_filtered, reduction = "pca", dims = 1:40)
  combined_filtered <- FindClusters(combined_filtered, resolution = 0.5)
  
  p1 <- DimPlot(combined_filtered, reduction = "umap", label = TRUE, split.by = 'samplename', pt.size = 2)
  p1 + ggtitle(i)

  p1pca <- DimPlot(combined_filtered, reduction = "pca", label = TRUE, split.by = 'samplename', pt.size = 2)
  p1pca + ggtitle(i)
  
  if("Pymt" %in% all.genes & "Gfpluc" %in% all.genes){
    
    p2 <- FeaturePlot(combined_filtered, features = c('Pymt','Gfpluc'), split.by = 'samplename', keep.scale = "feature", max.cutoff="q90", pt.size = 2)
    p2 + ggtitle(i)
    
    p3 <- FeatureScatter(object = combined_filtered, feature1 = 'Pymt', feature2 = 'Gfpluc', pt.size = 2)
    p3 + ggtitle(i)
  
    creatplot(p2,filename = paste(sub_dir,"/UMAP-PyMTAndGFP-",i,".jpeg",sep = ""))
    creatplot(p3,filename = paste(sub_dir,"/Scatter-PyMTAndGFP-",i,".jpeg",sep = ""))
    
  }
  
  p4 <- FeaturePlot(combined_filtered, features = c('percent.Ribosomal'), split.by = 'samplename', pt.size = 2)
  p4 + ggtitle(i)
  
  for (m in unlist(markerstoplot)) {
    if(m %in% all.genes){
      markerplot <- FeaturePlot(combined_filtered, features = c(m), split.by = 'samplename', pt.size = 2)
      markerplot + ggtitle(i)
      creatplot(markerplot,filename = paste(sub_dir,"/UMAP-",i,"_Marker-",m,".jpeg",sep = ""), w=1600)
      markerplot <- FeaturePlot(combined_filtered, features = c(m), reduction = "pca", split.by = 'samplename', pt.size = 2)
      markerplot + ggtitle(i)
      creatplot(markerplot,filename = paste(sub_dir,"/PCA-",i,"_Marker-",m,".jpeg",sep = ""), w=1600)
    }
  }
  
  umapbycluster <- DimPlot(combined_filtered, reduction = "umap", group.by = 'samplename', split.by = "seurat_clusters", pt.size = 2)
  pcabycluster <- DimPlot(combined_filtered, reduction = "pca", group.by = 'samplename', split.by = "seurat_clusters", pt.size = 2)
  
  
  cells_rankings <- AUCell_buildRankings(combined_filtered[['RNA']]@counts, plotStats=TRUE, splitByBlocks=TRUE)
  cells_AUC <- AUCell_calcAUC(allmarkers, cells_rankings)
  dataauc <- getAUC(cells_AUC)
  AUCannotations <- apply(dataauc, 2, function(x){names(which.max(x))})
  combined_filtered[["AUCannotations"]] <- AUCannotations
  clusterandanno <- paste(combined_filtered$seurat_clusters, "_", combined_filtered$AUCannotations, sep="")
  Idents(combined_filtered) <- combined_filtered@meta.data$AUCannotations
  clusterassign <- data.frame(table(clusterandanno))
  splitvalues <- data.frame(do.call('rbind', strsplit(as.character(clusterassign$clusterandanno), split = "_", fixed=TRUE)))
  clusterassign <- cbind(clusterassign,splitvalues)
  colnames(clusterassign) <- c("id","Freq","Cluster","CellType")
  clusterassign <-ddply(clusterassign, "Cluster", transform, TotalCells=sum(Freq))
  clusterassign$percentcells <- 100*(clusterassign$Freq/clusterassign$TotalCells)
  clusterassignsubset <- clusterassign[clusterassign$percentcells>10,]
  clusterassignsubset[order(as.numeric(clusterassignsubset$Cluster)),]
  clusterassignsubsetmax <- clusterassignsubset %>% group_by(Cluster) %>% filter(percentcells == max(percentcells, na.rm=TRUE))
  
  # to keep only 1 option for 1 cluster
  clusterassignsubsetmax <- clusterassignsubsetmax[!duplicated(clusterassignsubsetmax[,c('Cluster')]),]
  
  write.table(clusterassignsubsetmax, file = "clusternumber.txt")
  
  clusterlabel <- data.frame(clusterassignsubsetmax$Cluster,clusterassignsubsetmax$CellType)
  rownames(clusterlabel) <- clusterlabel$clusterassignsubsetmax.Cluster
  clusterlabel <- clusterlabel[order(as.numeric(clusterlabel$clusterassignsubsetmax.Cluster)),]
  
  combined_filtered$celllabel <- combined_filtered$seurat_clusters
  Idents(combined_filtered) <- "seurat_clusters"
  new.cluster.ids <- clusterlabel$clusterassignsubsetmax.CellType
  names(new.cluster.ids) <- levels(combined_filtered)
  Idents(combined_filtered) <- "celllabel"
  combined_filtered <- RenameIdents(combined_filtered, new.cluster.ids)
  
  umapplotwithcelltype <- DimPlot(combined_filtered, reduction = "umap", label.size = 6, label = TRUE, split.by = 'samplename', pt.size = 2) + NoLegend() + ggtitle(i)
  
  combined_filtered[["cluster_and_celltype"]] <- paste(combined_filtered$seurat_clusters,"_",new.cluster.ids[combined_filtered$celllabel], sep = "")
  
  celltypeplot <- data.frame(table(combined_filtered$samplename,new.cluster.ids[combined_filtered$celllabel]))
  write.table(celltypeplot, file = "clusternumber2.txt")
  celltypedist <- ggplot(celltypeplot, aes(x=Var1, y=Freq, fill=Var2)) + 
    geom_col(position = "fill") + 
    #scale_fill_brewer(palette = c25) +
    scale_fill_manual(values = allcolors[unique(unname(new.cluster.ids))]) +
    #scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Paired"))(17)) +
    theme_classic() + theme(legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=12)) +
    ggtitle(i)
  
  creatplot(p1,filename = paste(sub_dir,"/UMAP-Clusters-",i,".jpeg",sep = ""), w=1600)
  creatplot(p1pca,filename = paste(sub_dir,"/PCA-Clusters-",i,".jpeg",sep = ""), w=1600)
  creatplot(p4,filename = paste(sub_dir,"/UMAP-PercentRibo-",i,".jpeg",sep = ""), w=1600)
  creatplot(umapbycluster,filename = paste(sub_dir,"/UMAP-SplitByCluster-",i,".jpeg",sep = ""), w=1600)
  creatplot(pcabycluster,filename = paste(sub_dir,"/PCA-SplitByCluster-",i,".jpeg",sep = ""), w=1600)
  creatplot(umapplotwithcelltype,filename = paste(sub_dir,"/UMAP-With-CellType-AUCELL",i,".jpeg",sep = ""), w=1600)
  creatplot(celltypedist,filename = paste(sub_dir,"/CellTypeDistribution-AUCELL",i,".jpeg",sep = ""))
  
  table(combined_filtered$seurat_clusters,combined_filtered$samplename)
  
  celltypegseadotplot <- performgsea(combined_filtered)
  creatplot(celltypegseadotplot,filename = paste(sub_dir,"/CellTypeDistribution-GSEA",i,".jpeg",sep = ""), w=1600)
  
  allcells.markers <- FindAllMarkers(combined_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  allcells.markers <- allcells.markers[allcells.markers$p_val_adj<0.05,]
  write.table(allcells.markers,file = "AllMarkers.xls",quote = FALSE, sep = "\t", col.names=NA)
  
  combined_filtered[["cd45Pos"]] <- "Cd45Neg"
  combined_filtered$cd45Pos[WhichCells(combined_filtered, expression = Ptprc > 0)] <- "Cd45Pos"
  test <- as.data.frame.matrix(as.matrix(t(table(combined_filtered$cd45Pos, combined_filtered$AUCannotations))))
  test <- cbind(clusterassignsubsetmax,test[clusterassignsubsetmax$CellType,])
  test$Ratio <- test$Cd45Pos/test$Cd45Neg
  test$Cd45 <- "Neg"
  test[test$Ratio>1,]$Cd45 <- "Pos"
  test <- as.data.frame(test)
  rownames(test) <- test$id
  combined_filtered$Cd45status <- test[combined_filtered$cluster_and_celltype,]$Cd45
  forpie <- as.data.frame.table(table(combined_filtered$Cd45status,combined_filtered$samplename))
  for (z in sort(unique(forpie$Var2))) {
    forpiesubset <- forpie[forpie$Var2 == z,]
    forpiesubset$percent <- round((forpiesubset$Freq/sum(forpiesubset$Freq))*100)
    jpeg(paste("PieChart",z,"CD45.jpeg",sep = "_"), width=600, height=600)
      pie(forpiesubset$Freq, labels = paste(forpiesubset$Var1, sep = " ", forpiesubset$percent, "%"),  col = rainbow(length(forpiesubset$Var1)), main = paste("CD45",z,sep = "_"))
    dev.off()
  }
  
  combined_filtered[["celltype2"]] <- test[combined_filtered$cluster_and_celltype,]$CellType
  forpie <- as.data.frame.table(table(combined_filtered$celltype2,combined_filtered$Cd45status,combined_filtered$samplename))
  forpie <- forpie[forpie$Freq != 0,]
  for (conditionforpie in sort(unique(forpie$Var3))) {
    forpiecondition <- forpie[forpie$Var3 == conditionforpie,]
    for (z in sort(unique(forpiecondition$Var2))) {
      forpiesubset <- forpiecondition[forpiecondition$Var2 == z,]
      forpiesubset$percent <- round((forpiesubset$Freq/sum(forpiesubset$Freq))*100)
      jpeg(paste("PieChart",z,conditionforpie,"CD45.jpeg",sep = "_"), width=1000, height=1000)
        pie(forpiesubset$Freq, labels = paste(forpiesubset$Var1, sep = " ", forpiesubset$percent, "%"),  col = rainbow(length(forpiesubset$Var1)), main = paste("CD45",z,conditionforpie,sep = "_"))
      dev.off()
    }
  }
  
}

