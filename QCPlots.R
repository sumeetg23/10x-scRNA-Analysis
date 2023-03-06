qcplots <- function(seuratobject) {
  
  qc.metrics <- as_tibble(seuratobject[[]], rownames="Cell.Barcode") 
  
  currentdir <- getwd()
  
  pathtodir <- paste(currentdir,"QCPlots",sep = "/")
  
  creatdir(pathtodir)
  
  QCPlot1 <- VlnPlot(data, features=c("nCount_RNA","nFeature_RNA","percent.mt", "percent.Ribosomal","percent.Largest.Gene","percent.pymt"), group.by = "Experiment", split.by = "hash.ID", split.plot = TRUE)
  QCPlot2 <- VlnPlot(data, features=c("nCount_RNA","nFeature_RNA","percent.mt", "percent.Ribosomal","percent.Largest.Gene","percent.pymt"), group.by = "Experiment", split.by = "hash.ID", split.plot = TRUE, pt.size = 0)
  QCPlot3 <- VlnPlot(data, features=c("percent.gfp","percent.pymt"), group.by = "Experiment", split.by = "hash.ID", split.plot = TRUE)
  QCPlot4 <- qc.metrics %>% 
    ggplot(aes(x=Experiment, fill=samplename)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")
  
  QCPlot5 <- qc.metrics %>%
    ggplot(aes(color=Experiment, x=percent.mt, fill= samplename)) + 
    geom_density(alpha = 0.2) + 
    ggtitle("Distribution of Percent Mitocondrial Reads") +
    theme_classic()

  QCPlot6 <- qc.metrics %>%
    ggplot(aes(color=Experiment, x=percent.Largest.Gene, fill= samplename)) + 
    geom_density(alpha = 0.2) + 
    ggtitle("Distribution of Percentage Largest Gene") +
    theme_classic()
  
  QCPlot7 <- qc.metrics %>% 
    ggplot(aes(color=Experiment, x=nCount_RNA, fill= samplename)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    expand_limits(x = 1)
  
  # Visualize the distribution of genes detected per cell via histogram
  QCPlot8 <- qc.metrics %>% 
    ggplot(aes(color=Experiment, x=nFeature_RNA, fill= samplename)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10()
  
  # Visualize the distribution of genes detected per cell via boxplot
  QCPlot9 <- qc.metrics %>% 
    ggplot(aes(x=samplename, y=log10(nFeature_RNA), fill=Experiment)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NGenes")
  
  creatplot(QCPlot1, filename = "InitialQC-VinPlot-AllFeatures.jpeg", w=1200, h=600)
  creatplot(QCPlot2, filename = "InitialQC-VinPlot-NoPoints.jpeg", w=1200, h=600)
  creatplot(QCPlot3, filename = "InitialQC-VinPlot-GFPAndPyMT.jpeg", w=600, h=300)
  creatplot(QCPlot4, filename = "InitialQC-NumberOfCells.jpeg", w=600, h=600)
  creatplot(QCPlot5, filename = "InitialQC-PercentMito.jpeg", w=600, h=600)
  creatplot(QCPlot6, filename = "InitialQC-PercentLargestGene.jpeg", w=600, h=600)
  creatplot(QCPlot7, filename = "InitialQC-CountofRNAMolecules.jpeg", w=600, h=600)
  creatplot(QCPlot8, filename = "InitialQC-Hist-GenesPerCell.jpeg", w=600, h=600)
  creatplot(QCPlot9, filename = "InitialQC-BoxPlot-GenesPerCell.jpeg", w=600, h=600)

  for(i in unique(qc.metrics$Experiment)){
    titlename = paste("InitialQC","Scatter","GenesvsUMIPerCellWithPercentMT",i,".jpeg",sep="-")
    creatplot(qc.metrics[qc.metrics$Experiment == i,] %>% 
      ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
      geom_point() + 
      scale_colour_gradient(low = "gray90", high = "black") +
      theme_classic() +
      scale_x_log10() + 
      scale_y_log10() +
      ggtitle(i),filename = titlename, w=600, h=600)
  }
  
  setwd(currentdir)
}

creatdir <- function(sub_dir){

  # check if sub directory exists 
  if (!dir.exists(sub_dir)){
    
    dir.create(sub_dir)
    
  }
  setwd(sub_dir)
  
}

creatplot <- function(plot, filename="plot.jpeg", type="jpeg", w=600, h=800){
  
  if(type == "jpeg"){
    jpeg(filename, width=w, height=h)
      print(plot)
    dev.off()
  }
  
}