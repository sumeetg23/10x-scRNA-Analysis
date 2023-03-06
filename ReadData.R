readalldata <- function(sampleinfo){

  for (i in 1:nrow(sampleinfo)){
    
    data_dir <- sampleinfo$data_dir[i]
    data <- Read10X(data.dir = data_dir)
    allcells <- CreateSeuratObject(counts = data$`Gene Expression`)
    allcells[["HTO"]] <- CreateAssayObject(counts = data$`Antibody Capture`)
    allcells[["percent.mt"]] <- PercentageFeatureSet(object = allcells, pattern = "^mt-*")
    allcells[["percent.Ribosomal"]] <- PercentageFeatureSet(object = allcells, pattern = "^Rp[ls]")
    allcells[["percent.pymt"]] <- PercentageFeatureSet(object = allcells, pattern = "Pymt")
    allcells[["percent.gfp"]] <- PercentageFeatureSet(object = allcells, pattern = "Gfpluc")
    
    allcells[rownames(allcells) != "Malat1",] -> data.nomalat
    apply(data.nomalat@assays$RNA@counts,2,max) -> data.nomalat$largest_count
    apply(data.nomalat@assays$RNA@counts,2,which.max) -> data.nomalat$largest_index
    rownames(data.nomalat)[data.nomalat$largest_index] -> data.nomalat$largest_gene
    100 * data.nomalat$largest_count / data.nomalat$nCount_RNA -> data.nomalat$percent.Largest.Gene
    data.nomalat$largest_gene -> allcells$largest_gene
    data.nomalat$percent.Largest.Gene -> allcells$percent.Largest.Gene
    rm(data.nomalat)
    
    
    if(sampleinfo$sample[i] == "Day35") {
      DefaultAssay(allcells) <- "HTO"
      allcells[rownames(allcells) != "Hashtag-3",] -> newdata
      newdata[rownames(newdata) != "Hashtag-4",] -> newdata
      DefaultAssay(allcells) <- "RNA"
      allcells[["HTO"]] <- newdata[["HTO"]]
    }
    
    allcells <- NormalizeData(allcells, assay = "HTO", normalization.method = "CLR")
    allcells <- HTODemux(allcells, assay = "HTO", positive.quantile = 0.99)
    print(table(allcells$HTO_classification.global))
    Idents(allcells) <- "HTO_maxID"
    Idents(allcells) <- "HTO_classification.global"
    
    if(sampleinfo$sample[i] == "Day35") {
      neg.combined <- subset(allcells, subset = HTO_classification.global == "Negative" & percent.pymt >0 )
      neg.combined$hash.ID <- "Hashtag-2"
    }
    
    allcells <- subset(allcells, idents = "Singlet")
    print(allcells)

    if(sampleinfo$sample[i] == "Day35") {
      allcells <- merge(allcells, neg.combined)
    }
    
    allcells <- subset(allcells, subset = HTO_maxID != "Hashtag-3")
    print(allcells)
    allcells <- subset(allcells, subset = HTO_maxID != "Hashtag-4")
    print(allcells)
    
    allcells$Experiment <- sampleinfo$sample[i]
    allcells <- RenameCells(allcells, add.cell.id = i)
    
    if(i > 1) {
      combined <- merge(combined, allcells)
    } else {
      combined <- allcells
    }
    
  }
  
  combined$samplename <- combined$hash.ID
  combined$samplename[(which(str_detect(combined$samplename, "Hashtag-1")))] <- "WT"
  combined$samplename[(which(str_detect(combined$samplename, "Hashtag-2")))] <- "YS-OE"
  
  return(combined)
  
}