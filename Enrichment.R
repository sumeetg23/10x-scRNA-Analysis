library(msigdbr)
library(presto)
library(fgsea)
library(stats)

getgeneset <- function(x){
  setofgenes <- msigdbr(species = "Mus musculus", category = x)
  return(setofgenes)
}

#hallmarksets <- msigdbr(species = "Mus musculus", category = "H")

performgsea <- function(combined_filtered){
  
  celltypesig <- getgeneset("C8")
  genesettouse<- celltypesig %>% split(x = .$gene_symbol, f = .$gs_name)
  
  all.genes <- wilcoxauc(combined_filtered, 'seurat_clusters')
  dplyr::count(all.genes, group)
  
  rm(tophits)
  rm(allhits)
  rm(allcbind)
  
  for (i in sort(unique(combined_filtered$seurat_clusters))) {
    print(paste("for cluster ",i))
    subset.genes<- all.genes %>% dplyr::filter(group == i) %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
    ranks<- deframe(subset.genes)
    fgseaRes <- fgsea(genesettouse, stats = ranks, nperm = 1000)
    fgseaRes$cluster <- i
    if(exists("tophits") && is.data.frame(tophits)){
      tophits <- rbind(tophits, (head(fgseaRes %>% arrange(desc(NES)) %>% arrange(padj) %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme))))
      allcbind <- rbind(allcbind, fgseaRes %>% arrange(desc(NES)) %>% arrange(padj) %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme))
    }
    else {
      tophits <- (head(fgseaRes %>% arrange(desc(NES)) %>% arrange(padj) %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme)))
      allcbind  <- fgseaRes %>% arrange(desc(NES)) %>% arrange(padj) %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme)
    }
    
    fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
    
    colnames(fgseaRes)[2:(ncol(fgseaRes))] <- paste(colnames(fgseaRes)[2:(ncol(fgseaRes))],"_Cluster",i, sep = "")
    
    if(exists("allhits") && is.data.frame(allhits)){
      allhits <- merge(fgseaRes, allhits, by = "pathway")
    }
    else {
      allhits <- fgseaRes
    }
    
  }
  
  allcbindsubset <- subset(allcbind, pathway %in% unique(tophits$pathway))
  allcbindsubset <- allcbindsubset[allcbindsubset[,padj <0.05],]
  
  splittophits <- split(tophits, by="cluster")
  
  rm(forhc)
  for(i in 1:length(splittophits)){
    x <- data.frame(NES=splittophits[[i]]$NES, row.names=splittophits[[i]]$pathway)
    colnames(x) <- paste("NES_Cluster",i,sep = "")
    if(exists("forhc") && is.data.frame(forhc)){
      forhc <- combine(forhc,x)
    }
    else {
      forhc <- x
    }
  }
  
  forhc[is.na(forhc)] <- 0
  clustering <- hclust(dist(t(forhc)))
  #plot(clustering)
  
  rownames(allhits) <- allhits$pathway
  allhitsub <- allhits[unique(allcbindsubset$pathway),]
  rownames(allhitsub) <- allhitsub$pathway
  #allhitsub <- allhitsub[,grepl( "NES" , names( allhitsub ) )]
  allhitsub <- allhitsub %>% dplyr:: select(grep("NES", colnames(allhitsub)))
  clustering <- hclust(dist(t(allhitsub)))
  
  dotplot <- ggplot(tophits, aes(x=factor(cluster,level=str_replace(colnames(allhitsub)[clustering$order],"NES_Cluster", "")), 
                      y=pathway, 
                      colour=padj,
                      size=NES)) + geom_point()
  
  return(dotplot)
  
}


