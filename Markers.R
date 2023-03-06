library(GSEABase)
library(Polychrome)

mousify <- function(a){
  return(paste0(substr(a,1,1), tolower(substr(a,2,nchar(a)))))
}

getPanglaoDBmarkers <- function(){
  
  markers <- getallmarkers()
  genes <- markers$official.gene.symbol
  genes <- sapply(genes, mousify)
  markers$official.gene.symbol <- genes[markers$official.gene.symbol]
  celllist <- unique(markers$cell.type)
  for (i in 1:length(celllist)){
    if(i == 1){
      listofgenesets <- GeneSet(markers[markers$cell.type == celllist[i],]$official.gene.symbol, setName=celllist[i])
    }
    else {
      listofgenesets <- append(listofgenesets,GeneSet(markers[markers$cell.type == celllist[i],]$official.gene.symbol, setName=celllist[i]))
    }
  }
  listofgenesets <- append(listofgenesets,GeneSet(c("Epcam","Krt8","Krt18","Pymt","Gfpluc","Cdh1"), setName="Cancer cells"))
  fullset <- GeneSetCollection(listofgenesets)
  
  return(fullset)
  
}

getallcelltypes <- function(){
  
  markers <- getallmarkers()
  
  unqcelltypes <- unique(markers$cell.type)
  
  return(unqcelltypes)
  
}

getallmarkers <- function(){
  
  markers <- read.csv("PanglaoDB_markers_27_Mar_2020.tsv.gz", sep="\t")
  markers <- markers[markers$species != "Hs",]
  markers <- data.frame(markers)
  markers <- markers[grepl("Liver|Immune system|Vasculature", markers$organ),]
  
  return(markers)
  
}

getallcolors <- function(){
  
  markers <- getallmarkers()
  
  colpal <- createPalette(length(unique(markers$cell.type))+1,  c("#ff0000", "#00ff00", "#0000ff"))
  
  colpal <- setNames(colpal, c(unique(markers$cell.type),"Cancer cells"))
  
  return(colpal)
  
}