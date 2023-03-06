# 10x-scRNA-Analysis
10X Single Cell Analysis

This set of scripts analyzes 10x genomics 3' RNA-Seq data (that has been hashtagged with totalseq product from biolegends to allow for sample multiplex on 10x machine) that has been processed through cellranger.

Currently the scripts use cell markers defined in pangalo db to identify cell type for each cluster.

Input: 
- Location of Gene Count Matrix generated by 10x (2nd column is a unique sample name) - samplelist7.txt
- Genes that need to visualized on UMAP - markers-UMAP.txt
- paramters.yaml - Input file for conitions, working directory (i.e. output directory for results), Input text file (file with the location of gene count matrix) and threshold to us for Minimum Number of UMI (nCount_RNA), Minimum Number of Genes (nFeature_RNA), max percent mitocondrial fraction (percent.mt) and Max Number of UMI (max_nCount_RNA)

          Example parameters.yaml
          workingdir: /home/sgupta/RClass/ESR1-Analysis/
          Input: samplelist4.txt
          Condition:
           - Day1
           - Day35
          params:
           Day1:
            nCount_RNA: 1000
            nFeature_RNA: 600
            percent.mt: 20
           Day35:
            nCount_RNA: 600
            nFeature_RNA: 400
            percent.mt: 20
            max_nCount_RNA: 100000


