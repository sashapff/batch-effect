# batch-effect

### Goal

Reanalyse RNA-seq data by [Gilad & Mizrahi-Man, 2015](https://f1000research.com/articles/4-121) to reveal and correct batch effect.

### Objectives

* using already calculated counts, obtain genes differential expression (FPKM) RNA-seq data in Mus musculus and Homo sapiens
* analyze differential expression data to reveal 
* correct batch effect with different methods and compare the results

### Software used

* Cufflinks 2.2.1
* R 4.0.3
* ComBat 3.38.0
* Harmony 0.1

### Results

All methods and results can be found in R notebooks *BeforeCorrection.html*, *Combat.html*, *Harmony.html* and R scripts *BeforeCorrection.Rmd*, *Combat.Rmd*, *Harmony.Rmd*.

Example:

load input

```
datasets = as.data.frame(scan("Stanford_datasets.txt",list(setname="",seqBatch="",species="",tissue=""),sep="\t"))

fpkmMat <- as.matrix(read.table('Stanford_datasets_fpkmMat.txt',header=FALSE,sep='\t'))

rawCounts <- as.matrix(read.table('Stanford_datasets_rawCountsMat.txt',header=FALSE,sep='\t'))
colnames(rawCounts) <- datasets$setname
geneDetails <- as.data.frame(scan("ortholog_GC_table.txt",skip=1,list(mouse_name="",mouse_GC = 0.0,human_name = "",human_GC=0.0)))

rownames(rawCounts) <- geneDetails$human_name 
rownames(datasets) <- datasets$setname
rownames(geneDetails) <- geneDetails$human_name 
rownames(fpkmMat) <- geneDetails$human_name 
```

filter out lower 30% + mitochondrial genes

```
rowSums = apply(rawCounts,1,function(x) sum(x))
quantile(rowSums,probs = 0.3)

filter <- apply(rawCounts,1,function(x) sum(x)>2947.9 )
mt <- grep("mt-",geneDetails$mouse_name)
filteredNames <- setdiff(rownames(rawCounts[filter,]),rownames(rawCounts[mt,])) 
filteredFpkm <- fpkmMat[filteredNames,]
```

perform log transformation

```
logTransformed_fpkmMat = log2(filteredFpkm+1)
colnames(logTransformed_fpkmMat) <- datasets$setname
```

plot correlation heatmaps, pearson correlation, complete linkage

```
library(pheatmap)
pheatmap(cor(logTransformed_fpkmMat))
```

![heatmap](https://github.com/sashapff/batch-effect/blob/main/heatmap_example.png)

### References

1.	Gilad, Y. & Mizrahi-Man, O. A reanalysis of mouse ENCODE comparative gene expression data. F1000Research 4, 121 (2015).
2.	Lin, S. et al. Comparison of the transcriptional landscapes between human and mouse tissues. Proc. Natl. Acad. Sci. U. S. A. 111, 17224–17229 (2014).
3.	Johnson, W., Li, C. & Rabinovic, A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostat. Oxf. Engl. 8, 118–27 (2007).
4.	Abstract. https://www.bu.edu/jlab/wp-assets/ComBat/Abstract.html.
5.	immunogenomics/harmony. (Raychaudhuri Lab, 2020).
6.	Trapnell, C. et al. Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks. Nat. Protoc. 7, 562–578 (2012).
