#HPC_GSEA
#Day20 Prosensory GSEA
library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(zinbwave, quietly = TRUE)
library(scRNAseq, quietly = TRUE)
library(matrixStats, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(biomaRt, quietly = TRUE)
library(Seurat, quietly = TRUE)
library(edgeR, quietly = TRUE)
library(iDEA, quietly = TRUE)
library(zingeR, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(BiocParallel, quietly = TRUE)

#The following script was run on Indiana University's SLATE HPC

register(MulticoreParam(8))
print("Cores = Registered")
# Read in the data
d20.ps <- readRDS(file = "/N/slate/.../Prosensory_noPur.rds")
print("Seurat object read in")

#Create SCE(?) object
d20.sce <- SummarizedExperiment(as.matrix(d20.ps@assays$RNA@counts),
                                colData = d20.ps@meta.data)
rm(d20.ps)
print("SCE created, Seurat removed")

#Filter to only look at highly expressed genes
filter <- rowSums(assay(d20.sce)) > 5
d20.sce <- d20.sce[filter,]
assay(d20.sce) %>% log1p %>% rowVars -> vars #This line is cursed
names(vars) <- rownames(d20.sce)
vars <- sort(vars, decreasing = T)

#Subset to only look at the top 15000 genes
d20.sce <- d20.sce[names(vars[1:15000]),]
assayNames(d20.sce)[1] <- "counts"


#Construct Inputs for DESeq2 Test
d20.colData <- data.frame(Condition = as.factor(colData(d20.sce)$Condition))
iCounts <- assay(d20.sce, i = "counts")

d20.colData$Condition <- relevel(d20.colData$Condition, ref = "IWP2")

design <- model.matrix(~d20.colData$Condition)

dse <- DESeqDataSetFromMatrix(countData = iCounts, colData = d20.colData, design = ~ Condition)
weights <- zingeR::zeroWeightsLS(counts = iCounts, 
                                 design = design,
                                 maxit = 500, normalization = "DESeq2_poscounts",
                                 colData = d20.colData,
                                 designFormula = ~ Condition, 
                                 verbose = TRUE)

assays(dse)[["weights"]] <- weights

dse <- DESeq2::DESeq(dse,
                     test = "Wald", 
                     sfType = "poscounts", 
                     useT = TRUE, 
                     betaPrior = TRUE, 
                     modelMatrixType="expanded",
                     parallel = TRUE)
# dse <- DESeq2::estimateSizeFactors(dse, type="poscounts")
# dse <- estimateDispersions(dse)
# dse <- nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-(length(unique(d20.colData$Condition))-1) )

res <- results(dse, cooksCutoff = F, parallel = T)
res$log2FoldChange %>% summary()
res$lfcSE %>% summary()

#Save DESeq results
full.results <- results(dse, parallel = T)
full.results.df <- data.frame(gene = row.names(full.results),
                              log2FoldChange = full.results$log2FoldChange,
                              pvalue = full.results$pvalue,
                              padj = full.results$padj,
                              baseMean = full.results$baseMean,
                              logFC_SE = full.results$lfcSE)
write.csv(full.results.df, file = "/N/slate/.../Prosensory_DESeq2_First.csv")

#GSEA Testing
completeGS <- readRDS(file = "/N/slate/.../CompleteSets.rds")

de.summary <- data.frame(log2FoldChange = res$log2FoldChange,
                         lfcSE2 = res$logFC_SE^2,
                         row.names = res$gene)
rm(res)

d20.idea <- CreateiDEAObject(de.summary[de.summary$log2FoldChange > 0, ], completeGS, num_core = 8)
d20.idea <- iDEA.fit(d20.idea, modelVariant = F)
d20.idea <- iDEA.louis(d20.idea)
write.csv(d20.idea@gsea, file = "/N/slate/.../Upreg_NB_PS_CustomGS_GSEA.csv")

rm(d20.idea)

d20.idea <- CreateiDEAObject(de.summary[de.summary$log2FoldChange < 0, ], completeGS, num_core = 8)
d20.idea <- iDEA.fit(d20.idea, modelVariant = F)
d20.idea <- iDEA.louis(d20.idea)
write.csv(d20.idea@gsea, file = "/N/slate/.../Downreg_NB_PS_CustomGS_GSEA.csv")
