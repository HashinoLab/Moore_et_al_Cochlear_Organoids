#HPC_GSEA
#CDH7 Bi-allelic knockout Day20 Prosensory GSEA
library(tidyverse, quietly = TRUE)
library(scRNAseq, quietly = TRUE)
library(matrixStats, quietly = TRUE)
library(iDEA, quietly = TRUE)
library(zingeR, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(BiocParallel, quietly = TRUE)
library(data.table)
register(MulticoreParam(8))
print("Cores = Registered")

#The following script was run on Indiana University's SLATE HPC

day109 <- readRDS(file = "C:/Users/.../day109AutoClustRd.rds")

d109.hc.se <- SummarizedExperiment(as.matrix(day109.hc@assays$RNA@counts),
                                   colData = day109.hc@meta.data[c("Condition", "Day", "seurat_clusters")])

saveRDS(d109.hc.se, file = "C:/Users/.../Day109/GSEA/D109_HC_SE.rds")

print("SCE read in")

#Filter to only look at highly expressed genes
filter <- rowSums(assay(d109.hc.se)) > 5
d109.hc.se <- d109.hc.se[filter,]
assay(d109.hc.se) %>% log1p %>% rowVars -> vars #This line is cursed
names(vars) <- rownames(d109.hc.se)
vars <- sort(vars, decreasing = T)

#Subset to only look at the top 15000 genes
d109.hc.se <- d109.hc.se[names(vars[1:15000]),]
assayNames(d109.hc.se)[1] <- "counts"

#Construct Inputs for DESeq2 Test
d109.colData <- data.frame(Condition = colData(d109.hc.se)$Condition)
iCounts <- assay(d109.hc.se, i = "counts")

#Set factor levels for variables of interest
d109.colData$Condition <- as.factor(d109.colData$Condition)
#Ensure that the control group is the first factor level
d109.colData$Condition <- relevel(d109.colData$Condition, ref = "CHIR") #Remember to change this line to appropriate factor level!!!

design <- model.matrix(~d109.colData$Condition)

dse <- DESeqDataSetFromMatrix(countData = iCounts, colData = d109.colData, design = ~ Condition)
weights <- zingeR::zeroWeightsLS(counts = iCounts, 
                                 design = design,
                                 maxit = 500, normalization = "DESeq2_poscounts",
                                 colData = d109.colData,
                                 designFormula = ~ Condition, 
                                 verbose = TRUE)

assays(dse)[["weights"]] <- weights

rm(d109.hc.se) #Remove extraneous objects before parallelization begins
dse <- DESeq2::DESeq(dse,
                     test = "Wald", 
                     sfType = "poscounts", 
                     useT = TRUE, 
                     betaPrior = TRUE, 
                     modelMatrixType="expanded",
                     parallel = TRUE)

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
write.csv(full.results.df, file = "/N/slate/.../D109_HC_DESeq2.csv")

#Read in the genesets of interest
completeGS <- readRDS(file = "/N/slate/.../CompleteSets.rds")
de.summary <- data.frame(log2FoldChange = res$log2FoldChange,
                         lfcSE2 = res$lfcSE^2,
                         row.names = rownames(res))
print("Removing extraneous variables from the environment")
rm(res, full.results.df, full.results, d109.colData, design, weights, dse)

#GSEA for upregulated genes
d70.idea <- CreateiDEAObject(de.summary[de.summary$log2FoldChange > 0, ], completeGS, num_core = 8)
d70.idea <- iDEA.fit(d70.idea, modelVariant = F)
d70.idea <- iDEA.louis(d70.idea)
write.csv(d70.idea@gsea, file = "/N/slate/.../D109_HC_Upreg_GSEA.csv")
rm(d70.idea)

res <- results(dse, cooksCutoff = F, parallel = T)
res$log2FoldChange %>% summary()
res$lfcSE %>% summary()

completeGS <- readRDS(file = "/N/slate/.../CompleteSets.rds")
full.results <- fread(file = "/N/slate/.../D109_HC_DESeq2.csv")[-1,]

de.summary <- data.frame(log2FoldChange = full.results$log2FoldChange,
                         lfcSE2 = full.results$logFC_SE^2,
                         row.names = full.results$gene)

d70.idea <- CreateiDEAObject(de.summary, completeGS, num_core = 8)
rm(completeGS, full.results, de.summary)
d70.idea <- iDEA.fit(d70.idea, modelVariant = F)
d70.idea <- iDEA.louis(d70.idea)
write.csv(d70.idea@gsea, file = "/N/slate/.../D109_HC_Bidirectional_GSEA.csv")
