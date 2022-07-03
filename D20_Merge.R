#JingIntegration
library(dplyr)
library(ggplot2)
library(Seurat)
library(sctransform)
library(tidyverse)


## ----------------------------- Load in Data  ---------------------------------
folder <- "C:/ExampleFolder"

jing.list <- list()
#Initialize list of Seurat objects where ConA = JN1, ConB = JN2, and ConC = JN3
ConA <- paste0(folder, "P2G_D20_JN1/outs/filtered_feature_bc_matrix")
jing.list[[1]] <- Read10X(ConA) %>%
  CreateSeuratObject()
ConB <- paste0(folder, "P2G_D20_JN2")
jing.list[[2]] <- Read10X(ConB) %>%
  CreateSeuratObject()
ConC <- paste0(folder, "P2G_D20_JN3/outs/filtered_feature_bc_matrix")
jing.list[[3]] <- Read10X(ConC) %>%
  CreateSeuratObject()


#Annotate conditions
jing.list[[1]][["Condition"]] <- "CHIR"
jing.list[[2]][["Condition"]] <- "CHIR+PUR"
jing.list[[3]][["Condition"]] <- "IWP2"

## ----------------------------- Filtering  -----------------------------------

#Filters for (log-transformed) mitochondrial percentage genes
filterPlot <- function(obj) {
  obj[["mt.per"]] <- log1p(PercentageFeatureSet(obj, pattern = "^MT-"))
  mt.per.vec <- obj$mt.per
  obj <- subset(obj, mt.per < mean(mt.per.vec) + 3*sd(mt.per.vec) &
                  mt.per > mean(mt.per.vec) - 3*sd(mt.per.vec))
  return(obj)
}

#Filtering for housekeeping gene expression
#Appends metadata column "hk" that stores the log normalized expression
#of a selected housekeeping gene
#if subset = TRUE, will subset the seurat object discarding cells ~std.dev~
#standard deviations away from the mean 
housekeeping.filter <- function(object, housekeeping.gene = NULL, 
                                subset = F, std.dev = 2) {
  
  #Sanitize input
  housekeeping.gene <- as.character(housekeeping.gene)
  
  #creates a vector containing the reads found in each cell
  expression.vector <- object@assays$RNA@counts[housekeeping.gene, ] %>%
    
    #log normalize the expression.vector
    #log plus one is used to avoid negative infinities
    log1p()
  
  #Add metadata column to the Seurat object to facilitate subsetting
  object <- AddMetaData(object, metadata = expression.vector,
                        col.name = "hk")
  
  # Determines whether or not the object should be subset such that
  # cells greater or less than ~std.dev~ standard deviations are discarded
  if (subset){
    object <- subset(object, hk > mean(object$hk) - std.dev*sd(object$hk) &
                       hk < mean(object$hk) + std.dev*sd(object$hk))
  }
  
  return(object)
}

#Filters each dataset in jing.list according to mitochondrial percentage,
#gives number of cells filtered for each dataset
#cell count to counts.list. Each group is downsampled to sample.size <- min(counts.list)
counts.list.post <- counts.list.pre <- numeric()
for (i in 1:3) {
  counts.list.pre[i] <- length(Cells(jing.list[[i]]))
  jing.list[[i]] <- filterPlot(jing.list[[i]])
  counts.list.post[i] <- length(Cells(jing.list[[i]]))
}
diff <- counts.list.pre - counts.list.post
diff

#same thing but for RPL27
#Need to maintain consistent pre-processing throughout all datasets
counts.list.post.2 <- numeric()
for (i in 1:3) {
  jing.list[[i]] <- housekeeping.filter(jing.list[[i]], housekeeping.gene = "RPL27", subset = T, std.dev = 2)
  counts.list.post.2[i] <- length(Cells(jing.list[[i]]))
}
diff.2 <- counts.list.post - counts.list.post.2
diff.2

#Merge 
day20 <- merge(jing.list[[1]], jing.list[[2]]) %>%
  merge(jing.list[[3]])

## ----------------------------- Seurat workflow -------------------------------

#SCTransform
day20 <- SCTransform(day20, vars.to.regress = "mt.per") %>%
  RunPCA()

#Clustering
day20.auto <- FindNeighbors(day20, dims = 1:6, k.param = 1290) %>%
  FindClusters(resolution = 0.12) %>%
  RunUMAP(dims = 1:6, n.neighbors = 1290)

#Calculating ICVI score
clusters <- as.integer(Idents(day20.auto))
sil.score <- intCriteria(day20@tools$LoadEmbeddings,
                         part = clusters,
                         crit = "Sil")
sil.score


## ------------------------------ Plotting ------------------------------
sparse <- function() {
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank()
  )
}

save.plot <- function(plot, name){
  folder <- "C:/Users/.../Results/"
  ggsave(filename = paste0(folder, "_", name, ".jpg"),
         plot = plot,
         device = jpeg(),
         units = "in",
         height = 8,
         width = 8,
         dpi = 200)
  dev.off()
}

save.plot(DimPlot(day20.auto, pt.size = 1.2, order = rev(c(0,1,2,3))),
          "DimTest")

save.plot(DimPlot(day20.auto, cells = sample(Cells(day20.auto), 15000), pt.size = 1.2, group.by = "Condition"),
          "ConTest")

gene.list <- c("SULF1", "EDN3", "JAG1",
               "PAX8", "EPCAM", "FBXO2", "PAX2", "CDH1", "CDH2", "ESPN", "COL9A2",
               "NEUROG1", "NEUROD1",
               "UBE2C", "CDK1",
               "OTX1", "OTX2", "HES1", "PTCH1", "OTOL1", "LFNG", "GLI1",
               "DLX5", "DLX6", "GBX2", "MSX1", "HMX1", "HMX2", "OC90")

feature.cols = c("lightgrey" , "#002385")

for (i in seq_len(length(gene.list))) {
  save.plot(FeaturePlot(day20.auto,
                        cells = sample(Cells(day20.auto), 15000),
                        features = gene.list[i],
                        pt.size = 1.2, cols = feature.cols) +
              NoLegend() +
              sparse(),
            name = paste0(gene.list[i], "_FeaturePlot"))
}


## ------------------------------ Stacked Histogram ----------------------------
day20.auto$cluster <- day20.auto@active.ident

IWP2 <- "#B856D7"
CHIR <- "#55A0FB"
PUR <- "#3ce7a6"

new.idents <- c("Prosensory", "Prosensory (IWP2)", "Neuroblasts", "Cycling")
names(new.idents) <- levels(day20.auto@active.ident)
day20.t <- RenameIdents(day20.auto, new.idents)
day20.auto <- day20.t
rm(day20.t)

plot.data <- day20.auto@meta.data %>%
  group_by(cluster, Condition) %>%
  summarise(Count = n()) %>%
  mutate(Prop = Count / sum(Count))

plot.data$C2 <- con.switch[plot.data$Condition] %>% unlist()

con.switch <- list(
  CHIR = "CTRL",
  `CHIR+PUR` = "PUR",
  IWP2 = "PUR+IWP2"
)

f2E <- ggplot(plot.data, aes(x = cluster, y = Prop, fill = C2)) +
  geom_col() +
  labs(y = "Cluster Composition") +
  scale_color_manual(values = c(IWP2, CHIR, PUR)) +
  scale_fill_manual(values = c( CHIR, PUR, IWP2)) +
  guides(fill = guide_legend(direction = "vertical", ncol = 3)) +
  theme_cowplot(font_size = 22) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle =45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center")

f2E
save.plot(f2E, "Prop_Bar_Recolor")
save.eps(f2E, "Prop_Bar_Recolor")


save.plot <- function(plot, name){
  folder <- "C:/Users/.../Results/"
  ggsave(filename = paste0(folder, "_", name, ".jpg"),
         plot = plot,
         device = jpeg(),
         units = "in",
         height = 8,
         width = 6,
         dpi = 300)
}

save.eps <- function(plot, name){
  folder <- "C:/Users/.../Results/"
  ggsave(filename = paste0(folder, "_", name, ".eps"),
         plot = plot,
         device = "eps",
         units = "in",
         height = 8,
         width = 6,
         dpi = 300)
  dev.off()
}

## ------------------------------ Prosensory Analysis  -------------------------

prosense <- subset(day20, idents = levels(day20)[1:2])
prosense <- SCTransform(prosense, vars.to.regress = "mt.per") %>%
  FindNeighbors(dims = 1:6, k.param = 1360) %>%
  FindClusters(resolution = 0.24) %>%
  RunUMAP(dims = 1:6, n.neighbors = 680)

umapPlot <- DimPlot(prosense, group.by = "Condition", pt.size = 1.6, cols = c(CHIR, PUR, IWP2), cells = sample(Cells(prosense), 16000)) +
  NoLegend() + sparse()


## ------------------------------ Split Violin Plots  -------------------------

split.prosense <- prosense
split.prosense@active.ident <- factor(as.factor(split.prosense$Condition), levels = c("IWP2", "CHIR", "CHIR+PUR"))
split.prosense <- subset(split.prosense, idents = c("CHIR", "IWP2"))
split.prosense@active.ident <- as.factor(rep("Day 20", length(Cells(split.prosense))))
split.prosense$Condition <- factor(as.factor(split.prosense$Condition), levels = c("IWP2", "CHIR"))

gene.list.combined <- c("DLX5", "MSX1", "GPR155", 'ACSL4', "CXCR4", 
                        "TAGLN", "UBE2C", "JAG1", "OTX1", "NR2F2", "EDN3", 
                        "RSPO3", "GAS1", "PTCH1", "SULF1", "LRP2", 'EPCAM', 'FBXO2')


for(i in 1:length(gene.list.combined)){
  v.plot <- VlnPlot(split.prosense,
                    features = gene.list.combined[i],
                    split.by = "Condition",
                    split.plot = T,
                    cols = c(IWP2, CHIR),
                    pt.size = 0,
                    log = F,
                    y.max = y.iter)
  v.plot$layers[[1]]$aes_params$size <- 1.25
  v.plot <- v.plot +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    xlab("") +
    sparse() +
    NoLegend()
  name <- paste0(gene.list.combined[i], "_Day20_Vln")
  save.png(v.plot, name, y.iter)
  save.eps(v.plot, name, y.iter)
}

## ------------------------------ Feature Plots  ------------------------------

for (i in seq_len(length(gene.list.combined))) {
  save.plot(FeaturePlot(prosense,
                        cells = sample(Cells(prosense), 15000),
                        features = gene.list[i],
                        pt.size = 1.2, cols = feature.cols) +
              NoLegend() +
              sparse(),
            name = paste0(gene.list[i], "_FeaturePlot"))
}


## ------------------------------ Volcano Plot  -------------------------

library(EnhancedVolcano)
library(data.table)

#Define Colors
IWP2 <- "#B856D7"
CHIR <- "#55A0FB"

deseq.20 <- fread(file = "/Users/....Prosensory_DESeq2_First.csv")
lab.gen.20<- c("NR2F1", "NR2F2", "GAS1", "RSPO3", "EDN3", "LRP2", "SULF1", "PTCH1", "OTX1",
               "JAG1", "ASCL4", "DLX5", "TAGLN", "MSX1", "GPR15", "CXCR4")
rownames(deseq.20) <- deseq.20$gene
logfc.threshold <- 1
p.threshold = 1e-10

#This creates new entries for the genes of interest. This means that they're plotted twice,
#which prevents them from being buried underneath genes we're not interested in
#Also, rownames can't be duplicated, so a 1 is automatically appended to the end of each gene
deseq.20 <- rbind(deseq.20, deseq.20[(rownames(deseq.20) %in% lab.gen.109), ])

#This creates a long character vector composed of the color for each gene. Each entry corresponds to one gene. 
color.key <- ifelse(deseq.20$log2FoldChange > logfc.threshold & deseq.20$padj < p.threshold, CHIR,
                    ifelse(deseq.20$log2FoldChange < -logfc.threshold & deseq.20$padj < p.threshold, IWP2,
                           "lightgrey"))

#This sets the colors for the specific genes you want to highlight
color.key[(rownames(deseq.20) %in% c(lab.gen.109, paste0(lab.gen.109, "1"))) & deseq.20$log2FoldChange < -0.5] <- "darkred"
color.key[(rownames(deseq.20) %in% c(lab.gen.109, paste0(lab.gen.109, "1"))) & deseq.20$log2FoldChange > -1] <- "#3E3C9A"

#Name your color.key. Enhanced volcano doesn't work unless the color.key is named.
#Names are displayed as a plot legend.
names(color.key)[color.key == "#B856D7"] <- "IWP2"
names(color.key)[color.key == "#55A0FB"] <- "CHIR"
names(color.key)[color.key == "lightgrey"] <- "Below_Threshold"
names(color.key)[color.key == "darkred"] <- "Selected"
names(color.key)[color.key == "#3E3C9A"] <- "SelectedCHIR"


DESeq2.volcano.selectGenes <- EnhancedVolcano(deseq.20,
                                              lab = rownames(deseq.20),
                                              x = "log2FoldChange",
                                              y = "padj",
                                              title = "Day 20 Prosensory: IWP2 vs CHIR",
                                              pCutoff = 1e-10,
                                              FCcutoff = 1, #This sets the cutoff outside the bounds of the plot, so we don't see them
                                              pointSize = 4,
                                              colCustom = color.key,
                                              selectLab = "lab.gen.20",
                                              legendPosition = "right",
                                              labSize = 5,
                                              drawConnectors = T,
                                              arrowheads = FALSE,
                                              maxoverlapsConnectors = Inf) + 
  xlim(-4.5, 4.5) + 
  NoLegend()
DESeq2.volcano.selectGenes

## ------------------------------ Bubble Plot  -------------------------
library(data.table)
library(iDEA)
library(tidyverse)
library(stringr)
data("humanGeneSets")
data("humanGeneSetsInfo")


## Assemble GeneSetInfo for custom gene sets
custom.sets <- readRDS("/Users/alexander/Dropbox/Mac/Documents/HashinoLab/GSEA_Sets/CustomSets.rds")
cs.names <- names(custom.sets)

#Get relative positions for each gene set as they were built
idx.hm <- str_extract(cs.names, "[^_]+") %>% grep(pattern = "^HALLMARK") #Gets location of hallmark sets
idx.pos <- str_extract(cs.names, "[^_]+") %>% grep(pattern = "^chr") #Gets location of positional sets
idx.pos <- c(idx.pos, 312) #Ad Hoc correction for the Mitochondrial positional gene set 
idx.tf <- seq_len(length(cs.names))[!seq_len(length(cs.names)) %in% c(idx.hm, idx.pos)] #This is dogshit programming, but it works

#This function is specific to the CustomSets ordering
constructGeneInfo <- function(catLabels = c("hm", "c1", "c3")) {
  catNameMap <- vector(mode = "character", length = 851)
  catNameMap[idx.hm] <- catLabels[1]
  catNameMap[idx.pos] <- catLabels[2]
  catNameMap[idx.tf] <- catLabels[3]
  return(catNameMap)
}

customGeneSetsInfo <- data.frame(catName = constructGeneInfo(c("hm", "c1", "c3")),
                                 subcatName = constructGeneInfo(c("hm", "positional", "gtrd")),
                                 gsetName = constructGeneInfo(c("hm.hm", "c1.positional", "c3.gtrd")),
                                 gsetBioName = constructGeneInfo(c("Hallmark", "Positional", "Transcription factors")),
                                 gset = cs.names)

humanGeneSetsInfo <- rbind(humanGeneSetsInfo, customGeneSetsInfo)

## Create data frame for bubble plot, Combining Up+Downregulated genes into one chart


includedCats <- c("GO BIOLOGICAL PROCESS",
                  "GO MOLECULAR FUNCTION",
                  "GO CELLULAR COMPONENT",
                  "REACTOME",
                  "KEGG",
                  "PID",
                  "HALLMARK",
                  "POSITIONAL",
                  "TRANSCRIPTION FACTORS")

caseCats <- c("GO biological process",
              "GO molecular function",
              "GO cellular component",
              "Reactome",
              "KEGG",
              "PID",
              "Hallmark",
              "Positional",
              "Transcription factors")

plotdata<- fread("/Users/alexander/Dropbox/Mac/Documents/HashinoLab/NatBiotech/Day_20/Upreg_NB_PS_CustomGS_GSEA.csv")
bp.up <- data.table::merge.data.table(plotdata, humanGeneSetsInfo, by.x = "annot_id", by.y = "gset") #I think this is easier in SQL
bp.up$Category <- droplevels(bp.up$gsetBioName)

bp.up <- bp.up[toupper(bp.up$Category) %in% includedCats, ] #subset data, removing all geneset categories not in includedCats
bp.up$Category <- droplevels(bp.up$Category) %>% factor(caseCats) # Drop unused levels, reorder factor so that plot legend matches x position
bp.up <- bp.up[order(match(toupper(bp.up$Category), includedCats)), ] #Organizes the data so that categories are grouped on the x axis
bp.up$IDNum <- row.names(bp.up) %>% as.numeric() #IDNum will define x position on bubble plot
bp.up$Log10_Pvalue_Louis <- -1*log10(bp.up$pvalue_louis)

plotdata<- fread("/Users/alexander/Dropbox/Mac/Documents/HashinoLab/NatBiotech/Day_20/Downreg_NB_PS_CustomGS_GSEA.csv")
bp.down <- data.table::merge.data.table(plotdata, humanGeneSetsInfo, by.x = "annot_id", by.y = "gset") #I think this is easier in SQL
bp.down$Category <- droplevels(bp.down$gsetBioName)

bp.down <- bp.down[toupper(bp.down$Category) %in% includedCats, ] #subset data, removing all geneset categories not in includedCats
bp.down$Category <- droplevels(bp.down$Category) %>% factor(caseCats) # Drop unused levels, reorder factor so that plot legend matches x position
bp.down <- bp.down[order(match(toupper(bp.down$Category), includedCats)), ] #Organizes the data so that categories are grouped on the x axis
bp.down$IDNum <- row.names(bp.down) %>% as.numeric() #IDNum will define x position on bubble plot
bp.down$Log10_Pvalue_Louis <- +1*log10(bp.down$pvalue_louis)

bp.up <- bp.up[, -(2:7)]
bp.down <- bp.down[ , -(2:7)]

#Ensure that different gene sets scale the same on the x-axis for both up and down 
# regulated genes
floor <- 0
for(i in seq_len(length(caseCats))) {
  up.rows <- which(bp.up$Category == caseCats[i])
  down.rows <- which(bp.down$Category == caseCats[i])
  if (length(up.rows) > length(down.rows)) {
    bp.down[down.rows, "IDNum"] <- (bp.down[down.rows, "IDNum"] - floor) * length(up.rows)/length(down.rows) + floor
    next.rows <- seq_len(nrow(bp.down))[seq_len(nrow(bp.down)) > max(down.rows)]
    bp.down[next.rows, "IDNum"] <- bp.down[next.rows, "IDNum"] + length(up.rows) - length(down.rows) 
    floor <- min(bp.down[next.rows, "IDNum"])
  } else if (length(up.rows) < length(down.rows)) {
    bp.up[up.rows, "IDNum"] <- (bp.up[up.rows, "IDNum"] - floor) * length(down.rows)/length(up.rows) + floor
    next.rows <- seq_len(nrow(bp.up))[seq_len(nrow(bp.up)) > max(up.rows)]
    bp.up[next.rows, "IDNum"] <- bp.up[next.rows, "IDNum"] + length(down.rows) - length(up.rows) 
    floor <- min(bp.up[next.rows, "IDNum"])
  }
}

#Combine Up and Down regulated data
bp.data <- rbind(bp.up, bp.down)

SigSelectUp <- c("KDM7A_TARGET_GENES",
                 "GO_CHROMATIN_BINDING",
                 "MORC2_TARGET_GENES",
                 "RAG1_TARGET_GENES",
                 "GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION",
                 "ZFHX3_TARGET_GENES",
                 "GCM2_TARGET_GENES",
                 "GO_CHROMATIN",
                 "SETD7_TARGET_GENES",
                 "ZNF711_TARGET_GENES",
                 "MZF1_TARGET_GENES",
                 "GO_CHROMATIN_MODIFICATION",
                 "GO_TRANSCRIPTION_COACTIVATOR_ACTIVITY",
                 "PHF2_TARGET_GENES",
                 "GTF2E2_TARGET_GENES",
                 "ZNF423_TARGET_GENES",
                 "BARHL1_TARGET_GENES"
)

SigSelectDown <- c(
  "SETD7_TARGET_GENES",
  "DLX6_TARGET_GENES",
  "GO_REGULATION_OF_WNT_SIGNALING_PATHWAY",
  "GO_POSITIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY",
  "GO_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY",
  "GO_POSITIVE_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY",
  "GO_WNT_SIGNALING_PATHWAY",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "WNT_SIGNALING",
  "KEGG_WNT_SIGNALING_PATHWAY")



#%%%%%%%% String Manipulation Functions for bubble plot labels %%%%%%%%%%%%%%%%%

#To ensure genes in GTRD set remain capitalized
target.genes <- sapply(customGeneSetsInfo$gset[customGeneSetsInfo$catName == "c3"],
                       FUN = function(x) str_split(x, pattern = "_")[[1]][1])

CamelCaseConverter <- function(char) {
  new.char <- str_split(char, pattern = "_") %>% unlist()
  if (new.char[1] %in% target.genes){
    new.char[-1] <- sapply(new.char[-1], capitalize)
  } else if (new.char[1] %in% c("GO", "PID", "KEGG", "REACTOME", "HALLMARK")){
    new.char <- sapply(new.char[-1], capitalize)
  } else {
    new.char <- sapply(new.char, capitalize)
  } 
  return(paste(new.char, collapse = " "))
}

capitalize <- function(char) {
  char <- str_replace(tolower(char), "[[:alpha:]]{1}", toupper(str_extract(char, "[[:alpha:]]{1}")))
  if (char %in% c("Tgf", "Gli", "Myc", "Wnt")) { char <- toupper(char) }
  if (char == "Mtorc1") { char <- "mTORC1" } 
  if (char == "activ") { char <- "active" }
  if(char %in% c("In", "To", "By", "Of")) { return( tolower(char) ) } else { return( char ) }
}


SigUp <- bp.data[which((bp.data$annot_id %in% SigSelectUp) & bp.data$Log10_Pvalue_Louis > 0),]
SigUp$Term <- sapply(SigSelectUp[order(match(SigSelectUp, SigUp$annot_id))],
                     CamelCaseConverter)


SigDown <- bp.data[which( (bp.data$annot_id %in% SigSelectDown) & bp.data$Log10_Pvalue_Louis < 0),]
SigDown$Term <- sapply(SigSelectDown[order(match(SigSelectDown, SigDown$annot_id))],
                       CamelCaseConverter)


## ggplot function for bubble plot

library(ggplot2)
library("ggrepel")
library(RColorBrewer)
options(ggrepel.max.overlaps = Inf)


nine.colors <- brewer.pal(12, "Set3")[c(1,3,12,7,4,5,6,8, 10)]


p4 <- ggplot(bp.data, aes(x = IDNum, y = Log10_Pvalue_Louis, color = Category)) +
  geom_point(shape = 19, alpha=1, size = 6) + #aes(fill = Category, shape = 19),  +
  labs(x = "",
       y = expression(paste(bold(-log[10]), bold("("), bolditalic(p), bold("-value)")))) +
  ggtitle(label = expression(paste("Differentially expressed genesets: IWP-2 vs. CHIR prosensory cells at day 20"))) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        #axis.text.x=element_blank(),
        plot.title = element_text(size = 30, face="bold"),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'grey80'),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 40, face = 'bold'),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.title=element_text(size=24,face = 'bold'),
        legend.text=element_text(size=24),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'white'),
        plot.background = element_blank(),
        legend.key = element_rect(color = "transparent", fill = "transparent")) +
  geom_hline(yintercept = -1*log10(0.05), col = 'black',linetype = 2,size=2) +
  geom_hline(yintercept = log10(0.05), col = 'black',linetype = 2,size=2) +
  scale_y_continuous(breaks = c(-20, -10, 0, 10, 20), labels = c(20, 10, 0, 10, 20)) +
  scale_color_manual(values= nine.colors) +
  scale_fill_manual(values = nine.colors)+
  theme(legend.direction = "vertical")+
  theme(legend.position = c(0.2, 0.85))+
  theme(legend.box = "horizontal")+
  theme(legend.title.align = 0) +
  geom_text_repel(
    data = SigUp,
    aes(label = Term),
    col = 'black',
    size = 8,
    nudge_y = 3,
    label.size = 0.75,
    box.padding = 0.005,
    force = 25,
    force_pull = .2,
    nudge_x = 2,
    min.segment.length = 0,
    segment.size = 1,
    max.time = 2) +
  geom_text_repel(
    data = SigDown,
    aes(label = Term),
    col = 'black',
    size = 8,
    nudge_y = -3,
    label.size = 0.75,
    box.padding = 0.005,
    force = 50,
    force_pull = 1,
    nudge_x = 2,
    min.segment.length = 0,
    segment.size = 1)
