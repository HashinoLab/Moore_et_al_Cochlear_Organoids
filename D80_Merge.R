# Day 80 Merge (identical workflow to Day109_Merge)
library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)
library(EnhancedVolcano)

## ----------------------------- Load in Data  ---------------------------------

chir <- Read10X("C:/Users/.../filtered_feature_bc_matrix") %>%
  CreateSeuratObject()
iwp2 <- Read10X("C:/Users/.../filtered_feature_bc_matrix") %>%
  CreateSeuratObject()

#Add metadata
chir$Day <- "80"
iwp2$Day <- "80"

chir$Condition <- "CHIR"
iwp2$Condition <- "IWP2"
chir$mt.per <- PercentageFeatureSet(chir, pattern = "^MT-")
VlnPlot(chir, features = "mt.per")


## -----------------------------  Filter & Merge  ------------------------------
#Filters for (log-transformed) mitochondrial percentage genes
filterPlot <- function(obj) {
  obj[["mt.per"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  mt.per.vec <- obj$mt.per
  obj <- subset(obj, mt.per < 12.5)
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


#GroupDatasets into a list
list.80 <- list(chir, iwp2)

#Filters each dataset in list.80 according to mitochondrial percentage,
#gives number of cells removed for each dataset

#I use apply functions to process the data quickly/cleanly. Here's an explainer in case you're not familiar: https://www.guru99.com/r-apply-sapply-tapply.html
#each column in a seurat object represents one cell, so calling ncol is equivalent to length(Cells(seurat.object))
counts.list.pre <- sapply(list.80, ncol)
list.80 <- lapply(list.80, filterPlot)
counts.list.post <- sapply(list.80, ncol)
diff <- counts.list.pre - counts.list.post
print(diff)


#same thing but filtering based on housekeeping gene expression, I use RPL27
list.80 <- lapply(list.80, housekeeping.filter, housekeeping.gene = "RPL27", subset = T)
counts.list.post.2 <- sapply(list.80, ncol)
diff.2 <- counts.list.post - counts.list.post.2
print(diff.2)

day80 <- merge(list.80[[1]], list.80[[2]])
rm(iwp2, list.80, chir)


## -----------------------------  Seurat workflow  ------------------------------
day80 <- SCTransform(day80, vars.to.regress = "mt.per")
day80 <- RunPCA(day80)


day80 <- FindNeighbors(day80, dims = 1:10, k.param = 510) %>%
  FindClusters(resolution = 1.1) %>%
  RunUMAP(n.neighbors = 510, dims = 1:10)

#Calculating ICVI score
clusters <- as.integer(Idents(day109.hires))
sil.score <- intCriteria(day109@tools$LoadEmbeddings,
                         part = clusters,
                         crit = "Sil")
sil.score

saveRDS(day80, "C:/.../Day80_AutoClustrD.rds")

## -----------------------  Build Gene Table -----------------------------------
library(readxl)
library(data.table)

HCs <- c("LFNG","GRXCR1", "IRX2", "PVALB", "MYO7A", "GFI1",
        "PCP4", "EPCAM", "PCP4", "SOX2", 
        "OTOF", "BCL11B", "ENO2", "LHX3" )

SCs <- c("OTOG", "LGR5")

`Cochlear HCs` <- c("PTGIR", "GRP", "TESC",  "CORO2A", "HES6",
                 "GATA3", "NR2F1", "GNG8", "LMOD3", "B3GALT5",
                  "INSM1", "LMOD3", "SMTN", "EFNA5")

`Vestibular HCs` <- c("ZBBX","DNAH6", "DNAH5", "TEKT1", "TCTEX1D1", "DNAJC5B")

OHCs <- c("STRIP2", "NEUROD6", "IKZF2", "SLC26A5")

IHCs <- c("SLC17A8")

Conflicting <- c("CALB1", "OCM", "KCNJ13", "ENDOD1")


#Create data frame (by making list, ensuring each elemnt has the same length, then
#recombining. There's almost certainly an easier way to do this.)
hc.marker.list <- list(HCs = HCs, `Cochlear HCs` = `Cochlear HCs`,
                       `Vestibular HCs` = `Vestibular HCs`,
                       OHCs = OHCs, IHCs = IHCs,
                       Conflicting = Conflicting, SCs = SCs)

col.max <- lapply(hc.marker.list, length) %>% 
  unlist() %>%
  max()

resize <- function(vec, size) {
  length(vec) <- size
  return(vec)
}

hc.marker.list <- lapply(hc.marker.list, resize, size = col.max)
hc.markers.1 <- data.frame(matrix(unlist(hc.marker.list),
                                  nrow = col.max))
names(hc.markers.1) <- names(hc.marker.list)


#Read in data from excel
hc.markers <- readxl::read_xlsx("C:/Users/.../HC_Markers2.xlsx")
missing.cols <- names(hc.markers.1)[!names(hc.markers.1) %in% names(hc.markers)]
hc.markers[missing.cols] <- NA

missing.cols <- names(hc.markers)[!names(hc.markers) %in% names(hc.markers.1)]
hc.markers.1[missing.cols] <- NA

hc.markers.2 <- rbind(hc.markers, hc.markers.1)
hc.markers <- hc.markers.2

hc.markers[5:6, 1] <- c("EPCAM", "FBXO2")
write.csv(hc.markers, "C:/Users/.../MarkerMasterList.csv")





#Vln Genes
gene.str <- "NR2F1

NR2F2

GATA3

INSM1

ZNF503

FGF8

GNG8

LFNG

FGFR3

LGR5

RPRM

CD164L2

ZBBX

TEKT1

SKOR1

AMPD3

VEPH1

TCTEX1D1

PCDH20

NEUROD6"

genes <- strsplit(gene.str, split = "\n") %>% unlist()
genes <- genes[1+2*(0:20)]

## ------------------------- Violin Plots ------------------------------------
#Define Colors
IWP2 <- "#B856D7"
CHIR <- "#55A0FB"

#Set Folder
folder <- "C:/Users/.../Day80/"

hc.markers <- read.csv("C:/Users/.../MarkerMasterList.csv", header = T)[-1]

#save.plot for Vln
save.plot <- function(plot, name, folder = "C:/Users/.../Day80_New/"){
  ggsave(filename = paste0(folder, name, ".jpg"),
         plot = plot,
         device = jpeg(),
         units = "in",
         height = 8,
         width = 8,
         dpi = 300)
  dev.off()
}

#Subset HairCells for Violin Plots 

split.day80 <- subset(day80, idents = c(5,6))
split.day80@active.ident <- factor(as.factor(split.day80$Condition), levels = c("IWP2", "CHIR"))
split.day80@active.ident <- as.factor(rep("Day 80", length(Cells(split.day80))))
split.day80$Condition <- factor(as.factor(split.day80$Condition), levels = c("IWP2", "CHIR"))

for (y.iter in c(2, 2.25, 2.5, 2.75)){
  for(i in 1:length(genes)){
    v.plot <- VlnPlot(split.day80, 
                      features = genes[i], 
                      split.by = "Condition", 
                      split.plot = T, 
                      cols = c(IWP2, CHIR),
                      pt.size = 0,
                      log = F) 
    v.plot$layers[[1]]$aes_params$size <- 1
    v.plot <- v.plot +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.line = element_line(size = 1), axis.ticks.y = element_line(size=1),
            axis.text.y = element_text(size = 14)) + 
      xlab("") +
      #sparse() +
      NoLegend()
    name <- paste0("WithAxis", genes[i], "_Day80_Vln")
    save.png(v.plot, name, y.iter)
    save.eps(v.plot, name, y.iter)
  }
}

save.png <- function(plot, name, y.iter){
  plot.path <- paste0("/Users/alexander/Dropbox/Mac/Documents/HashinoLab/NatBiotech/Plots/D80_requests/",
                      as.character(y.iter),
                      "_", name,
                      ".png")
  ggsave(filename = plot.path,
         plot = plot,
         device = "png",
         units = "in",
         height = 4,
         width = 4,
         dpi = 300)
}

save.eps <- function(plot, name, y.iter){
  plot.path <- paste0("/Users/alexander/Dropbox/Mac/Documents/HashinoLab/NatBiotech/Plots/D80_requests/",
                      as.character(y.iter),
                      "_", name,
                      ".eps")
  ggsave(filename = plot.path,
         plot = plot,
         device = "eps",
         units = "in",
         height = 4,
         width = 4,
         dpi = 300)
}


#Iterate through columns in the data frame
for (col in colnames(hc.markers)){
  #Set and create folder to save files to
  marker.folder <- paste0(folder, "VlnPlots/", col, "/")
  dir.create(marker.folder)
  #Read in genes in dataframe, remove NA valuse
  genes <- hc.markers[[col]] %>%
    discard(is.na)
  #Save genes as VlnPlot
  for (gene in genes) {
    #If the gene passed isn't found within the dataset, an error is thrown which
    #stops the loop
    #So, the plotting function has to be wrapped in a tryCatch Statement
    plot <- tryCatch( 
      { VlnPlot(split.day80,
                features = gene,
                split.by = "Condition",
                split.plot = T,
                cols = c(IWP2, CHIR),
                pt.size = 0) +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()) + 
          xlab("Hair Cells") +
          sparse() +
          NoLegend()
      },
      error = function(cond) {
        message("Hey bud, couldn't find the gene you were looking for")
        message(cond)
        message("\n")
        return(NULL)
      }
    )
    if (!is.null(plot)){
      save.plot(plot, name = paste0(gene, "_Day80_Vln"), folder = marker.folder)
    }
  }
}


## --------------------------   Feature Plots   -------------------------------

#Slightly different dimensions
save.plot <- function(plot, name, folder = "/Users/.../D80/"){
  ggsave(filename = paste0(folder,name, ".eps"),
         plot = plot,
         device = "eps",
         units = "in",
         height = 8,
         width = 8,
         dpi = 300)
}

feature.cols = c("lightgrey" , "#002385")


for (col in colnames(hc.markers)){
  #Set and create folder to save files to
  marker.folder <- paste0(folder, "FeaturePlots/", col, "/")
  dir.create(marker.folder)
  #Read in genes in dataframe, remove NA valuse
  genes <- hc.markers[[col]] %>%
    discard(is.na)
  #Save genes as VlnPlot
  for (gene in genes) {
    plot <- tryCatch( 
      { FeaturePlot(day80,
                    features = gene, 
                    cols = feature.cols,
                    max.cutoff = 2,
                    pt.size = 1.8,
                    order = T) + 
          NoLegend() + 
          sparse()
      },
      error = function(cond) {
        message("Hey bud, couldn't find the gene you were looking for")
        message(cond)
        message("\n")
        return(NULL)
      }
    )
    if (!is.null(plot)){
      save.plot(plot, name = paste0(gene, "_Day80_Feature"))
    }
  }
}


## --------------------------   Blend Plots  -------------------------------  

#Define Colors
IWP2 <- "#B856D7"
CHIR <- "#55A0FB"

blend.plot <- FeaturePlot(day80,
                features = c("DLX5", "GATA3"),
                pt.size = 1.4,
                blend = T,
                order = T,
                blend.threshold = 0.1,
                cols = c("lightgrey",  "#5d55fb", "#B856D7")) + 
  theme(axis.text = element_text(size = 22),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(0.25, "cm"),
        axis.title = element_text(size = 22)
  )

save.plot(blend.plot, name = "blended_BigTxtLegend")

save.plot <- function(plot, name, folder = "C:/Users/...l/Day80/Supplemental/"){
  ggsave(filename = paste0(folder, name, ".jpg"),
         plot = plot,
         device = jpeg(),
         units = "in",
         height = 8,
         width = 32,
         dpi = 300)
  dev.off()
}


## ------------------------------- UMAP Projection ---------------------------
save.eps <- function(plot, name){
  ggsave(filename = paste0("C:/Users/.../Day80_New/UMAP/", name, ".eps"),
         plot = plot,
         device = "eps",
         units = "in",
         height = 8,
         width = 8,
         dpi = 300)
}

save.png <- function(plot, name){
  ggsave(filename = paste0("C:/Users/.../Day80_New/UMAP/", name, ".png"),
         plot = plot,
         device = "png",
         units = "in",
         height = 8,
         width = 8,
         dpi = 300)
}

levels(Idents(day80.reorder)) <- c("\nCHIR SCs\n", "\nIWP2 SCs\n", "\nCHIR Transitional\n", "\nIWP2 Transitional\n", "\nCHIR HCs\n", "\nIWP2 HCs\n",
                                   "\nNC / SCPs\n", "\nEarly DRG\n", "\nMaturing DRG\n", "\nNC Mesenchyme\n", "\nCycling\n")

palette.coo <- c(
  "#53d497", #0
  "#597D8C", #1
  "#6FA397", #2
  "#FEDD72", #3
  "#7fc763", #4
  "#584fdb", #5 
  "#d76cf5", #6
  "#FFC2B4", #7
  "#D4AD80", #8
  "#99D6EA", #9
  "#6798C0") #10

uplot <- DimPlot(day80, cols = palette.coo, label = F, pt.size = 1.6) + NoLegend()


save.eps(uplot, "D80_UMAP_V6")
save.png(uplot, "D80_UMAP_V6")
uplot.sparse <- uplot + sparse()
save.eps(uplot.sparse, "D80_UMAP_V6_sparse")
save.png(uplot.sparse, "D80_UMAP_V6_sparse")

## ----------- Bubble Plot ----------------

## Create data frame for bubble plot, Combining Up+Downregulated genes into one chart

library(data.table)
library(iDEA)
library(tidyverse)
library(stringr)
data("humanGeneSets")
data("humanGeneSetsInfo")

custom.sets <- readRDS("C:/Users/.../CustomSets.rds")
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

plotdata<- fread("C:/Users/.../D80_HC_Upreg_GSEA.csv")
bp.up <- data.table::merge.data.table(plotdata, humanGeneSetsInfo, by.x = "annot_id", by.y = "gset") #I think this is easier in SQL
bp.up$Category <- droplevels(bp.up$gsetBioName)

bp.up <- bp.up[toupper(bp.up$Category) %in% includedCats, ] #subset data, removing all geneset categories not in includedCats
bp.up$Category <- droplevels(bp.up$Category) %>% factor(caseCats) # Drop unused levels, reorder factor so that plot legend matches x position
bp.up <- bp.up[order(match(toupper(bp.up$Category), includedCats)), ] #Organizes the data so that categories are grouped on the x axis
bp.up$IDNum <- row.names(bp.up) %>% as.numeric() #IDNum will define x position on bubble plot
bp.up$Log10_Pvalue_Louis <- -1*log10(bp.up$pvalue_louis)


plotdata<- fread("C:/Users/.../d80_HC_Downreg_GSEA.csv")
bp.down <- data.table::merge.data.table(plotdata, humanGeneSetsInfo, by.x = "annot_id", by.y = "gset") #I think this is easier in SQL
bp.down$Category <- droplevels(bp.down$gsetBioName)

bp.down <- bp.down[toupper(bp.down$Category) %in% includedCats, ] #subset data, removing all geneset categories not in includedCats
bp.down$Category <- droplevels(bp.down$Category) %>% factor(caseCats) # Drop unused levels, reorder factor so that plot legend matches x position
bp.down <- bp.down[order(match(toupper(bp.down$Category), includedCats)), ] #Organizes the data so that categories are grouped on the x axis
bp.down$IDNum <- row.names(bp.down) %>% as.numeric() #IDNum will define x position on bubble plot
bp.down$Log10_Pvalue_Louis <- +1*log10(bp.down$pvalue_louis)

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

bp.data <- rbind(bp.up, bp.down)

SigSelectDown <- c("GO_CILIUM_MORPHOGENESIS",
                   "GO_MICROTUBULE_BASED_PROCESS",
                   "GO_CILIUM",
                   "FOXJ2_TARGET_GENES",
                   "GO_MICROTUBULE_BASED_MOVEMENT",
                   "GO_MOTILE_CILIUM",
                   "GO_MICROTUBULE_BUNDLE_FORMATION",
                   "GO_MICROTUBULE_CYTOSKELETON_ORGANIZATION",
                   "GO_CILIARY_TIP",
                   "GO_PROTEIN_TRANSPORT_ALONG_MICROTUBULE"
)

SigSelectUp <- c("ZNF711_TARGET_GENES",
                 "GO_VOLTAGE_GATED_CATION_CHANNEL_ACTIVITY",
                 "MORC2_TARGET_GENES",
                 "GO_POSTTRANSCRIPTIONAL_REGULATION_OF_GENE_EXPRESSION",
                 "GO_VOLTAGE_GATED_ION_CHANNEL_ACTIVITY",
                 "KDM7A_TARGET_GENES",
                 "GO_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_MEMBRANE",
                 "GO_PROTEIN_LOCALIZATION_TO_MEMBRANE",
                 "GO_IONOTROPIC_GLUTAMATE_RECEPTOR_COMPLEX",
                 "GO_GATED_CHANNEL_ACTIVITY",
                 "GO_VOLTAGE_GATED_POTASSIUM_CHANNEL_ACTIVITY")


SigUp <- bp.data[which((bp.data$annot_id %in% SigSelectUp) & bp.data$Log10_Pvalue_Louis > 0),]
SigUp$Term <- sapply(SigSelectUp[order(match(SigSelectUp, SigUp$annot_id))],
                     CamelCaseConverter)

SigDown <- bp.data[which( (bp.data$annot_id %in% SigSelectDown) & bp.data$Log10_Pvalue_Louis < 0),]
SigDown$Term <- sapply(SigSelectDown[order(match(SigSelectDown, SigDown$annot_id))],
                       CamelCaseConverter)

target.genes <- sapply(customGeneSetsInfo$gset[customGeneSetsInfo$catName == "c3"],
                       FUN = function(x) str_split(x, pattern = "_")[[1]][1])

CamelCaseConverter <- function(char) {
  new.char <- str_split(char, pattern = "_") %>% unlist()
  if (length(grep("^GO", new.char)) > 0) {
    new.char[-grep("^GO", new.char)] <- sapply(new.char[-grep("^GO", new.char)], capitalize)
  } else if (new.char[1] %in% target.genes){
    new.char[-1] <- sapply(new.char[-1], capitalize)
  }
  return(paste(new.char, collapse = " "))
}

capitalize <- function(char) {
  str_replace(tolower(char), "[[:alpha:]]{1}", toupper(str_extract(char, "[[:alpha:]]{1}")))
}


library(ggplot2)
library("ggrepel")
library(RColorBrewer)

nine.colors <- brewer.pal(12, "Set3")[c(1,3,12,7,4,5,6,8, 10)]

p4 <- ggplot(bp.data, aes(x = IDNum, y = Log10_Pvalue_Louis, color = Category)) +
  geom_point(shape = 19, alpha=1, size = 6) + #aes(fill = Category, shape = 19),  +
  labs(x = "",
       y = expression(paste(bold(-log[10]), bold("("), bolditalic(p), bold("-value)")))) +
  ggtitle(label = expression(paste("Differentially expressed genesets: IWP-2 vs. CHIR hair cells at day 109"))) +
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
  geom_hline(yintercept = 1, col = 'black',linetype = 2,size=2) +
  geom_hline(yintercept = -1, col = 'black',linetype = 2,size=2) +
  # guides(color = guide_legend(order = 1,override.aes = list(alpha = 1,size=7)),
  #        size = guide_legend(order = 2,override.aes = list(alpha = 1,shape=21)),
  #        fill = FALSE) + labs(size = "Gene set size")+
  scale_color_manual(values= nine.colors) +
  scale_fill_manual(values = nine.colors)+
  theme(legend.direction = "vertical")+
  theme(legend.position = c(0.2, 0.85))+
  theme(legend.box = "horizontal")+
  theme(legend.title.align = 0) +
  geom_label_repel(
    data = rbind(SigDown, SigUp),
    aes(label = Term),
    col = 'black',
    size = 8,
    nudge_y = 0.6,
    max.overlaps = 1000) +
  ylim(c(-15, 9))


## ---------------------  Volcano Plot -----------------------------------------

library(EnhancedVolcano)
library(dplyr)

#Define Colors
IWP2 <- "#B856D7"
CHIR <- "#55A0FB"

#Define Specific Genes to Highlight
chir.v.iwp2.80 <- fread(file = "/Users/.../D80_HC_DESeq2.csv") %>% as.data.frame()
chir.v.iwp2.80$log2FoldChange <- chir.v.iwp2.80$log2FoldChange * -1
lab.gen.80 <- c("NR2F1", "NR2F2", "GATA3", "INSM1", "ZNF503", "GNG8", "HES6", "TMPRSS3",
                 "CD164L2", "TEKT1", "SKOR1", "VEPH1", "TCTEX1D1", "PCDH20", "NEUROD6", # "ZBBX", "TEKT2",
                 "FGF8", "MEIS2", "NDRG1")
rownames(chir.v.iwp2.80) <- chir.v.iwp2.80$gene
logfc.threshold <- 1
p.threshold = 1e-10

#This creates new entries for the genes of interest. This means that they're plotted twice,
#which prevents them from being buried underneath genes we're not interested in
#Also, rownames can't be duplicated, so a 1 is automatically appended to the end of each gene
chir.v.iwp2.80 <- rbind(chir.v.iwp2.80, chir.v.iwp2.80[(rownames(chir.v.iwp2.80) %in% lab.gen.80), ])


#This creates a long character vector composed of the color for each gene. Each entry corresponds to one gene. 
color.key <- ifelse(chir.v.iwp2.80$log2FoldChange > logfc.threshold & chir.v.iwp2.80$padj < p.threshold, CHIR,
                    ifelse(chir.v.iwp2.80$log2FoldChange < -logfc.threshold & chir.v.iwp2.80$padj < p.threshold, IWP2,
                           "lightgrey"))

#This sets the colors for the specific genes you want to highlight
color.key[(rownames(chir.v.iwp2.80) %in% c(lab.gen.80, paste0(lab.gen.80, "1"))) & chir.v.iwp2.80$log2FoldChange < -0.5] <- "darkred"
color.key[(rownames(chir.v.iwp2.80) %in% c(lab.gen.80, paste0(lab.gen.80, "1"))) & chir.v.iwp2.80$log2FoldChange > -1] <- "#3E3C9A"

#Name your color.key. Enhanced volcano doesn't work unless the color.key is named.
#Names are displayed as a plot legend.
names(color.key)[color.key == "#B856D7"] <- "IWP2"
names(color.key)[color.key == "#55A0FB"] <- "CHIR"
names(color.key)[color.key == "lightgrey"] <- "Below_Threshold"
names(color.key)[color.key == "darkred"] <- "Selected"
names(color.key)[color.key == "#3E3C9A"] <- "SelectedCHIR"


DESeq2.volcano.selectGenes <- EnhancedVolcano(chir.v.iwp2.80,
                                              lab = rownames(chir.v.iwp2.80),
                                              x = "log2FoldChange",
                                              y = "padj",
                                              title = "Day 80 HCs: IWP2 vs CHIR",
                                              pCutoff = 1e-10,
                                              FCcutoff = 1, #This sets the cutoff outside the bounds of the plot, so we don't see them
                                              pointSize = 4,
                                              colCustom = color.key,
                                              selectLab = "lab.gen.80",
                                              legendPosition = "right",
                                              labSize = 5,
                                              drawConnectors = T,
                                              arrowheads = FALSE,
                                              maxoverlapsConnectors = Inf) + 
  xlim(-4.5, 4.5) + 
  NoLegend()

DESeq2.volcano.selectGenes
