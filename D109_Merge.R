#Day 109 - Submission Version
library(Seurat)
library(tidyverse)
library(ggplot2)
library(purrr)
library(clusterCrit)

## ----------------------------- Load in Data  ---------------------------------

chir <- Read10X("C:/Users/.../filtered_feature_bc_matrix") %>%
  CreateSeuratObject()
iwp2 <- Read10X("C:/Users/.../filtered_feature_bc_matrix") %>%
  CreateSeuratObject()
iwp2.small <- Read10X("C:/Users/.../filtered_feature_bc_matrix") %>%
  CreateSeuratObject()


#Add metadata
chir$Day <- "109"
iwp2$Day <- "109"
iwp2.small$Day <- "109"

chir$Condition <- "CHIR"
iwp2$Condition <- "IWP2"
iwp2.small$Condition <- "IWP2"
chir$mt.per <- PercentageFeatureSet(chir, pattern = "^MT-")


## ----------------------------- Filtering  -----------------------------------

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
list.109 <- list(chir, iwp2, iwp2.small)

#Filters each dataset in list.109 according to mitochondrial percentage,
#gives number of cells filtered for each dataset
#cell count to counts.list. Each group is downsampled to sample.size <- min(counts.list)
counts.list.post <- counts.list.pre <- numeric()
for (i in 1:3) {
  counts.list.pre[i] <- length(Cells(list.109[[i]]))
  list.109[[i]] <- filterPlot(list.109[[i]])
  counts.list.post[i] <- length(Cells(list.109[[i]]))
}
diff <- counts.list.pre - counts.list.post
diff

#same thing but for RPL27
#Need to maintain consistent pre-processing throughout all datasets
counts.list.post.2 <- numeric()
for (i in 1:3) {
  list.109[[i]] <- housekeeping.filter(list.109[[i]], housekeeping.gene = "RPL27", subset = T, std.dev = 2)
  counts.list.post.2[i] <- length(Cells(list.109[[i]]))
}
diff.2 <- counts.list.post - counts.list.post.2
diff.2


#Preprocessing complete
# Sequentially merge the samples (w/o down sampling)

length(Cells(list.109[[1]])) #5971
length(Cells(list.109[[2]])) #9561
length(Cells(list.109[[3]])) #512

list.109[[3]]$size <- "small.dataset"
list.109[[2]]$size <- "large.dataset"
list.109[[1]]$size <- "large.dataset"

merge.1 <- merge(list.109[[3]], list.109[[2]])
day109 <- merge(list.109[[1]], merge.1)
length(Cells(day109))
rm(iwp2, iwp2.small, list.109, chir)

## ----------------------------- Seurat workflow -------------------------------

day109 <- SCTransform(day109, vars.to.regress = "mt.per")
day109 <- RunPCA(day109)
day109 <- FindNeighbors(day109, k.param = 410, dims = 1:9) %>%
  FindClusters(resolution = 1.0) %>% 
  RunUMAP(n.neighbors = 410, dims = 1:9)

#Calculating ICVI score
clusters <- as.integer(Idents(day109.hires))
sil.score <- intCriteria(day109@tools$LoadEmbeddings,
                         part = clusters,
                         crit = "Sil")
sil.score

saveRDS(day109, "C:/.../day109AutoClustRd.rds")

## ----------------------------- DEA --------------------------------------

coch.v.ut <- FindMarkers(day109, ident.1 = 2, ident.2 = 4, logfc.threshold = 0.5)
coch.v.ut[6] <- rownames(coch.v.ut)
frame.names <- names(coch.v.ut)
frame.names[3:4] <- c("pct.CHIR", "pct.IWP2")
frame.names[6] <- "Gene"
names(coch.v.ut) <- frame.names
coch.v.ut <- arrange(coch.v.ut, desc(avg_logFC))
write.csv(coch.v.ut, file = "C:/Users/.../HC_Comparison_d109.csv")


coch.v.ut.sc <- FindMarkers(day109, ident.1 = 1, idents.2 = 5, logfc.threshold = 0.5)
coch.v.ut.sc[6] <- rownames(coch.v.ut.sc)
frame.names <- names(coch.v.ut.sc)
frame.names[3:4] <- c("pct.CHIR", "pct.IWP2")
frame.names[6] <- "Gene"
names(coch.v.ut.sc) <- frame.names
coch.v.ut <- arrange(coch.v.ut, desc(avg_logFC))
write.csv(coch.v.ut, file = "C:/Users/.../SC_Comparison_d109.csv")

d109.mark <- FindAllMarkers(day109, logfc.threshold = 0.69, min.pct = 0.3, only.pos = T)

top.marks <- d109.mark %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice(1:25)

## ----------------------------- Plotting  ---------------------------------

#Define condition colors
IWP2 <- "#B856D7"
CHIR <- "#55A0FB"

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

#Basically just need to define DimPlot, feature plot, and Violin plot functions, then 
#pass in a gene list


con.plot <- DimPlot(day109, group.by = "Condition", pt.size = 1.8,
                    cols = c(CHIR, IWP2), cells = sample(Cells(day109), 14000)) +
  NoLegend() + sparse()
save.plot(con.plot, name = "Day109_By_Condition_NoLegend_sparse")


## ----------------------------- Violin & Feature Plots -------

split.day109 <- subset(day109, idents = c(2,4))
split.day109@active.ident <- factor(as.factor(split.day109$Condition), levels = c("IWP2", "CHIR"))
split.day109@active.ident <- as.factor(rep("Day 20", length(Cells(split.day109))))
split.day109$Condition <- factor(as.factor(split.day109$Condition), levels = c("IWP2", "CHIR"))


#Set Folder
folder <- "C:/Users/Examplefolder"

save.png <- function(plot, name, y.iter){
  plot.path <- paste0("/Users/alexander/Dropbox/Mac/Documents/HashinoLab/NatBiotech/Plots/D109_requests/",
                      as.character(y.iter),
                      "/", name,
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
  plot.path <- paste0("/Users/alexander/Dropbox/Mac/Documents/HashinoLab/NatBiotech/Plots/D109_requests/",
                      as.character(y.iter),
                      "/", name,
                      ".eps")
  ggsave(filename = plot.path,
         plot = plot,
         device = "eps",
         units = "in",
         height = 4,
         width = 4,
         dpi = 300)
}

gene.list.combined <- c("NR2F1", "GATA3", "INSM1", "HES6", "TMPRSS3", "ZNF503", "FGF8",
                         "GNG8", "LFNG", "CD164L2", "ZBBX", "TEKT1", "SKOR1", "AMPD3", "VEPH1", 
                         "TCTEX1D1", "PCDH20", "NEUROD6", "NDRG1", "MEIS2")

for (y.iter in c(2, 2.25, 2.5, 2.75)){
  for(i in 1:length(gene.list.combined)){
    v.plot <- VlnPlot(split.day109, 
                      features = gene.list.combined[i], 
                      split.by = "Condition", 
                      split.plot = T, 
                      cols = c(IWP2, CHIR),
                      pt.size = 0,
                      log = F,
                      y.max = y.iter) 
    v.plot$layers[[1]]$aes_params$size <- 0.75
    v.plot <- v.plot +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
      xlab("") +
      #sparse() +
      NoLegend()
    name <- paste0("WithAxis", gene.list.combined[i], "_Day109_Vln")
    save.png(v.plot, name, y.iter)
    save.eps(v.plot, name, y.iter)
    v.plot <- v.plot + sparse()
    name <- paste0(gene.list.combined[i], "_Day109_Vln")
    save.png(v.plot, name, y.iter)
    save.eps(v.plot, name, y.iter)
  }
}

#Read in data
hc.markers <- readxl::read_xlsx("C:/Users/.../HC_Markers2.xlsx")

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
      { VlnPlot(split.day109,
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
        return(NULL)
      }
    )
    if (!is.null(plot)){
      save.plot(plot, name = paste0(gene, "_Day109_Vln"), folder = marker.folder)
    }
  }
}

#  Feature Plots to OneDrive

#Slightly different dimensions
save.plot <- function(plot, name, folder = "C:/Users/.../FeaturePlots/"){
  ggsave(filename = paste0(folder,name, ".jpg"),
         plot = plot,
         device = jpeg(),
         units = "in",
         height = 8,
         width = 8,
         dpi = 300)
  dev.off()
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
      { FeaturePlot(day109,
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
        return(NULL)
      }
    )
    if (!is.null(plot)){
      save.plot(plot, name = paste0(gene, "_Day109_Feature"))
    }
  }
}


## ----------------------------- Blend Plot   ---------------------------------

blend.109 <- FeaturePlot(day109,
            features = c("FGF20", "FGFR3"), #c("MEIS2", "GATA3"),
            pt.size = 1.6,
            blend = T,
            order = T,
            blend.threshold = 0.25,
            cols = c("lightgrey",  "#5d55fb", "#B856D7"), #c("lightgrey",  "#0CE994", "#E83C7E"),
            max.cutoff = 2.0)+ 
  theme(axis.text = element_text(size = 22),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(0.25, "cm"),
        axis.title = element_text(size = 22)
  )

save.plot <- function(plot, name, folder = marker.folder){
  ggsave(filename = paste0(folder, name, ".jpg"),
         plot = plot,
         device = jpeg(),
         units = "in",
         height = 8,
         width = 32,
         dpi = 300)
  dev.off()
}
save.eps <- function(plot, name, folder = marker.folder){
  ggsave(filename = paste0(folder, name, ".eps"),
         plot = plot,
         device = "eps",
         units = "in",
         height = 8,
         width = 32,
         dpi = 300)
}

save.plot(blend.109, name = "FGFR3_FGF20_BlendPlot_DV_Colors")
save.eps(blend.109, name = "FGFR3_FGF20_BlendPlot_DV_Colors")



## ----------------------------- Bar Chart --------------------------------
library(cowplot)
library(forcats)

IWP2 <- "#B856D7"
CHIR <- "#55A0FB"

# Stacked Histogram

plot.data <- day109@meta.data %>%
  select(Condition, seurat_clusters) %>%
  mutate(Group = fct_collapse(seurat_clusters, #This is not the most efficient way to do this, but...
                              `Mature DRGs` = "0",
                              `CTRL Transitional` = "1",
                              `CTRL HCs` = "2",
                              `Immature DRGs` = "3",
                              `PUR + IWP2 HCs` = "4",
                              `PUR + IWP2 Transitional` = "5",
                              `Mixed Transitional` = "6",
                              `NC/Mesenchyme` = "7",
                              `FP Neurons` = "8") ) %>%
  mutate(Group = fct_relevel(Group,
                             "CTRL Transitional",
                             "Mixed Transitional",
                             "PUR + IWP2 Transitional",
                             "PUR + IWP2 HCs",
                             "CTRL HCs",
                             "FP Neurons",
                             "Immature DRGs",
                             "Mature DRGs",
                             "NC/Mesenchyme"
                             )) %>%
  group_by(Group, Condition) %>%
  summarise(Count = n()) %>%
  mutate(Prop = Count / sum(Count))

f2E <- ggplot(plot.data, aes(x = Group, y = Prop, fill = Condition)) +
  geom_col() +
  labs(y = "Cluster Composition") +
  guides(fill = guide_legend(direction = "vertical", ncol = 3)) +
  scale_fill_manual(values = c(CHIR, IWP2) ) +
  theme_cowplot(font_size = 22) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle =45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center") + NoLegend()
f2E

ggsave(filename = "/Users/alexander/Dropbox/Mac/Documents/HashinoLab/NatBiotech/Figures/D109_StackedHist.png",
       plot = f2E,
       width = 6,
       height = 10,
       units = "in",
       dpi = 300,
       device = "png")
ggsave(filename = "/Users/alexander/Dropbox/Mac/Documents/HashinoLab/NatBiotech/Figures/D109_StackedHist.eps",
       plot = f2E,
       width = 8,
       height = 10,
       units = "in",
       dpi = 300,
       device = "eps")

save.plot <- function(plot, name){
  ggsave(filename = paste0(folder, name, ".tiff"),
         plot = plot,
         device = tiff(),
         units = "in",
         height = 10,
         width = 6,
         dpi = 300)
  dev.off()
}

plot.data <- day80@meta.data %>%
  select(Condition, seurat_clusters) %>%
  mutate(Group = fct_collapse(seurat_clusters,
                              `Neural` = c("1", "2", "8", "9", "10"),
                              `Supporting Cells` = c("0", "3", "4", "7"),
                              `Hair Cells` = c("5", "6"))) %>%
  mutate(Group = fct_relevel(Group,
                             "Neural", "Supporting Cells", "Hair Cells")) %>%
  group_by(Condition, Group) %>%
  summarise(Count = n()) %>%
  mutate(Prop = Count / sum(Count))

f2E <- ggplot(plot.data, aes(x = Condition, y = Prop, fill = Group)) +
  geom_col() +
  labs(y = "Clusters Composition") +
  guides(fill = guide_legend(direction = "vertical", ncol = 3)) +
  theme_cowplot(font_size = 22) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle =45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center")
f2E

save.plot <- function(plot, name){
  ggsave(filename = paste0(folder, name, ".tiff"),
         plot = plot,
         device = tiff(),
         units = "in",
         height = 10,
         width = 6,
         dpi = 300)
  dev.off()
}


## ----------------------------- Percent Expressing Testing   -----------------

#CALB1, MYO7A, and PCP4 are all ubiquitous in Hair cells
day109.hc <- subset(day109, idents = c(2,4))
day109.hc.iwp2 <- subset(day109.hc, Condition == "IWP2")
day109.hc.chir <- subset(day109.hc, Condition == "CHIR")


ncol(day109.hc.iwp2) #2029 Total Cells

sum(GetAssayData(day109.hc, "counts")["CALB1", ] > 0) # = 1778 cells, or 38%
sum(GetAssayData(day109.hc, "counts")["MYO7A", ] > 0) # = 1927 cells, or 42%

sum(GetAssayData(day109.hc.iwp2, "counts")["GATA3", ] > 0) # = 1370 cells, or 79%

iwp2.exp <- list(SLC26A5 = sum(GetAssayData(day109.hc.iwp2, "counts")["SLC26A5", ] > 0), # = 33 cells, or 1.6%
                 SCL17A8 = sum(GetAssayData(day109.hc.iwp2, "counts")["SLC17A8", ] > 0), # = 130 cells, or 6.4%
                 SLC26A5.pct = sum(GetAssayData(day109.hc.iwp2, "counts")["SLC26A5", ] > 0) * 100 / ncol(day109.hc.iwp2),
                 SCL17A8.pct = sum(GetAssayData(day109.hc.iwp2, "counts")["SLC17A8", ] > 0) * 100 / ncol(day109.hc.iwp2))

sum(GetAssayData(day109.hc.iwp2, "counts")["SLC17A8", ] > 0) # = 130 cells, or 6.4%

chir.exp <- list(SLC26A5 = sum(GetAssayData(day109.hc.chir, "counts")["SLC26A5", ] > 0),
                 SCL17A8 = sum(GetAssayData(day109.hc.chir, "counts")["SLC17A8", ] > 0),
                 SLC26A5.pct = sum(GetAssayData(day109.hc.chir, "counts")["SLC26A5", ] > 0) * 100 / ncol(day109.hc.chir),
                 SCL17A8.pct = sum(GetAssayData(day109.hc.chir, "counts")["SLC17A8", ] > 0) * 100 / ncol(day109.hc.chir))
res <- cbind(iwp2.exp,
             chir.exp)

sum(GetAssayData(day109.hc.chir, "counts")["SLC26A5", ] > 0)*100 #/ ncol(day109.hc.chir) # = 9 cells, or 0.35%
sum(GetAssayData(day109.hc.chir, "counts")["SLC17A8", ] > 0)*100 #/ ncol(day109.hc.chir) # = 216 cells, or 8.5%


## ----------------------------- UMAP HEX Colors ----
d109 <- readRDS(file = "/Users/.../day109AutoClustRd.rds")

save.eps <- function(plot, name){
  ggsave(filename = paste0("C:/Users", name, ".eps"),
         plot = plot,
         device = "eps",
         units = "in",
         height = 8,
         width = 8,
         dpi = 300)
}

save.png <- function(plot, name){
  ggsave(filename = paste0("C:/Users", name, ".png"),
         plot = plot,
         device = "png",
         units = "in",
         height = 8,
         width = 8,
         dpi = 300)
}

#Bright Pastels for relevant groups
palette.coo <- c(
  "#626E99",
  "#7fc763", 
  "#584fdb",   
  "#597D8C", 
  "#d76cf5", 
  "#FFC2B4", 
  "#FEDD72", 
  "#6FA397",
  "#D4AD80")

uplot <- DimPlot(d109, cols = palette.coo, label = F, pt.size = 1.6) + NoLegend()
uplot

save.eps(uplot, "D109_UMAP_V6")
save.png(uplot, "D109_UMAP_V6")
uplot.sparse <- uplot + sparse()
save.eps(uplot.sparse, "D109_UMAP_V6_sparse")
save.png(uplot.sparse, "D109_UMAP_V6_sparse")


## ----------------------------- Bubble Plot ----------------

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

plotdata<- fread("C:/Users/.../D109_HC_Upreg_GSEA.csv")
bp.up <- data.table::merge.data.table(plotdata, humanGeneSetsInfo, by.x = "annot_id", by.y = "gset") #I think this is easier in SQL
bp.up$Category <- droplevels(bp.up$gsetBioName)

bp.up <- bp.up[toupper(bp.up$Category) %in% includedCats, ] #subset data, removing all geneset categories not in includedCats
bp.up$Category <- droplevels(bp.up$Category) %>% factor(caseCats) # Drop unused levels, reorder factor so that plot legend matches x position
bp.up <- bp.up[order(match(toupper(bp.up$Category), includedCats)), ] #Organizes the data so that categories are grouped on the x axis
bp.up$IDNum <- row.names(bp.up) %>% as.numeric() #IDNum will define x position on bubble plot
bp.up$Log10_Pvalue_Louis <- -1*log10(bp.up$pvalue_louis)


plotdata<- fread("C:/Users/.../D109_HC_Downreg_GSEA.csv")
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


## ----------------------------- Volcano Plot -----------------------------------------

library(EnhancedVolcano)
library(dplyr)

#Define Colors
IWP2 <- "#B856D7"
CHIR <- "#55A0FB"



#Define Specific Genes to Highlight
chir.v.iwp2.109 <- fread(file = "/Users/.../D109_HC_DESeq2.csv") %>% as.data.frame()
chir.v.iwp2.109$log2FoldChange <- chir.v.iwp2.109$log2FoldChange * -1
lab.gen.109 <- c("NR2F1", "NR2F2", "GATA3", "INSM1", "ZNF503", "GNG8", "HES6", "TMPRSS3",
                 "CD164L2", "TEKT1", "SKOR1", "VEPH1", "TCTEX1D1", "PCDH20", "NEUROD6", # "ZBBX", "TEKT2",
                 "FGF8", "MEIS2", "NDRG1")
rownames(chir.v.iwp2.109) <- chir.v.iwp2.109$gene
logfc.threshold <- 1
p.threshold = 1e-10

#This creates new entries for the genes of interest. This means that they're plotted twice,
#which prevents them from being buried underneath genes we're not interested in
#Also, rownames can't be duplicated, so a 1 is automatically appended to the end of each gene
chir.v.iwp2.109 <- rbind(chir.v.iwp2.109, chir.v.iwp2.109[(rownames(chir.v.iwp2.109) %in% lab.gen.109), ])


#This creates a long character vector composed of the color for each gene. Each entry corresponds to one gene. 
color.key <- ifelse(chir.v.iwp2.109$log2FoldChange > logfc.threshold & chir.v.iwp2.109$padj < p.threshold, CHIR,
                    ifelse(chir.v.iwp2.109$log2FoldChange < -logfc.threshold & chir.v.iwp2.109$padj < p.threshold, IWP2,
                           "lightgrey"))

#This sets the colors for the specific genes you want to highlight
color.key[(rownames(chir.v.iwp2.109) %in% c(lab.gen.109, paste0(lab.gen.109, "1"))) & chir.v.iwp2.109$log2FoldChange < -0.5] <- "darkred"
color.key[(rownames(chir.v.iwp2.109) %in% c(lab.gen.109, paste0(lab.gen.109, "1"))) & chir.v.iwp2.109$log2FoldChange > -1] <- "#3E3C9A"

#Name your color.key. Enhanced volcano doesn't work unless the color.key is named.
#Names are displayed as a plot legend.
names(color.key)[color.key == "#B856D7"] <- "IWP2"
names(color.key)[color.key == "#55A0FB"] <- "CHIR"
names(color.key)[color.key == "lightgrey"] <- "Below_Threshold"
names(color.key)[color.key == "darkred"] <- "Selected"
names(color.key)[color.key == "#3E3C9A"] <- "SelectedCHIR"


DESeq2.volcano.selectGenes <- EnhancedVolcano(chir.v.iwp2.109,
                                              lab = rownames(chir.v.iwp2.109),
                                              x = "log2FoldChange",
                                              y = "padj",
                                              title = "Day 109 HCs: IWP2 vs CHIR",
                                              pCutoff = 1e-10,
                                              FCcutoff = 1, #This sets the cutoff outside the bounds of the plot, so we don't see them
                                              pointSize = 4,
                                              colCustom = color.key,
                                              selectLab = "lab.gen.109",
                                              legendPosition = "right",
                                              labSize = 5,
                                              drawConnectors = T,
                                              arrowheads = FALSE,
                                              maxoverlapsConnectors = Inf) + 
  xlim(-4.5, 4.5) + 
  NoLegend()

DESeq2.volcano.selectGenes
