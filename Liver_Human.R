#Livar and Gut Data Mice

# Import libraries
library(phyloseq);packageVersion("phyloseq")	# v1.27.0
library(ggplot2);packageVersion("ggplot2")	# v3.1.0
library(extrafont)
library(RColorBrewer)
library(plyr)
library(scales)
library(vegan)	# v2.5-3
library(gplots); packageVersion("gplots")	#v3.0.1
library(reshape2)

# create plot folder
path <- "./plot"
dir.create(path)
setwd("~/Dropbox (Mason Lab)/Maria_SaxenaLab/Liver_QIIME/Box6")
source("~/Dropbox (Mason Lab)/Maria_SaxenaLab/Liver_QIIME/Josh_Human/microbiome_functions.R")


# Define default theme for figures
default.theme <- theme(axis.title = element_text(size = 12, family = "Arial", 
                                                 face = "bold"),
                       axis.line= element_line(color = "black"),
                       axis.text = element_text(size = 10, family = "Arial"),
                       legend.text = element_text(family = "Arial", size = 10),
                       legend.title = element_text(family = "Arial", size = 12, face = "bold"),
                       strip.text = element_text(family = "Arial", size = 8),
                       panel.grid = element_blank(),
                       panel.background = element_blank())

# Define function to import data
import_data <- function(otu = "otu_table.biom",
                        tree = NULL,
                        mapping = "mapping_file.txt",
                        parseFunction = parse_taxonomy_greengenes){
  # import the biom table to a phylosqe object
  # otu is the otu table with the OTU ID and Abundance
  # mapping is the mapping file with all the information about samples
  
  otutable <- import_biom(otu,
                          treefilename = tree,
                          parseFunction = parse_taxonomy_greengenes)
  
  meta <- import_qiime_sample_data(mapping)
  colnames(meta)[1] <- "SampleID"
  phylo <- merge_phyloseq(otutable, meta)
  return(phylo)
}

pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}

# import fecal data
phylo <- import_data(otu="otu_table_abx.biom", 
                     mapping = "mapping_file_abx.txt", 
                     tree = "rep_set_filtered_abx.tre")
levels(sample_data(phylo)$SampleType) <- list(Gut = "Fecal", Liver = "Liver")
tax_table(phylo) <- tax_table(phylo)[, 1:7]

# Root the tree
require(ape)
new.outgroup <- pick_new_outgroup(phy_tree(phylo))
phy_tree(phylo) <- ape::root(phy_tree(phylo), 
                             outgroup=new.outgroup, resolve.root=TRUE)


# Subset data
control <- subset_samples(phylo, Description == "Control")
control <- prune_taxa(taxa_sums(control) > 0, control)

gut <- subset_samples(phylo, SampleType == "Gut")
gut <- prune_taxa(taxa_sums(gut) > 0, gut)

liver <- subset_samples(phylo, SampleType == "Liver")
liver <- prune_taxa(taxa_sums(liver) > 0, liver)

# Alpha Diversity
measures <- c("Observed", "ACE", "Chao1", "Shannon", "Simpson")
plot_richness(control, x = "SampleType", color = "SampleType", measures = measures) +
  geom_boxplot(width = .5, position = "dodge", size = .8) +
  default.theme +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 8, family = "Arial")) + 
  labs(x = "", y = "Alpha Diversity")

# PD Diversity
library(picante)

control.otu <- as(otu_table(control), "matrix")
control.otu <- t(control.otu)
control.meta <- as(sample_data(control), "data.frame")
control.pd <- pd(control.otu, phy_tree(control))
control.pd$SampleType <- control.meta[rownames(control.pd), "SampleType"]

ggplot(control.pd, aes(x = SampleType, y = PD, color = SampleType)) +
  geom_point() +
  geom_boxplot(width = .3, position = "dodge", size = .8) + 
  default.theme +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 8, family = "Arial"))

#p value alpha diversity
alp <- estimate_richness(control, measures = measures)
alp$PD <- control.pd$PD
wil <- t(sapply(alp, 
                function(x) unlist(
                  wilcox.test(x~sample_data(control)$SampleType))))

#Rel. Abundance Gut vs. Liver

# Phylum
glomrank <- "Phylum"
control.glom <- tax_glom(control, glomrank)
control.glom.norm <- transform_sample_counts(control.glom,
                                             function(x) x/sum(x))
control.glom.md <- psmelt(control.glom.norm)
control.glom.dat <- ddply(control.glom.md, c("SampleType", "Phylum"),
                          summarise,
                          mean = mean(Abundance), sd = sd(Abundance),
                          N = length(Abundance), se = sd/sqrt(N))

ggplot(control.glom.dat, aes(x= reorder(Phylum, -mean), y = mean, group = SampleType, fill = SampleType)) +
  geom_bar(stat = "identity", position = "dodge", width = .5, colour = "black") +
  geom_errorbar(aes(ymin = mean, ymax = mean + sd), 
                stat = "identity", 
                position = position_dodge(.5), width = .1) +
  scale_y_continuous(labels = percent) +
  labs(x = "Phylum",y = "Relative Abundance") +
  default.theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
#p value
control.pvalue <- data.frame(Phylum = unique(as.character(control.glom.md$Phylum)),
                             row.names = unique(as.character(control.glom.md$Phylum)),
                             pvalue = 0)
for (p in rownames(control.pvalue)){
  wil <- wilcox.test(Abundance~SampleType, data = control.glom.md, subset = Phylum == p)
  control.pvalue[p, "pvalue"] <- wil$p.value
}

#Heatmap
control.otu <- as.data.frame(otu_table(control.glom.norm))
control.tax <- as(tax_table(control.glom.norm), "matrix")
control.tax <- as.data.frame(control.tax)
control.meta <- as(sample_data(control.glom.norm), "data.frame")

control.dat <- apply(control.otu, 1, scale)
rownames(control.dat) <- colnames(control.otu)
merged <- merge(control.meta, control.dat, by = "row.names", all.x = FALSE)
reorder <- merged[order(merged$SampleType, merged$Row.names), ]

cols <- c(Gut="#940000", PDA="#002ba1")
sidebar <- cols[reorder$SampleType]

control.htm <- reorder[, 8:ncol(reorder)]
rownames(control.htm) <- reorder$Row.names
control.htm <- t(control.htm)
rownames(control.htm) <- control.tax[rownames(control.htm), "Genus"]
control.htm <- as.matrix(control.htm)

scaleredblue <- colorRampPalette(colors = c("red","black","green"))(200)


heatmap.2(control.htm, Rowv = TRUE, Colv = FALSE, dendrogram = "row",scale = "none", 
          trace = "none", col = scaleredblue, xlab = "Sample", ylab = "Genus", 
          key.xlab = "z-score",  margins = c(10, 30), ColSideColors = sidebar, 
          density.info= "none", key.title = NULL, cexRow = 1.5, cexCol = 1.5, 
          main = NULL, keysize = .7, hclustfun = function(x) hclust(x, method="ward.D2"), 
          distfun = function(x) dist(x, method="euclidean"))
legend("top", bty = "n",legend = names(cols), fill = cols, pch = ".", horiz = TRUE)


# Lefse
## Export and Convert to Biom table
control.otus <- as(otu_table(control), "matrix")
control.tax <- as(tax_table(control), "matrix")
control.otus <- as.data.frame(control.otus)
control.tax <- as.data.frame(control.tax)


for (t in colnames(control.tax)){
  if(!is.null(control.tax$taxonomy)){
    control.tax$taxonomy <- paste0(control.tax$taxonomy, "; ", control.tax[, t])
  }else {control.tax$taxonomy <- control.tax[, t]}
}
control.output <- merge(control.otus, control.tax, by = "row.names")
rownames(control.output) <- control.output$Row.names
colnames(control.output)[1] <- "OTUs"
control.output <- control.output[, -c(12:22)]
write.table(control.output,
            file.path("./lefse/Control_OTU.txt"),
            row.names = F, quote = F, sep = "\t")
control.meta <- as(sample_data(control), "data.frame")
write.table(control.meta[, -c(2:4,6)],
            file.path("./lefse/Control_meta.txt"),
            row.names = F, quote = F, sep = "\t")
