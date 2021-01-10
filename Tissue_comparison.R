#Tissue Comparison
setwd("~/Dropbox/Maria_SaxenaLab/Liver_Microbiome_2nd_review/tissue/")
source("~/Dropbox/Maria_SaxenaLab/Liver_QIIME/Josh_Human/microbiome_functions.R")
library(phyloseq) #V1.27.0
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(ggpubr)
library(tidyr)
library(reshape2)
library(gplots)
library(data.table)
library(FSA)
library(dplyr)
library(vegan)

#Import files 
tissue_phylo <- import_data("otu_table_nochimera_tissue.biom", mapping = "metadata_tissue.txt", tree = "rep_set_chimerafree_tissue.tre")
metadata<-data.frame(sample_data(tissue_phylo))
head(metadata)
tax<-data.frame(tax_table(tissue_phylo))

tissue_phylo_trim=subset_samples(tissue_phylo, SampleID !="GFLiver" & SampleType != "BufferControl" &  Description != "CCLK5" &  Description != "LymphNode")


#Beta Diversity

#p value
dist <- phyloseq::distance(tissue_phylo_trim, method = "bray")
perma <- adonis(dist~Description, data = as(sample_data(tissue_phylo_trim), "data.frame"), permutations = 1000)

GP.ord <- ordinate(tissue_phylo_trim, "PCoA", "bray")
plot_ordination(tissue_phylo, GP.ord,  color="Description", axes =1:2)+ 
  stat_ellipse(type = "t",linetype = 1,alpha=0.5)+
  theme_bw(base_size = 20, base_line_size = 0.5)+ geom_point(size=5)+labs(color="Treatment")+ scale_color_brewer(palette = c("Dark2")) +
  annotate("text", x = -0.3, y = 0.6, 
           label = "paste(italic(p), \" = 0.00099\")", parse = TRUE, size=6)

#Heatmap Genus
tissue_phylo_glom <- tax_glom(tissue_phylo_trim, "Genus")
top5 = prune_taxa(names(sort(taxa_sums(tissue_phylo_glom), TRUE))[1:40], tissue_phylo_glom)
tissue.trans <- microbiome::transform(top5, "log10")
otu_trans=as.data.frame(otu_table(tissue.trans))
taxa_names=as.data.frame(tax_table(tissue_phylo_trim))
taxa_names=taxa_names[rownames(otu_trans),]

metadata=metadata[order(metadata$Description),]
ann=metadata %>% select("Description")
col_order <- row.names(ann)
otu_trans <- otu_trans[, col_order] #Rorder columns
brewer.pal(8, name = "Dark2")
Description=c( "#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E")
names(Description)=c("Duodenum",  "GallBladder", "LiverWT",  "Pancreas", "Spleen")
ann_colors= list(Description=Description)


pheatmap(as.matrix(otu_trans), border_color = "black",fontsize_row=15, color = rev(brewer.pal(n=11,name = "RdYlBu")),
         fontsize_col = 10, scale = "none", 
         cluster_cols = F, cluster_rows = T,
         labels_col = NULL,labels_row=paste0(taxa_names$Genus),  annotation_col = ann, annotation_colors = ann_colors)

