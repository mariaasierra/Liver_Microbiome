## Liver and Gut human data Josh 
# Script Maria A Sierra 
setwd("~/Desktop/Maria/QIIME/Josh_Human/R_inputs/")
setwd("~/Dropbox (Mason Lab)/Maria_SaxenaLab/Liver_QIIME/Josh_Human/R_inputs")

source("~/Dropbox (Mason Lab)/Maria_SaxenaLab/Liver_QIIME/Josh_Human/microbiome_functions.R")
library(phyloseq) #V1.27.0
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(ggpubr)
library(tidyr)
library(reshape2)
library(gplots)
library(data.table)
library(metafor)
library(vegan)
devtools::install_github("vmikk/metagMisc", force = TRUE)
library(metagMisc)
library(PhyloMeasures)

#Import files 
human_phylo <- import_data("otu_table_nochimera.json.biom", mapping = "mapping_file.txt", tree = "rep_set_chimerafree.tre")
metadata<-data.frame(sample_data(human_phylo))
tax<-data.frame(tax_table(human_phylo))
newhuman_phylo = subset_samples(human_phylo, SampleID != "PL16" & SampleID != "PF16")
metadata<-data.frame(sample_data(newhuman_phylo))


#Rarefaction
calculate_rarefaction_curves <- function(human_phylo, measures, depths) {
  estimate_rarified_richness <- function(human_phylo, measures, depth) {
    if(max(sample_sums(human_phylo)) < depth) return()
    physeq.tree <- prune_samples(sample_sums(human_phylo) >= depth, human_phylo)
    
    rarified_psdata <- rarefy_even_depth(human_phylo, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }  
  
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, human_phylo = human_phylo, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(human_phylo, c('Observed'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(human_phylo)), by.x = 'Sample', by.y = 'row.names')

pdf("Rarefaction_curve.pdf", height = 10, width = 10)
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Description,
    group = Sample)) + geom_line() +
  xlab("No. of sequences") + ylab("No. of OTUs") + scale_colour_manual(values =c( "red", "darkblue")) +
 theme(plot.title = element_text(size=40,lineheight=1, vjust=1, hjust = 0.5,
                                  family="Arial", face = "bold",margin = margin(10, 0, 10, 0)),
        axis.text.y = element_text(color="black", size=20), 
        axis.text.x = element_text(color="black", size=20 ),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 20),
       axis.title.x = element_text(size = 20),
       legend.text = element_text(size = 18),
       legend.title = element_text(face = "bold",size = 18))
dev.off()

#Diversity indices

measures <- c("Observed", "ACE", "Chao1", "Shannon", "Simpson", "FaithPD")
rich_values<-estimate_richness(newhuman_phylo, measures = measures)
#write.csv(rich_values, "rich_values.csv")
#human_phylogeny <- import_data("otu_table_nochimera.json.biom", mapping = "mapping_file.txt", tree = "rep_set_chimerafree.tre")
#rich_values<-estimate_richness(human_phylogeny, measures = measures)
PD=phyloseq_phylo_div(newhuman_phylo, measures = c("PD"))

rich_values_total= cbind(PD, rich_values)

#write.csv(rich_values, "rich_values_wo_phylo.csv")

pdf("plot_richness.pdf", height = 10, width = 10)
plot_richness(newhuman_phylo, x = "Description", color = "Description", measures = measures)+
  geom_boxplot() +
  default.theme + scale_color_manual(values=c( "red", "darkblue")) +
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 10, family = "Arial"),
        legend.text = element_text(size = 18),
        legend.title = element_text(face = "bold",size = 18)) +
  labs(x = "", y = "Alpha Diversity")
dev.off()

# calculate t-test
t(sapply(rich_values, function(x) unlist(t.test(x~sample_data(newhuman_phylo)$Description)[c("estimate","p.value","statistic","conf.int")])))

#Barplot diversity
newhuman_phylo.norm <- transform_sample_counts(newhuman_phylo, function(x) x/sum(x))

pdf("Diversity_phylum_wo_norm.pdf", height = 10, width = 10)
plot_bar(newhuman_phylo.norm,  x="Sample", y="Abundance", fill = "Phylum") + geom_bar(stat="identity", position="stack") + 
  labs(x = "Human Samples", y = "Relative abundance") + theme(plot.title = element_text(size=40,lineheight=1, vjust=1, hjust = 0.5,
                                                                family="Arial", face = "bold",margin = margin(10, 0, 10, 0)),
                                      axis.text.y = element_text(color="black", size=30), 
                                      axis.text.x = element_text(face = "bold"),
                                      axis.ticks.x=element_blank(),
                                      axis.title.y = element_text(size = 20))
dev.off()

ps <- tax_glom(newhuman_phylo.norm, "Phylum")
ps1 <- merge_samples(ps, "Description")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Phylum")

##Top phyla
#Relative Abundance plots
cols <- c("Fecal" = "red", "Liver" = "blue")

#Phylum level
pdf("abundance_plot_phylum.pdf", height = 6, width = 8)
abundance_plot_data=plot_abund(newhuman_phylo, glomrank = "Phylum", group = "Description", cols = cols )
abundance_plot_data
dev.off()
#write.csv(abundance_plot_data$data, "abundance_plot_data.csv")

#Class level
abundance_plot_class=plot_abund(human_phylo, glomrank = "Class", group = "Description", cols = cols )
pdf("abundance_plot_class.pdf", height = 6, width = 8)
abundance_plot_class
dev.off()
#write.csv(abundance_plot_class$data, "abundance_plot_class.csv")

# Genus Heatmap
gm <- tax_glom(newhuman_phylo, "Phylum")
gm.norm <- transform_sample_counts(gm, function(x) x/sum(x))
otu <- as(otu_table(gm.norm), "matrix")
meta <- as(sample_data(gm), "data.frame")
meta <- meta[order(meta$Description, meta$SampleID), ]
tax <- data.frame(as(tax_table(gm), "matrix"))
otu.center <- apply(otu, 1, scale)
rownames(otu.center) <- colnames(otu)
otu.center <- otu.center[match(rownames(meta), rownames(otu.center)), ]
otu.center <- t(otu.center)

top.otus <- names(sort(taxa_sums(gm.norm), TRUE))[1:5]
hmp.dat <- otu.center[top.otus, ]
rownames(hmp.dat) <- tax[rownames(hmp.dat), "Phylum"]
colbar <- cols[meta$Description]
row.hc <- hclust(dist(hmp.dat))
row.dd <- as.dendrogram(row.hc)

labels(row.dd)  # dendrogram leaves
rownames(hmp.dat)   # original data labels

row.dd.ordered<-order.dendrogram(row.dd) # the index of dendrogram 

heatmap.2(hmp.dat, Rowv = row.dd.ordered, 
          Colv = FALSE, dendrogram = "row", scale = "none", trace = "none", 
          col = colorRampPalette(colors = c("red","black","darkblue"))(300), labCol = NA,
          labRow=as.expression(lapply(rownames(hmp.dat), function(a) bquote(italic(.(a))))),
          key.xlab = "z-score", density.info = "none", margins = c(10, 15), 
          ColSideColors = colbar, key.title = NA, cexRow = 1.5, main = NULL, keysize = 0.9)  + theme(legend.position = "left")

#Phylum level
library(microbiome)

#Normalize samples
human_phylo_norm = transform_sample_counts(human_phylo,  function(x) x / sum(x))
plot_bar(human_phylo_norm,  x="Sample", y="Abundance", fill = "Phylum") + 
  geom_bar(stat="identity", position="stack")+  
  labs(x = "", y = "Abundance") + 
  theme(plot.title = element_text(size=40,lineheight=1, 
                                  vjust=1, hjust = 0.5,
                                  family="Arial", face = "bold",margin = margin(10, 0, 10, 0)),
                                        axis.text.y = element_text(color="black", size=30), 
                                        axis.text.x = element_text(face = "bold"),
                                        axis.ticks.x=element_blank(),
                                        axis.title.y = element_text(size = 25))


# PCoA Gut vs Liver
ord = ordinate(newhuman_phylo, "PCoA", "wUnifrac")
(ordplot <- plot_ordination(newhuman_phylo, ord, "Samples", color="Description", axes = 1:2))
pdf("PCoA_GutvsLiver.pdf", height = 10, width = 14)
ordplot +  stat_ellipse(type = "t",linetype = 2,alpha=0.5) + geom_point(size=4, alpha=0.5) +
  geom_point(size=4, alpha=0.5) + theme_test() + 
  scale_color_manual(values =c( "red", "darkblue")) +
  annotate("text", x = -0.35, y = 0.45, 
           label = "paste(italic(p), \" = 0.00099\")", parse = TRUE, size=10)  + 
  theme(axis.text.y = element_text(color="black",size=20), 
        axis.text.x = element_text(color="black", size=20),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        legend.title = element_text( size=20, face = "bold"), legend.text = element_text(size=15, face = "bold"))

png(filename = "GutvsLiver_Human.png",
    width = 600, height = 600, units = "px")
plot_ordination(newhuman_phylo, ord, type = "samples", color = "Description") +
  geom_point(size = 5,shape = 20) + stat_ellipse() + default.theme + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 17),
                                                                           legend.position = "none") +labs( x ="PC1 [31.5%]", y = "PC2 [21.1%]")+
  scale_color_manual(values = c("red", "blue")) +
  annotate("text", x = -0.30, y = -0.39 ,
           label = "paste(italic(p), \" = 0.001\")", parse = TRUE, size=10)

dev.off()

dist <- phyloseq::distance(newhuman_phylo, method = "wUnifrac")
perma <- adonis(dist~Description, data = as(sample_data(newhuman_phylo), "data.frame"), permutations = 1000)
dist <- permutest(betadisper(dist, sample_data(newhuman_phylo)$Description), permutations = 1000)

# PCoA Female vs Male
ord = ordinate(newhuman_phylo, "PCoA", "bray")
(ordplot <- plot_ordination(newhuman_phylo, ord, "Samples", color="Sex", axes = 1:2))
pdf("PCoA_Sex.pdf", height = 8, width = 8)
ordplot +  stat_ellipse(type = "t",linetype = 2,alpha=0.5) + geom_point(size=4, alpha=0.5) +
  geom_point(size=4, alpha=0.5) + theme_test() + 
  scale_color_manual(values =c( "hotpink3", "olivedrab4"))
dev.off()

dist <- phyloseq::distance(newhuman_phylo, method = "bray")
perma <- adonis(dist~Description, data = as(sample_data(newhuman_phylo), "data.frame"), permutations = 1000)
dist <- permutest(betadisper(dist, sample_data(newhuman_phylo)$Sex), permutations = 1000)

# PCoA Female vs Male Total

ord = ordinate(human_phylo, "PCoA", "bray")
(ordplot <- plot_ordination(human_phylo, ord, "Samples", color="Sex", axes = 1:2))
pdf("PCoA_Sex.pdf", height = 8, width = 8)
ordplot +  stat_ellipse(type = "t",linetype = 2,alpha=0.5) + geom_point(size=4, alpha=0.5) +
  geom_point(size=4, alpha=0.5) + theme_test() + 
  scale_color_manual(values =c( "hotpink3", "olivedrab4"))
dev.off()

dist <- phyloseq::distance(human_phylo, method = "bray")
perma <- adonis(dist~Description, data = as(sample_data(human_phylo), "data.frame"), permutations = 1000)
dist <- permutest(betadisper(dist, sample_data(human_phylo)$Sex), permutations = 1000)

# PCoA Female vs Male Liver

liver_samples <- subset_samples(human_phylo, Description == "Liver")
ord = ordinate(liver_samples, "PCoA", "bray")
(ordplot <- plot_ordination(liver_samples, ord, "Samples", color="Sex", axes = 1:2))
pdf("PCoA_Sex_liver.pdf", height = 8, width = 8)
ordplot +  stat_ellipse(type = "t",linetype = 2,alpha=0.5) + geom_point(size=4, alpha=0.5) +
  geom_point(size=4, alpha=0.5) + theme_test() + 
  scale_color_manual(values =c( "hotpink3", "olivedrab4"))
dev.off()

dist <- phyloseq::distance(human_phylo, method = "bray")
perma <- adonis(dist~Sex, data = as(sample_data(human_phylo), "data.frame"), permutations = 1000)
dist <- permutest(betadisper(dist, sample_data(human_phylo)$Sex), permutations = 1000)

# PCoA HCC Liver

metadata<-data.frame(sample_data(newhuman_phylo))
metadata$Diagnosis <- as.character(metadata$Diagnosis)
metadata$Diagnosis[metadata$Diagnosis == 'Cancer'] <- 'Other'
metadata$Diagnosis[metadata$Diagnosis == 'Cholecystitis'] <- 'Other'
metadata$Diagnosis[metadata$Diagnosis == 'Cyst'] <- 'Other'
metadata$Diagnosis[metadata$Diagnosis == 'Metaplasia'] <- 'Other'
OTU<-otu_table(newhuman_phylo)
metadata<-sample_data(metadata)
tax<-tax_table(newhuman_phylo)
tree<- phy_tree(newhuman_phylo)
human_phylo2 = merge_phyloseq(OTU, metadata, tax,tree) 
liver_HCC <- subset_samples(human_phylo2, Description == "Liver")

ord = ordinate(liver_HCC, "PCoA", "wUnifrac")
(ordplot <- plot_ordination(liver_HCC, ord, "Samples", color="Diagnosis", axes = 1:2))

pdf("PCoA_HCC_liver.pdf",  height = 8, width = 12)

ordplot +  stat_ellipse(type = "t",linetype = 2,alpha=0.5) + geom_point(size=4, alpha=0.5) +
  geom_point(size=4, alpha=0.5) + theme_test() +
  scale_color_manual(values =c( "indianred3", "orange3")) +
  annotate("text", x = -0.45, y = 0.53, 
           label = "paste(italic(p), \" = 0.4136\")", parse = TRUE, size=10)  + 
  theme(axis.text.y = element_text(color="black",size=20), 
        axis.text.x = element_text(color="black", size=20),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        legend.title = element_text( size=20, face = "bold"), legend.text = element_text(size=15, face = "bold"))


png(filename = "HCCvsOther_LiverHuman.png",
    width = 600, height = 600, units = "px")
plot_ordination(liver_HCC, ord, type = "samples", color = "Diagnosis") +
  geom_point(size = 5,shape = 20) + stat_ellipse() + default.theme + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 17),
                                                                           legend.position = "right",
                                                                           legend.text = element_text(size=15),
                                                                           legend.title = element_text(size=20)) +labs( x ="PC1 [31.5%]", y = "PC2 [21.1%]")+
  scale_color_manual(values = c("indianred3", "orange3")) +
  annotate("text", x = -0.3, y = 0.65 ,
           label = "paste(italic(p), \" = 0.4136\")", parse = TRUE, size=10)

dev.off()


dist <- phyloseq::distance(liver_HCC, method = "bray")
perma <- adonis(dist~Diagnosis, data = as(sample_data(liver_HCC), "data.frame"), permutations = 1000)
dist <- permutest(betadisper(dist, sample_data(liver_HCC)$Diagnosis), permutations = 1000)

# PCoA 60 age people liver
metadata<-data.frame(sample_data(human_phylo))
#ifelse (metadata$Age>60) metadata$Age = "Elder" else metadata$Age = "Young"
human_phylo <- subset_samples(human_phylo, Description == "Liver")
metadata$Age[metadata$Age<= 60] <- 1#'Young'
metadata$Age[metadata$Age >= 60] <- 2 #'Elder'
metadata=as.data.frame(metadata)
metadata$Age[metadata$Age == 1] = "Young"
metadata$Age[metadata$Age == 2] = "Elder"
metadata<-sam_data(metadata)
OTU<-otu_table(human_phylo)
tax<-tax_table(human_phylo)
human_phylo3 = merge_phyloseq(OTU, metadata, tax) 

ord = ordinate(human_phylo3, "PCoA", "bray")
(ordplot <- plot_ordination(human_phylo3, ord, "Samples", color="Age", axes = 1:2))
pdf("PCoA_60_liver.pdf", height = 8, width = 8)
ordplot +  stat_ellipse(type = "t",linetype = 2,alpha=0.5) + geom_point(size=4, alpha=0.5) +
  geom_point(size=4, alpha=0.5) + theme_test() 
dev.off()

