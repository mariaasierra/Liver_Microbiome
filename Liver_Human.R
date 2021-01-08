## Liver and Gut Human Data 
# Script Maria A Sierra - 2019
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
#Diversity indices
measures <- c("Observed", "ACE", "Chao1", "Shannon", "Simpson", "FaithPD")
rich_values<-estimate_richness(newhuman_phylo, measures = measures)
PD=phyloseq_phylo_div(newhuman_phylo, measures = c("PD"))
rich_values_total= cbind(PD, rich_values)

plot_richness(newhuman_phylo, x = "Description", color = "Description", measures = measures)+
  geom_boxplot() +
  default.theme + scale_color_manual(values=c( "red", "darkblue")) +
  theme(axis.text.x = element_text(angle = 90, size = 15),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = 10, family = "Arial"),
        legend.text = element_text(size = 18),
        legend.title = element_text(face = "bold",size = 18)) +
  labs(x = "", y = "Alpha Diversity")
# calculate t-test
t(sapply(rich_values, function(x) unlist(t.test(x~sample_data(newhuman_phylo)$Description)[c("estimate","p.value","statistic","conf.int")])))

##Top phyla
#Relative Abundance plots
cols <- c("Fecal" = "red", "Liver" = "blue")

#Phylum level
abundance_plot_data=plot_abund(newhuman_phylo, glomrank = "Phylum", group = "Description", cols = cols )
abundance_plot_data

newhuman_phylo.norm <- transform_sample_counts(newhuman_phylo, function(x) x/sum(x))

ord = ordinate(newhuman_phylo, "PCoA", "wUnifrac")
(ordplot <- plot_ordination(newhuman_phylo, ord, "Samples", color="Description", axes = 1:2))
plot_ordination(newhuman_phylo, ord, type = "samples", color = "Description") +
  geom_point(size = 5,shape = 20) + stat_ellipse() + default.theme + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 17),
                                                                           legend.position = "none") +labs( x ="PC1 [31.5%]", y = "PC2 [21.1%]")+
  scale_color_manual(values = c("red", "blue")) +
  annotate("text", x = -0.30, y = -0.39 ,
           label = "paste(italic(p), \" = 0.001\")", parse = TRUE, size=10)


