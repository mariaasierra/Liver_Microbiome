#Abx Experiment in Mice

source("~/microbiome_functions.R")
#Beta Diverstiy Antibiotics Suppl Fig S3
phylo <- import_data("Data/otu_table_abx.biom", mapping = "Data/mapping_file_abx.txt", tree = "Data/rep_set_filtered_abx.tre")
tax_table(phylo) <- tax_table(phylo)[, ranks]
metadata<-sample_data(phylo)

#Metronidazole
#Metro <- subset_samples(phylo, Description == "Metro" | Description == "Control")
Metro <- subset_samples(phylo, Description == "Metronidazole")
Metro <- prune_taxa(taxa_sums(Metro) > 0, Metro)
new.outgroup <- pick_new_outgroup(phy_tree(Metro))
phy_tree(Metro) <- ape::root(phy_tree(Metro), 
                             outgroup=new.outgroup, resolve.root=TRUE)
dist <- phyloseq::distance(Metro, method = "wUnifrac")
perma <- adonis(dist~SampleType, data = as(sample_data(Metro), "data.frame"), permutations = 1000)
disp <- permutest(betadisper(dist, sample_data(Metro)$SampleType), permutations = 1000)
ord <- ordinate(Metro, method = "PCoA", distance = dist)
evals <- ord$values$Eigenvalues
#pdf("Metronidazole_LivervsGut.pdf", height = 5, width = 7)
plot_ordination(Metro, ord, type = "samples", color = "SampleType") +
  geom_point(size = 5,shape = 20) + stat_ellipse() + default.theme + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 17)) +
  scale_color_manual(values = c("red", "blue"))
#dev.off()

#Vancomycin
Van <- subset_samples(phylo, Description == "Vancomycin")
Van <- prune_taxa(taxa_sums(Van) > 0, Van)
new.outgroup <- pick_new_outgroup(phy_tree(Van))
phy_tree(Van) <- ape::root(phy_tree(Van), 
                           outgroup=new.outgroup, resolve.root=TRUE)
dist <- phyloseq::distance(Van, method = "wUnifrac")
perma <- adonis(dist~SampleType, data = as(sample_data(Van), "data.frame"), permutations = 1000)
disp <- permutest(betadisper(dist, sample_data(Van)$SampleType), permutations = 1000)
ord <- ordinate(Van, method = "PCoA", distance = dist)
evals <- ord$values$Eigenvalues
#pdf("Vancomycin_LivervsGut.pdf", height = 5, width = 7)
plot_ordination(Van, ord, type = "samples", color = "SampleType") +
  geom_point(size = 5,shape = 20) + stat_ellipse() + default.theme + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 17)) +
  scale_color_manual(values = c("red", "blue"))
dev.off()
str(ord)

#Neomycin
Neomycin <- subset_samples(phylo, Description == "Neomycin")
Neomycin <- prune_taxa(taxa_sums(Neomycin) > 0, Neomycin)
new.outgroup <- pick_new_outgroup(phy_tree(Neomycin))
phy_tree(Neomycin) <- ape::root(phy_tree(Neomycin), 
                                outgroup=new.outgroup, resolve.root=TRUE)
dist <- phyloseq::distance(Neomycin, method = "wUnifrac")
perma <- adonis(dist~SampleType, data = as(sample_data(Neomycin), "data.frame"), permutations = 1000)
disp <- permutest(betadisper(dist, sample_data(Neomycin)$SampleType), permutations = 1000)
ord <- ordinate(Neomycin, method = "PCoA", distance = dist)
evals <- ord$values$Eigenvalues
#pdf("Neomycin_LivervsGut.pdf", height = 5, width = 7)
plot_ordination(Neomycin, ord, type = "samples", color = "SampleType") +
  geom_point(size = 5,shape = 20) + stat_ellipse() + default.theme + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 17)) +
  scale_color_manual(values = c("red", "blue"))
dev.off()

#Multi_Anti
Multi_Anti <- subset_samples(phylo, Description == "Multi_Anti")
Multi_Anti <- prune_taxa(taxa_sums(Multi_Anti) > 0, Multi_Anti)
new.outgroup <- pick_new_outgroup(phy_tree(Multi_Anti))
phy_tree(Multi_Anti) <- ape::root(phy_tree(Multi_Anti), 
                                  outgroup=new.outgroup, resolve.root=TRUE)
dist <- phyloseq::distance(Multi_Anti, method = "wUnifrac")
perma <- adonis(dist~SampleType, data = as(sample_data(Multi_Anti), "data.frame"), permutations = 1000)
disp <- permutest(betadisper(dist, sample_data(Multi_Anti)$SampleType), permutations = 1000)
ord <- ordinate(Multi_Anti, method = "PCoA", distance = dist)
evals <- ord$values$Eigenvalues
#pdf("Multi_Anti_LivervsGut.pdf", height = 5, width = 7)
plot_ordination(Multi_Anti, ord, type = "samples", color = "SampleType") +
  geom_point(size = 5,shape = 20) + stat_ellipse() + default.theme + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 17)) +
  scale_color_manual(values = c("red", "blue"))
dev.off()

#Rel Abundance Antibiotics groups
# Subset data
control <- subset_samples(phylo, Description == "Control")
control <- prune_taxa(taxa_sums(control) > 0, control)
sample_data(control)

gut <- subset_samples(phylo, SampleType == "Fecal")
gut <- prune_taxa(taxa_sums(gut) > 0, gut)

liver <- subset_samples(phylo, SampleType == "Liver")
liver <- prune_taxa(taxa_sums(liver) > 0, liver)


control <- subset_samples(control, Description != "Multi_Anti")
control <- prune_taxa(taxa_sums(control) > 0, control)

#Barplot Rel. Abundance
glomrank <- "Phylum"
gut.glom <- tax_glom(gut, glomrank)
gut.glom.norm <- transform_sample_counts(gut.glom, function(x) x/sum(x))
gut.glom.md <- psmelt(gut.glom.norm)
gut.glom.dat <- ddply(gut.glom.md, c("Description", "Phylum"),
                      summarise,
                      mean = mean(Abundance), sd = sd(Abundance),
                      N = length(Abundance), se = sd/sqrt(N))
ggplot(gut.glom.dat, aes(x= reorder(Phylum, -mean), y = mean, group = Description, fill = Description)) +
  geom_bar(stat = "identity", position = "dodge", width = .5, colour = "black") +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), 
                stat = "identity", 
                position = position_dodge(.5), width = .1) +
  scale_y_continuous(labels = percent) +
  labs(x = "Phylum",y = "Relative Abundance") +
  scale_fill_brewer(palette = "Set1", "") +
  default.theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 8))
groups <- unique(as.character(gut.glom.md$Description))
gut.phyla <- unique(as.character(gut.glom.md$Phylum))

gut.pvalue <- data.frame(t(combn(groups, 2)), stringsAsFactors = FALSE)
colnames(gut.pvalue) <- c("Group1", "Group2")
gut.pvalue <- data.frame(gut.pvalue,Phylum=rep(gut.phyla,each=nrow(gut.pvalue)),
                         pvalue = 0)


for(i in rownames(gut.pvalue)){
  group1 <- gut.pvalue[i, "Group1"]
  group2 <- gut.pvalue[i, "Group2"]
  wil <- wilcox.test(Abundance~Description, 
                     data = gut.glom.md,
                     subset = Description %in% c(group1, group2) &
                       Phylum == gut.pvalue[i, "Phylum"])
  gut.pvalue[i, "pvalue"] <- wil$p.value
} 


liver.glom <- tax_glom(liver, glomrank)
liver.glom.norm <- transform_sample_counts(liver.glom, function(x) x/sum(x))
liver.glom.md <- psmelt(liver.glom.norm)
liver.glom.dat <- ddply(liver.glom.md, c("Description", "Phylum"),
                        summarise,
                        mean = mean(Abundance), sd = sd(Abundance),
                        N = length(Abundance), se = sd/sqrt(N))
ggplot(liver.glom.dat, aes(x= reorder(Phylum, -mean), y = mean, group = Description, fill = Description)) +
  geom_bar(stat = "identity", position = "dodge", width = .5, colour = "black") +
  geom_errorbar(aes(ymin = mean, ymax = mean + se), 
                stat = "identity", 
                position = position_dodge(.5), width = .1) +
  scale_y_continuous(labels = percent) +
  labs(x = "Phylum",y = "Relative Abundance") +
  scale_fill_brewer(palette = "Set1", "") +
  default.theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.text = element_text(size = 8))

liver.phyla <- unique(as.character(liver.glom.md$Phylum))

liver.pvalue <- data.frame(t(combn(groups, 2)), stringsAsFactors = FALSE)
colnames(liver.pvalue) <- c("Group1", "Group2")
liver.pvalue <- data.frame(liver.pvalue,Phylum=rep(liver.phyla,each=nrow(liver.pvalue)),
                           pvalue = 0)

for(i in rownames(liver.pvalue)){
  group1 <- liver.pvalue[i, "Group1"]
  group2 <- liver.pvalue[i, "Group2"]
  wil <- wilcox.test(Abundance~Description, 
                     data = liver.glom.md,
                     subset = Description %in% c(group1, group2) &
                       Phylum == liver.pvalue[i, "Phylum"])
  liver.pvalue[i, "pvalue"] <- wil$p.value
} 

# Remove "Multi_Anti" Group
liver.alp <- subset_samples(liver, Description != "Multi_Anti")
liver.alp <- prune_taxa(taxa_sums(liver.alp) > 0, liver.alp)

#Alpha diversity indices

measures <- c("Observed", "ACE", "Chao1", "Shannon", "Simpson")
plot_richness(liver.alp, x = "Description", color = "Description", measures = measures) +
  geom_boxplot(width = .5, position = "dodge", size = .8) +
  default.theme +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 8, family = "Arial")) + 
  labs(x = "", y = "Alpha Diversity") +
  scale_color_manual(values = brewer.pal(5, "Set1")[-3])

liver.alp.otu <- as(otu_table(liver.alp), "matrix")
liver.alp.otu <- t(liver.alp.otu)
liver.alp.meta <- as(sample_data(liver.alp), "data.frame")
liver.alp.pd <- pd(liver.alp.otu, phy_tree(liver.alp))
liver.alp.pd$Description <- liver.alp.meta[rownames(liver.alp.pd), "Description"]
ggplot(liver.alp.pd, aes(x = Description, y = PD, color = Description)) +
  geom_point() +
  geom_boxplot(width = .3, position = "dodge", size = .8) + 
  default.theme +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 8, family = "Arial"))+
  scale_color_manual(values = brewer.pal(5, "Set1")[-3])

