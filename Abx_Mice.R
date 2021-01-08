
source("~/Data/microbiome_functions.R")

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
