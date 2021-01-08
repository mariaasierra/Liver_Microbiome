# Author: Qianhao Li
# Date: 04/25/2019


library(phyloseq);packageVersion("phyloseq")  #V1.27.0
library(ggplot2);packageVersion("ggplot2")  #V3.1.0
library(plyr); packageVersion("plyr") #V1.8.4
library(extrafont)
library(scales)
library(RColorBrewer)
library(ggpubr)

# Directory for figures
plot_path <- "./Plots"
if (!file.exists(plot_path)){
	dir.create(plot_path)
}


# Default Theme for Figures
default.theme <- theme(text = element_text(family = "Arial", color = "black"),
                       axis.text = element_text(size = 10),
                       axis.title = element_text(size = 12, face = "bold"),
                       legend.text = element_text(size = 8),
                       legend.key = element_rect(fill = "white"),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), 
                       panel.background = element_blank(),
                       axis.line = element_line(color = "black"))

# Define function to import data
import_data <- function(otu = "otu_table.biom",
                        tree = "tree",
                        mapping = "mapping_file.txt",
                        parseFunction = parse_taxonomy_greengenes){
  # import the biom table to a phylosqe object
  # otu is the otu table with the OTU ID and Abundance
  # mapping is the mapping file with all the information about samples
  
  otutable <- import_biom(otu,
                          treefilename = tree,
                          parseFunction = parseFunction)
  
  meta <- import_qiime_sample_data(mapping)
  colnames(meta)[1] <- "SampleID"
  phylo <- merge_phyloseq(otutable, meta)
  return(phylo)
}

reads_sum <- function(phyloseq){
  # To Show the total number of taxas and samples in the project
  
  # Create a data frame includes the counts of each taxa
  readsum <- data.frame(nreads = sort(taxa_sums(phyloseq), TRUE),
                        sorted = 1:ntaxa(phyloseq), 
                        type = "OTUs")
  
  # Add the count of taxas for each sample
  readsum <- rbind(readsum,
                   data.frame(nreads = sort(sample_sums(phyloseq), TRUE),
                              sorted = 1:nsamples(phyloseq),
                              type = "Samples"))
  
  readsum_map <- aes_string(x = "sorted", y = "nreads")
  
  p <- ggplot(readsum, readsum_map) +
    geom_bar(stat = "identity") +
    scale_y_log10() +
    facet_wrap(~type, 1, scales = "free") +
    ggtitle("Total number of reads")
  return (p)
}


taxa_prev <- function(phy_obj, s = 1){
  # "s" means at least how many samples should have the specific OTUs
  dat <- as(otu_table(phy_obj), "matrix")
  res <- rowSums(dat > 0)
  to_be_kept <- names(res[res >= s])
  return(to_be_kept)
}

filter_by_nsample <- function(phy_obj, max = 10){
  dat <- data.frame(n = 0, otus = ntaxa(phy_obj), seqs = sum(taxa_sums(phy_obj)))
  for (s in seq(max)){
    to_be_kept <- taxa_prev(phy_obj, s = s)
    sub <- prune_taxa(to_be_kept, phy_obj)
    dat <- rbind(dat, 
                 data.frame(n = s, otus = ntaxa(sub), seqs = sum(taxa_sums(sub))))
  }
  par(mar = c(5,5,3,5))
  plot.new()
  plot(dat$n, dat$otus, type = "b", ylab = "Num of OTUs", xlab = ">= Num of Samples", col = "red")
  par(new = TRUE)
  plot(dat$n, dat$seqs, type = "b", axes = FALSE, xlab = "", ylab = "", col = "blue")
  axis(side=4, at = pretty(range(dat$seqs)), col = "blue")
  mtext("Num of Seqs", side = 4, col = "blue", line = 4)
  return(dat)
}

filter_by_nreads <- function(phy_obj, ratio = 0.00005){
	dat <- data.frame(n = 0, otus = ntaxa(phy_obj), seqs = sum(taxa_sums(phy_obj)))
	for (s in seq(sum(taxa_sums(phy_obj))*ratio)){
		sub <- prune_taxa(taxa_sums(phy_obj) >= s, phy_obj )
		dat <- rbind(dat, 
                 data.frame(n = s, otus = ntaxa(sub), seqs = sum(taxa_sums(sub))))
	}
	par(mar = c(5,5,3,5))
	plot.new()
  	plot(dat$n, dat$otus, type = "b", ylab = "Num of OTUs", xlab = ">= Num of Reads", col = "red")
  	par(new = TRUE)
  	plot(dat$n, dat$seqs, type = "b", axes = FALSE, xlab = "", ylab = "", col = "blue")
  	mtext("Num of Seqs", side = 4, col = "blue", line = 4)
  	axis(side=4, at = pretty(range(dat$seqs)), col = "blue")
  	return(dat)
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

lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

darken <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}

ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

plot_alpha <- function(phy_obj, x = NULL,
             measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "PD"),
             cols = NULL,
             d.theme = NULL){
  # Create Box Plot of Alpha Diversity
    # "x" would be the x-axis for the plot
    # "measures" are the alpha diveristy algorithms to be calculated
    
    if ("PD" %in% measures){
      library(picante)
      set.seed(1234)
      new.outgroup <- pick_new_outgroup(phy_tree(phy_obj))
      phy_tree(phy_obj) <- ape::root(phy_tree(phy_obj),
                                   outgroup = new.outgroup, resolve.root = TRUE)

      alp.data <- estimate_richness(phy_obj, measures = measures[-match("PD", measures)])
      otu <- t(as(otu_table(phy_obj), "matrix"))
      alp.pd <- pd(otu, phy_tree(phy_obj))
      alp.data$PD <- alp.pd$PD
    }else{
      alp.data <- estimate_richness(phy_obj, measures = measures)
    }
    rownames(alp.data) <- sample_data(phy_obj)$SampleID
    alp.data <- merge(alp.data, sample_data(phy_obj), by = "row.names")

    plot.list <- list()
    for (i in seq(length(measures))){
        p <- ggplot(alp.data, aes_string(x = x, y = measures[i], fill = x)) +
              geom_boxplot(width = .5, position = "dodge") +
              scale_fill_manual(values = cols) +
              d.theme +
              theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      plot.list[[i]] <- p
    }
    return(plot.list)
}

alpha_pvalue <- function(phy_obj,
             measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "PD"),
             group = NULL,
             test = "wilcox"){
  if ("PD" %in% measures){
      library(picante)
      set.seed(1234)
      new.outgroup <- pick_new_outgroup(phy_tree(phy_obj))
      phy_tree(phy_obj) <- ape::root(phy_tree(phy_obj),
                                   outgroup = new.outgroup, resolve.root = TRUE)

      alp.data <- estimate_richness(phy_obj, measures = measures[-match("PD", measures)])
      otu <- t(as(otu_table(phy_obj), "matrix"))
      alp.pd <- pd(otu, phy_tree(phy_obj))
      alp.data$PD <- alp.pd$PD
    }else{
      alp.data <- estimate_richness(phy_obj, measures = measures)
    }
    rownames(alp.data) <- sample_data(phy_obj)$SampleID
    alp.data <- merge(alp.data, sample_data(phy_obj), by = "row.names")

  if (test == "wilcox"){
    if (length(unique(alp.data[[group]])) < 2){
      print("There should be more than 1 factors to be compared. Please check the data.")
      break
  } else if(length(unique(alp.data[[group]])) == 2){
      alp.pvalue <- data.frame(Measures = measures, pvalue = 0, row.names = measures)
      for (m in measures){
        test <- wilcox.test(get(m)~get(group), data = alp.data)
        alp.pvalue[m, "pvalue"] <- test$p.value
      }
    } else if(length(unique(alp.data[[group]])) > 2){
      pvalue.list <- list()
      for (m in measures){
         pvalues <- pairwise.wilcox.test(alp.data[[m]], alp.data[[group]], p.adj = "none")$p.value
         pvalues <- data.frame(melt(pvalues), Measures = m)
         colnames(pvalues)[1:3] <- c("Group1", "Group2", "pvalues")
         pvalues <- pvalues[!is.na(pvalues$pvalues), ]
         pvalue.list[[m]] <- pvalues[, c("Group1", "Group2", "Measures", "pvalues")]
      }
      alp.pvalue <- ldply(pvalue.list, rbind)
      alp.pvalue <- alp.pvalue[, -1]
    }
  } else{
    if (length(unique(alp.data[[group]])) < 2){
      print("There should be more than 1 factors to be compared. Please check the data.")
      break
    } else if(length(unique(alp.data[[group]])) == 2){
      alp.pvalue <- data.frame(Measures = measures, pvalue = 0, row.names = measures)
      for (m in measures){
        test <- t.test(get(m)~get(group), data = alp.data)
        alp.pvalue[m, "pvalue"] <- test$p.value
      }
    } else if(length(unique(alp.data[[group]])) > 2){
      pvalue.list <- list()
      for (m in measures){
        pvalues <- pairwise.t.test(alp.data[[m]], alp.data[[group]], p.adj = "none")$p.value
        pvalues <- data.frame(melt(pvalues), Measures = m)
        colnames(pvalues)[1:3] <- c("Group1", "Group2", "pvalues")
        pvalues <- pvalues[!is.na(pvalues$pvalues), ]
        pvalue.list[[m]] <- pvalues[, c("Group1", "Group2", "Measures", "pvalues")]
      }
      alp.pvalue <- ldply(pvalue.list, rbind)
      alp.pvalue <- alp.pvalue[, -1]
    }
  }
  output <- list()
  output[[1]] <- alp.data
  output[[2]] <- alp.pvalue
  return(output)
}

plot_beta <- function(phy_obj,
            group = NULL,
            dis.method = "wUnifrac",
            beta.method = "PCoA",
            cols = NULL,
            ellipse = TRUE){
  require(vegan)
  # Root the tree
  if (dis.method %in% c("wUnifrac", "Unifrac")){
  	set.seed(1234)
  	new.outgroup <- pick_new_outgroup(phy_tree(phy_obj))
  	phy_tree(phy_obj) <- ape::root(phy_tree(phy_obj),
                                 outgroup = new.outgroup, resolve.root = TRUE)
  }

  phy.dist <- phyloseq::distance(phy_obj, method = dis.method)

  set.seed(1234)
  perma <- adonis(phy.dist~get(group), data = as(sample_data(phy_obj), "data.frame"), permutations = 1000)

  phy.ord <- ordinate(phy_obj, method = beta.method, distance = phy.dist)
  p <- plot_ordination(phy_obj, phy.ord, type = "samples", color = group) +
      geom_point(size = 5, shape = 20) +
      scale_color_manual(values = cols, "") +
      scale_shape(solid = FALSE)
  if(ellipse){p <- p + stat_ellipse()}

  beta_div <- list()
  beta_div[["pvalue"]] <- perma
  beta_div[["plot"]] <- p
  return(beta_div)
}

beta_3d <- function(phy_obj,
                group = NULL,
                dis.method = "wUnifrac",
                beta.method = "PCoA",
                cols = NULL,
                sphere.size = 1, radius = 1){
  # cols should be named vector with colors' names
  require(rgl)
  require(vegan)

  # Root the Tree
  new.outgroup <- pick_new_outgroup(phy_tree(phy_obj))
  phy_tree(phy_obj) <- ape::root(phy_tree(phy_obj),
                                   outgroup = new.outgroup, resolve.root = TRUE)

  phy.dist <- phyloseq::distance(phy_obj, method = dis.method)
  phy.ord <- ordinate(phy_obj, method = beta.method, distance = phy.dist)

  set.seed(1234)
  perma <- adonis(phy.dist~get(group), data = as(sample_data(phy_obj), "data.frame"), permutations = 1000)
  print(paste0("The pvalue is ", perma$aov.tab$`Pr(>F)`[1]))

  x <- phy.ord$vectors[,1]; y <- phy.ord$vectors[,2]; z <- phy.ord$vectors[,3]
  # Scale the data
  x1 <- (x - min(x))/(max(x) - min(x))
  y1 <- (y - min(y))/(max(y) - min(y))
  z1 <- (z - min(z))/(max(z) - min(z))
  # Build 3D plot
  next3d()
  par3d(windowRect = 50 + c(0, 0, 200, 200), fontname = "Arial", cex = 3)
  rgl.bg(color = "white", fogtype = "none")

  size <- sphere.size*((100/length(x1))^(1/3))*0.015
  radius <- radius/median(radius)
  rgl.spheres(x1, y1, z1, r = 0.02, color = cols[as.character(sample_data(phy_obj)[[group]])], radius = size * radius) 
  # Add x, y, and z Axes
  var1 <- percent(phy.ord$values$Relative_eig[1])
  var2 <- percent(phy.ord$values$Relative_eig[2])
  var3 <- percent(phy.ord$values$Relative_eig[3])
  xlab <- paste0("Axis.1 [", var1, "]")
  ylab <- paste0("Axis.2 [", var2, "]")
  zlab <- paste0("Axis.3 [", var3, "]")
  rgl.lines(c(0, 1), c(0, 0), c(0, 0), color = "black")
  rgl.lines(c(0, 0), c(0, 1), c(0, 0), color = "black")
  rgl.lines(c(0, 0), c(0, 0), c(0, 1), color = "black")
  rgl.texts(rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1)), text = c(xlab, ylab, zlab), color = "black", adj = c(0.5, -0.8))

  # Add ellipse
  groups <- sample_data(phy_obj)[[group]]
  levs <- levels(groups)
  for (i in 1:length(levs)){
    g <- levs[i]
    selected <- groups == g
    xx <- x1[selected]; yy <- y1[selected]; zz <- z1[selected]
    ellips <- ellipse3d(cov(cbind(xx, yy, zz)), centre = c(mean(xx), mean(yy), mean(zz)), level = 0.95)
    shade3d(ellips, col = cols[g], alpha = 0.1, lit = FALSE)
    wire3d(ellips, col = cols[g], alpha = 0.1, lit = FALSE)
  }
}

plot_abund <- function(phy_obj,
             glomrank = "Phylum",
             group = NULL,
             cols = NULL,
             facet = FALSE,
             merge_unclass = FALSE){
  # Used to plot relative abundance at Phylum, Genus and Species levels

  phy.glom <- tax_glom(phy_obj, glomrank)
  phy.glom.norm <- transform_sample_counts(phy.glom, function(x) x/sum(x))
  if (merge_unclass){
    phy.glom.norm <- merge_taxa(phy.glom.norm,
                  taxa_names(phy.glom.norm)[tax_table(phy.glom.norm)[, glomrank] %in% c("unidentified", "uncultured", "uncultured bacterium")])
    tax_table(phy.glom.norm)[is.na(tax_table(phy.glom.norm)[, glomrank]), glomrank] <- "UnClassified"
  }
  if (glomrank %in% c("Genus", "Species")){
    top20 <- names(sort(taxa_sums(phy.glom.norm), TRUE))[1:5]
    tax_table(phy.glom.norm)[!(rownames(tax_table(phy.glom.norm)) %in% top20), glomrank] <- "Others"
    phy.glom.norm <- merge_taxa(phy.glom.norm,
                  taxa_names(phy.glom.norm)[tax_table(phy.glom.norm)[, glomrank] == "Others"])
    tax_table(phy.glom.norm)[is.na(tax_table(phy.glom.norm)[, glomrank]), glomrank] <- "Others"
  }
  
  md <- psmelt(phy.glom.norm)
  dat <- ddply(md, c(group, glomrank), summarise,
       mean = mean(Abundance), sd = sd(Abundance),
       N = length(Abundance), se = sd/sqrt(N))

  if (glomrank %in% c("Genus", "Species")){
    
    index <- ddply(dat, c(glomrank), summarise, mean = mean(mean))
    index <- as.character(index[order(index$mean), glomrank])
    dat[[glomrank]] <- factor(dat[[glomrank]], levels = index)

    my_cols <- c("#8CD0FF","#201625","#04F757","#FFF69F","#C2FF99","#809693","#DDEFFF",
           "#6F0062","#FF913F","#3A2465","#FF4A46","#997D87","#452C2C","#FFFF00",
           "#A1C299","#FFB500","#7B4F4B","#CB7E98","#5A0007","#7900D7","#013349")
    p <- ggplot(dat, aes_string(x = group, y = "mean", fill = glomrank)) +
        geom_bar(stat = "identity", position = "stack", width = .3) +
        scale_fill_manual(values = my_cols) +
        theme(legend.position = "right") +
        scale_y_continuous(labels = percent) + 
        ylab("Relative Abundance")
  } else{
    if (facet){
      p <- ggplot(dat, aes_string(x = group, y = "mean", fill = group)) +
        geom_bar(position = "dodge", stat = "identity", width = .5, color = "black") +
        geom_errorbar(aes(ymin = mean, ymax = mean + se), position = position_dodge(.5), width = .1) +
        scale_fill_manual(values = cols, "") +
        scale_y_continuous(labels = percent) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title.x = element_blank(), 
              legend.position = "right",
              strip.text = element_text(size = 10),
              strip.background = element_rect(fill = "white")) +
        ylab("Relative Abundance") +
        facet_wrap(~get(glomrank), scales = "free_y")
    } else{
      p <- ggplot(dat, aes_string(x = paste0("reorder(", glomrank, ", -mean)"), y = "mean", fill = group)) +
        geom_bar(position = "dodge", stat = "identity", width = .5, color = "black") +
        geom_errorbar(aes(ymin = mean, ymax = mean + se), position = position_dodge(.5), width = .1) +
        scale_fill_manual(values = cols, "") +
        scale_y_continuous(labels = percent) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        ylab("Relative Abundance")
    } 
  }
  output <- list()
  output[["data"]] <- md
  output[["plot"]] <- p
  return(output)
}

taxa_pvalue <- function(phy_obj,
                        glomrank = "Phylum",
                        group = NULL){
  phy.glom <- tax_glom(phy_obj, glomrank)
  phy.glom.norm <- transform_sample_counts(phy.glom, function(x) x/sum(x))
  if (glomrank %in% c("Genus", "Species")){
  	phy.glom.norm <- merge_taxa(phy.glom.norm,
                  taxa_names(phy.glom.norm)[tax_table(phy.glom.norm)[, glomrank] == "unidentified"])
  	tax_table(phy.glom.norm)[is.na(tax_table(phy.glom.norm)[, glomrank]), glomrank] <- "UnClassified"
  }

  md <- psmelt(phy.glom.norm)

  taxas <- as.character(unique(md[[glomrank]]))
  if(length(unique(md[[group]])) == 2){
    p.table <- data.frame(Taxas = taxas, pvalue = 0, row.names = taxas)
    for (t in taxas){
      test <- wilcox.test(Abundance~get(group), data = md, subset = md[[glomrank]] == t)
      p.table[t, "pvalue"] <- test$p.value
    }
  } else if (length(unique(md[[group]])) > 2){
    p.table.list <- list()
    for (t in taxas){
      sub <- md[md[[glomrank]] == t, ]
      pvalues <- pairwise.wilcox.test(sub$Abundance, sub[[group]], p.adj = "none")$p.value
      pvalues <- data.frame(melt(pvalues), Taxas = t)
      colnames(pvalues)[1:3] <- c("Group1", "Group2", "pvalues")
      pvalues <- pvalues[!is.na(pvalues$pvalues), ]
      p.table.list[[t]] <- pvalues[, c("Group1", "Group2", "Taxas", "pvalues")]
    }
    p.table <- ldply(p.table.list, rbind)
    p.table <- p.table[, -1]
  }
  return(p.table)
}

phylum_box <- function(phy_obj, group = NULL, cols = NULL){
  phy.glom <- tax_glom(phy_obj, "Phylum")
  phy.glom.norm <- transform_sample_counts(phy.glom, function(x) x/sum(x))

  md <- psmelt(phy.glom.norm)
  p <- ggplot(md, aes_string(x = group, y = "Abundance", fill = group)) +
       geom_boxplot(position = "dodge",width = .5) +
       scale_fill_manual(values = cols, "") +
       scale_y_continuous(labels = percent) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
       ylab("Relative Abundance") +
       facet_wrap(~Phylum, scale = "free_y")
  return(p)
}

genus_heatmap <- function(phy_obj, 
                          group = "Group", 
                          sub.genus = 20, 
                          scale.color = colorRampPalette(colors = c("red","black","green"))(200), 
                          limit = 2, 
                          side.cols = NULL, 
                          save = TRUE, 
                          output = "genus_heatmap.tiff", 
                          width = 10, 
                          height = 10, 
                          legend.position = "topright", 
                          column.tree = FALSE, 
                          margins = c(10, 15), 
                          cexRow = 1.5, 
                          key.size = 0.6){
  
  gm <- tax_glom(phy_obj, "Genus")
  gm.norm <- transform_sample_counts(gm, function(x) x/sum(x))

  otu <- as(otu_table(gm.norm), "matrix")
  meta <- as(sample_data(gm), "data.frame")
  meta <- meta[order(meta[[group]], meta$SampleID), ]
  tax <- data.frame(as(tax_table(gm), "matrix"))

  otu.center <- apply(otu, 1, scale)
  rownames(otu.center) <- colnames(otu)
  otu.center <- otu.center[match(rownames(meta), rownames(otu.center)), ]
  otu.center <- t(otu.center)

  if (is.numeric(sub.genus)){
    top.otus <- names(sort(taxa_sums(gm.norm), TRUE))[1:sub.genus]
    hmp.dat <- otu.center[top.otus, ]
    rownames(hmp.dat) <- tax[rownames(hmp.dat), "Genus"]
  }else{
    hmp.dat <- otu.center
    rownames(hmp.dat) <- tax[rownames(hmp.dat), "Genus"]
    hmp.dat <- hmp.dat[sub.genus, ]}
  
  hmp.dat[hmp.dat > limit] <- limit
  hmp.dat[hmp.dat < -limit] <- -limit

  colbar <- side.cols[meta[[group]]]
  dd <- ifelse(column.tree, "both", "row")

  if (save){
    tiff(output, res = 300, units = "in", width = width, height = height)
  }
  heatmap.2(hmp.dat, Rowv = TRUE, Colv = column.tree, dendrogram = dd, scale = "none", trace = "none", col = scale.color, labCol = NA, key.xlab = "z-score", density.info = "none", margins = margins, ColSideColors = colbar, key.title = NA, cexRow = cexRow, main = NULL, keysize = key.size)
  if(is.character(legend.position)){
    legend(legend.position, xpd = TRUE, bty = "n", legend = names(side.cols), fill = side.cols, pch = ".", horiz = FALSE)
  }else{
    legend(legend.position[1], legend.position[2], xpd = TRUE, bty = "n", legend = names(side.cols), fill = side.cols, pch = ".", horiz = FALSE)
  }
  if(save){dev.off()}
}

# import taxonomy information
parse_taxonomy_default = function(char.vec){
  # Remove any leading empty space
  char.vec = gsub("^[[:space:]]{1,}", "", char.vec)
  # Remove any trailing space
  char.vec = gsub("[[:space:]]{1,}$", "", char.vec)
  if( length(char.vec) > 0 ){
    # Add dummy element (rank) name
    names(char.vec) = paste("Rank", 1:length(char.vec), sep="")
  } else {
    warning("Empty taxonomy vector encountered.")
  }
  return(char.vec)
}

#' @examples
#' taxvec3 = c("D_0__Bacteria", "D_1__Firmicutes", "D_2__Bacilli", "D_3__Staphylococcaceae")
#' parse_taxonomy_silva(taxvec3)
parse_taxonomy_silva_128 <- function(char.vec){
  # Use default to assign names to elements in case problem with greengenes prefix
  char.vec = parse_taxonomy_default(char.vec)
  # Check for unassigned taxa
  if (char.vec["Rank1"] == "Unassigned") {
    char.vec <- c(Rank1="D_0__Unassigned", Rank2="D_1__Unassigned", Rank3="D_2__Unassigned", Rank4="D_3__Unassigned",
                  Rank5="D_4__Unassigned", Rank6="D_5__Unassigned", Rank7="D_6__Unassigned")
  }
  # Define the meaning of each prefix according to SILVA taxonomy
  Tranks = c(D_0="Kingdom", D_1="Phylum", D_2="Class", D_3="Order", D_4="Family", D_5="Genus", D_6="Species")
  # Check for prefix using regexp, warn if there were none. trim indices, ti
  ti = grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
  if( length(ti) == 0L ){
    warning(
      "No silva prefixes were found. \n",
      "Consider using parse_taxonomy_delfault() instead if true for all OTUs. \n",
      "Dummy ranks may be included among taxonomic ranks now."
    )
    # Will want to return without further modifying char.vec
    taxvec = char.vec
    # Replace names of taxvec according to prefix, if any present...
  } else {
    # Format character vectors for Ambiguous taxa
    if( length(ti) < 7 ){
      for (key in names(char.vec)) {
        if ( char.vec[key] == "Ambiguous_taxa" ) {
          tax_no <- (as.numeric(substr(key, 5, 5)) - 1)
          char.vec[key] = sprintf("D_%s__Ambiguous_taxa", tax_no)
        }
      }
      # Reset the trimmed indicies if Ambiguous taxa
      ti = grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
    }
    # Remove prefix using sub-"" regexp, call result taxvec
    taxvec = gsub("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", "", char.vec)
    # Define the ranks that will be replaced
    repranks = Tranks[substr(char.vec[ti], 1, 3)]
    # Replace, being sure to avoid prefixes notK present in Tranks
    names(taxvec)[ti[!is.na(repranks)]] = repranks[!is.na(repranks)]
  }
  return(taxvec)
}

