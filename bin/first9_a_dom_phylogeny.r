## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite', 'ggtree', 'RColorBrewer')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the 34 residues
rd34 <- readAAStringSet("data/sp2_adomains_9range_extracted.faa")

# Make a phylogenetic tree of 34 residues
phylo_fin <- treeio::read.newick("data/sp2_adomains_9range_extracted.nwk")
phylo_fin$tip.label
# phylo_fin$node.label <- gsub(phylo_fin$node.label, "X", "")        
#phylo_fin$node.label <- paste0(word(phylo_fin$node.label, sep = "_", 1), ".", word(phylo_fin$node.label, sep = "_", 2))
head(phylo_fin$node.label)

# Make a phylogenetic tree
pl <- ggtree(phylo_fin) #layout = "circular", size = 2)
tofix <- pl$data$label[!pl$data$isTip]
# tofix
# tofix <- gsub("X", "", tofix)    
# tofix <- paste0(word(tofix, sep = "_", 1), ".", word(tofix, sep = "_", 2))
# pl$data$label[!pl$data$isTip] <- tofix
# table(word(pl$data$label[pl$data$isTip], -1, sep = "_")) # looks good
# pl$data$label <- gsub("\\.NA", "", pl$data$label)
# Append metadata
word(pl$data$label, 2, sep = "_")
dd <- data.frame(label = pl$data$label,
                 grp = ifelse(pl$data$isTip, word(pl$data$label, 2, sep = "_"), "NONE"), 
                 size = ifelse(as.numeric(pl$data$label) > 0.75, 0.1, 0))

p2 <- pl %<+% dd
tail(pl$data$label)
table(dd$size)
# Set the color palette
table(dd$grp)

# pal1 <- colorRampPalette(colors=brewer.pal(8, "Set1"))(length(dd$grp))
# pal1[pal1 == "#377EB8"] <- "#92D050"
# pal1[pal1 == "#A65628"] <- "gray68"
# pal1[pal1 == "#F781BF"] <- "#A65628"
# pal1[pal1 == "#4DAF4A"] <- "#377EB8"
# pal1[pal1 == "#FFFF33"] <- "goldenrod"

pal1 <- colorRampPalette(colors=brewer.pal(8, "Spectral"))(length(unique(dd$grp)))

# pal1 <- rainbow(length(dd$grp))

# pal2 <- c(pal1, "#F781BF", "blue1", "darkorchid1", "navy", #"black", 
#           "gray68", "plum1", "blue1",
#           "deepskyblue", "gold", "darkorchid1", 
#           "deeppink2", "lightslateblue",
#           "lightblue2", "darkseagreen1")
# palette(pal2)
# pal2

pdf("data/sp2_9extracted_FastTree_nodes_labeled_monomers.pdf", width = 10, height = 60)
par(mar=c(0.001,0.001,0.001,0.001))
ptree <- p2 +
  aes(color = grp) +
  scale_color_manual(values = pal1) +
  #ggplot2::xlim(-0.1, NA) +
  geom_nodepoint(size = p2$data$size[!p2$data$isTip], color = "gray40") +
  # geom_tiplab(aes(label = word(label, sep = "_", 2)), size = 1.5) +
  geom_tiplab(aes(label = label), size = 1.5) +
  theme(legend.position = "none") +
  geom_treescale(offset = 1, fontsize = 4, x=3, y=1)

ptree
dev.off()
