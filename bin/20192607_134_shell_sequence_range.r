## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite', 'ggtree', 'RColorBrewer')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the 34 residues
# rd34 <- readAAStringSet("data/sp2_adomains_9range_extracted.faa")
rd34 <- readAAStringSet("data/134_shell_sp2_bacterial_untailored_range.faa")

head(rd34)
table(duplicated(rd34))

# Remove duplicates?
# Leave in for now

# HMMAlign
# hmmalign --trim -o 134_shell_sp2_bacterial_untailored_range.sto AMP-binding.hmm 134_shell_sp2_bacterial_untailored_range.faa

# Convert stockholm to FASTA
# http://sequenceconversion.bugaco.com/converter/biology/sequences/stockholm_to_fasta.php



# FastTree
# FastTree -gamma <134_shell_sp2_bacterial_untailored_range_aligned.fasta> 134_shell_sp2_bacterial_untailored.nwk

# names(rd34) <- gsub("_\\[none\\]", "", names(rd34))
names(rd34)

# Make a phylogenetic tree of 34 residues
phylo_fin <- treeio::read.newick("data/134_shell_sp2_bacterial_untailored.nwk")
phylo_fin$tip.label
# phylo_fin$node.label <- gsub(phylo_fin$node.label, "X", "")        
#phylo_fin$node.label <- paste0(word(phylo_fin$node.label, sep = "_", 1), ".", word(phylo_fin$node.label, sep = "_", 2))
head(phylo_fin$node.label)

# Read in the tip labels
euk <- read_csv("data/taxonomy_of_sp2_doms.csv") %>%
  mutate(newnam = gsub("\\\t", "_", sp2_id )) %>%
  mutate(kingdom = case_when(phylum == "Ascomycota" ~ "Fungi",
                             phylum == "Cyanobacteria" ~ "Cyanobacteria",
                             TRUE ~ "Bacteria"))
table(euk$kingdom)

which_fung <- euk %>%
  dplyr::filter(kingdom == "Fungi") %>%
  dplyr::select(newnam) %>%
  pull()

# Make a phylogenetic tree
pl <- ggtree(phylo_fin, layout = "circular", size = 0.5)

pl$data$label
euk$newnam
# euk$newnam <- paste0(euk$newnam, "_")
euk$kingdom[pl$data$label %in% euk$newnam]

drf <- data.frame(euk$kingdom, euk$newnam, stringsAsFactors = F)
dtl <- data.frame(pl$data$label, stringsAsFactors = F)
colnames(drf) <- c("kingdom", "lbl")
colnames(dtl) <- c("lbl")
jnd <- left_join(dtl, drf, by = "lbl")
head(jnd)
nrow(jnd)
length(pl$data$label)

dd <- data.frame(label = pl$data$label,
                 grp = ifelse(pl$data$isTip, word(pl$data$label, 2, sep = "_"), "NONE"), 
                 size = ifelse(as.numeric(pl$data$label) > 0.75, 0.1, 0),
                 kingdom = jnd$kingdom)

jnd$kingdom

p2 <- pl %<+% dd
tail(pl$data$label)
table(dd$size)

# Set the color palette
pal1 <- colorRampPalette(colors=brewer.pal(8, "Set1"))(length(dd$grp))
pal1[pal1 == "#377EB8"] <- "#92D050"
pal1[pal1 == "#A65628"] <- "gray68"
pal1[pal1 == "#F781BF"] <- "#A65628"
pal1[pal1 == "#4DAF4A"] <- "#377EB8"
pal1[pal1 == "#FFFF33"] <- "goldenrod"

pal1 <- colorRampPalette(colors=brewer.pal(8, "Spectral"))(length(unique(dd$grp)))

# pal1 <- rainbow(length(dd$grp))

pal2 <- c(pal1, "#F781BF", "blue1", "darkorchid1", "navy", #"black",
          "gray68", "plum1", "blue1",
          "deepskyblue", "gold", "darkorchid1",
          "deeppink2", "lightslateblue",
          "lightblue2", "darkseagreen1")
palette(pal2)
pal2

pdf("data/134_shell_sp2_tree_circular.pdf", width = 15, height = 15)
par(mar=c(0.001,0.001,0.001,0.001))
ptree <- p2 +
  aes(color = grp) +
  scale_color_manual(values = pal2)
  # ggplot2::xlim(NA, 4) +
  # geom_tippoint(aes(x = 6, shape = kingdom),  size = 0.5)
  # scale_shape_manual(values=c(15, 16, 17, 18)) +
  # geom_nodepoint(size = p2$data$size[!p2$data$isTip], color = "gray40") +
  # geom_tiplab2(aes(label = label), size = 0.5) +
  # #aes(label = word(label, sep = "_", 2)), size = 1.5) +
  # theme(legend.position = "none") +
  # geom_treescale(offset = 1, fontsize = 4, x=3, y=1)

ptree
dev.off()

# Make a phylogenetic tree
pl <- ggtree(phylo_fin, size = 0.5)

dd <- data.frame(label = pl$data$label,
                 grp = ifelse(pl$data$isTip, word(pl$data$label, 2, sep = "_"), "NONE"), 
                 size = ifelse(as.numeric(pl$data$label) > 0.75, 0.1, 0),
                 kingdom = jnd$kingdom)



p2 <- pl %<+% dd
p2

pdf("data/134_shell_sp2_tree.pdf", width = 20, height = 90)
par(mar=c(0.001,0.001,0.001,0.001))
ptree <- p2 +
  aes(color = grp) +
  scale_color_manual(values = pal2) +
  ggplot2::xlim(NA, 7) + 
# geom_tippoint(aes(x = 6, shape = kingdom),  size = 0.5)
# scale_shape_manual(values=c(15, 16, 17, 18)) +
 geom_nodepoint(size = p2$data$size[!p2$data$isTip], color = "gray40") +
 geom_tiplab(aes(label = label), size = 2) +
# #aes(label = word(label, sep = "_", 2)), size = 1.5) +
 # theme(legend.position = "none") +
 geom_treescale(offset = 1, fontsize = 4, x=3, y=1)

ptree
dev.off()

