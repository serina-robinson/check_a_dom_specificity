## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite', 'ggtree', 'RColorBrewer', 'randomcoloR')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the sequences
rdall <- readAAStringSet("data/sp2.adomains.faa")
names(rdall) <- gsub("\\\t", "_", names(rdall))
head(names(rdall))

# Align sequences
# al <- AlignSeqs(rdall)
# a_adj <- AdjustAlignment(al)
# writeXStringSet(a_adj, "data/sp2.adomains.DECIPHER.aligned.adjusted.faa")

# Read in the phylogenetic tree
phylo_fin <- treeio::read.newick("data/sp2_full_length_bacterial_fungal_fasttree.nwk")
phylo_fin$tip.label
# phylo_fin$node.label <- gsub(phylo_fin$node.label, "X", "")        
#phylo_fin$node.label <- paste0(word(phylo_fin$node.label, sep = "_", 1), ".", word(phylo_fin$node.label, sep = "_", 2))


# Read in the tip labels
# colnames(euk)
euk <- read_csv("data/taxonomy_of_sp2_doms.csv") %>%
  mutate(newnam = gsub("\\\t", "_", sp2_id )) %>%
  mutate(kingdom = case_when(phylum == "Ascomycota" ~ "Fungi",
                             phylum == "Cyanobacteria" ~ "Cyanobacteria",
                             TRUE ~ "Bacteria"))
table(euk$kingdom)
rd34 <- rdall

which_fung <- euk %>%
  dplyr::filter(kingdom == "Fungi") %>%
  dplyr::select(newnam) %>%
  pull()


# 34 signatures
fungs <- rd34[names(rd34) %in% which_fung]
head(fungs)

nofungs <- rd34[!names(rd34) %in% which_fung]
length(nofungs)

# Full length sequences
nofungs2 <- rdall[!names(rd34) %in% which_fung]
fung2 <- rdall[names(rd34) %in% which_fung]
names(fung2)
# writeXStringSet(nofungs2, "data/sp2_1043_adomains_full_length_no_fungs.faa")

# Make a phylogenetic tree
pl <- ggtree(phylo_fin) #layout = "circular", size = 2)
tofix <- pl$data$label[!pl$data$isTip]
tofix <- gsub("X", "", tofix)    
tofix <- paste0(word(tofix, sep = "_", 1), ".", word(tofix, sep = "_", 2))
pl$data$label[!pl$data$isTip] <- tofix
tofix
table(word(pl$data$label[pl$data$isTip], -1, sep = "_")) # looks good
pl$data$label <- gsub("\\.NA", "", pl$data$label)

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
# pal1 <- colorRampPalette(colors=brewer.pal(8, "Set1"))(length(dd$grp))
# pal1[pal1 == "#377EB8"] <- "#92D050"
# pal1[pal1 == "#A65628"] <- "gray68"
# pal1[pal1 == "#F781BF"] <- "#A65628"
# pal1[pal1 == "#4DAF4A"] <- "#377EB8"
# pal1[pal1 == "#FFFF33"] <- "goldenrod"
# length(unique(dd$grp))
# # pal1 <- colorRampPalette(colors=brewer.pal(8, "Spectral"))(length(unique(dd$grp)))
# 
# # pal1 <- rainbow(length(dd$grp))
# 
# pal2 <- c(pal1, "#F781BF", "blue1", "darkorchid1", "navy", #"black",
#           "gray68", "plum1", "blue1",
#           "deepskyblue", "gold", "darkorchid1",
#           "deeppink2", "lightslateblue",
#           "lightblue2", "darkseagreen1")
pal <- colorRampPalette(brewer.pal(8,"Set1"))(8)
pal2 <- c("gray80", "dodgerblue", "goldenrod", "chartreuse", pal[c(1, 3:5, 8)], "blue", "gold1", "black", distinctColorPalette(length(unique(dd$grp))))

pdf("data/sp2_tree_with_taxonomy_legend.pdf", width = 20, height = 40)
par(mar=c(0.001,0.001,0.001,0.001))
ptree <- p2 +
  aes(color = grp) +
  scale_color_manual(values = pal2) +
  # ggplot2::xlim(NA, 4) +
  geom_tippoint(aes(x = 4, shape = kingdom)) +
  scale_shape_manual(values=c(15, 16, 17, 18)) +
  geom_nodepoint(size = p2$data$size[!p2$data$isTip], color = "gray40") +
  geom_tiplab(aes(label = label), size = 1) +
    #aes(label = word(label, sep = "_", 2)), size = 1.5) +
  theme(legend.position = "right") +
  geom_treescale(offset = 1, fontsize = 4, x=3, y=1)
  
ptree
dev.off()

