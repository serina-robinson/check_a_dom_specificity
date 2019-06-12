## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite', 'ggtree', 'RColorBrewer')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the 34 residues
# rd34 <- readAAStringSet("data/sp2_adomains_9range_extracted.faa")
rd34 <- readAAStringSet("~/Documents/Wageningen_UR/github/amplicon_pred/data/sp2_fixed_adomains_9range_extracted.faa")

head(rd34)
table(duplicated(rd34))
# names(rd34) <- gsub("_\\[none\\]", "", names(rd34))
names(rd34)

# Make a phylogenetic tree of 34 residues
phylo_fin <- treeio::read.newick("data/sp2_adomains_9range_extracted.nwk")
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


which_fung
names(rd34)
fungs <- rd34[names(rd34) %in% which_fung]
head(fungs)
length(fungs)
nofungs <- rd34[!names(rd34) %in% which_fung]
length(nofungs)
writeXStringSet(nofungs, "~/Documents/Wageningen_UR/github/amplicon_pred/data/sp2_adomains_9range_fungal_excluded.faa")
grep("D-Ser", names(nofungs))
grep("A", nofungs)
grep("-", nofungs)

checkvitt <- readAAStringSet("~/Documents/Wageningen_UR/github/adenylpred_dev/data/vittorio_amplicons_aligned_NRPS_Adomain.fasta")
checkvitt <- AAStringSet(gsub("X", "-", checkvitt))
checkvitt
grep("X", checkvitt)
length(checkvitt)
# Remove duplicates

writeXStringSet(checkvitt, "~/Documents/Wageningen_UR/github/adenylpred_dev/data/vittorio_amplicons_aligned_NRPS_Adomains_Xsremoved.fasta")
# writeXStringSet(fungs, "data/50_fungal_Adomains.faa")
# write_csv(euk, "datafungal_sequences.csv")

# Make a phylogenetic tree
pl <- ggtree(phylo_fin) #layout = "circular", size = 2)
#tofix <- pl$data$label[!pl$data$isTip]
#tofix <- gsub("X", "", tofix)    
# tofix <- paste0(word(tofix, sep = "_", 1), ".", word(tofix, sep = "_", 2))
#pl$data$label[!pl$data$isTip] <- tofix
#tofix
#table(word(pl$data$label[pl$data$isTip], -1, sep = "_")) # looks good
# pl$data$label <- gsub("\\.NA", "", pl$data$label)
pl$data$label
euk$newnam <- paste0(euk$newnam, "_")
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
pal1

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

pdf("data/sp2_9extracted_tree_with_taxonomy_legend.pdf", width = 20, height = 40)
par(mar=c(0.001,0.001,0.001,0.001))
ptree <- p2 +
  aes(color = grp) +
  scale_color_manual(values = pal1) +
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

