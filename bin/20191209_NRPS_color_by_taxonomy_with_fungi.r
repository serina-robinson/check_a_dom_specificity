# Install packages
pacman::p_load("Biostrings", "DECIPHER", "ggtree", "randomcoloR", "stringr", "tidyverse", "RColorBrewer", "readxl")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the taxonomy data
tax <-  read_csv("data/taxonomy_of_sp2_doms.csv", col_names = T) %>%
  #::mutate(label = gsub("\\\t", "_", sp2_id)) %>%
  dplyr::mutate(label_id = word(sp2_id, sep = "\\t", 1))
head(tax)
#tax$label %in% 

# Read in the sp2 data
sp2 <- readAAStringSet("data/sp2.adomains.faa")
names(sp2) <- gsub("\\\t", "_", names(sp2))

# Read in the tree
phylo_fin <- treeio::read.newick("data/7674_sp2_and_antismashdb_v2_60_seqs_154.nwk")
pl <- ggtree(phylo_fin, layout = "circular", size = 0.25)
length(pl$data$label[pl$data$isTip])

# Create annotation data frame
dd <- data.frame(label = pl$data$label,
                 grp = as.character(ifelse(pl$data$isTip, word(pl$data$label, 2, sep = "_"), "UNKNOWN")),
                 label_id = word(pl$data$label, sep = "_", 1), stringsAsFactors = F)
tax$label_id[which(tax$label_id %in% dd$label_id == F)] # all included
dim(tax)

# Change the substrate group to be unknown when not experimentally verified
dd$grp[grepl("[[:digit:]]+", dd$grp)] <- "UNKNOWN"
dd$grp[dd$grp == "AMP"] <- "UNKNOWN"

# Merg with taxonomy data
dd_merg <- dd %>%
  left_join(., tax, by = "label_id")
table(dd_merg$class)
table(!is.na(dd_merg$class))
head(dd_merg$class)

# Change so if the frequency of the substrate is less than 10
dd_annot <- dd_merg %>% 
  add_count(grp) %>%
  dplyr::mutate(broad_grp = case_when((n < 8 & grp != "His") ~ "rare_aa",
                                      TRUE ~ grp)) %>%
  dplyr::mutate(tipsize = case_when(broad_grp != "UNKNOWN" ~ 0.5,
                                    TRUE ~ 0))

dd_annot$tipsize[dd_annot$tipsize == 0] <- NA
dd_annot$class[is.na(dd_annot$class)] <- "Unknown_taxonomy_or_antiSMASHdb_v2_entry"
table(dd_annot$class)

sort(table(dd_annot$grp), decreasing = T)

dd_towrite <- 


length(table(dd_annot$class))
# unique(dd_annot$broad_grp[dd_annot$tipsize != 0.00001])

# Set the color palette
pal1 <- colorRampPalette(colors=brewer.pal(8, "Set1"))(8)
pal1[pal1=="#377EB8"] <- "#FF7F00"
pal2 <- c(pal1[1:4], "dodgerblue", "#F781BF", "darkseagreen1", "goldenrod",  "#A65628", "black", "blue1", "gray75")
pal2
# # pal1[pal1 == "#E41A1C"] <- 
# pal1[pal1 == "#F781BF"] <- "black"
# pal1[pal1 == "#FF7F00"] <- "violet"
# pal1[pal1 == "#FFFF33"] <- "blue1"
# pal1[pal1=="#377EB8"] <- "#FF7F00"
# pal2 <- c(pal1, "dodgerblue",  "gray75")
                    # "deepskyblue", "gold", "darkorchid1","plum1",
                    # "deeppink2", "lightslateblue", "darkorchid1",
                    # "lightblue2", "darkseagreen1", "black", "palevioletred4", "gray75", "steelblue4")
pal2

#pal2[3] <- 

#           "dodgerblue", "plum1", "darkgreen",
#           "deepskyblue", "gold", "darkorchid1", 
#           "deeppink2", "lightslateblue",
#           "lightblue2", "darkseagreen1", "black", "palevioletred4", "gray75", "steelblue4")
# length(pal2)

# Make a tree by coloring the branches of interest
p2 <- pl %<+% dd_annot

#pal2[levels(as.factor(p2$data$broad_grp)) == "UNKNOWN"] <- "gray75"
#pal2[levels(as.factor(p2$data$broad_grp)) == "Tyr"] <- "blue3"

ptree <- p2 +
  aes(color = class, alpha = 0.5) +
  scale_color_manual(values = pal2) +
  # eom_tippoint(aes(x = 8), size = ifelse(grepl(paste0(word(names(sp2), sep = "_", 1), collapse = "|"), p2$data$label[p2$data$isTip]), 0.5, NA), alpha = 0.5) +
  geom_point(aes(x = 8), size = ifelse(grepl(paste0(word(names(sp2), sep = "_", 1), collapse = "|"), p2$data$label), 0.5, NA), alpha = 0.5) +
  # geom_tiplab2(aes(label = ifelse(p2$data$broad_grp == "UNKNOWN", "", p2$data$label)), size = 0.5) +
  # geom_tiplab2(aes(label = label, subset = broad_grp != "UNKNOWN"), size = 0.1) +
  geom_tiplab2(aes(label = paste0(label, "_", class), subset = broad_grp != "UNKNOWN"), size = 0.2) +
  # geom_tiplab2(aes(label = label), size = 0.2) +
  theme(legend.position = "right")
# geom_treescale(offset = 1, fontsize = 4, x = 10, y = 10) 


pdf("data/20191209_FastTree_sp2_antismashdb_v2_specificity_colored_taxonomy_with_fungal.pdf", width = 10, height = 10)
par(mar=c(0.001,0.001,0.001,0.001))
ptree
dev.off()

# ggsave("data/7610_FastTree_sp2_antismashdb_v2_seqs_unlabeled_circular_nodes_labeled.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)
