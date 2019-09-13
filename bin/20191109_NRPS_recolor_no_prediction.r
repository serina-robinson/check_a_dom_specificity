# Install packages
pacman::p_load("Biostrings", "DECIPHER", "ggtree", "randomcoloR", "stringr", "tidyverse", "RColorBrewer")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")
  
# Read in the dataset
aa <- readAAStringSet("data/antismashdb_adoms_nostops_filtered_154extracted.faa")
length(aa)
names(aa) <- gsub(paste0(c("\\|", "\\-", "\\[none\\]", "\\(", "\\)",  "\\+", "\\.1"), collapse = "|"), "_", names(aa))
# names(aa) <- gsub(paste0(c("\\[none\\]", "\\(", "\\)",  "\\+", "\\.1"), collapse = "|"), "", names(aa))
names(aa)

# Remove the gaps
aadf <- as.character(aa)

nogaps <- gsub("-", "", aa)
nogaps_aa <- AAStringSet(nogaps)
# writeXStringSet(nogaps_aa, "data/antismashdb_adoms_nostops_filtered_154extracted_namsfixed_nogaps.faa")

# Remove the bars which are invalid identifiers

# CD-HIT
# source activate cd-hit
# cd-hit -i antismashdb_adoms_nostops_filtered_154extracted.faa -o antismashdb_60 -c 0.6 -n 4

# Check the length of clustered hits
clstr <- readAAStringSet("data/antismashdb_v2_60")
length(clstr)

clstr50 <- readAAStringSet("data/antismashdb_v2_50")
length(clstr50)

clstr40 <- readAAStringSet("data/antismashdb_v2_40")
length(clstr40)

# Find the clustered hits back in the original sequence alignment file
orig_seq <- AAStringSet(aa[names(aa) %in% names(clstr)])
length(orig_seq)
# writeXStringSet(orig_seq, "data/6581_antismashdb_v2_CDHIT60_alned_seqs.faa")
head(orig_seq)
# Alignment (optional) 
# aa.al <- AlignSeqs(clstr) # Note: takes significant time to run
# writeXStringSet(aa.al, "data/6581_seqs_aligned_DECIPHER.faa")

# Combine with sp2 training set
sp2 <- readAAStringSet("data/sp2_bacterial_adomains_154extracted_untailored_unsplit.faa")
length(sp2)

comb <- c(sp2, orig_seq)
# BrowseSeqs(comb, "data/sp2_and_antismashdb_seqs.html")
# writeXStringSet(comb, "data/sp2_and_antismashdb_v2_60_seqs_154.faa")
table(width(comb))
length(comb)

# Submit FastTree job
# fasttree -gamma <sp2_and_antismashdb_v2_60_seqs_154.faa> data/sp2_and_antismashdb_v2_60_seqs_154.nwk

# Submit RaxML job? 6581 sequences will take forever...

# Read in the tree
phylo_fin <- treeio::read.newick("data/sp2_and_antismashdb_v2_60_seqs_154.nwk")
pl <- ggtree(phylo_fin, layout = "circular", size = 0.25)
length(pl$data$label[pl$data$isTip])
# Create annotation data frame
dd <- data.frame(label = pl$data$label,
                 grp = as.character(ifelse(pl$data$isTip, word(pl$data$label, 2, sep = "_"), "UNKNOWN")), stringsAsFactors = F) 
#size = ifelse(as.numeric(pl$data$label) > 0.75, 0.1, 0), stringsAsFactors = F)

# Change the substrate group to be unknown when not experimentally verified
dd$grp[grepl("[[:digit:]]+", dd$grp)] <- "UNKNOWN"
dd$grp[dd$grp == "AMP"] <- "UNKNOWN"

# Change so if the frequency of the substrate is less than 10
dd_annot <- dd %>% 
  add_count(grp) %>%
  dplyr::mutate(broad_grp = case_when((n < 8 & grp != "His") ~ "rare_aa",
                                      TRUE ~ grp)) %>%
  dplyr::mutate(tipsize = case_when(broad_grp != "UNKNOWN" ~ 0.5,
                                    TRUE ~ 0))

dd_annot$tipsize[dd_annot$tipsize == 0] <- NA
sort(table(dd_annot$grp), decreasing = T)

# unique(dd_annot$broad_grp[dd_annot$tipsize != 0.00001])

# Set the color palette
pal1 <- colorRampPalette(colors=brewer.pal(8, "Set1"))(8)
pal1[pal1 == "#377EB8"] <- "#92D050"
pal1[pal1 == "#A65628"] <- "gray42"
pal1[pal1 == "#F781BF"] <- "#A65628"
pal1[pal1 == "#4DAF4A"] <- "#377EB8"
pal1[pal1 == "#FFFF33"] <- "goldenrod"
pal2 <- c(pal1, "#F781BF", "blue1", "darkorchid1", "navy", #"black", 
          "dodgerblue", "plum1", "darkgreen",
          "deepskyblue", "gold", "darkorchid1", 
          "deeppink2", "lightslateblue",
          "lightblue2", "darkseagreen1", "black", "palevioletred4", "gray75", "steelblue4")
length(pal2)

# Make a tree by coloring the branches of interest
p2 <- pl %<+% dd_annot
n <- length(table(p2$data$broad_grp))
set.seed(1234)
pal2 <- distinctColorPalette(n)
pal2[levels(as.factor(p2$data$broad_grp)) == "Gly"] <- "black"
pal2[levels(as.factor(p2$data$broad_grp)) == "UNKNOWN"] <- "gray75"
pal2[levels(as.factor(p2$data$broad_grp)) == "Tyr"] <- "navy"

ptree <- p2 +
  aes(color = broad_grp, alpha = 0.5) +
  scale_color_manual(values = pal2) +
  # eom_tippoint(aes(x = 8), size = ifelse(grepl(paste0(word(names(sp2), sep = "_", 1), collapse = "|"), p2$data$label[p2$data$isTip]), 0.5, NA), alpha = 0.5) +
  geom_point(aes(x = 8), size = ifelse(grepl(paste0(word(names(sp2), sep = "_", 1), collapse = "|"), p2$data$label), 0.5, NA), alpha = 0.5) +
  # geom_tiplab2(aes(label = ifelse(p2$data$broad_grp == "UNKNOWN", "", p2$data$label)), size = 0.5) +
  # geom_tiplab2(aes(label = label, subset = broad_grp != "UNKNOWN"), size = 0.1) +
  geom_tiplab2(aes(label = label, subset = broad_grp != "UNKNOWN"), size = 0.2) +
  # geom_tiplab2(aes(label = label), size = 0.2) +
  theme(legend.position = "right")
# geom_treescale(offset = 1, fontsize = 4, x = 10, y = 10) 


pdf("data/20191109_FastTree_sp2_antismashdb_v2_specificity_colored_no_preds.pdf", width = 10, height = 10)
par(mar=c(0.001,0.001,0.001,0.001))
ptree
dev.off()

# ggsave("data/7610_FastTree_sp2_antismashdb_v2_seqs_unlabeled_circular_nodes_labeled.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)
