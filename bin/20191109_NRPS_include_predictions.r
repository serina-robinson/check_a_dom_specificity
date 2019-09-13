# Install packages
pacman::p_load("Biostrings", "DECIPHER", "ggtree", "randomcoloR", "stringr", "tidyverse", "RColorBrewer")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the dataset
aa <- readAAStringSet("data/antismashdb_adoms_nostops_filtered_154extracted.faa")
names(aa) <- gsub(paste0(c("\\|", "\\-", "\\[none\\]", "\\(", "\\)",  "\\+", "\\.1"), collapse = "|"), "_", names(aa))

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
                 grp = as.character(ifelse(pl$data$isTip, word(pl$data$label, 2, sep = "_"), NA)), stringsAsFactors = F) 
#size = ifelse(as.numeric(pl$data$label) > 0.75, 0.1, 0), stringsAsFactors = F)

# Change the substrate group to be unknown when not experimentally verified
dd$grp[grepl("[[:digit:]]+", dd$grp)] <-  NA
dd$grp[dd$grp == "AMP"] <- NA

# Change so if the frequency of the substrate is less than 10
dd$grp[dd$grp == "di-OH-Bz"] <- "Dhb"
names(sp2)[grep("bOH-Tyr", names(sp2))]

dd$grp[grep("bOH-Tyr", pl$data$label)] <- "Bht"
head(dd)

dd_annot <- dd %>% 
  add_count(grp) %>%
  dplyr::mutate(broad_grp = case_when(n < 8 & !grp %in% c("His", "Pip", "Bht", "Dhb") ~ "rare_aa",
                                      grp == "Abu|Ival" ~ "rare_aa",
                                      # grp == "His" ~ "His",
                                      # grp == "Pip" ~ "Pip",
                                      # grp ==  "Bht" ~ "Bht",
                                      # grp == "Dhb" ~ "Dhb",
                                      TRUE ~ grp)) %>%
  dplyr::mutate(tipsize = case_when(broad_grp != "UNKNOWN" ~ 0.5,
                                    TRUE ~ 0)) %>%
  dplyr::mutate(substrate = broad_grp) %>%
  dplyr::mutate(test_label = gsub("_", "", label))


dd_annot$tipsize[dd_annot$tipsize == 0] <- NA
sort(table(dd_annot$broad_grp), decreasing = T)

# Read in Barbara's predictions from SANDPUMA1
sandpuma1 <- read_delim("data/sandpuma_output.tsv", delim = "\t", col_names = F)
conflicts <- read_delim("data/proteinogenic_conflict_ASM_SVM.txt", delim = "\t", col_names = F) # conflicts...maybe color black??
no_call_ASM <- read_delim("data/proteinogenic_novel_ASM.txt", delim = "\t", col_names = F) # should be left gray?
rare_maxcount15 <- read_delim("data/proteinogenic_rare_maxcount_15.txt", delim = "\t", col_names = F) # add colorful circles??


# Filter out predictions we don't want to color
tocol <- sandpuma1 %>%
  tidyr::separate_rows(X1, sep = "\\|\\|\\|") %>%
  dplyr::filter(!X1 %in% conflicts$X1) %>%
  dplyr::filter(!X1 %in% no_call_ASM$X1) %>%
  dplyr::filter(X2 == "SVM") %>%
  dplyr::mutate(substrate = tools::toTitleCase(X3)) %>%
  dplyr::mutate(fixnam = gsub(paste0(c("\\|", "\\-", "\\[none\\]", "\\(", "\\)",  "\\+", "\\.1"), collapse = "|"), "_", X1)) %>%
  dplyr::mutate(test_label = gsub("_", "", fixnam)) %>%
  dplyr::mutate(rare_maxcount = ifelse(X1 %in% rare_maxcount15$X1, X1, NA))
table(tocol$X1 %in% rare_maxcount15$X1)
table(tocol$rare_maxcount)

tocol$substrate[tocol$substrate == "N/a"] <- NA
table(tocol$substrate)

# Intersection
table(tocol$test_label %in% dd_annot$test_label) # 3373 intersect 

# Combine everything
dd_comb <- dd_annot %>%
  left_join(., tocol, by = c("test_label")) %>%
  dplyr::mutate(substrate = ifelse(is.na(substrate.x), substrate.y, substrate.x)) %>%
  dplyr::mutate(substrate = case_when(is.na(substrate) ~ "UNKNOWN",
                                           TRUE ~ substrate))

# Annotate the tree
p2 <- pl %<+% dd_comb
table(dd_comb$substrate)

# Set color palette
n <- length(table(p2$data$substrate))
n
set.seed(1234)
pal2 <- distinctColorPalette(n)
pal2[levels(as.factor(p2$data$substrate)) == "Glu"] <- "maroon"
pal2[levels(as.factor(p2$data$substrate)) == "UNKNOWN"] <- "gray75"
pal2[levels(as.factor(p2$data$substrate)) == "rare_aa"] <- "blue3"
pal2[levels(as.factor(p2$data$substrate)) == "Ala"] <- "orange4"
# pal2[levels(as.factor(p2$data$substrate)) == "Dab"] <- "forestgreen"
pal2

ptree <- p2 +
  aes(color = substrate, alpha = 0.5) +
  scale_color_manual(values = pal2) +
  geom_point(aes(x = 8), size = ifelse(p2$data$substrate != "UNKNOWN", 0.5, NA), alpha = 0.5) +
  geom_point(aes(x = 8.25), size = ifelse(is.na(p2$data$rare_maxcount), NA, 1), alpha = 0.5, shape = 15, color = "black") +
  # geom_tippoint(aes(x = 8), size = ifelse(grepl(paste0(word(names(sp2), sep = "_", 1), collapse = "|"), p2$data$label[p2$data$isTip]), 0.5, NA), alpha = 0.5) +
  # geom_tiplab2(aes(label = label, subset = substrate != "UNKNOWN"), size = 0.2) +
  geom_tiplab2(aes(label = label, subset = !is.na(rare_maxcount)), size = 0.2, color = "black") +
  geom_tiplab2(aes(label = label, subset = grepl(paste0(word(names(sp2), sep = "_", 1), collapse = "|"), p2$data$label[p2$data$isTip])), size = 0.2) +
  theme(legend.position = "right")
# geom_treescale(offset = 1, fontsize = 4, x = 10, y = 10) 

pdf("data/20191109_FastTree_sp2_antismashdb_v2_specificity_colored_SANDPUMA_SVM_predictions_labeled.pdf", width = 10, height = 10)
par(mar=c(0.001,0.001,0.001,0.001))
ptree
dev.off()

# ggsave("data/7610_FastTree_sp2_antismashdb_v2_seqs_unlabeled_circular_nodes_labeled.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

### Now using the same coloring, but ignoring the SVM predictions
# ptree <- p2 +
#   aes(color = substrate, alpha = 0.5) +
#   scale_color_manual(values = pal2) +
#  # geom_point(aes(x = 8), size = ifelse(p2$data$substrate != "UNKNOWN", 0.5, NA), alpha = 0.5) +
#   geom_tippoint(aes(x = 8), size = ifelse(grepl(paste0(word(names(sp2), sep = "_", 1), collapse = "|"), p2$data$label[p2$data$isTip]), 0.5, NA), alpha = 0.5) +
#   # geom_tiplab2(aes(label = label, subset = substrate != "UNKNOWN"), size = 0.2) +
#   geom_tiplab2(aes(label = label, subset = grepl(paste0(word(names(sp2), sep = "_", 1), collapse = "|"), p2$data$label[p2$data$isTip])), size = 0.2) +
#   theme(legend.position = "right")
# # geom_treescale(offset = 1, fontsize = 4, x = 10, y = 10)
# 
# 
# pdf("data/20191109_FastTree_sp2_antismashdb_v2_specificity_colored_no_preds.pdf", width = 10, height = 10)
# par(mar=c(0.001,0.001,0.001,0.001))
# ptree
# dev.off()