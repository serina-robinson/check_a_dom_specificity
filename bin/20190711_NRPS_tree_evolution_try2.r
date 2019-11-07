# Install packages
pacman::p_load("Biostrings", "DECIPHER", "ggtree", "randomcoloR", "stringr", 
               "tidyverse", "RColorBrewer", "gplots", "ape", "dendextend", "phylogram")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the dataset
aa <- readAAStringSet("data/antismashdb_adoms_nostops_filtered_154extracted.faa")
oldnams <- names(aa)
names(aa) <- gsub(paste0(c("\\|", "\\-", "\\[none\\]", "\\(", "\\)",  "\\+", "\\.1"), collapse = "|"), "_", names(aa))
newnams <- names(aa)
# write_tsv(data.frame(cbind(oldnams, newnams)), "output/antismashv2_NRPS_names_mapping.txt")
length(newnams)

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
# 
# clstr50 <- readAAStringSet("data/antismashdb_v2_50")
# length(clstr50)
# 
# clstr40 <- readAAStringSet("data/antismashdb_v2_40")
# length(clstr40)

# Find the clustered hits back in the original sequence alignment file
orig_seq <- AAStringSet(aa[names(aa) %in% names(clstr)])
# writeXStringSet(orig_seq, "data/6581_antismashdb_v2_CDHIT60_alned_seqs.faa")

# Alignment (optional) 
# aa.al <- AlignSeqs(clstr) # Note: takes significant time to run
# writeXStringSet(aa.al, "data/6581_seqs_aligned_DECIPHER.faa")

# Combine with sp2 training set
sp2 <- readAAStringSet("data/sp2_bacterial_adomains_154extracted_untailored_unsplit.faa")
sp2_ids <- word(names(sp2), sep = "_", 1)

# Combine with fungal sequences
sp2_full <- readAAStringSet("data/sp2.adomains.faa")
names(sp2_full) <- gsub("\\\t","_", names(sp2_full))
sp2_full_ids <- word(names(sp2_full), sep = "_", 1)

fung_ids <- setdiff(sp2_full_ids, sp2_ids) # fungal seqs
fungs <- sp2_full[grep(paste0(fung_ids, collapse = "|"), names(sp2_full))] # 64 sequences
ids_tofix <- c("a0098", "a0102", "a0108", "b0077")
subs_toreplace <- c("Orn", "Gln", "Thr|Ser", "Pro")
mustfix <- names(fungs)[grep(paste0(ids_tofix, collapse = "|"), names(fungs))]
newnams <- paste0(word(mustfix, sep = "_", 1), "_", subs_toreplace, "_", word(mustfix, sep = "_", 3), word(mustfix, sep = "\\.1", 2))
names(fungs)[grep(paste0(ids_tofix, collapse = "|"), names(fungs))] <- newnams
#writeXStringSet(fungs, "data/sp2_fungal_only_full_length.faa")
fung_154 <- readAAStringSet("data/sp2_fungal_only_154aa.faa")

# Fix the fungal discrepancies (thanks to Barbara)
sp2_comb <- AAStringSet(c(fung_154, sp2))
table(width(sp2_comb))

comb <- AAStringSet(c(sp2_comb, orig_seq))
# writeXStringSet(comb, "data/7674_sp2_and_antismashdb_v2_60_seqs_154.faa")

# Submit FastTree job
# fasttree -gamma <7674_sp2_and_antismashdb_v2_60_seqs_154.faa> 7674_sp2_and_antismashdb_v2_60_seqs_154.nwk

# Submit RaxML job? 6581 sequences will take forever...

# Read in the tree
phylo_fin <- treeio::read.newick("data/7674_sp2_and_antismashdb_v2_60_seqs_154.nwk")
head(phylo_fin$tip.label)
pl <- ggtree(phylo_fin, layout = "circular", size = 0.25)

# Create annotation data frame
dd <- data.frame(label = pl$data$label,
                 grp = as.character(ifelse(pl$data$isTip, word(pl$data$label, 2, sep = "_"), NA)), stringsAsFactors = F) 
#size = ifelse(as.numeric(pl$data$label) > 0.75, 0.1, 0), stringsAsFactors = F)

# Change the substrate group to be unknown when not experimentally verified
dd$grp[grepl("[[:digit:]]+", dd$grp)] <-  NA
dd$grp[dd$grp == "AMP"] <- NA

# Change so if the frequency of the substrate is less than 10
dd$grp[dd$grp == "di-OH-Bz"] <- "Dhb"
dd$grp[grep("bOH-Tyr", pl$data$label)] <- "Bht"

dd_annot <- dd %>% 
  add_count(grp) %>%
  dplyr::mutate(broad_grp = case_when(n < 8 & !grp %in% c("His", "Pip", "Bht", "Dhb") ~ "rare_aa",
                                      TRUE ~ grp)) %>%
  dplyr::mutate(tipsize = case_when(broad_grp != "UNKNOWN" ~ 0.5,
                                    TRUE ~ 0)) %>%
  dplyr::mutate(substrate = tolower(broad_grp)) %>%
  dplyr::mutate(test_label = gsub("_", "", label))

dd_annot$tipsize[dd_annot$tipsize == 0] <- NA
sort(table(dd_annot$substrate), decreasing = T)

# Read in Barbara's predictions from SANDPUMA1
sandpuma1 <- read_delim("data/sandpuma_output.tsv", delim = "\t", col_names = F)
conflicts <- read_delim("data/proteinogenic_conflict_ASM_SVM.txt", delim = "\t", col_names = F) # conflicts...maybe color black??
no_call_ASM <- read_delim("data/proteinogenic_novel_ASM.txt", delim = "\t", col_names = F) # should be left gray?
rare_maxcount15 <- read_delim("data/proteinogenic_rare_maxcount_15.txt", delim = "\t", col_names = F) # add colorful circles??
tail(rare_maxcount15)

table(tocol$X1 %in% rare_maxcount15$X1)
tocol$X1
# Read in Barbara's PARAS predictions
paras <- read_delim("data/antismashdb_paras_predictions.txt", delim = "\t", col_names = T) %>%
  #dplyr::mutate(Prediction = tools::toTitleCase(Prediction)) %>%
  dplyr::rename(paras_prediction = Prediction) %>%
  dplyr::mutate(ID = gsub("_\\[none\\]", "", ID))


# Filter out predictions we don't want to color
tocol <- sandpuma1 %>%
  tidyr::separate_rows(X1, sep = "\\|\\|\\|") %>%
  dplyr::filter(!X1 %in% conflicts$X1) %>%
  dplyr::filter(!X1 %in% no_call_ASM$X1) %>%
  dplyr::filter(X2 == "SVM") %>%
  dplyr::mutate(ID = gsub("_\\[none\\]", "", X1)) %>%
  #dplyr::mutate(substrate = tools::toTitleCase(X3)) %>%
  dplyr::mutate(SVM_pred = X3) %>%
  dplyr::mutate(fixnam = gsub(paste0(c("\\|", "\\-", "\\[none\\]", "\\(", "\\)",  "\\+", "\\.1"), collapse = "|"), "_", X1)) %>%
  dplyr::mutate(test_label = gsub("_", "", fixnam)) %>%
  dplyr::mutate(rare_maxcount = ifelse(X1 %in% word(rare_maxcount15$X1, sep = "\\|\\|\\|", 1), X1, NA))

word(rare_maxcount15$X1, sep = "\\|\\|\\|", 1)[283]
rare_maxcount15$X1

# Combination
sandpuma_paras_comb <- left_join(tocol, paras, by = "ID")
sandpuma_paras_comb$SVM_pred[sandpuma_paras_comb$SVM_pred == "N/A"] <- NA

# Intersection
table(sandpuma_paras_comb$test_label %in% dd_annot$test_label) # 3373 intersect 

# Combine everything
dd_comb <- dd_annot %>%
  left_join(., sandpuma_paras_comb, by = c("test_label")) %>%
  dplyr::select(-X1, -X2, -X3) %>%
  dplyr::mutate(SVM_pred_comb = ifelse(is.na(substrate), SVM_pred, substrate)) %>%
  dplyr::mutate(SVM_pred_comb = case_when(is.na(SVM_pred_comb) ~ "UNKNOWN",
                                          TRUE ~ SVM_pred_comb)) %>%
  dplyr::mutate(paras_pred_comb = ifelse(is.na(substrate), paras_prediction, substrate)) %>%
  dplyr::mutate(paras_pred_comb = case_when(is.na(paras_pred_comb) ~ "UNKNOWN",
                                            TRUE ~ paras_pred_comb)) %>%
  dplyr::mutate(true_substrate = case_when(is.na(substrate) ~ "UNKNOWN",
                                           TRUE ~ substrate)) %>%
  dplyr::select(-SVM_pred, -paras_prediction, -substrate) 

head(dd_comb)

# Annotate the tree
p2 <- pl %<+% dd_comb

# Set color palette
n1 <- length(unique(p2$data$true_substrate))
n1
n2 <- length(unique(p2$data$paras_pred_comb))
n2
n3 <- length(unique(p2$data$SVM_pred_comb))
n3

setdiff(unique(p2$data$true_substrate), unique(p2$data$paras_pred_comb))
setdiff(unique(p2$data$paras_pred_comb), unique(p2$data$true_substrate))
setdiff(unique(p2$data$SVM_pred), unique(p2$data$true_substrate))

set.seed(1234)
pal2 <- distinctColorPalette(n1)
p2$data$paras_pred_comb
pal2[levels(as.factor(p2$data$paras_pred_comb)) == "ser"] <- "maroon"
pal2[levels(as.factor(p2$data$paras_pred_comb)) == "UNKNOWN"] <- "gray75"
pal2[levels(as.factor(p2$data$paras_pred_comb)) == "rare_aa"] <- "blue3"
pal2[levels(as.factor(p2$data$paras_pred_comb)) == "asn"] <- "orange4"
pal2[levels(as.factor(p2$data$paras_pred_comb)) == "pro"] <- "forestgreen"

# table(p2$data$paras_pred_comb != "UNKNOWN")
# cbind(p2$data$SVM_pred, p2$data$paras_pred_comb)

table(p2$data$paras_pred_comb)
table(p2$data$SVM_pred_comb)


ptree <- p2 +
  aes(color = paras_pred_comb, alpha = 0.5) +
  scale_color_manual(values = pal2, name = "") +
  geom_point(aes(x = 7.0), size = ifelse(p2$data$SVM_pred_comb != "UNKNOWN", 0.5, NA), alpha = 0.5) +
  geom_point(aes(x = 7.25), size = ifelse(p2$data$paras_pred_comb != "UNKNOWN", 0.5, NA), shape = 17) + 
  geom_point(aes(x = 7.5), size = ifelse(is.na(p2$data$rare_maxcount), NA, 1), alpha = 0.5, shape = 15, color = "black") +
  # geom_tippoint(aes(x = 8), size = ifelse(grepl(paste0(word(names(sp2), sep = "_", 1), collapse = "|"), p2$data$label[p2$data$isTip]), 0.5, NA), alpha = 0.5) +
  # geom_tiplab2(aes(label = label, subset = substrate != "UNKNOWN"), size = 0.2) +
  geom_tiplab2(aes(label = label, subset = !is.na(rare_maxcount)), size = 0.2, color = "black") +
  # geom_tiplab2(aes(label = label, subset = grepl(paste0(word(names(sp2), sep = "_", 1), collapse = "|"), p2$data$label[p2$data$isTip])), size = 0.2) +
  theme(legend.position = "right", legend.title = NULL)
# geom_treescale(offset = 1, fontsize = 4, x = 10, y = 10) 


angles <- ptree$data$angle[!is.na(ptree$data$rare_maxcount)]


rounded_angles <- plyr::round_any(angles, 10)
maxcount_names <- ptree$data$label[!is.na(ptree$data$rare_maxcount)]
angle_factors <- factor(rounded_angles, labels = 1:11)
length(rounded_angles)
names_mapping <- data.frame(cbind(oldnams, newnams), stringsAsFactors = F)


maxcount_df <- data.frame(ptree$data$label[!is.na(ptree$data$rare_maxcount)],
                          ptree$data$SVM_pred_comb[!is.na(ptree$data$rare_maxcount)],
                          ptree$data$paras_pred_comb[!is.na(ptree$data$rare_maxcount)],
                          ptree$data$true_substrate[!is.na(ptree$data$rare_maxcount)],
                          ptree$data$angle[!is.na(ptree$data$rare_maxcount)],
                          rounded_angles, angle_factors, stringsAsFactors = F)

colnames(maxcount_df) <- c("label", "svm_pred", "paras_pred", "true_substrate", "exact_angle", "rounded_angle_grp", "grp_no")
mapr <- names_mapping[names_mapping$newnams %in% maxcount_df$label,]

fin_maxcount_df <- maxcount_df %>%
  inner_join(., mapr, by = c("label" = "newnams")) %>%
  dplyr::mutate(old_label = word(oldnams, 1, sep = "_\\[none\\]")) %>%
  dplyr::select(old_label, label, svm_pred, paras_pred, exact_angle, rounded_angle_grp, grp_no)
head(fin_maxcount_df)

# names_mapping$oldnams[maxcount_df$label %in% names_mapping$newnams]

head(maxcount_df)
write_csv(fin_maxcount_df, "output/NRPS_Adomains_for_John_Chu.csv")
table(maxcount_df$paras_pred)


pal3 <- col2hex(pal2)
pal3

dd_all <- data.frame(cbind(ptree$data$label, ptree$data$isTip, ptree$data$rare_maxcount, ptree$data$true_substrate, ptree$data$SVM_pred_comb, ptree$data$paras_pred_comb)) 
colnames(dd_all) <- c("label", "isTip", "rare_maxcount", "true_substrate", "svm_pred_comb", "paras_pred_comb")

dd_fin <- dd_all %>%
  dplyr::filter(isTip == T) %>%
  dplyr::mutate(rare_maxcount_15 = ifelse(is.na(rare_maxcount), "False", "True")) %>%
  dplyr::mutate(color = pal3[as.factor(paras_pred_comb)]) %>%
  dplyr::select(-isTip, -rare_maxcount)

head(dd_fin)

dd_lab <- dd_fin %>%
  dplyr::mutate(bool = ifelse(rare_maxcount_15 == "True", 1, 0)) %>%
  dplyr::filter(bool == 1) %>%
  dplyr::select(label, bool) %>%
  write_csv(., "data/rare_maxcount_labs.txt", col_names = F)
head(dd_lab)


# write_csv(dd_fin, "data/7674_sp2_and_antismashdb_v2_60seqs_substrate_prediction_annotation_file.csv", col_names = T)
# write_tsv(dd_fin, "data/7674_sp2_and_antismashdb_v2_60seqs_substrate_prediction_annotation_file.txt", col_names = T)

writeLines(paste0(pal3, collapse = " "))
rep(1, length(pal3))
writeLines(paste0(levels(dd_fin$paras_pred_comb), collapse = " "))

dd_trim <- dd_fin %>%
  # dplyr::mutate(hashlab = paste0('#', label)) %>%
  dplyr::select(label, color, paras_pred_comb) %>%
  write.table(., "data/7674_sp2_and_antismashdb_v2_60seqs_substrate_prediction_annotation_file_space_delim.txt", row.names = F, sep = " ", quote = F, col.names = F)

pdf("data/20190911_FastTree_sp2_antismashdb_v2_specificity_colored_SANDPUMA_SVM_PARAS_predictions_combined.pdf", width = 10, height = 10)
par(mar=c(0.001,0.001,0.001,0.001))
ptree
dev.off()

pdf(file="output/SANDPUMA_SVM_PARAS_predictions_shape_legend.pdf",bg="white", height = 5, width = 5)
plot.new()
par(mar=c(0.001,0.001,0.001,0.001))
legend("center", pch=c(16, 17, 15), col = "gray75",
       legend=c("SVM predictions", "PARAS predictions", "Good (rare aa) candidates for \n experimental verification in \n collaboration with John Chu"), bty="n")
dev.off()

# dev.off()

## Now cut the tree at different depths
ptree$data$paras_pred_comb
ptree$data$true_substrate
ptree$data$branch.length


# dendro <- ReadDendrogram("data/7674_sp2_and_antismashdb_v2_60_seqs_154.nwk")
dend1 <- ReadDendrogram(textConnection(write.tree(phylo_fin)))
write.dendrogram(dend1, "output/NRPS_dendrogram_written.txt")
hcdend <- as.hclust(dend1)
table(ptree$data$paras_pred_comb)

ptree$data$label <- paste0(ptree$data$paras_pred_comb, "_", ptree$data$label)
ggtree_phylo <- treeio::as.phylo(ptree)
# treeio::write.tree(ggtree_phylo, "output/NRPS_tree_with_substrates_out.nwk")

# Now try using the phylogram package
# dend2 <- as.dendrogram.phylo(ggtree_phylo)
# length(table(word(labels(dend2), sep = "_", 1)))
# length(labels(dend2))

# Dendcut into 3 groups
# dendcut <- cutree(dend2, k = 3)

# dendcut10 <- cutree(dend2, k = 100)


# dendcut_ll <- list()
# for(i in 4:7673) {
#   dendcut_ll[[i]] <- cutree(dend2, k = i)
# }

# Now try converting to 