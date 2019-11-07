# # Install packages
pacman::p_load("Biostrings", "DECIPHER", "ggtree", "randomcoloR", "stringr", 
               "tidyverse", "RColorBrewer", "gplots", "ape", "dendextend", "phylogram")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the sp2 domains
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

# Make phylogenetic tree
sp2_comb <- AAStringSet(c(fung_154, sp2))
writeXStringSet(sp2_comb, "data/sp2_bacterial_and_fungal_combined_154aa_extracted.faa")
head(sp2_comb)

# Read in the tree
tr <- treeio::read.tree("data/sp2_bacterial_fungal_154aa_fasttree.nwk")



# Convert to dendrogram
dend2 <- as.dendrogram.phylo(tr)

# Cut at different heights
dendcut <- cutree(dend2, k = 10)
dend5 <- cutree(dend2, k = 5)
dend15 <- cutree(dend2, k = 15)
dend20 <- cutree(dend2, k = 20)
dend25 <- cutree(dend2, k = 25)
dend3 <- cutree(dend2, k = 3)

# Combine all 
dend_comb <- data.frame(c(dend5, dendcut, dend15, dend20, dend25), stringsAsFactors = F)




dendh3 <- cutree(dend2, h = 3)
