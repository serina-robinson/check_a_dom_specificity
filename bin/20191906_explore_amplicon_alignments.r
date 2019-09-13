## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite', 'ggtree', 'RColorBrewer', 'data.table')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the data 
aa2 <- readAAStringSet("~/Documents/Wageningen_UR/github/amplicon_pred/data/aa-sequences.fasta")
aa_sm <- aa2[1:100]
aa_sm

# Read in the aligned a-domain seqs
# aa2_al <- readAAStringSet("~/Documents/Wageningen_UR/github/amplicon_pred/data/aa-sequences_37_extracted.fasta")
aa2_al <- readAAStringSet("~/Documents/Wageningen_UR/github/amplicon_pred/data/vittorio_sequences_34range_extracted.faa")
head(aa2_al)
head(aa2_al)
# Set random seed
set.seed(1234)

# Remove all sequences with gaps
aa2_nogaps <- aa2_al[-grep("-----", aa2_al)]
length(aa2_nogaps) # 49,580

noX <- gsub("X", "-", aa2_nogaps)
grep("X", noX)
length(noX)
tmp <- AAStringSet(noX)
names(tmp)
writeXStringSet(AAStringSet(noX), "data/49580_vittorio_amplicons_range_extracted.faa")

writeXStringSet(AAStringSet(noX), "~/Documents/Wageningen_UR/github/adenylpred_scratch/data/49580_vittorio_amplicons_range_extracted.faa")

writeXStringSet(aa2_nogaps[1:10], "~/Documents/Wageningen_UR/github/adenylpred_scratch/data/10_test_vittorio_amplicons_range_extracted.faa")

noX <- gsub("X", "-", aa2_al)
tmp <- AAStringSet(noX)
length(names(tmp))
writeXStringSet(tmp, "data/51914_vittorio_amplicons_range_extracted.faa")

writeXStringSet(tmp, "~/Documents/Wageningen_UR/github/adenylpred_scratch/data/51914_vittorio_amplicons_range_extracted.faa")


full <- aa2[-grep("-", aa2_al)]

length(aa2_al)
rand <- full[sample(size = 100, x = length(full), replace = F)]
length(rand)

# Align with the PheA sequence
phea <- readAAStringSet("~/Documents/Wageningen_UR/github/amplicon_pred/data/1AMU:A|PDBID|CHAIN|SEQUENCE")

# Try alignment 
comb <- c(phea, rand)
aa.al <- AlignSeqs(comb)
BrowseSeqs(aa.al)
head(aa2_al[1:100])
