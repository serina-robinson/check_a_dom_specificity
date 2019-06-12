## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in gapped seqs
adoms <- readAAStringSet("data/adoms_with_gaps.faa")
adom_nms <- names(adoms)
adom_seqs <- as.character(adoms)
dtf <- data.frame(cbind(names(adoms), adom_seqs))
colnames(dtf) <- c("A-domain names", "34 active site residue")
write_csv(dtf, "data/adoms_gapped_with_taxonomy.csv")

