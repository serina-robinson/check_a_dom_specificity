# Install packages
pacman::p_load("gridExtra", "ggseqlogo", "cowplot", "Biostrings", "DECIPHER", "ggplot2", "tidyverse")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Some sample data
data(ggseqlogo_sample)
head(pfms_dna)
head(seqs_aa)

ggseqlogo:::probability_method(seqs_aa[[1]], decreasing = F, seq_type = 'auto', namespace = NULL)
tmp <- data.frame(ggseqlogo:::probability_method(seqs_aa[[1]], decreasing = F, seq_type = 'auto', namespace = NULL)) %>%
  dplyr::filter(position == 1)

mat <- ggseqlogo:::makePFM(seqs_aa[[1]], seq_type = 'auto', namespace = NULL, keep_letter_mat = F)
head(mat)

# Calculate make