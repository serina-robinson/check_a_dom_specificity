## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite', 'ggtree', 'RColorBrewer', 'data.table')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in training set
ts <- readAAStringSet("~/Downloads/sp2_bacterial_adomains_untailored_unsplit.faa")
length(ts)
head(ts)
