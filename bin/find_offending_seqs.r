Install packages
pacman::p_load("DECIPHER", "plot3D", "plot3Drgl", "ade4", "data.table",
               "phangorn", "cowplot", "RColorBrewer", "phylobase", "treeio", "ggtree", "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# pacman::p_load("tidyverse", "Biostrings")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/amplicon_pred/")

# Read in the sp2 adomains
sp2 <- readAAStringSet("data/sp2_adomains_tabbed.faa")
sp_rem <- sp2[-grep("X", sp2)]
length(sp_rem)
writeXStringSet(sp_rem, "~/Documents/Wageningen_UR/github/amplicon_pred/data/sp2_adomains_tabbed.faa")


# Read in the 34 trimmed amplicons
sptrim <- readAAStringSet("data/sp2_34.faa")
sp9 <- AAStringSet(substr(sptrim, 1, 9))
writeXStringSet(sp9, "~/Documents/Wageningen_UR/github/amplicon_pred/data/sp9_trimmed.faa")
