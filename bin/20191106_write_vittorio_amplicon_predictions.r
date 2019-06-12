## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite', 'ggtree', 'RColorBrewer', 'data.table')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the data 
aa2 <- readAAStringSet("~/Documents/Wageningen_UR/github/adenylpred_dev/data/vittorio_amplicons_aligned_NRPS_Adomains_Xsremoved.fasta")
names(aa2)
aa2[1]

# Read in the HMM Aligned
aa_al <- readAAStringSet("~/Documents/Wageningen_UR/github/adenylpred_dev/data/vittorio_mod_hmm_aligned.fasta")
aa_df <- data.frame(cbind(names(aa_al), as.character(aa_al)))
head(aa_df)

# Read in the predictions
monomers <- fread("~/Documents/Wageningen_UR/github/adenylpred_dev/data/20191106_vittorio_monomer_predictions.txt", data.table = F)

# Read in the group predictions
grps <- fread("~/Documents/Wageningen_UR/github/adenylpred_dev/data/vittorio_amplicons_groups_predictions.txt", data.table = F)

aa_df
# Combine
datcomb <- monomers %>%
  inner_join(grps, by = "Query_name")

datcomb$Query_name <- gsub("_\\[none\\]", "", datcomb$Query_name)
colnames(datcomb)[2:3] <- gsub("\\.x", "_monomer", colnames(datcomb)[2:3])
colnames(datcomb)[4:5] <- gsub("\\.y", "_groups", colnames(datcomb)[4:5])
colnames(datcomb)

datcomb_ord <- datcomb[order(datcomb$Probability_score_monomer, decreasing = T),]
head(datcomb_ord)

write_csv(data.frame(table(datcomb_ord$Prediction_monomer)), "data/monomer_prediction_table.csv")
write_csv(data.frame(table(datcomb_ord$Prediction_groups)), "data/groups_prediction_table.csv")


summary(datcomb_ord$Probability_score_monomer)
table(datcomb_ord$Probability_score_monomer < 0.5)
table(datcomb_ord$Probability_score_groups < 0.5)

write_csv(datcomb_ord, "data/Amplicon_specificity_predictions_20191106.csv")
