## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite', 'ggtree', 'RColorBrewer', 'data.table')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the data 
aa2 <- readAAStringSet("~/Documents/Wageningen_UR/github/adenylpred_dev/data/vittorio_amplicons_aligned_NRPS_Adomains_Xsremoved.fasta")

# Read in the HMM Aligned
aa_al <- readAAStringSet("~/Documents/Wageningen_UR/github/adenylpred_dev/data/vittorio_mod_hmm_aligned.fasta")
names(aa_al) <- gsub(" ", "", names(aa_al))
aa_df <- data.frame(cbind(names(aa_al), as.character(aa_al)))
colnames(aa_df) <- c("Query_name", "Seq")
aa_df$Query_name <- as.character(aa_df$Query_name)

# Read in the predictions
monomers <- fread("~/Documents/Wageningen_UR/github/adenylpred_dev/data/20191106_vittorio_monomer_predictions.txt", data.table = F)

# Read in the group predictions
grps <- fread("~/Documents/Wageningen_UR/github/adenylpred_dev/data/vittorio_amplicons_groups_predictions.txt", data.table = F)

# Combine
datcomb <- monomers %>%
  inner_join(grps, by = "Query_name") 
head(datcomb)

datcomb$Query_name <- gsub("_\\[none\\]", "", datcomb$Query_name)
colnames(datcomb)[2:3] <- gsub("\\.x", "_monomer", colnames(datcomb)[2:3])
colnames(datcomb)[4:5] <- gsub("\\.y", "_groups", colnames(datcomb)[4:5])
datcomb$Query_name

aa_df$Query_name
datcomb2 <- aa_df %>%
  inner_join(datcomb, by = "Query_name")

head(datcomb2)
datcomb_ord <- datcomb2[order(datcomb2$Probability_score_monomer, decreasing = T),] 
  
write_csv(data.frame(table(datcomb_ord$Prediction_monomer)), "data/monomer_prediction_table.csv")
write_csv(data.frame(table(datcomb_ord$Prediction_groups)), "data/groups_prediction_table.csv")

summary(datcomb_ord$Probability_score_monomer)
table(datcomb_ord$Probability_score_monomer < 0.5)
table(datcomb_ord$Probability_score_groups < 0.5)

write_csv(datcomb_ord, "data/Amplicon_specificity_predictions_with_seqs_20191106.csv")
