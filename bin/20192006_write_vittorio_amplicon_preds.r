## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite', 'ggtree', 'RColorBrewer', 'data.table')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the data 
aa2 <- readAAStringSet("~/Documents/Wageningen_UR/github/adenylpred_scratch/data/51914_vittorio_amplicons_range_extracted.faa")
names(aa2)
length(aa2)

# Read in the HMM Aligned
aa_al <- readAAStringSet("~/Documents/Wageningen_UR/github/adenylpred_dev/data/vittorio_mod_hmm_aligned.fasta")
aa_df <- data.frame(cbind(names(aa_al), as.character(aa_al)))
head(aa_df)

# Read in the predictions
monomers <- fread("~/Documents/Wageningen_UR/github/adenylpred_scratch/data/51914_vittorio_amplicons_monomers_predictions.txt", data.table = F)

# Read in the group predictions
grps <- fread("~/Documents/Wageningen_UR/github/adenylpred_scratch/data/51914_vittorio_amplicons_groups_predictions.txt", data.table = F)

aa_df
# Combine
datcomb <- monomers %>%
  inner_join(grps, by = "Query_name")

datcomb$Query_name <- gsub("_\\[none\\]", "", datcomb$Query_name)
colnames(datcomb)[2:3] <- gsub("\\.x", "_monomer", colnames(datcomb)[2:3])
colnames(datcomb)[4:5] <- gsub("\\.y", "_groups", colnames(datcomb)[4:5])
colnames(datcomb)
head(datcomb)

datcomb_ord <- datcomb[order(datcomb$Probability_score_monomer, decreasing = T),]

head(datcomb_ord)

png("data/probability_distribution.png", width = 4, height = 3, units = 'in', res = 600)
ggplot(data = datcomb_ord) +
  geom_density(aes(x = Probability_score_monomer), alpha = 0.3, fill = "pink") +#, stat = "identity") +
  scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2.2), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_density(aes(x = Probability_score_groups), alpha = 0.3, fill = "dodgerblue") + #stat = "identity") +
  xlab("Prediction probability score") +
  geom_vline(aes(xintercept = 0.5),
             color="gray", linetype="dashed", size=1)
dev.off()

ggplot(data = datcomb_ord, aes(x = Probability_score_groups)) +
  geom_density(alpha = 0.2, fill = "dodgerblue") +
  scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2.2), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
 

monomer_tab <- ifelse(as.numeric(datcomb_ord$Probability_score_monomer) > 0.5, datcomb_ord$Prediction_monomer, "No confident prediction")
dtb <- data.frame(table(monomer_tab))
res <- dtb[order(dtb$Freq, decreasing = T),]

grp_tab <- ifelse(as.numeric(datcomb_ord$Probability_score_groups) > 0.5, datcomb_ord$Prediction_groups, "No confident prediction")
table(grp_tab)
dtb <- data.frame(table(grp_tab))
dtb
res2 <- dtb[order(dtb$Freq, decreasing = T),]

write_csv(res, "data/20192006_monomer_prediction_table.csv")
write_csv(res2, "data/20192006_groups_prediction_table.csv")
# 0.5 cutoff
datcomb_tr <- datcomb_ord[datcomb_ord$Probability_score_groups >= 0.5,]
dim(datcomb_tr)[1]/(dim(datcomb)[1])

write_csv(datcomb_tr, "data/Amplicon_predictions_prob_groups_greater_0.6_20192006.csv")
datcomb_tr_mon <- datcomb_ord[datcomb_ord$Probability_score_monomer >= 0.5,]
dim(datcomb_tr_mon)[1]/(dim(datcomb)[1])

datcomb_poor <- datcomb_ord[datcomb_ord$Probability_score_monomer <= 0.5,]
dim(datcomb_poor)
datcomb_poor$Query_name %in% 

write_csv(datcomb_tr_mon, "data/Amplicon_predictions_prob_monomers_greater_0.5_20192006.csv")
summary(datcomb_tr_mon$Probability_score_groups)






# 0.6 cutoff 
datcomb_tr <- datcomb_ord[datcomb_ord$Probability_score_groups > 0.6,]
write_csv(datcomb_tr, "data/0.6_Amplicon_predictions_prob_groups_greater_20192006.csv")
datcomb_tr_mon <- datcomb_ord[datcomb_ord$Probability_score_monomer > 0.6,]
write_csv(datcomb_tr_mon, "data/0.6 Amplicon_predictions_prob_monomers_greater_0.6_20192006.csv")
summary(datcomb_tr_mon$Probability_score_groups)

# 0.7 cutoff 
datcomb_tr <- datcomb_ord[datcomb_ord$Probability_score_groups > 0.7,]
write_csv(datcomb_tr, "data/0.7_Amplicon_predictions_prob_groups_greater_20192006.csv")
dim(datcomb_tr)
datcomb_tr_mon <- datcomb_ord[datcomb_ord$Probability_score_monomer > 0.7,]
write_csv(datcomb_tr_mon, "data/0.7_Amplicon_predictions_prob_monomers_greater_20192006.csv")
dim(datcomb_tr_mon)
summary(datcomb_tr_mon$Probability_score_groups)

# 0.8 cutoff 
datcomb_tr <- datcomb_ord[datcomb_ord$Probability_score_groups > 0.8,]
write_csv(datcomb_tr, "data/0.8_Amplicon_predictions_prob_groups_greater_0.8_20192006.csv")
dim(datcomb_tr)
datcomb_tr_mon <- datcomb_ord[datcomb_ord$Probability_score_monomer > 0.8,]
write_csv(datcomb_tr_mon, "data/0.8_Amplicon_predictions_prob_monomers_greater_0.8_20192006.csv")
dim(datcomb_tr_mon)
summary(datcomb_tr_mon$Probability_score_groups)


colnames(aa_df) <- c("Query_name", "Seq")
aa_df$Query_name <- trimws(as.character(aa_df$Query_name))
head(aa_df$Query_name)
datcomb$Query_name <- trimws(datcomb$Query_name)
datcomb$Query_name[1:10]
aa_df$Query_name[1:10]
fin_merg <- aa_df %>%
  inner_join(datcomb, by = "Query_name")
head(fin_merg)


#write_csv(fin_merg, "data/Amplicon_specificity_predictions_20192006.csv")


# Read in old predictions
old_pred <- read_csv("data/Amplicon_specificity_predictions_20191106.csv")


# Read in the data 
aa2 <- readAAStringSet("~/Documents/Wageningen_UR/github/adenylpred_dev/data/vittorio_amplicons_aligned_NRPS_Adomains_Xsremoved.fasta")
names(aa2)
length(aa2)

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
colnames(datcomb)[2:3] <- gsub("\\.x", "_monomer_old", colnames(datcomb)[2:3])
colnames(datcomb)[4:5] <- gsub("\\.y", "_groups_old", colnames(datcomb)[4:5])
colnames(datcomb)


alljoin <- fin_merg %>%
  inner_join(datcomb, by = "Query_name")
dim(alljoin)  
head(alljoin)

sp2 <- readAAStringSet("data/51914_vittorio_amplicons_range_extracted.faa")
which_poor <- sp2[gsub("_\\[none\\]", "", names(sp2)) %in% datcomb_poor$Query_name]
head(which_poor)
length(which_poor)
writeXStringSet(which_poor, "data/26510_amplicons_with_low_confidence_predictions.faa")
length(grep("--", which_poor))

3114/26510
# write_csv(alljoin, "data/Amplicon_predictions_old_new_comparison.csv")

     