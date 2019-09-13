# Install packages
pacman::p_load("gridExtra", "ggseqlogo", "cowplot", "Biostrings", "DECIPHER", "ggplot2", "tidyverse")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Create data frame of signatures
trainset <- readAAStringSet("data/sp2_bacterial_adomains_untailored_unsplit.faa")

# Extract the 34 sequences
trset34 <- readAAStringSet("data/sp2_bacterial_adomains_untailored_unsplit_34_aa_xtracted.faa")

# Replace \\[none]
names(trset34) <- gsub("_\\[none\\]", "", names(trset34))
head(names(trset34))
table(word(names(trset34), 2, sep = "_"))

# 
# comb_dedup <- trset34[!duplicated(trset34)]
comb_dedup <- trset34
# Create a data frame
dat <- data.frame(comb_dedup) %>%
  rownames_to_column() %>%
  mutate(spec = word(rowname, sep = "_", 2)) %>%
  add_count(spec) %>%
  arrange(desc(n))
  # dplyr::filter(spec == 1)
dim(dat)
# write_csv(dat, "data/sp2_untailored_sigs_for_HMMlogos.csv")

colnames(dat) <- c("domain_names", "aa34", "spec", "n")



# Function to pull aa spec
PullCode <- function(aa, dat) {
  aa2 <- dat[dat$spec == aa,]
  aa3 <- as.character(aa2$aa34)
  # aa3 <- gsub("-", "", as.character(aa2$stachelhaus.code))
  p <- ggplot() + geom_logo(aa3, method = "p", col_scheme = 'chemistry') + theme_logo() + 
    theme(legend.position = 'none', axis.text.x = element_blank()) +
    ggtitle(paste0(aa, " (", dat$n[dat$spec == aa][1], ")" ))
  return(p)
}

PullMat <- function(aa, dat) {
  aa2 <- dat[dat$spec == aa,]
  aa3 <- as.character(aa2$aa34)
  mat <- ggseqlogo:::makePFM(aa3, seq_type = 'auto', namespace = NULL, keep_letter_mat = F)
  return(mat)
}

# Amino acids

### Negatively charged
# Asp aspartic acid (asp, solubility -55)
# Glu glutamic acid (glu, )
### Positively charged
# Arginine (lys, solubility -23)
# Asparagine (asn)
### Neutral
# Serine (ser, solubility, -10)
# Glycine (gly solubility, 0)

# leucine
# aas <- c("asp", "glu", "arg", "lys", "ser", "gly")
aas <- unique(dat$spec)
pl <- list()
for(i in 1:length(aas)){
  pl[[i]] <- PullCode(aas[i], dat)
} 

mat <- list()
for(i in 1:length(aas)) {
  mat[[i]] <- PullMat(aas[i], dat)
}
names(mat) <- aas



matbind <- bind_rows(mat, .id = NULL)
dim(matbind) # 680 by 75
head(matbind)



dtf <- data.frame(matbind)

head(dtf)
dim(dtf)
rownames(dtf)
dim(dtf)
dtf[,1]
aas_for_df <- rep(attr(mat[1]$Leu, "namespace"), 34)
pos_for_df <- rep(1:34, each = 20)

findf <- data.frame(cbind(aas_for_df, pos_for_df, dtf))

write_csv(findf, "data/probability_matrices_for_each_residue.csv")

# Output
# pdf("data/HMMlogos_34_extracted_sp2_trainingset.pdf", width = 150, height = 130)
# do.call("grid.arrange", c(pl, ncol = 1))
# dev.off()

