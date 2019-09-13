# Install packages
pacman::p_load("gridExtra", "ggseqlogo", "cowplot", "Biostrings", "DECIPHER", "ggplot2", "tidyverse")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")


# Create data frame of signatures
trainset <- readAAStringSet("data/sp2_bacterial_adomains_untailored_unsplit.faa")


traindf <- data.frame(trainset) %>%
  rownames_to_column() %>%
  mutate(id = word(rowname, sep = "_", 1)) %>%
  mutate(spec = word(rowname, sep = "_", 2))
head(traindf)


colnames(traindf)[1:2] <- c("domain_names", "fullseq")

# Extract the 34 sequences
trset34 <- readAAStringSet("data/sp2_bacterial_adomains_untailored_unsplit_34_aa_xtracted.faa")
trset_range <- readAAStringSet("data/sp2_bacterial_untailored_full_range_extracted.faa")
width(trset_range)


# Replace \\[none]
names(trset_range) <- gsub("_\\[none\\]", "", names(trset_range))

table(word(names(trset_range), 2, sep = "_"))


cdat <- data.frame(trset_range) %>%
  rownames_to_column() %>%
  mutate(spec = word(rowname, sep = "_", 2)) %>%
  add_count(spec) %>%
  arrange(desc(n)) %>%
  dplyr::select(spec, n) %>%
  distinct()


# write_csv(cdat, "data/monomer_counts.csv")

comb_dedup <- trset_range[!duplicated(trset_range)]

# Create a data frame
dat <- data.frame(comb_dedup) %>%
  rownames_to_column() %>%
  mutate(spec = word(rowname, sep = "_", 2)) %>%
  add_count(spec) %>%
  arrange(desc(n))
# dplyr::filter(spec == 1)
dim(dat)
head(dat)


cdat# write_csv(dat, "data/sp2_untailored_sigs_for_HMMlogos.csv")

colnames(dat) <- c("domain_names", "aa_range", "spec", "n")

# Function to pull aa spec
PullCode <- function(aa, dat) {
  aa2 <- dat[dat$spec == aa,]
  aa3 <- as.character(aa2$aa_range)
  # aa3 <- gsub("-", "", as.character(aa2$stachelhaus.code))
  p <- ggplot() + geom_logo(aa3, method = "p", col_scheme = 'chemistry') + theme_logo() + 
    theme(plot.title = element_text(size=100),
          legend.position = 'none', axis.text.x = element_blank()) +
    ggtitle(paste0(aa, " (", dat$n[dat$spec == aa][1], ")" ))
  return(p)
}

PullMat <- function(aa, dat) {
  aa2 <- dat[dat$spec == aa,]
  aa3 <- as.character(aa2$aa_range)
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
attr(mat[30]$`Ile|Val`, "namespace")

aas_for_df <- rep(attr(mat[1]$Leu, "namespace"), (dim(dtf)[1]/(20)))
pos_for_df <- rep(1:(dim(dtf)[1]/(20)), each = 20)

findf <- data.frame(cbind(aas_for_df, pos_for_df, dtf))

write_csv(findf, "data/probability_matrices_for_each_residue_fullrange.csv")


# Output
# pdf("data/HMMlogos_full_range_trainingset.pdf", width = 150, height = 350)
# do.call("grid.arrange", c(pl, ncol = 1))
# dev.off()


# Combine the dat with the aa34
# aa34 <- read_csv("data/sp2_untailored_sigs_for_HMMlogos.csv")
# dim(aa34)
# comb2 <- left_join(dat, aa34, by = c("spec", "domain_names")) %>%
#   dplyr::rename(n_range = n.x,
#                 n_34 = n.y) %>%
#   mutate(id = word(domain_names, sep = "_", 1)) %>%
#   dplyr::select(id, spec, domain_names, aa_range, aa34)
# head(comb2)
# 
# comb3 <- right_join(traindf, comb2, by = c("id", "spec")) %>%
#   dplyr::select(-domain_names.x)
# head(comb3)
# # write_csv(comb3, "data/Trainingset_for_HMMlogos.csv")
