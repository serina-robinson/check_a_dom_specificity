# Install packages
pacman::p_load("gridExtra", "reshape2", "plyr", "data.table", "ggseqlogo", "cowplot", "Biostrings", "DECIPHER", "ggplot2", "tidyverse")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the original
sp2 <- readAAStringSet("data/sp2_bacterial_adomains_untailored_unsplit.faa")
sort(table(word(names(sp2), sep = "_", 2)), decreasing = T)

# Read in the CD-HIT sets
cdhit <- list.files('data/cd_hit_a_doms/', pattern = "sp2_\\d{2}$", full.names = T) %>%
  rev() %>%
  map(~readAAStringSet(.))

head(cdhit)
ss <-word(list.files('data/cd_hit_a_doms/', pattern = "sp2_\\d{2}$", full.names = T), sep = "_", -1) %>%
             rev()

lens <- unlist(lapply(cdhit, length))
nams <- lapply(cdhit, function(x) { names(x)})
nams
specs <- lapply(nams, function(x) word(x,  sep = "_", 2))
# specs <- lapply(nams, function(x) word(x,  sep = "_", -1))
spectabs <- lapply(specs, function(x) sort(table(x), decreasing = T))
head(spectabs)
names(spectabs) <- ss
df <- ldply(spectabs, data.frame)
head(df)
df_wide <- dcast(df, .id ~ x)
colnames(df_wide)[1] <- "cd_hit_perc_id_cutoff"

df2 <- df_wide[order(df_wide$cd_hit_perc_id_cutoff, decreasing = TRUE),]
head(df2)
nams
df2[is.na(df2)] <- 0
write_csv(df2, "data/sp2_cd_hit_untailored_aas_cutoffs.csv")
# %>% 

# cdhit 
# map_df(~readAAStringSet(.))
