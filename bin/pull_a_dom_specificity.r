## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'taxize')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

## Read in the sequences
sp2 <- readAAStringSet('data/sp2_adomains.faa')
names(sp2)

# Pull from BGC the source organism
bgc2 <- read_delim('data/20171028_MIBIG_scrape.txt', delim = "\t") %>%
  janitor::clean_names() %>%
  mutate(Genus = word(organism, 1, sep = " "))
colnames(bgc2)[1] <- "bgc"

sp2df <- data.frame(names(sp2))
colnames(sp2df) <- "sp2_id"
findf <- sp2df %>%
  mutate(bgc = word(sp2_id, sep = "\\\t", -2)) %>%
  left_join(., bgc2, by = "bgc") %>%
  dplyr::select(sp2_id, bgc, main_product, biosynthetic_class, Genus) %>%
  janitor::clean_names()

genun <- unique(findf$genus)
length(genun)

# Pull the taxid for each genus
# taxid <- classification(genun, db = 'ncbi')
tax_trim <- taxid[unlist(lapply(taxid, length)) == 3]
head(tax_trim)
superkingdom <- unlist(lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^superkingdom$", tax_trim[[x]]$rank),1] }))
class <-  lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^class$", tax_trim[[x]]$rank),1] })
class_l <- unlist(lapply(1:length(class), function(x) ifelse(length(class[[x]]) == 0, "Unknown", class[[x]])))
order2 <-  lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^order$", tax_trim[[x]]$rank),1] })
order_l <- unlist(lapply(1:length(order2), function(x) ifelse(length(order2[[x]]) == 0, "Unknown", order2[[x]])))
phylum <-  lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^phylum$", tax_trim[[x]]$rank),1] })
phylum_l <-  unlist(lapply(1:length(phylum), function(x) ifelse(length(phylum[[x]]) == 0, "Unknown", phylum[[x]])))
genus <- lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^genus$", tax_trim[[x]]$rank),1] })
genus_l <- unlist(lapply(1:length(genus), function(x) ifelse(length(genus[[x]])== 0, "Unknown", genus[[x]])))
species <- lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^species$", tax_trim[[x]]$rank),1] })
species_l <- unlist(lapply(1:length(species), function(x) ifelse(length(species[[x]]) == 0, "Unknown", species[[x]])))
family <- lapply(1:length(tax_trim), function(x) { tax_trim[[x]][grep("^family$", tax_trim[[x]]$rank),1] })
family_l <- unlist(lapply(1:length(family), function(x) ifelse(length(family[[x]]) == 0, "Unknown", family[[x]])))

taxdf <- data.frame(superkingdom, phylum_l, class_l, order_l, family_l, genus_l, species_l, stringsAsFactors = F)
colnames(taxdf) <- gsub("_l", "", colnames(taxdf))
head(taxdf)

# Merge with the other table
findf_gen <- dplyr::inner_join(findf, taxdf, by = "genus")
head(findf_gen)
write_csv(findf_gen, "data/taxonomy_of_sp2_doms.csv")
