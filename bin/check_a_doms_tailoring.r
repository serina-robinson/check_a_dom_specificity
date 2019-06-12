## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

## Read in the sequences
sp2 <- readAAStringSet('data/sp2_adomains.faa')
names(sp2)
sp2df <- data.frame(names(sp2))
head(sp2df)

# Find all specificities
# names(sp2) <- gsub("\\\t", "_", names(sp2))
sp2_sub <- data.frame(names(sp2), word(names(sp2), sep = "\\\t", 2))
colnames(sp2_sub)
colnames(sp2_sub) <- c("sp2_names", "substrate")

sp2_sub_key <- unique(sp2_sub[,1])
# grep("Cl", sp2_sub[,1])

# Find the methylated sequences
sub_to_keep <- sp2_sub %>%
  dplyr::filter(grepl("Me|OH|Cl", substrate))
nrow(sub_to_keep) # 97 sequences
head(sub_to_keep)

# write_csv(sub_to_keep, "data/tailored_aas.csv")

bgc_all <- word(sub_to_keep$sp2_names, sep = "\\\t", 5)
comb_df <- data.frame(cbind(sub_to_keep, bgc_all, stringsAsFactors = F))

bgc_ids <- bgc_all[!grepl("not_in_mibig", bgc_all)]

# Read in the JSON files for all proteins
allfils <- list.files( "~/Documents/Wageningen_UR/github/mibig_training_set_build_test/data/mibig_json_1.4/", full.names = T)
allfils
bgcs <- allfils[grep(paste0(bgc_ids, collapse = "|"), allfils)]
bgcs

# Find the publications from the json files
json_list <- list()
pub_list <- list()

for(i in 1:length(bgcs)) {
  json_list[[i]] <- fromJSON(bgcs[i], simplifyDataFrame = TRUE)
  names(json_list)[i] <- bgcs[i]
  tmp <- json_list[[i]]$general_params$publications
  tmp_spl <- trimws(unlist(str_split(tmp, pattern = ",")))
  pub_list[[i]] <- tmp_spl
}
names(pub_list) <- gsub(".json", "", word(bgcs, 2, sep = "//"))
to_pull <- unlist(pub_list)
to_pull_un <- unique(to_pull)

to_pull_no_sp_chars <- to_pull_un[-grep(paste0(c("-","unpublished", "10.1", "\\?"), collapse = "|"), to_pull_un)]
to_pull_final <- to_pull_no_sp_chars[to_pull_no_sp_chars!=""]
to_pull_final
writeLines(to_pull_final, sep = ",")
finall <- pub_list %>%
  map(~paste0(., collapse = ";"))
finall

unl <- unlist(finall)


dtf <- data.frame(cbind(unl, names(unl)), stringsAsFactors = F) %>%
  mutate(pmid_last = word(unl, -1, sep = ";"))

colnames(dtf) <- c("pmid_all", "bgc_all", "pmid_last")
nodoi <- dtf$pmid_last[!grepl("anie|ncomms", dtf$pmid_last)]
nodoi <- nodoi[nchar(nodoi) > 1]
write_csv(data.frame(nodoi), "data/only_pmids.csv", col_names = F)

mrged <- left_join(comb_df, dtf)
head(mrged)

# Read in pubmed info
pubmed_info <- read_csv("data/tailoring_rxn_pubmed_info.csv") %>%
  # dplyr::mutate(pmid_last = gsub("\\/pubmed/", "", pubmed_info$URL)) %>%
  dplyr::mutate(pmid_last = as.character(EntrezUID)) %>%
  dplyr::select(Title, Description, Details, Identifiers, pmid_last)

  


final_mrg <- left_join(mrged, pubmed_info, by = "pmid_last")
head(final_mrg)
write_csv(final_mrg, "data/tailored_aas_with_refs.csv")


