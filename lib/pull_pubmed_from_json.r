# Parse the names
prot_names <- word(names(no_nrps), 6, sep = "\\|")
bgc_names <- word(names(no_nrps), 1, sep = "\\|")
bgc_acc <- word(names(no_nrps), 7, sep = "\\|")
names(no_nrps)
# sort(table(prot_names), decreasing = T)[1:20]
prot_df <- data.frame(cbind(as.character(no_nrps), bgc_acc, prot_names, bgc_names))
colnames(prot_df) <- c("seqs", "acc", "prot_names", "bgcs")
head(prot_df)
# write.csv(cbind(as.character(no_nrps), prot_names), "output/mibig_prot_names.csv")

bgc_ids <- word(names(no_nrps), 1, sep = "\\|")
length(bgc_ids) # 355

# Read in the JSON files for all proteins
allfils <- list.files("~/Documents/Wageningen_UR/github/mibig_training_set_build/data/",full.names = T)
bgcs <- allfils[grep(paste0(bgc_ids, collapse = "|"), allfils)]
bgcs #only 297 had .json files

json_list <- list()
pub_list <- list()
for(i in 1:length(bgcs)) {
  json_list[[i]] <- fromJSON(bgcs[i], simplifyDataFrame = TRUE)
  names(json_list)[i] <- bgcs[i]
  tmp <- json_list[[i]]$general_params$publications
  tmp_spl <- trimws(unlist(str_split(tmp, pattern = ",")))
  pub_list[[i]] <- tmp_spl
}
# str(json_list[[2]])

names(pub_list) <- gsub(".json", "", word(bgcs, 2, sep = "//"))
names(json_list) <- names(pub_list)

to_pull <- unlist(pub_list)
to_pull_un <- unique(to_pull)

to_pull_no_sp_chars <- to_pull_un[-grep(paste0(c("-","unpublished", "10\\.", "\\?"), collapse = "|"), to_pull_un)]
to_pull_final <- to_pull_no_sp_chars[to_pull_no_sp_chars!=""]

finall <- pub_list %>%
  map(~paste0(., collapse = ";"))
finall

unl <- unlist(finall)

# write.table(to_pull_final, quote = F, row.names = F, col.names = F, "output/mibig_pmids_to_pull.txt")

dtf <- data.frame(cbind(unl, names(unl)), stringsAsFactors = F) %>%
  mutate(pmid_last = word(unl, -1, sep = ";"))