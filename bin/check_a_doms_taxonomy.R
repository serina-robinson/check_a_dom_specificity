## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'taxize')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

## Read in the sequences
sp2 <- readAAStringSet('data/sp2_adomains.faa')
names(sp2)

# Find all specificities
# names(sp2) <- gsub("\\\t", "_", names(sp2))
sp2_acc <- word(names(sp2), sep = "\\\t", 3)
sp2_acc_key <- unique(sp2_acc)
names(sp2)
sp2_acc_rem <- sp2_acc_key[!sp2_acc_key %in% c("simA_not_in_genbank",
                                               "chaP_not_in_genbank")]
                                               # "AAK1902.1",
                                               # "CAA11795.1",
                                               # "WP_028678148.1",
                                               # "WP_053065270.1",
                                               # "CFJ82030.1",
                                               # "AAF17280.1")]

# Taxize
uids <- genbank2uid(sp2_acc_rem, api_key="826357f5ff17c7ec62e583909071e94f9d08")
length(uids)
length(sp2_acc_rem) # missing 80 queries
head(uids)
unl_uids <- as.numeric(unlist(lapply(1:length(uids), function(x) { uids[[x]][1] })))
# [36] simA_not_in_genbank, [45] AAK1902.1, [288] CAA11795.1, [36] WP_028678148.1, [45] WP_053065270.1


# Make uid/acc key
uid_key <- cbind(sp2_acc_key, unl_uids)

# Classify
taxs <- classification(unl_uids, db = "ncbi",
                 api_key="826357f5ff17c7ec62e583909071e94f9d08")
head(taxs)
tax_trim <- taxs[unlist(lapply(taxs, length)) == 3]
tax_trim
superkingdom <- unlist(lapply(1:length(tax_trim), function(x) { tax_trim[[x]][2,1] }))
phylum <-  unlist(lapply(1:length(tax_trim), function(x) { tax_trim[[x]][4, 1] }))
taxdf <- data.frame(cbind(names(tax_trim), superkingdom, phylum), stringsAsFactors = F)
colnames(taxdf)[1] <- "uid"


# Create a dataframe
sp2_acc_df <- data.frame(sp2_acc, stringsAsFactors = F)
colnames(sp2_acc_df) <- "acc"

# Merge back into the big dataset
mrged <- dplyr::left_join(sp2_acc_df, taxdf, by = "acc")
# write_csv(mrged, "data/sp2_adomains_with_taxonomy.csv")

# Read in  the fungal A domains
taxdat <- read_csv("data/sp2_adomains_with_taxonomy.csv")
head(taxdat)

# Read
table(taxdat$superkingdom)

# Separate 