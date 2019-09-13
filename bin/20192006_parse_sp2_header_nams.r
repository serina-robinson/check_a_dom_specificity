## Install packages
pacman::p_load('Biostrings', 'DECIPHER', 'tidyverse', 'jsonlite', 'ggtree', 'RColorBrewer', 'data.table')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/check_a_dom_specificity/")

# Read in the predictions
preds <- read_csv("data/Amplicon_specificity_predictions_20192006.csv")
myun <- unique(preds$Prediction_monomer)
finun <- myun[myun != "abu|ival"]
length(finun)
length(myun)

# Remove abu|ival
abu|ival

# Read in the sp2
sp2 <- readAAStringSet("~/Downloads/sp2.adomains_George.faa")


# Read in Barbara's fixed
barb <- readAAStringSet("data/sp2_bacterial_adomains_untailored_unsplit.faa")
names(barb) <- gsub("_", "\t", names(barb))
accs <- word(names(barb), 1, sep = "\t")
dtf <- data.frame(cbind(accs, names(barb)), stringsAsFactors = F)
accs <- word(names(sp2), sep = "\t", 1)
dtf2 <- data.frame(cbind(accs , names(sp2)), stringsAsFactors = F)

dtf3 <- dtf2 %>%
  left_join(dtf, by = "accs")
colnames(dtf3)[2:3] <- c("tailored", "untailored")
1093 - 1029
table(is.na(dtf3$untailored))

dtf3$untailored[is.na(dtf3$untailored)] <- as.character(dtf3$tailored[is.na(dtf3$untailored)])
head(dtf3)

# orig <- readAAStringSet("data/sp2_1043_adomains_full_length_no_fungs.faa")
# length(barb)
# length(orig)

# Tailoring 
specs <- word(dtf3$untailored, sep = "\t", 2)
grep("OH-Orn",specs)

specs <- tolower(specs)


# Filter out
specnew <- specs %>%
  tibble::enframe(name = NULL) %>%
  mutate(new = case_when(!specs %in% finun ~ "other",
         TRUE ~ specs))
specnew
names(sp2)
newnams <- paste0(word(names(sp2), sep = "\t", 1), "\t",
                       specnew$new, "\t",
                       word(names(sp2), sep = "\t", 3), "\t",   
                            word(names(sp2), sep = "\t", 4), "\t",
                                 word(names(sp2), sep = "\t", 5), "\t",
                  word(names(sp2), sep = "\t", 6))

newnams
specs[grep("oh-orn", specs)]

names(sp2) <- newnams
table(word(names(sp2), sep = "\t", 2))
writeXStringSet(sp2, "data/sp2_adomains_names_fixed.faa")
