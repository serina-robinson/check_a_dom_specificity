# Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the 34aa HMM aligned sequences
aa34 <- readAAStringSet("~/Documents/Wageningen_UR/github/check_a_dom_specificity/data/sp2_adoms_34extract_hmmalign.faa")
name_key <- readAAStringSet("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/data/sp2_34extract_names_fixed_large_grps.faa")
names(aa34) <- names(name_key)
table(word(names(aa34), sep = "_", -1))

# Remove all duplicates
dedup <- aa34[!duplicated(aa34)]
length(dedup) # 843 
subst <- word(names(dedup), sep = "_", 2)
subst <- gsub("\\|", "\\.", subst)
subst <- gsub("\\-", "\\.", subst)
names(dedup) <- gsub("\\|", "\\.", names(dedup))
names(dedup) <- gsub("\\-", "\\.", names(dedup))
table(subst)
length(subst)
writeXStringSet(dedup, "data/843_34extract_hmmalign.faa")

dtf <- data.frame(names(dedup), dedup, subst, stringsAsFactors = F) %>%
  dplyr::filter(!grepl("reject", names.dedup.)) %>%
  dplyr::filter(!grepl("-", dedup)) %>%
  group_by(subst) %>%
  add_count(subst) %>%
  dplyr::filter(n > 8) 

dim(dtf) # 652
table(dtf$subst)

# Remove all sequences with specificity less than 8
dedup2 <- AAStringSet(dtf$dedup)
names(dedup2) <- dtf$names.dedup.
writeXStringSet(dedup2, "data/652_monomers_deduplicated.faa")

# First do large substrate groups
rdaln <- read.alignment("data/652_monomers_deduplicated.faa", format = "fasta")
str(rdaln)

# Convert to a vector of physicochemical properties
source("src/convert_aln_15aap.r")
aa <- convert_aln_15aap(rdaln) #15 physicochemical properties
aadf <- data.frame(aa,stringsAsFactors = F)
dim(aadf) # 843 x 510
rownames(aadf)
# write.csv(aadf, paste0("data/652_aa_seqs_510_feats_for_supervised_monomers_20190406.csv"), row.names=rownames(aadf), quote = F)

# Read in the data
rawdat <- read_csv("data/652_aa_seqs_510_feats_for_supervised_monomers_20190406.csv", skip_empty_rows = T)

rawdat <- data.frame(cbind(rawdat$X1), scale(rawdat[,2:ncol(rawdat)]), stringsAsFactors = F)
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, 2, sep = "_")
table(rawdat$clf)

# Remove invariable columns


# Remove the holdout test predictions
# dat <- rawdat[-grep("reject", rawdat$nms),] # 658 observations
# dim(dat)

# Set seed 
set.seed(20190506)
dim(rawdat)
# dat <- rawdat[!duplicated(rawdat),]
dat <- rawdat %>%
  dplyr::select(-RADA880108_27)
dim(dat)
table(dat$clf)


# Split into test and training data
dat_split <- initial_split(dat, strata = "clf")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("nms", "clf")]
x_test <- dat_test[,!colnames(dat_test) %in% c("nms", "clf")]
y_train <- as.factor(dat_train$clf)
y_test <- as.factor(dat_test$clf)
table(y_test)


# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$nms)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nms)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nms)

###==================================================================###

####### Random Forest with feature selection ########
# Random Forest using the ranger package for prediction
# Fitting mtry = 2, splitrule = gini, min.node.size = 1 on full training set

rf <- train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 5,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  # tuneGrid = tgrid,
  num.trees = 1000,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  importance = "permutation")
# preProcess = c("center", "scale"))

# saveRDS(rf, "data/model_comparisons/rf_xvalidated_small_subst_grp_repeatedCV_20190406.rds")
rfcv$results
rfcv$bestTune
getTrainPerf(rfcv)
rf <- readRDS("data/model_comparisons/rf_xvalidated_small_subst_grp_repeatedCV_20190406.rds")


# ROC curve
approx_roc_curve <- function(x, label) {
  x %>%
    pluck("pred") %>%
    roc_curve(obs, y_train) %>%
    mutate(model = label)
}

# ROC curve is not looking good...
approx_roc_curve(rf, "Random Forest") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

# Confusion matrix
getTrainPerf(rf)

# ROC curve
approx_roc_curve <- function(x, label) {
  x %>%
    pluck("pred") %>%
    roc_curve(obs, y_train) %>%
    mutate(model = label)
}
# saveRDS(rf, "data/model_comparisons/rf_xvalidated_small_subst_grp_20190406.rds")

# ROC curve is not looking good...
approx_roc_curve(rf, "Random Forest") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

# Confusion matrix
getTrainPerf(rf)

# Try prediction
rf_ml <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = as.character(rf$bestTune$splitrule), 
                mtry = rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
                importance = "permutation")
saveRDS(rf_ml, "data/model_comparisons/rf_NRPS_broad_monomers_repeatedCV_20190406.rds")
rf_ml <- readRDS("data/model_comparisons/rf_NRPS_broad_monomers_repeatedCV_20190406.rds")

rf_pred <- predict(rf, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$clf))
cm_rf
0.7975 - 0.7276
sink("data/model_comparisons/cm_rf_15aa_functional_class_repeatedCV_20190406.txt")
cm_rf
sink()

vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T), 
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:25)
vimp

pdf("data/model_comparisons/rf_feat_select_var_imp_NRPS_broad_monomers_repeatedCV_20190406.pdf", width = 6, height = 6)
ggplot(data = vimp, 
       aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
dev.off()


dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "F1", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "data/model_comparisons/rf_AA34_NRPS_Precision_Recall_Accuracy_Substrate_monomers.csv", row.names = T)

## Make a heatmap of confusion matrix results
cm_list <- list(cm_rf$table)
names(cm_list) <- c("rf_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Reference, Prediction)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "royalblue") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("data/model_comparisons/aa34_rf_NRPS_conf_matrices_subst_grp_20190406.jpeg", height=7, width=7, units='in')

###Naive Bayes method
nb_grid <- expand.grid(usekernel = TRUE, fL = 0, adjust = 1)

nb <- train(
  x = df_train, 
  y = y_train,
  method = "nb",
  tuneGrid = nb_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, 
                           repeats = 5,
                           verboseIter = T, classProbs = F,
                           savePredictions = "final"))

# Confusion matrix
getTrainPerf(nb_ml)
nb_ml$results

# saveRDS(nb, "data/model_comparisons/nb_NRPS_broad_monomers_repeatedCV_20190406.rds")
nb_ml <- readRDS("data/model_comparisons/nb_NRPS_broad_monomers_repeatedCV_20190406.rds")
nb_pred <- predict(nb_ml, newdata = form_test)

cm_nb <- confusionMatrix(nb_pred, as.factor(dat_test$clf))
cm_nb
0.7178 - 0.6421
dtl_feat_select <- data.frame(round(cm_nb$byClass[,colnames(cm_nb$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "data/model_comparisons/nb_AA34_Precision_Recall_Accuracy_Substrate_monomers.csv", row.names = T)

## Make a heatmap of confusion matrix results
cm_list <- list(cm_nb$table)
names(cm_list) <- c("nb_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "royalblue") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("data/model_comparisons/aa34_nb_conf_matrices_subst_grp_repeatedCV_20190406.jpeg", height=7, width=7, units='in')
