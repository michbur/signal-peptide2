library(seqinr)
library(hmeasure)
library(pbapply)
library(cvTools)
library(biogram)
library(signal.hsmm)
library(randomForest)
library(ggplot2)
library(hmeasure)
library(kernlab)
library(snow)

# working directory ------------------------------

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"

# functions ------------------------------------


train_models <- function(pos_train, neg_train, all_data, signif_ngrams) {
  train_dat <- all_data[c(pos_train, neg_train), ]
  train_df <- data.frame(as.matrix(count_specified(train_dat, signif_ngrams)), 
                         tar = as.factor(c(rep(1, length(pos_train)), 
                                           rep(0, length(neg_train)))
                         ))
  
  list(model_RF = randomForest(tar ~ ., train_df),
       model_rbf = ksvm(tar ~ ., train_df, kernel = "rbfdot"),
       model_lin = ksvm(tar ~ ., train_df, kernel = "vanilladot"))
}

test_models <- function(pos_test, neg_test, all_data, signif_ngrams, models) {
  test_dat <- all_data[c(pos_test, neg_test), ]
  test_df <- data.frame(as.matrix(count_specified(test_dat, signif_ngrams)))
  
  cbind(model = c("rf", "rbf", "linear"), 
        do.call(rbind, lapply(models, function(single_model)
          HMeasure(c(rep(1, length(pos_test)), 
                     rep(0, length(neg_test))), 
                   as.numeric(predict(single_model, test_df)) - 1)[["metrics"]]
        ))
  )
}