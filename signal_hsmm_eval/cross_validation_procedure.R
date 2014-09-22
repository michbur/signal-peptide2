pos_seqs <- read_uniprot("sept_signal.txt", euk = TRUE)
neg_seqs <- read.fasta("sept_neg.fasta", seqtype = "AA")
#remove sequences with atypical aminoacids
atyp_aa <- which(sapply(neg_seqs, function(i) any(i %in% c("X", "J", "Z", "B"))))
too_short <- which(sapply(neg_seqs, length) < 50)
neg_seqs_pure <- neg_seqs[-unique(c(atyp_aa, too_short))]



#jackknife --------------------------------
#proteins with signal peptides

jack_pos <- pblapply(1L:length(pos_seqs), function(protein_id) try({
  model_jack <- hsmm(pos_seqs[-protein_id], aaaggregation)
  predict.signal.hsmm(model_jack, pos_seqs[[protein_id]])
}, silent = TRUE))

model_jack_full <- hsmm(pos_seqs, aaaggregation)
jack_neg <- pblapply(1L:length(neg_seqs_pure), function(protein_id) try({
  predict.signal.hsmm(model_jack_full, neg_seqs_pure[[protein_id]])
}, silent = TRUE))

save(jack_pos, jack_neg, file = "jackknife.RData")

#cross-validation --------------------------------
#proteins with signal peptides

cl <- makeCluster(4, type = "SOCK")
clusterExport(cl, c("predict.signal.hsmm", "signal.hsmm_decision"))

multifolds_cl_work <- pblapply(1L:100, function(dummy_variable) { 
  pos_ids <- cvFolds(length(pos_seqs), K = 5)
  cv_neg <- neg_seqs_pure[sample(1L:length(neg_seqs_pure), length(pos_seqs))]
  
  fold_res <- lapply(1L:5, function(fold) {
    model_cv <- hsmm(pos_seqs[pos_ids[[4]][,][pos_ids[[5]] == fold]], aaaggregation)
    test_dat <- c(pos_seqs[pos_ids[[4]][,][pos_ids[[5]] != fold]],
                  cv_neg[pos_ids[[4]][,][pos_ids[[5]] != fold]])
    parLapply(cl, 1L:length(test_dat), function(protein_id) try({
      library(signal.hsmm)
      predict.signal.hsmm(model_cv, test_dat[[protein_id]])
    }, silent = TRUE))
  })
  lapply(fold_res, function(fold)
    sapply(fold, function(i)
      if(class(i) != "try-error") {
        i[[1]][["sp_probability"]]
      } else {
        NA
      }))
})

stopCluster(cl)

save(multifolds_cl_work, file = "crossval_part_work.RData")
save(multifolds_cl, file = "crossval_part.RData")



pos_ids <- cvFolds(length(pos_seqs), K = 5)
true_labels <- lapply(1L:5, function(fold) 
  c(rep(1, sum(pos_ids[[5]] != fold)),
    rep(0, sum(pos_ids[[5]] != fold))))

cv_res <- do.call(rbind, pblapply(c(multifolds_cl, multifolds_cl_work), 
                                  function(random_split) 
                                    do.call(rbind, lapply(1L:5, function(fold) {
                                      single_pred <- random_split[[fold]]
                                      NAs <- is.na(single_pred)
                                      res <- HMeasure(true_labels[[fold]][!NAs], single_pred[!NAs])[["metrics"]]
                                      TP <- as.numeric(res[["TP"]])
                                      FP <- as.numeric(res[["FP"]])
                                      TN <- as.numeric(res[["TN"]])
                                      FN <- as.numeric(res[["FN"]])
                                      cbind(res, MCC = (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)))
                                    }))))

save(cv_res, file = "crossval_full.RData")

