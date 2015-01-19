#evaluation of cross-validation
load(paste0(pathway, "cleave_pred100.RData"))

single_cv <- multifolds_cl_work[[1]][[1]]

#only predicted as positive
rowMeans(sapply(multifolds_cl_work, function(five_cv)
  rowMeans(sapply(five_cv, function(single_cv) {
    real_cs <- single_cv[single_cv[, 1] > 0.5 & !is.na(single_cv[, 4]), 4]
    hsmm_cs <- single_cv[single_cv[, 1] > 0.5 & !is.na(single_cv[, 4]), 2]
    rf_cs <- single_cv[single_cv[, 1] > 0.5 & !is.na(single_cv[, 4]), 3]
    c(mean(sqrt((real_cs - hsmm_cs)^2)), mean(sqrt((real_cs - rf_cs)^2)))
  }))))

#all really positive
rowMeans(sapply(multifolds_cl_work, function(five_cv)
  rowMeans(sapply(five_cv, function(single_cv) {
    real_cs <- single_cv[!is.na(single_cv[, 4]), 4]
    hsmm_cs <- single_cv[!is.na(single_cv[, 4]), 2]
    rf_cs <- single_cv[!is.na(single_cv[, 4]), 3]
    c(mean(sqrt((real_cs - hsmm_cs)^2)), mean(sqrt((real_cs - rf_cs)^2)))
  }))))
