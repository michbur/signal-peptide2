#evaluation of cross-validation
load(paste0(pathway, "cleave_pred100.RData"))

single_cv <- multifolds_cl_work[[1]][[1]]

#only predicted as positive
sapply(multifolds_cl_work, function(five_cv)
  rowMeans(sapply(five_cv, function(single_cv) {
    real_cs <- single_cv[single_cv[, 1] > 0.5 & !is.na(single_cv[, 4]), 4]
    hsmm_cs <- single_cv[single_cv[, 1] > 0.5 & !is.na(single_cv[, 4]), 2]
    rf_cs <- single_cv[single_cv[, 1] > 0.5 & !is.na(single_cv[, 4]), 3]
    c(mean(sqrt((real_cs - hsmm_cs)^2)), mean(sqrt((real_cs - rf_cs)^2)))
  })))

#all really positive
rowMeans(sapply(multifolds_cl_work, function(five_cv)
  rowMeans(sapply(five_cv, function(single_cv) {
    real_cs <- single_cv[!is.na(single_cv[, 4]), 4]
    hsmm_cs <- single_cv[!is.na(single_cv[, 4]), 2]
    rf_cs <- single_cv[!is.na(single_cv[, 4]), 3]
    c(mean(sqrt((real_cs - hsmm_cs)^2)), mean(sqrt((real_cs - rf_cs)^2)))
  }))))

cs_tab <- do.call(rbind, lapply(multifolds_cl_work, function(five_cv)
  do.call(rbind, lapply(five_cv, function(single_cv) 
    data.frame(single_cv[, -3]
  )))))

#predicted signal peptides which are really signal peptides
#pred_cs <- cs_tab[cs_tab[["sp_probability"]] > 0.5 & !is.na(cs_tab[["real"]]), ]

csdf <- data.frame(cs_tab, sp = as.factor(as.numeric(!is.na(cs_tab[["real"]]))), 
           predb = as.factor(round(cs_tab[["sp_probability"]], 0)))

ggplot(csdf, aes(x = sp_probability))
