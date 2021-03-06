#old version no protein names here
#evaluation of cross-validation
load(paste0(pathway, "cleave_pred100.RData"))

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
#hsmm + rf is subpar to only hsmm

#1L:10, because of the slow machine
cs_tab <- do.call(rbind, lapply(multifolds_cl_work[1L:10], function(five_cv)
  do.call(rbind, lapply(five_cv, function(single_cv) 
    data.frame(single_cv[, -3]
  )))))

#predicted signal peptides which are really signal peptides ---------------
#pred_cs <- cs_tab[cs_tab[["sp_probability"]] > 0.5 & !is.na(cs_tab[["real"]]), ]

csdf <- data.frame(cs_tab, sp = as.factor(as.numeric(!is.na(cs_tab[["real"]]))), 
           predb = as.factor(round(cs_tab[["sp_probability"]], 0)))

ggplot(csdf, aes(x = sp_probability, fill = sp)) + geom_density()

#sp length versus prediction succes -----------------------------
#true cleavage site data frame
tcsdf <- csdf[!is.na(cs_tab[["real"]]), c("sp_probability", "sp_end", "real")]
ggplot(tcsdf, aes(sp_end)) + geom_bar()

data.frame(table(tcsdf[["real"]]))
tcsdf <- data.frame(tcsdf, 
                    len_gr = cut(tcsdf[["real"]], c(6, 15, 24, 35, 89)))


sapply(levels(tcsdf[["len_gr"]]), function(i)
  mean(tcsdf[tcsdf[["len_gr"]] == i, "sp_probability"]))

#unreadable
ggplot(tcsdf, aes(x = abs(sp_end - real), fill = len_gr)) + 
  geom_density(alpha = 0.25)

ggplot(tcsdf, aes(x = sp_probability, fill = len_gr)) + 
  geom_density(alpha = 0.25)

#bad prediction for long signal peptides
ggplot(tcsdf[tcsdf[["len_gr"]] == "(35,89]", ], aes(x = sp_probability, y = real)) +
  geom_point()
