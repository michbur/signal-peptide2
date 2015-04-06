load(paste0(pathway, "cleave_pred80names.RData")) -> tmp
names80 <- multifolds_cl_work

load(paste0(pathway, "cleave_pred170names.RData")) -> tmp
names170 <- multifolds_cl_work

multifolds_cl_work <- c(names80, names170)
rm(names80, names170)



#calculate mean error of position prediction

mean(sapply(multifolds_cl_work, function(single_fold) {
  all_cv <- do.call(rbind, lapply(single_fold, function(single_cv) {
    sub_dat <- single_cv[!is.na(single_cv[, "real"]),
                         c("sp_end", "real")]
  }))
  mean(abs(all_cv[, "sp_end"] - all_cv[, "real"] - 1))
}))
#average position misplacement - 4.14 aa
