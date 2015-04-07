load(paste0(pathway, "cleave_pred80names.RData")) -> tmp
names80 <- multifolds_cl_work

load(paste0(pathway, "cleave_pred170names.RData")) -> tmp
names170 <- multifolds_cl_work

multifolds_cl_work <- c(names80, names170)
rm(names80, names170)

prot_names <- aggregate(Freq ~ all_cv, 
                        do.call(rbind, lapply(multifolds_cl_work, function(single_fold) {
                          all_cv <- unlist(lapply(single_fold, function(single_cv) {
                            sub_dat <- single_cv[!is.na(single_cv[, "real"]),
                                                 c("sp_probability", "real")]
                            rownames(sub_dat[sub_dat[, "sp_probability"] < 0.5, ])  
                          }))
                          data.frame(table(all_cv))
                        })), sum)
rm(multifolds_cl_work)

prot_names <- prot_names[order(prot_names[, 2], decreasing = TRUE), ]
#821 from 3897

prot_names[prot_names[, "Freq"] == 1000, ]

#commented to not overwrite
#write.csv2(prot_names, paste0(pathway, "FN.csv"))
