#position error for recomb poster

load(paste0(pathway, "cleave_pred80names.RData")) -> tmp
names80 <- multifolds_cl_work

load(paste0(pathway, "cleave_pred170names.RData")) -> tmp
names170 <- multifolds_cl_work

multifolds_cl_work <- c(names80, names170)
rm(names80, names170)


cs_tab <- do.call(rbind, lapply(multifolds_cl_work, function(five_cv)
  do.call(rbind, lapply(five_cv, function(single_cv) 
    #lets omit random forest predition
    data.frame(single_cv[, -3]
    )))))

#predicted signal peptides which are really signal peptides ---------------
#pred_cs <- cs_tab[cs_tab[["sp_probability"]] > 0.5 & !is.na(cs_tab[["real"]]), ]

csdf <- data.frame(cs_tab, sp = as.factor(as.numeric(!is.na(cs_tab[["real"]]))), 
                   predb = as.factor(round(cs_tab[["sp_probability"]], 0)))


#sp length versus prediction succes -----------------------------
#true cleavage site data frame
tcsdf <- csdf[!is.na(cs_tab[["real"]]), c("sp_probability", "sp_end", "real")]
#write.csv2(tcsdf, file = "tcsdf_poster.csv")


poster_data <- list(mean = mean(abs(tcsdf[["real"]] - tcsdf[["sp_end"]] - 1)),
                    median = mean(abs(tcsdf[["real"]] - tcsdf[["sp_end"]] - 1)),
                    position_table = data.frame(table(abs(tcsdf[["real"]] - tcsdf[["sp_end"]] - 1))),
                    metrics = do.call(cbind, lapply(multifolds_cl_work[1L:5], function(five_cv) {
                      rowMeans(sapply(five_cv, function(single_cv) 
                        #lets omit random forest predition
                        unlist(HMeasure(!is.na(single_cv[, "real"]), single_cv[, "sp_probability"])[["metrics"]])
                      ))})))
