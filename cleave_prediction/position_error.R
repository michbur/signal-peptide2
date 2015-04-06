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

round(median(abs(tcsdf[["real"]] - tcsdf[["sp_end"]] + 1)), 4)

ggplot(tcsdf, aes(x = abs(sp_end - real + 1))) + 
  geom_histogram() + 
  scale_y_continuous(name = "Density") +
  scale_x_continuous("Cleavage site error") + 
  xlim(0, 30)