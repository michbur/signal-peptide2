jack_pos_prob <- sapply(jack_pos, function(i) i[[1]][["sp_probability"]])
jack_neg_class <- sapply(jack_neg, function(i) class(i))
#algorithm does not work all the time
which(jack_neg_class == "try-error")

jack_neg_prob <- sapply(jack_neg[jack_neg_class != "try-error"], function(i) 
  i[[1]][["sp_probability"]])

library(hmeasure)

real_labels <- c(rep(1, length(jack_pos_prob)),
                 rep(0, length(jack_neg_prob)))
jack_res <- HMeasure(real_labels, c(jack_pos_prob, jack_neg_prob))[["metrics"]]
TP <- as.numeric(jack_res[["TP"]])
FP <- as.numeric(jack_res[["FP"]])
TN <- as.numeric(jack_res[["TN"]])
FN <- as.numeric(jack_res[["FN"]])

jack_res <- t(jack_res)
jack_res <- rbind(jack_res, MCC = (TP*TN - FP*FN)/((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)))
write.table(jack_res, file = "jackknife_table.txt", sep = "\t")

real_labels2 <- c(rep(1, length(jack_pos_prob)),
                  rep(0, length(jack_pos_prob)))
jack_res2 <- HMeasure(real_labels2, 
                      c(jack_pos_prob, 
                        jack_neg_prob[sample(length(jack_neg_prob), 
                                             length(jack_pos_prob))]))[["metrics"]]
TP <- as.numeric(jack_res2[["TP"]])
FP <- as.numeric(jack_res2[["FP"]])
TN <- as.numeric(jack_res2[["TN"]])
FN <- as.numeric(jack_res2[["FN"]])

jack_res2 <- t(jack_res2)
jack_res2 <- rbind(jack_res2, 
                   MCC = (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)))
write.table(jack_res2, file = "jackknife_table2.txt", sep = "\t")

TP <- 100
FP <- 1
TN <- 100
FN <- 1
(TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
