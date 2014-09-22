model_self <- hsmm(pos_seqs, aaaggregation)
preds_self <- predict.signal.hsmm(model_self, pos_seqs)
pred_self_df <- pred2df(preds_self)
save(pred_self_df, file = "selfcons_df.RData")
HMeasure(c(rep(1, length(preds_self)), 0), c(pred_self_df[, 1], 0))[["metrics"]]
