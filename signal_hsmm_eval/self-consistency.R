model_self <- hsmm(pos_seqs, aaaggregation)
preds_self <- predict.signal.hsmm(model_self, pos_seqs)
