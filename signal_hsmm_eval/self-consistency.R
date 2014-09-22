model_self <- hsmm(pos_seqs, aaaggregation)
preds_self <- predict.signal.hsmm(model_self, pos_seqs)

preds <- pred2df(preds_self)

strange_seqs <- pos_seqs[rownames(preds)[preds[, 1] < 0.5]]

write.fasta(strange_seqs, names(strange_seqs), "false_negatives.fasta")

content <- readLines("sept_signal.txt")
prot_ids <- grep("\\<ID   ", content)
prot_names <- content[grep("\\<ID   ", content)]

res <- sapply(names(strange_seqs), function(i) grepl(i, prot_names))

total <- res[, 1]
for(i in 2L:ncol(res))
  total <- total + res[, i]

fn_records <- which(as.logical(total))

sink(file = "false_negatives.txt")
lapply(fn_records, function(i) {
  cat(content[prot_ids[i]:(prot_ids[i + 1] - 1)], "\n", sep = "\n")
})
sink()

