pos_seqs <- read_uniprot(paste0(pathway, "sept_signal.txt"), euk = TRUE)
neg_seqs <- read.fasta(paste0(pathway, "sept_neg.fasta"), seqtype = "AA")cleaves <- sapply(pos_seqs, function(i) attr(i, "sig")[2])
atyp_aa <- which(sapply(neg_seqs, function(i) any(i %in% c("X", "J", "Z", "B"))))
too_short <- which(sapply(neg_seqs, length) < 50)
neg_seqs_pure <- neg_seqs[-unique(c(atyp_aa, too_short))]

pos_freqs <- t(sapply(pos_seqs, function(i) {
  sig_length <- attr(i, "sig")[2]
  table(factor(i[1L:sig_length], levels = a()[-1]))/sig_length
}))

negs <- sample(1L:length(neg_seqs_pure), length(pos_seqs))
neg_seqs_pure_small <- neg_seqs_pure[negs]

neg_freqs <- t(sapply(1L:length(pos_seqs), function(i) {
  sig_length <- attr(pos_seqs[[i]], "sig")[2]
  table(factor(neg_seqs_pure_small[[i]][1L:sig_length], levels = a()[-1]))/sig_length
}))

whole_table <- rbind(cbind(SP = rep("positive", length(pos_seqs)), data.frame(pos_freqs)),
      cbind(SP = rep("negative", length(pos_seqs)), data.frame(neg_freqs)))
write.csv2(whole_table, file = "aa_freqs.csv")
