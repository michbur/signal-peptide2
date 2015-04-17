pos_seqs <- read_uniprot(paste0(pathway, "sept_signal.txt"), euk = TRUE)
neg_seqs <- read.fasta(paste0(pathway, "sept_neg.fasta"), seqtype = "AA")
#remove sequences with atypical aminoacids
atyp_aa <- which(sapply(neg_seqs, function(i) any(i %in% c("X", "J", "Z", "B", "U"))))
too_short <- which(sapply(neg_seqs, length) < 80)
neg_seqs <- neg_seqs[-unique(c(atyp_aa, too_short))]

too_short <- which(sapply(pos_seqs, length) < 80)
pos_seqs <- pos_seqs[-c(too_short)]

cleaves <- tolower(t(sapply(pos_seqs, function(seq) {
  cleave_site <- attr(seq, "sig")[2]
  seq[(cleave_site-4):(cleave_site+3)]
})))

post_cleaves <- tolower(t(sapply(pos_seqs, function(seq) {
  cleave_site <- attr(seq, "sig")[2]
  seq[(cleave_site):(cleave_site+7)]
})))

pre_cleaves <- tolower(t(sapply(pos_seqs, function(seq) {
  cleave_site <- attr(seq, "sig")[2]
  seq[(cleave_site-6):(cleave_site+1)]
})))

tar <- c(rep(1, nrow(cleaves)), rep(0, 2*nrow(cleaves)))

all_data <- degenerate(matrix(l2n(rbind(cleaves, post_cleaves, pre_cleaves), 
                                  seq_type = "prot"), ncol = 8),
                       list(`1` = c(1, 6, 8, 10, 11, 18),
                            `2` = c(2, 13, 14, 16, 17),
                            `3` = c(5, 19, 20),
                            `4` = c(7, 9, 12, 15),
                            `4` = c(3, 4)))

pos_id <- which(tar == 1)
neg_id <- which(tar == 0)



#here fold lapply

res <- lapply(1L:20, function(repetition) {
  folds_pos <- cvFolds(length(pos_id), type = "random")
  folds_neg <- cvFolds(length(neg_id), type = "random")
  pblapply(1L:5, function(fold) {
    pos_train <- folds_pos[["subsets"]][folds_pos[["which"]] == fold]
    pos_test <- folds_pos[["subsets"]][folds_pos[["which"]] != fold]
    
    neg_train <- folds_neg[["subsets"]][folds_neg[["which"]] == fold]
    neg_test <- folds_neg[["subsets"]][folds_neg[["which"]] != fold]
    
    constr <- construct_ngrams(c(rep(1, length(pos_train)), 
                                 rep(0, length(neg_train))), 
                               all_data[c(pos_train, neg_train), ], 
                               1L:5, 4, gap = FALSE)
    constr_gap <- construct_ngrams(c(rep(1, length(pos_train)), 
                                     rep(0, length(neg_train))), 
                                   all_data[c(pos_train, neg_train), ],
                                   1L:5, 4, gap = TRUE)
    
    lapply(1L:4, function(n) {
      no_gap <- train_models(pos_train, neg_train, all_data, constr[[n]])
      gap <- train_models(pos_train, neg_train, all_data, constr_gap[[n]])
      cbind(fold = rep(fold, 6), n = rep(n, 6), 
            gap = c(rep(FALSE, 3), rep(TRUE, 3)),
            test_models(pos_test, neg_train, all_data, constr[[n]], no_gap),
            test_models(pos_test, neg_test, all_data, constr_gap[[n]], gap))
    })
  })
})

save(res, file = paste0(pathway, "build_cleave_20.RData"))
