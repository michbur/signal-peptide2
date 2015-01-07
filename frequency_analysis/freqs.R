pos_seqs <- read_uniprot(paste0(pathway, "sept_signal.txt"), euk = TRUE)
neg_seqs <- read.fasta(paste0(pathway, "sept_neg.fasta"), seqtype = "AA")
#remove sequences with atypical aminoacids
atyp_aa <- which(sapply(neg_seqs, function(i) any(i %in% c("X", "J", "Z", "B", "U"))))
too_short <- which(sapply(neg_seqs, length) < 50)
neg_seqs_pure <- neg_seqs[-unique(c(atyp_aa, too_short))]

library(seqinr)
library(signal.hsmm)
library(biogram)
library(tree)
library(randomForest)
library(kernlab)
library(e1071)
library(klaR)
library(pbapply)
library(cvTools)
library(nnet)
library(party)
library(ipred)
#here start loop
#shuffle pos

multifolds <- lapply(1L:2, function(n_gram) {
  
  pos_freqs <- t(pbsapply(pos_seqs, function(i) {
    sig_length <- attr(i, "sig")[2]
    as.matrix(count_ngrams(i[1L:sig_length], n_gram, a()[-1]))/sig_length
  }))
  
  pblapply(1L:150, function(dummy_variable) { 
    pos_ids <- cvFolds(nrow(pos_freqs), K = 5)
    
    cv_neg <- neg_seqs_pure[sample(1L:length(neg_seqs_pure), nrow(pos_freqs))]
    neg_freqs <- t(sapply(1L:nrow(pos_freqs), function(i) {
      sig_length <- attr(pos_seqs[[i]], "sig")[2]
      as.matrix(count_ngrams(cv_neg[[i]][1L:sig_length], n_gram, a()[-1]))/sig_length
    }))
    
    fold_res <- lapply(1L:5, function(fold) try({
      train_data <- data.frame(cbind(data.frame(rbind(pos_freqs[pos_ids[[5]] == fold, ],
                                                      neg_freqs[pos_ids[[5]] == fold, ])),
                                     target = c(rep("pos", sum(pos_ids[[5]] == fold)), 
                                                rep("neg", sum(pos_ids[[5]] == fold)))))
      
      test_data <- data.frame(cbind(data.frame(rbind(pos_freqs[pos_ids[[5]] != fold, ],
                                                     neg_freqs[pos_ids[[5]] != fold, ])),
                                    target = c(rep("pos", sum(pos_ids[[5]] != fold)), 
                                               rep("neg", sum(pos_ids[[5]] != fold)))))
      
      modelRF <- randomForest(target ~ ., data = train_data)
      modelSVMlin <- ksvm(target ~ ., data = train_data, prob.model = TRUE,
                          kernel = "vanilladot")
      modelSVMrbf <- ksvm(target ~ ., data = train_data, prob.model = TRUE,
                          kernel = "rbfdot")
      modelNB <- naiveBayes(target ~ ., data = train_data)
      modelNN <- nnet(target ~ ., data = train_data, size=10, MaxNWts = 2e6)
      modelKNN <-ipredknn(target ~ ., data = train_data)
      
      RF_pos <- predict(modelRF, test_data, type = "prob")[, 2]
      SVMlin_pos <- predict(modelSVMlin, test_data, type = "prob")[, 2]
      SVMrbf_pos <- predict(modelSVMrbf, test_data, type = "prob")[, 2]
      NB_pos <- predict(modelNB, test_data, type = "raw")[, 2]
      NN_pos <- predict(modelNN, test_data, type="raw")
      KNN_pos <- predict(modelKNN, test_data)
      
      data.frame(SVMlin_pos = SVMlin_pos,
                 SVMrbf_pos = SVMrbf_pos,
                 NB_pos = NB_pos,
                 RF_pos = RF_pos,
                 NN_pos = NN_pos,
                 KNN_pos = KNN_pos)
    }, silent = TRUE))
  })
})
save(multifolds, file = paste0(pathway, "freq_analysis.RData"))


