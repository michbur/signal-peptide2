pos_seqs <- read_uniprot(paste0(pathway, "sept_signal.txt"), euk = TRUE)
neg_seqs <- read.fasta(paste0(pathway, "sept_neg.fasta"), seqtype = "AA")
#remove sequences with atypical aminoacids
atyp_aa <- which(sapply(neg_seqs, function(i) any(i %in% c("X", "J", "Z", "B", "U"))))
too_short <- which(sapply(neg_seqs, length) < 50)
neg_seqs <- neg_seqs[-unique(c(atyp_aa, too_short))]




pos_ids <- cvFolds(length(pos_seqs), K = 5)
neg_ids <- cvFolds(length(neg_seqs), K = 5)

i = 1
train_dat <- pos_ids[["which"]] == i


# #cross-validation --------------------------------
# #proteins with signal peptides
# 
# cl <- makeCluster(4, type = "SOCK")
# clusterExport(cl, c("predict.signal.hsmm", "signal.hsmm_decision"))
# 
# multifolds_cl_work <- pblapply(1L:100, function(dummy_variable) { 
#   pos_ids <- cvFolds(length(pos_seqs), K = 5)
#   cv_neg <- neg_seqs_pure[sample(1L:length(neg_seqs_pure), length(pos_seqs))]
#   
#   fold_res <- lapply(1L:5, function(fold) {
#     model_cv <- hsmm(pos_seqs[pos_ids[[4]][,][pos_ids[[5]] == fold]], aaaggregation)
#     test_dat <- c(pos_seqs[pos_ids[[4]][,][pos_ids[[5]] != fold]],
#                   cv_neg[pos_ids[[4]][,][pos_ids[[5]] != fold]])
#     parLapply(cl, 1L:length(test_dat), function(protein_id) try({
#       library(signal.hsmm)
#       predict.signal.hsmm(model_cv, test_dat[[protein_id]])
#     }, silent = TRUE))
#   })
#   lapply(fold_res, function(fold)
#     sapply(fold, function(i)
#       if(class(i) != "try-error") {
#         i[[1]][["sp_probability"]]
#       } else {
#         NA
#       }))
# })
# 
# stopCluster(cl)