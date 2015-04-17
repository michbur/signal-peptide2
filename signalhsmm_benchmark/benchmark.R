
#real cleavage sites
cleave_sites <- sapply(benchmark_dat, function(protein) 
  ifelse(is.null(attr(protein, "sig")[2]), NA, attr(protein, "sig")[2]))
is_signal <- !(is.na(cleave_sites))

#train two instances of signal.hsmm


signal.hsmm2010 <- train_hsmm(read_uniprot(paste0(pathway, "pub_pos_train.txt"), euk = TRUE),
                              aaaggregation)

signal.hsmm1987 <- train_hsmm(read_uniprot(paste0(pathway, "pos_ultrahard_data.txt"), euk = TRUE),
                              aaaggregation)

signal.hsmm2010_preds <- pred2df(predict(signal.hsmm2010, benchmark_dat))

signal.hsmm1987_preds <- pred2df(predict(signal.hsmm1987, benchmark_dat))

colnames(signal.hsmm2010_preds) <- c("signal.peptide", "sp.start", "sig.end")
colnames(signal.hsmm1987_preds) <- c("signal.peptide", "sp.start", "sig.end")

bench_list <- list(phobius = read_phobius("./signalhsmm_benchmark/eval_phobius.txt"),
                   predsi = read_predsi("./signalhsmm_benchmark/eval_predsi.txt"),
                   philius = read_philius("./signalhsmm_benchmark/eval_philius.xml"),
                   spnotm = read_signalp41("./signalhsmm_benchmark/eval_signalp41notm.txt"),
                   sptm = read_signalp41("./signalhsmm_benchmark/eval_signalp41tm.txt"),
                   signalhsmm2010 = signal.hsmm2010_preds,
                   signalhsmm1987 = signal.hsmm1987_preds)

bench_list[["signalhsmm2010"]][["sig.end"]] <- bench_list[["signalhsmm2010"]][["sig.end"]] - 1
bench_list[["signalhsmm1987"]][["sig.end"]] <- bench_list[["signalhsmm1987"]][["sig.end"]] - 1

positions <- sapply(bench_list, function(predictor)
  mean(abs(cleave_sites[is_signal & predictor[["signal.peptide"]] > 0.5] - 
             predictor[is_signal & predictor[["signal.peptide"]] > 0.5, "sig.end"]),
       na.rm = TRUE))

metrics <- do.call(rbind, lapply(bench_list, function(predictor)
  HMeasure(is_signal, predictor[["signal.peptide"]])[["metrics"]]))[, c("AUC", "H", "Gini", ), ]

benchmark_dat <- data.frame(metrics, positions = positions)
save(benchmark_dat, file = "benchmark_dat.RData")


lapply(bench_list[6L:7], function(predictor)
  cleave_sites[is_signal & predictor[["signal.peptide"]] > 0.5] - 2 -
             predictor[is_signal & predictor[["signal.peptide"]] > 0.5, "sig.end"]
  )

