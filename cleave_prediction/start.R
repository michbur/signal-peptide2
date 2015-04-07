library(seqinr)
library(hmeasure)
library(pbapply)
library(cvTools)
library(parallel)
library(snow)
library(biogram)
library(signal.hsmm)
library(randomForest)
library(ggplot2)
library(hmeasure)

# working directory ------------------------------

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"

# SIGNAL-HSMM ------------------------------------


calc_t <- function(list_prots, aa_list) {
  nhc <- t(vapply(list_prots, find_nhc, rep(0, 4)))
  
  n_region <- NULL
  h_region <- NULL
  c_region <- NULL
  rest <- NULL
  
  for(i in 1L:length(list_prots)){
    region_starts <- nhc[i, ]
    n_region <- c(n_region, list_prots[[i]][1:(region_starts[2] - 1)])
    h_region <- c(h_region, list_prots[[i]][region_starts[2]:(region_starts[3] - 1)])
    c_region <- c(c_region, list_prots[[i]][region_starts[3]:(region_starts[4] - 1)])
    rest <- c(rest, list_prots[[i]][region_starts[4]:length(list_prots[[i]])])
  }
  
  t1 <- rep(0, length(aa_list))
  temp <- table(degenerate(n_region, aa_list))
  t1[as.numeric(names(temp))] <- temp
  names(t1) <- 1:length(aa_list)
  
  t2 <- rep(0, length(aa_list))
  temp <- table(degenerate(h_region, aa_list))
  t2[as.numeric(names(temp))] <- temp
  names(t2) <- 1:length(aa_list)
  
  t3 <- rep(0, length(aa_list))
  temp <- table(degenerate(c_region, aa_list))
  t3[as.numeric(names(temp))] <- temp
  names(t3) <- 1:length(aa_list)
  
  t4 <- rep(0, length(aa_list))
  temp <- table(degenerate(rest, aa_list))
  t4[as.numeric(names(temp))] <- temp
  names(t4) <- 1:length(aa_list)
  
  len_c <- nhc[, "cs"] - nhc[, "start_c"]
  len_h <- nhc[, "start_c"] - nhc[, "start_h"]
  len_n <- nhc[, "start_h"] - nhc[, "start_n"]
  lengths <- matrix(c(len_n, len_h, len_c), ncol = 3)
  colnames(lengths) <- c("n", "h", "c")
  
  list(mean_cs = mean(nhc[, 4]), sd_cs = sd(nhc[, 4]), t1 = t1, t2 = t2, t3 = t3, t4 = t4, 
       lengths = lengths)
}

measure_region <- function(region, max.length = 32) {
  lengths <- table(region)
  res <- rep(0, max.length)
  lengths <- lengths[as.numeric(names(lengths))>0] #removing lengths smaller than 1
  
  start_l <- min(as.numeric(names(lengths)))
  end_l <- max(as.numeric(names(lengths)))
  if(prod(start_l:end_l %in% as.numeric(names(lengths)))){
    max_length <- length(lengths) #if all lengths are present in training set
  } else{
    max_length <- 1
    sl <- sum(lengths)
    while(sum(lengths[1:max_length])/sl <= 0.51) {
      max_length <- which.min(start_l:end_l %in% as.numeric(names(lengths))) - 1
      start_l <- start_l + 1  #to assure that max_length wouldn't be too small
      max_length <- ifelse(max_length == 0, length(lengths), max_length)
    }
  }
  max_length <- min(max_length, max.length)
  
  prop_lengths <- lengths[1:max_length]/sum(lengths[1:max_length])
  res[as.numeric(names(prop_lengths))] <- prop_lengths
  res
}


hsmm <- function(train_data, aa_group, max.length = 32) {
  train_data <- lapply(train_data, toupper)
  ts <- calc_t(train_data, aa_group)
  
  t1 <- ts[["t1"]]
  t2 <- ts[["t2"]]
  t3 <- ts[["t3"]]
  t4 <- ts[["t4"]]
  
  overall <- t4 #table(degenerate(unlist(analized_sequences), aa5))
  overall.probs <- overall/sum(overall)          
  overall.probs.log = log(overall.probs) #for viterbi
  
  lengths <- ts[["lengths"]]
  params <- apply(lengths, 2, measure_region, max.length = max.length)
  params <- cbind(params, rep(1/max.length, max.length))
  
  #setting params for hmm -------
  additional_margin = 10
  pipar <- c(1,0,0,0)
  tpmpar <- matrix(c(0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1,
                     0, 0, 0, 0), 4, byrow = TRUE)
  od <- matrix(c((t1/sum(t1))[1:4],
                 (t2/sum(t2))[1:4],
                 (t3/sum(t3))[1:4],
                 (t4/sum(t4))[1:4]), 4, byrow = TRUE)
  
  list(pipar = pipar, tpmpar = tpmpar, od = od, 
       overall_probs_log = overall.probs.log, params = params)
}

#make it later into a method
predict.signal.hsmm <- function(signal.hsmm_model, test_data) {
  
  if (class(test_data) == "numeric" || class(test_data) == "factor" || 
        class(test_data) == "data.frame" || class(test_data) == "matrix")
    stop("Input data must have class 'SeqFastaAA', 'character' or 'list'.")
  
  
  if(class(test_data) == "SeqFastaAA" || 
       class(test_data) == "character") {
    #single input
    decisions <- signal.hsmm_decision(test_data, aa_group = aaaggregation, 
                                      pipar = signal.hsmm_model[["pipar"]], 
                                      tpmpar = signal.hsmm_model[["tpmpar"]], 
                                      od = signal.hsmm_model[["od"]], 
                                      overall_probs_log = signal.hsmm_model[["overall_probs_log"]], 
                                      params = signal.hsmm_model[["params"]])
    decisions <- list(decisions)
    names(decisions) <- attr(test_data, "name")
  } else {
    #list input
    decisions <- lapply(test_data, function(prot)
      signal.hsmm_decision(prot, aa_group = aaaggregation, 
                           pipar = signal.hsmm_model[["pipar"]], 
                           tpmpar = signal.hsmm_model[["tpmpar"]], 
                           od = signal.hsmm_model[["od"]], 
                           overall_probs_log = signal.hsmm_model[["overall_probs_log"]], 
                           params = signal.hsmm_model[["params"]]))
  }
  class(decisions) <- "hsmm_pred_list"
  decisions
}

signal.hsmm_decision <- function(prot, aa_group, pipar, tpmpar, 
                                 od, overall_probs_log, params) {
  if (length(prot) == 1) {
    prot <- strsplit(prot, "")[[1]]
    if ("name" %in% names(attributes(prot)))
      attr(prot, "name") <- "undefined_name"
    if (length(prot) == 1)
      stop("Input sequence is too short.")
  }
  if(!is_protein(prot))
    stop("Atypical aminoacids detected, analysis cannot be performed.")
  
  deg_sample <- as.numeric(degenerate(toupper(prot)[1L:50], aa_group))
  #remove atypical amino acids
  deg_sample <- na.omit(deg_sample)
  viterbi_res <- duration_viterbi(deg_sample, pipar, tpmpar, od, params)
  viterbi_path <- viterbi_res[["path"]]
  c_site <- ifelse(any(viterbi_path == 3), 
                   max(which(viterbi_path == 3)), 
                   length(deg_sample))
  #get probabilities of signal peptide model
  prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
  #get probabilities of no signal peptide model
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
  prob.total <- exp(prob.signal - prob.non)
  res <- list(sp_probability = unname(1 - 1/(1 + prob.total)),
              sp_start = 1,
              sp_end = c_site,
              struc = viterbi_path,
              prot = toupper(prot[1L:70]),
              name = attr(prot, "name"),
              str_approx = 0)
  class(res) <- "hsmm_pred"
  
  #structure approximation - if atypical (normally negative signal peptide)
  while(!all(1L:4 %in% res[["struc"]])) {
    res[["struc"]] <- c(res[["struc"]], which.min(1L:4 %in% res[["struc"]]))
    res[["str_approx"]] <- res[["str_approx"]] + 1
  }
  
  res
}

#Selection of Functional Signal Peptide Cleavage Sites from a Library of Random Sequences
cleave_aggregation <- list(`1` = c("k", "r", "h"),
                           `2` = c("a", "g", "s", "v", "i"), 
                           `3` = c("l", "m", "f", "w", "c", "u"), 
                           `4` = c("t", "n", "q", "d", "e", "p", "y"))


