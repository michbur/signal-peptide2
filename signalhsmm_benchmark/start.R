read_predtat <- function(connection) {
  all_lines <- readLines(connection)
  #get decisions
  sig_pep <- grepl("Sec signal peptide predicted", all_lines)
  #get cleavage sites
  cleave_sites <- sapply(sapply(all_lines[sig_pep], function(i) 
    strsplit(i, "cleavage site: ", fixed = TRUE)[[1]][2], USE.NAMES = FALSE), function(j)
      strsplit(j, " ")[[1]][c(1, 3)], USE.NAMES = FALSE)
  sig_start <- rep(NA, length(all_lines))
  sig_end <- sig_start
  sig_start[sig_pep] <- as.numeric(cleave_sites[1, ])
  sig_end[sig_pep] <- as.numeric(cleave_sites[2, ])
  
  res <- data.frame(signal.peptide = sig_pep, 
                    sig.start = sig_start,
                    sig.end = sig_end) 
  rownames(res) <- sapply(strsplit(all_lines, " - "), function(i) i[1])
  res
}


read_signalp41 <- function(connection) {
  all_lines <- readLines(connection)
  do.call(rbind, lapply(all_lines, function(i) {
    line <- strsplit(i, " ")[[1]]
    line <- line[!line == ""]
    res <- data.frame(signal.peptide = line[10] == "Y",
                      sig.start = ifelse(line[10] == "Y", 1, NA),
                      sig.end = ifelse(line[10] == "Y", as.numeric(line[5]) - 1, NA))
    rownames(res) <- line[1]
    res
  }))
}

read_signalp3_nn <- function(connection) {
  all_lines <- readLines(connection)
  do.call(rbind, lapply(all_lines, function(i) {
    line <- strsplit(i, " ")[[1]]
    line <- line[!line == ""]
    res <- data.frame(signal.peptide = line[14] == "Y",
                      sig.start = ifelse(line[14] == "Y", 1, NA),
                      sig.end = ifelse(line[14] == "Y", as.numeric(line[6]) - 1, NA))
    rownames(res) <- line[1]
    res
  }))
}

read_signalp3_hmm <- function(connection) {
  all_lines <- readLines(connection)
  do.call(rbind, lapply(all_lines, function(i) {
    line <- strsplit(i, " ")[[1]]
    line <- line[!line == ""]
    res <- data.frame(signal.peptide = line[7] == "Y",
                      sig.start = ifelse(line[7] == "Y", 1, NA),
                      sig.end = ifelse(line[7] == "Y", as.numeric(line[4]) - 1, NA))
    rownames(res) <- line[1]
    res
  }))
}


read_predsi <- function(connection) {
  dat <- read.table(connection, sep = "\t")
  data.frame(signal.peptide = dat[, 4] == "Y",
             sig.start = ifelse(dat[, 4] == "Y", 1, NA),
             sig.end = ifelse(dat[, 4] == "Y", as.numeric(dat[, 3]), NA),
             row.names = dat[, 1])
}

read_phobius <- function(connection) {
  all_lines <- readLines(connection)
  all_lines <- all_lines[-1]
  splited <- strsplit(all_lines, " ")
  #remove "" characters
  purged <- t(sapply(splited, function(i) i[i != ""]))
  cl_sites <- sapply(which(purged[, 3] == "Y"), function(i)
    as.numeric(strsplit(strsplit(purged[i,4], "/")[[1]][1], "c")[[1]][[2]]))
  res <- data.frame(signal.peptide = purged[, 3] == "Y",
                    sig.start = ifelse(purged[, 3] == "Y", 1, NA),
                    sig.end = rep(NA, nrow(purged)), 
                    row.names = purged[, 1])
  res[purged[, 3] == "Y", "sig.end"] <- cl_sites
  res
}

read_philius <- function(connection) {
  require(XML)
  all_dat <- xmlToList(xmlTreeParse(connection, asTree = TRUE))
  seq_dat_id <- 1L:(length(all_dat)/2)*2
  #data for table
  table_dat <- sapply(seq_dat_id, function(i) 
    unlist(all_dat[i][[1]][[1]][c(24, 22)]))
  cleaved <- sapply(table_dat, function(i)
    !(is.null(i[1]) || is.na(i[1])))
  res <- data.frame(signal.peptide = cleaved,
                    sig.start = ifelse(cleaved, 1, NA),
                    sig.end = rep(NA, length(seq_dat_id)),
                    row.names = unlist(all_dat[1L:(length(all_dat)/2)*2 - 1]))
  res[cleaved, "sig.end"] <- as.numeric(sapply(table_dat[cleaved], function(i) i[2]))
  res
}