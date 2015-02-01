pos_seqs <- read_uniprot(paste0(pathway, "sept_signal.txt"), euk = TRUE)

#cleavage site sequence
css <- tolower(t(sapply(pos_seqs, function(seq) {
  cleave_site <- attr(seq, "sig")[2]
  seq[(cleave_site-4):(cleave_site+3)]
})))

#cleavage site length
csl <- sapply(pos_seqs, function(seq) attr(seq, "sig")[2])

#csdf <- data.frame(apply(css, 2, l2n, seq_type = "prot"), csl)

csdf <- data.frame(names(csl), 
                   do.call(cbind, lapply(1L:ncol(css), function(i)
                     data.frame(factor(css[, i], level = tolower(a())[-1])))), 
                   csl,
                   cut(csl, c(6, 15, 24, 35, 89)))

colnames(csdf) <- c("prot", paste0("P", c(paste0(".", 4:1), 1L:4)), "position", "positionf")

ggplot(csdf, aes(x = position)) + geom_density() +
  scale_x_continuous("Cleavage site position") +
  geom_vline(xintercept = 15, colour = "red") +
  geom_vline(xintercept = 34, colour = "red")
  




summary(csdf[["position"]])
table(csdf[["positionf"]])


ggplot(csdf, aes(x = P.1, fill = positionf)) + geom_bar(position = "dodge")
