library(seqinr)
library(hmeasure)
library(pbapply)
library(cvTools)
library(parallel)
library(snow)
library(biogram)
library(signal.hsmm)

# working directory ------------------------------

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"


