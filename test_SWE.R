rm(list = ls())
setwd("~/EnKF_GTE")
library(tidyverse)
library(foreach)
library(doParallel)
source('Code/test_SWE_single.R')

args = commandArgs(trailingOnly = TRUE)
mt = args[1]
b = as.numeric(args[2])
print(paste(mt, b))
test_SWE(mt, b)