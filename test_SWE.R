rm(list = ls())
setwd("~/EnKF_GTE")
library(tidyverse)
library(foreach)
library(doParallel)
source('Code/test_SWE_single.R')

args = commandArgs(trailingOnly = TRUE)
mt = args[1]
k.set = as.numeric(args[2])
b = as.numeric(args[3])
print(paste(mt, k.set, b))

test_SWE(mt, k.set, b)