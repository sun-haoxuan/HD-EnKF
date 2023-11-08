rm(list = ls())
setwd("E:/Project/EnKF_GTE")
# setwd("~/EnKF_GTE")
library(tidyverse)
library(foreach)
library(doParallel)
source('Code/test_SWE_single.R')

args = commandArgs(trailingOnly = TRUE)
mt = args[1]
b = as.numeric(args[2])
print(paste(mt, b))

analyse = test_SWE('band-iterwithinfl', 2)
