rm(list = ls())
library(tidyverse)
library(foreach)
library(doParallel)

dir.root = '~/L96'
setwd(dir.root)
'./Code/Assimilation' %>% list.files(pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)
'./Code/Calculation' %>% list.files(pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)
'./Code/Module' %>% list.files(pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)

dir.save.xt = 'Input/TrueState/'
dir.save.dM = 'Input/distMat/'
dir.save.H = 'Input/ObsMat/'
dir.save.R = 'Input/ObsErr/'
dir.save.yo = 'Input/Observation/'

p = 40
q = 40
n = 30
f = 8
b = 1

load(paste0(dir.save.xt, 'L96_p', p, '.Rdata'))
load(paste0(dir.save.dM, 'L96_p', p, '.Rdata'))
load(paste0(dir.save.yo, 'L96_p', p, '_q', q, '_boot', b, '.Rdata'))
load(paste0(dir.save.R, 'L96_p', p, '_q', q, '_boot', b,  '.Rdata'))
load(paste0(dir.save.H, 'L96_p', p, '_q', q, '_boot', b,  '.Rdata'))

method = 'HD-EnKF'
source(paste0('Code/Assimilation/Methods/', method, '.R'))
DAinfo = get_DAmethod(method)

option = list(
  boot = b,
  p = p,
  q = q,
  S = 2000,
  S.seq = 4,
  F.true = 8,
  h = 0.05,
  sigma = 1,
  rho = 0.5,
  state = list(
    x.t = x.t,
    y.o = y.o,
    H = H,
    R = Rlist$R,
    R.inv = Rlist$R.inv,
    R.inv.root = Rlist$R.inv.root
  ),
  DA = list(
    n = n,
    F.set = f,
    sigma.in = 0.1,
    DAfun = DAfun,
    DAinfo = DAinfo
  )
)

analysis = L96_Analysis(option)
