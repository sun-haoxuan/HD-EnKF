source('Code/Functions/cal_HPHR.R')
source('Code/Functions/cal_RMSE.R')
source('Code/Functions/SWE_State.R')
source('Code/Functions/SWE_Forecast.R')
source('Code/Functions/SWE_Analyse.R')
source('Code/Functions/cal_taper_cov_SWE.R')

test_SWE = function(mt, k.set = 5e4, b = 0){
  opt = list(
    ## Generate true state
    nx = 50, # grid size of x
    ny = 31, # grid size of ys
    L = 5e+5,
    D = 3e+5,
    g = 9.8, # gravitational constant
    f = 1e-4,
    k = NA,
    dt = 30, # hardwired timestep
    h0 = 50,
    h1 = 5.5,
    h2 = 3.325,
    ## Generate observations
    id = sort(sample(50, 10)),
    sigma.obs = c(sqrt(0.5), sqrt(0.5), 1), # observation error
    rho = 0.5, ## observation
    vf = 1, # for Wrong R
    ## Assimilation Step
    S.o.seq = 360, # assimilate every 3h
    S.o.start = 60, ## Start of the observation
    S.o.end = 5820, # assimilate 48h, then forecast 24h
    S = 8700,
    # S.o.end = 11580, # assimilate 48h, then forecast 24h
    # S = 14460,
    ## Assimilation
    n = 100,
    k.true = 1e+4,
    k.set = k.set,
    sigma.i = 1,
    method = mt,
    boot = b
  )
  
  set.seed(b)
  state = SWE_State(opt)
  analyse = SWE_Analyse(state, opt)
  
  save(analyse, file = paste0('Output/Result_SWE_max10/', mt, '_K', k.set, '_B', b, '.Rdata'))
  
  print(
    data.frame(
      method = mt,
      boot = b,
      cv = analyse$converge,
      RMSE_pre = sqrt(mean(analyse$error[1, 2941:5820] ^ 2)),
      RMSE_post = sqrt(mean(analyse$error[1, 5821:8700] ^ 2))
      # RMSE_d1 = sqrt(mean(analyse$error[1, 1 * 2880 - 2819:60] ^ 2)),
      # RMSE_d2 = sqrt(mean(analyse$error[1, 2 * 2880 - 2819:60] ^ 2)),
      # RMSE_d3 = sqrt(mean(analyse$error[1, 3 * 2880 - 2819:60] ^ 2)),
      # RMSE_d4 = sqrt(mean(analyse$error[1, 4 * 2880 - 2819:60] ^ 2)),
      # RMSE_d5 = sqrt(mean(analyse$error[1, 5 * 2880 - 2819:60] ^ 2))
      # RMSE_u_pre = sqrt(mean(analyse$error[2, 2941:5820] ^ 2)),
      # RMSE_u_post = sqrt(mean(analyse$error[2, 5821:8700] ^ 2)),
      # RMSE_v_pre = sqrt(mean(analyse$error[3, 2941:5820] ^ 2)),
      # RMSE_v_post = sqrt(mean(analyse$error[3, 5821:8700] ^ 2)),
      # RMSE_h_pre = sqrt(mean(analyse$error[4, 2941:5820] ^ 2)),
      # RMSE_h_post = sqrt(mean(analyse$error[4, 5821:8700] ^ 2))
    )
  )
  
}