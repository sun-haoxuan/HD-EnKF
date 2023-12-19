rm(list = ls())
setwd("~/EnKF_GTE/Output")
library(tidyverse)
library(foreach)
library(doParallel)

resultall = function(dir.open, p, q, n, S.seq, mt, f, B){
  for(b in 1:B){
    load(paste0(dir.open, 'L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/', mt, '_F', f, '_B', b, '.Rdata'))
    print(paste('B', b, 'RMSE', round(sqrt(mean(analyse$error[1001:2000] ^ 2)), 3),
                'obj', round(sqrt(mean(analyse$likelihood[251:500] ^ 2)), 3)))
  }
}

for(p in c(40, 100, 200)){
# for(p in c(500)){
  for(n in c(20, 30, 40)){
    q = p; S.seq = 4
    dir.save = paste0('L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
    dir.create(dir.save)
    
    for(mt in c('standard', 'oracle', 'inflation')){
      dir.open = paste0('L96_test/L96_param5/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
      for(f in 4:12){
        rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
          {
            load(file = paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))

            data.frame(
              boot = b,
              RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
            )
          }
        b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
        file.copy(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'),
                  paste0(dir.save, mt, '_F', f, '.Rdata'))
      }
    }

    for(mt in c('banding', 'tapering', 'GC')){
      if(n %in% c(20, 30)){
        dir.open = paste0('L96_final_1121/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
        for(f in 4:12){
          file.copy(paste0(dir.open, mt, '_F', f, '.Rdata'),
                    paste0(dir.save, mt, '_F', f, '.Rdata'))
        }
      }else{
        dir.open = paste0('L96_test/L96_param5/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
        for(f in 4:12){
          rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
            {
              load(file = paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))

              data.frame(
                boot = b,
                RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
              )
            }
          b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
          file.copy(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'),
                    paste0(dir.save, mt, '_F', f, '.Rdata'))
        }
      }
    }

    for(mt in c('banding-inflation', 'tapering-inflation', 'GC-inflation')){
      dir.open = paste0('L96_test/L96_set11_x5_combinemethods/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
      for(f in 4:12){
        rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
          {
            load(file = paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))

            data.frame(
              boot = b,
              RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
            )
          }
        b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
        file.copy(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'),
                  paste0(dir.save, mt, '_F', f, '.Rdata'))
      }
    }
    
    # for(mt in c('standard', 'oracle', 'inflation', 'banding-inflation', 'tapering-inflation', 'GC-inflation')){
    #   dir.open = paste0('L96_test/L96_set12_x5_p500/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
    #   for(f in 4:12){
    #     rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
    #       {
    #         load(file = paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
    #         
    #         data.frame(
    #           boot = b,
    #           RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
    #         )
    #       }
    #     b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
    #     file.copy(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'),
    #               paste0(dir.save, mt, '_F', f, '.Rdata'))
    #   }
    # }
    
  }
}

file.remove('Output/L96_final/L96_p200_q200_n30_seq4/standard_F12.Rdata')
file.copy('Output/L96_test/L96_set13_p200_n30/L96_p200_q200_n30_seq4/standard_F12_B11.Rdata',
          'Output/L96_final_1125/L96_p200_q200_n30_seq4/standard_F12.Rdata')

## Summary ####
methods = c('standard', 'oracle', 'inflation','banding', 'tapering', 'GC')
param = data.frame()
for(p in c(40, 100, 200, 500)){
  for(n in c(20, 30, 40)){
    for(mt in methods){
      for(f in 4:12){
        param = rbind(param, data.frame(
          p = p, n = n, mt = mt, f = f
        ))
      }
    }
  }
}

result = foreach(i = 1:nrow(param), .combine = 'rbind') %do%
  {
    mt = param$mt[i]
    f = param$f[i]
    p = q = param$p[i]
    n = param$n[i]
    S.seq = 4
    
    if(mt %in% c('banding', 'tapering', 'GC')){
      load(file = paste0('L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/',
                         mt, '-inflation_F', f, '.Rdata'))
    }else{
      load(file = paste0('L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/',
                        mt, '_F', f, '.Rdata'))
    }

    data.frame(
      p = p, n = n,
      method = mt,
      F.set = f,
      RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2)),
      obj = mean(analyse$likelihood[251:500]),
      iteration = mean(analyse$iteration[251:500]),
      bandwidth = mean(analyse$bandwidth[251:500])
    )
  }
result$method[which(result$method == 'banding')] = 'HD-EnKF(BL)'
result$method[which(result$method == 'tapering')] = 'HD-EnKF(CZZ)'
result$method[which(result$method == 'GC')] = 'HD-EnKF(GC)'
write.csv(result, 'L96_final_1125/Result/result-all_1218.csv', row.names = F)

result = read.csv('L96_final_1125/Result/result-all_1129.csv')
# result$RMSE[which(is.na(result$RMSE))] = Inf
g = result %>%
  mutate(method = factor(method, c('standard', 'oracle', 'inflation', 
                                   'HD-EnKF(BL)', 'HD-EnKF(CZZ)', 'HD-EnKF(GC)'))) %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, color = method))+
  facet_wrap(vars(p, n), ncol = 3,
             labeller = label_bquote(p==.(p)~n==.(n)))+
  labs(x = 'Force term', y = 'RMSE')+
  # ggsci::scale_color_npg()+
  scale_color_manual(values = c('#8B0000', '#CC5500', '#006400', '#008080', '#00008B', '#800080'))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
print(g)
ggsave('L96_final_1125/Result/result-all_RMSE_3.png', g, width = 9, height = 10, dpi = 600)

## obj ####
result = read.csv('L96_final_1125/Result/result-all_1126.csv')
result$method[which(result$method == 'banding')] = 'HD-EnKF(BL)'
result$method[which(result$method == 'tapering')] = 'HD-EnKF(CZZ)'
result$method[which(result$method == 'GC')] = 'HD-EnKF(GC)'
foreach(i = c(40, 100, 200, 500), .combine = 'rbind') %do%
  {
    tab1 = result %>%
      filter(p == i & n == 20 & F.set == 12) %>%
      dplyr::select(method, RMSE, obj, bandwidth)
    tab2 = result %>%
      filter(p == i & n == 30 & F.set == 12) %>%
      dplyr::select(method, RMSE, obj, bandwidth)
    tab3 = result %>%
      filter(p == i & n == 40 & F.set == 12) %>%
      dplyr::select(method, RMSE, obj, bandwidth)
    tab1[, 1] %>%
      cbind(data.frame(NA)) %>%
      cbind(round(tab1[, 2:3], 2)) %>%
      cbind(round(tab1[, 4], 1)) %>%
      cbind(data.frame(NA)) %>%
      cbind(round(tab2[, 2:3], 2)) %>%
      cbind(round(tab2[, 4], 1)) %>%
      cbind(data.frame(NA)) %>%
      cbind(round(tab3[, 2:3], 2)) %>%
      cbind(round(tab3[, 4], 1))
  } %>%
  xtable::xtable() %>%
  print(include.rownames = F)


## bandwidth ####
p = q = 40
n = 30
f = 12
g = foreach(mt = c('banding', 'tapering', 'GC'), .combine = 'rbind') %do%
  {
    load(paste0('L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq4/', mt, '-inflation_F', f, '.Rdata'))
    if(mt == 'banding'){
      mt = 'HD-EnKF(BL)'; y = floor(analyse$bandwidth); bw = 10
    }else if(mt == 'tapering'){
      mt = 'HD-EnKF(CZZ)'; y = round(analyse$bandwidth, 1); bw = 10
    }else if(mt == 'GC'){
      mt = 'HD-EnKF(GC)'; y = analyse$bandwidth; bw = 10
    }
    data.frame(
      method = mt,
      x = 4 * 1:500,
      y = zoo::rollmean(analyse$bandwidth, bw, fill = 'extend')
    )
  } %>% ggplot()+
  geom_line(aes(x = x, y = y))+
  facet_wrap(vars(method), ncol = 1, scale = 'free_y')+
  labs(x = 'Time step', y = 'Selected bandwidth')+
  theme_bw()
ggsave('L96_final_1125/Result/result_bandwidth.png', g, width = 8, height = 4)  
  
g = foreach(mt = c('standard', 'oracle', 'inflation', 'banding', 'tapering', 'GC'), .combine = 'rbind') %do%
  {
    if(mt %in% c('banding', 'tapering', 'GC')){
      mt = paste0(mt, '-inflation')
    }
    load(paste0('L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq4/', mt, '_F', f, '.Rdata'))
    if(mt == 'banding-inflation'){
      mt = 'HD-EnKF(BL)'
    }else if(mt == 'tapering-inflation'){
      mt = 'HD-EnKF(CZZ)'
    }else if(mt == 'GC-inflation'){
      mt = 'HD-EnKF(GC)'
    }
    data.frame(
      method = mt,
      x = 1:2000,
      y = zoo::rollmean(analyse$error, 20, fill = 'extend')
    )
  } %>% 
  mutate(method = factor(method, c('standard', 'oracle', 'inflation', 'HD-EnKF(BL)', 'HD-EnKF(CZZ)', 'HD-EnKF(GC)'))) %>%
  ggplot()+
  geom_line(aes(x = x, y = y, col = method))+
  labs(x = 'Time step', y = latex2exp::TeX('RMSE$_i$'))+
  # ggsci::scale_color_npg()+
  scale_color_manual(values = c('#8B0000', '#CC5500', '#006400', '#008080', '#00008B', '#800080'))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
ggsave('L96_final_1125/Result/result_RMSE_t.png', g, width = 5, height = 5) 

## bandwidth table ####
result.bw = data.frame()
for(p in c(40, 100, 200)){
  for(n in c(20, 30, 40)){
    q = p; S.seq = 4
    for(f in 4:12){
      tmp = foreach(mt = c('banding', 'tapering', 'GC'), .combine = 'rbind') %do%
        {
          load(paste0('L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq4/',
                      mt, '-inflation_F', f, '.Rdata'))
          if(mt == 'banding'){
            mt = 'HD-EnKF(BL)'; y = analyse$bandwidth
          }else if(mt == 'tapering'){
            mt = 'HD-EnKF(CZZ)'; y = analyse$bandwidth
          }else if(mt == 'GC'){
            mt = 'HD-EnKF(GC)'; y = analyse$bandwidth
          }
          # if(mt == 'banding'){
          #   mt = 'HD-EnKF(BL)'; y = floor(analyse$bandwidth); bw = 10
          # }else if(mt == 'tapering'){
          #   mt = 'HD-EnKF(CZZ)'; y = round(analyse$bandwidth, 1); bw = 10
          # }else if(mt == 'GC'){
          #   mt = 'HD-EnKF(GC)'; y = analyse$bandwidth; bw = 10
          # }
          data.frame(
            p = p, n = n, F.set = f,
            method = mt,
            x = 4 * 1:500,
            y = analyse$bandwidth
            # y = zoo::rollmean(analyse$bandwidth, bw, fill = 'extend')
          )
        }
      result.bw = rbind(result.bw, tmp)
    }
  }
}

result.bw %>%
  filter(x > 1000) %>%
  group_by(p, n, F.set, method) %>%
  summarise(bw = paste0(round(mean(y), 1), ' (',
                        round(sd(y), 2), ')')) %>%
  # summarise(bw = paste0(round(mean(y), 1), ' (',
  #                       round(min(y), 1), ',',
  #                       round(max(y), 1), ')')) %>%
  ungroup() %>%
  filter(F.set == 12) %>% 
  dplyr::select(-F.set) %>% 
  spread(method, bw) %>%
  xtable::xtable(digits = 0) %>%
  print(include.rownames = F)

## taper only ####
setwd("~/EnKF_GTE/Output")
dir.create('L96_taper-only')
for(p in c(40, 100, 200)){
  for(n in c(20, 30, 40)){
    q = p; S.seq = 4
    dir.save = paste0('Output/L96_taper-only/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
    dir.create(dir.save)
    
    for(mt in c('banding', 'tapering', 'GC')){
      for(f in 4:12){
        dir.open = paste0('L96_test/L96_set15_taperonly/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
        if(f != 8){
          dir.open = paste0('L96_test/L96_set15_taperonly/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
        }else{
          dir.open = paste0('L96_test/L96_set16_taperonly_param2/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
        }
        rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
          {
            load(file = paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
            
            data.frame(
              boot = b,
              RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
            )
          }
        b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
        file.copy(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'),
                  paste0(dir.save, mt, '_F', f, '.Rdata'))
      }
    }
  }
}

## adjust on F = 8
read.csv("Output/L96_final_1125/Result/result-all_1129.csv") %>% filter(p == 40 & n == 30 & F.set == 8)
p = q = 40; n = 30; S.seq = 4
dir.open = paste0('Output/L96_test/L96_set16_taperonly_param2/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
dir.save = paste0('Output/L96_taper-only/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
resultall('Output/L96_test/L96_set16_taperonly_param2/',  40, 40, 30, 4, 'GC', 8, 10)
file.remove(paste0(dir.save, 'banding_F8.Rdata'))
file.copy(paste0(dir.open, 'banding_F8_B1.Rdata'), paste0(dir.save, 'banding_F8.Rdata'))
file.remove(paste0(dir.save, 'tapering_F8.Rdata'))
file.copy(paste0(dir.open, 'tapering_F8_B9.Rdata'), paste0(dir.save, 'tapering_F8.Rdata'))
file.remove(paste0(dir.save, 'GC_F8.Rdata'))
file.copy(paste0(dir.open, 'GC_F8_B7.Rdata'), paste0(dir.save, 'GC_F8.Rdata'))

## plot
p = q = 40; n = 30; S.seq = 4
dir1 = paste0('Output/L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
dir2 = paste0('Output/L96_taper-only/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
result = data.frame()
for(mt in c('banding', 'tapering', 'GC')){
  if(mt == 'banding'){taper = 'HD-EnKF(BL)'}
  if(mt == 'tapering'){taper = 'HD-EnKF(CZZ)'}
  if(mt == 'GC'){taper = 'HD-EnKF(GC)'}
  
  for(f in 4:12){
    load(paste0(dir1, mt, '-inflation_F', f, '.Rdata'))
    tmp = data.frame(
      taper = taper,
      method = 'tapering + inflation + iteration',
      F.set = f,
      RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
    )
    result = rbind(result, tmp)
    
    load(paste0(dir2, mt, '_F', f, '.Rdata'))
    tmp = data.frame(
      taper = taper,
      method = 'tapering only',
      F.set = f,
      RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
    )
    result = rbind(result, tmp)
  }
}
g = result %>%
  filter(taper == 'HD-EnKF(GC)') %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, linetype = method))+
  # facet_wrap(vars(taper), ncol = 3)+
  labs(x = 'Force term', y = 'RMSE')+
  ylim(c(0, 5.5))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
print(g)
ggsave('Output/L96_taper-only/result-taper-only_2.png', g, width = 4, height = 3.5, dpi = 600)

## eigen adjust ####
setwd("~/EnKF_GTE/Output")
dir.create('L96_eigen_adjust')
p = q = 100
n = 30
S.seq = 4

for(ev in c(30, 60, 90)){
  dir.save = paste0('L96_eigen_adjust/ev', ev, '/')
  dir.create(dir.save)
  for(mt in c('banding', 'tapering', 'GC')){
    for(f in 4:12){
      dir.open = paste0('L96_test/L96_set18_ev_num/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')

      rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
        {
          load(paste0(dir.open, mt, '-ev', ev, '_F', f, '_B', b, '.Rdata'))
          data.frame(
            boot = b,
            RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
          )
        }
      b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
      file.copy(paste0(dir.open, mt, '-ev', ev, '_F', f, '_B', b, '.Rdata'),
                paste0(dir.save, mt, '_F', f, '.Rdata'))
    }
  }
}

## plot
result = data.frame()
for(ev in c(30, 60, 90)){
  dir.save = paste0('L96_eigen_adjust/ev', ev, '/')
  for(mt in c('banding', 'tapering', 'GC')){
    if(mt == 'banding'){taper = 'HD-EnKF(BL)'}
    if(mt == 'tapering'){taper = 'HD-EnKF(CZZ)'}
    if(mt == 'GC'){taper = 'HD-EnKF(GC)'}
    for(f in 4:12){
      load(paste0(dir.save, mt, '_F', f, '.Rdata'))
      tmp = data.frame(
        taper = taper,
        ev = ev,
        F.set = f,
        RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
      )
      result = rbind(result, tmp)
    }
  }
}
result$ev = factor(result$ev)
g = result %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, linetype = ev))+
  facet_wrap(vars(taper), ncol = 3)+
  labs(x = 'Force term', y = 'RMSE')+
  # ylim(c(0, 5.5))+
  scale_linetype_discrete(name = 'number of eigenvalues \n and eigenvectors')+
  theme_bw()+
  theme(legend.position = 'bottom',
        legend.title = element_text(hjust = 0.5))
print(g)
ggsave('L96_eigen_adjust/result-eigen-adj.png', g, width = 8, height = 3, dpi = 600)

## smoothBW ####
for(p in c(40, 100, 200)){
  q = p
  n = 30
  S.seq = 4
  dir.open = paste0('L96_test/L96_set20_smoothBW/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
  dir.save = paste0('L96_smoothBW/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
  dir.create(dir.save)
  for(mt in c('banding-smoothBW', 'tapering-smoothBW', 'GC-smoothBW')){
    for(f in 4:12){
      rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
        {
          load(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
          data.frame(
            boot = b,
            RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
          )
        }
      b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
      file.remove(paste0(dir.save, mt, '_F', f, '.Rdata'))
      file.copy(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'),
                paste0(dir.save, mt, '_F', f, '.Rdata'))
    }
  }
}

## plot
result = data.frame()
for(p in c(40, 100, 200)){
  q = p
  n = 30
  S.seq = 4
  dir.open1 = paste0('Output/L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
  dir.open2 = paste0('Output/L96_smoothBW/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
  # for(mt in c('banding', 'tapering', 'GC')){
  #   if(mt == 'banding'){method1 = 'HD-EnKF(BL)'; method2 = 'HD-EnKF(BL, constant)'}
  #   if(mt == 'tapering'){method1 = 'HD-EnKF(CZZ)'; method2 = 'HD-EnKF(CZZ, constant)'}
  #   if(mt == 'GC'){method1 = 'HD-EnKF(GC)'; method2 = 'HD-EnKF(GC, constant)'}
  for(mt in c('GC')){
    if(mt == 'GC'){method1 = 'HD-EnKF'; method2 = 'HD-EnKF(constant)'}
    for(f in 4:12){
      load(paste0(dir.open1, mt, '-inflation_F', f, '.Rdata'))
      tmp = data.frame(
        p = p, n = n,
        method = method1,
        F.set = f,
        RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
      )
      result = rbind(result, tmp)
      
      load(paste0(dir.open2, mt, '-smoothBW_F', f, '.Rdata'))
      tmp = data.frame(
        p = p, n = n,
        method = method2,
        F.set = f,
        RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
      )
      result = rbind(result, tmp)
    }
  }
}
g = result %>%
  # mutate(method = factor(method, c('HD-EnKF(BL)', 'HD-EnKF(BL, constant)',
  #                                  'HD-EnKF(CZZ)', 'HD-EnKF(CZZ, constant)',
  #                                  'HD-EnKF(GC)', 'HD-EnKF(GC, constant)'))) %>%
  mutate(method = factor(method, c('HD-EnKF', 'HD-EnKF(constant)'))) %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, color = method))+
  facet_wrap(vars(p, n), ncol = 3,
             labeller = label_bquote(p==.(p)~n==.(n)))+
  # scale_color_manual(values = c( '#008080', '#FF8C00','#00008B', '#8B4500', '#800080', '#8B2300'))+
  scale_color_manual(values = c('#00008B', '#8B4500'))+
  labs(x = 'Force term', y = 'RMSE')+
  theme_bw()+
  theme(legend.position = 'bottom',
        legend.title = element_blank())
print(g)
ggsave('Output/L96_smoothBW/result-smoothBW_2.png', g, width = 8, height = 3.5, dpi = 600)

## All ev40 ####
methods = c('banding-ev40', 'tapering-ev40', 'GC-ev40')
param = data.frame()
for(p in c(40, 100, 200)){
  for(n in c(20, 30, 40)){
    for(mt in methods){
      for(f in 4:12){
        param = rbind(param, data.frame(
          p = p, n = n, mt = mt, f = f
        ))
      }
    }
  }
}

result = foreach(i = 1:nrow(param), .combine = 'rbind') %do%
  {
    mt = param$mt[i]
    f = param$f[i]
    p = q = param$p[i]
    n = param$n[i]
    S.seq = 4
    
    rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
      {
        load(file = paste0('Output/L96_test/L96_set22_ev40/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/',
                           mt, '_F', f, '_B', b, '.Rdata'))
        
        data.frame(
          boot = b,
          RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
        )
      }
    b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
    
    
    
    data.frame(
      p = p, n = n,
      method = mt,
      F.set = f,
      RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2)),
      obj = mean(analyse$likelihood[251:500]),
      iteration = mean(analyse$iteration[251:500]),
      bandwidth = mean(analyse$bandwidth[251:500])
    )
  }
result$method[which(result$method == 'banding')] = 'HD-EnKF(BL)'
result$method[which(result$method == 'tapering')] = 'HD-EnKF(CZZ)'
result$method[which(result$method == 'GC')] = 'HD-EnKF(GC)'

## All ev-p ####
methods = c('banding-ev-p', 'tapering-ev-p', 'GC-ev-p')
param = data.frame()
for(p in c(40, 100, 200)){
  for(n in c(20, 30, 40)){
    for(mt in methods){
      for(f in 4:12){
        param = rbind(param, data.frame(
          p = p, n = n, mt = mt, f = f
        ))
      }
    }
  }
}

result = foreach(i = 1:nrow(param), .combine = 'rbind') %do%
  {
    mt = param$mt[i]
    f = param$f[i]
    p = q = param$p[i]
    n = param$n[i]
    S.seq = 4
    
    rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
      {
        load(file = paste0('Output/L96_test/L96_set17_ev-p_param5/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/',
                           mt, '_F', f, '_B', b, '.Rdata'))
        
        data.frame(
          boot = b,
          RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
        )
      }
    b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
    
    file.copy(paste0('Output/L96_test/L96_set17_ev-p_param5/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/',
                     mt, '_F', f, '_B', b, '.Rdata'),
              paste0('Output/L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/',
                     mt, '_F', f, '.Rdata'))
    
    load(file = paste0('Output/L96_test/L96_set17_ev-p_param5/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/',
                       mt, '_F', f, '_B', b, '.Rdata'))
    data.frame(
      p = p, n = n,
      method = mt,
      F.set = f,
      RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2)),
      obj = mean(analyse$likelihood[251:500]),
      iteration = mean(analyse$iteration[251:500]),
      bandwidth = mean(analyse$bandwidth[251:500])
    )
  }
result$method[which(result$method == 'banding-ev-p')] = 'HD-EnKF(BL)'
result$method[which(result$method == 'tapering-ev-p')] = 'HD-EnKF(CZZ)'
result$method[which(result$method == 'GC-ev-p')] = 'HD-EnKF(GC)'
write.csv(result, 'Output/L96_final_1125/Result/result-ev-p_231219.csv', row.names = F)

result1 = read.csv('Output/L96_final_1125/Result/result-ev-p_231219.csv')
result2 = read.csv('Output/L96_final_1125/Result/result-all_1129.csv')
result = rbind(result1, result2 %>% filter(method %in% c('standard', 'oracle', 'inflation')))

g = result %>%
  filter(method %in% c('standard', 'oracle', 'inflation','HD-EnKF(GC)')) %>%
  mutate(method = {
    x = method
    x[which(x == 'HD-EnKF(GC)')] = 'HD-EnKF'
    x = factor(x, c('standard', 'oracle', 'inflation', 'HD-EnKF'))
    x
  }) %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, color = method))+
  facet_wrap(vars(p, n), ncol = 3,
             labeller = label_bquote(p==.(p)~n==.(n)))+
  labs(x = 'Force term', y = 'RMSE')+
  # ggsci::scale_color_npg()+
  scale_color_manual(values = c('#8B0000', '#CC5500', '#006400', '#00008B'))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
print(g)
ggsave('Output/L96_final_1125/Result/result-all_RMSE_4.png', g, width = 9, height = 10, dpi = 600)

result$method = factor(result$method, c('standard', 'oracle', 'inflation', 'HD-EnKF(BL)', 'HD-EnKF(CZZ)', 'HD-EnKF(GC)'))
foreach(i = c(40, 100, 200), .combine = 'rbind') %do%
  {
    tab1 = result %>%
      filter(p == i & n == 20 & F.set == 12) %>%
      dplyr::select(method, RMSE, obj, bandwidth) %>%
      arrange(method)
    tab2 = result %>%
      filter(p == i & n == 30 & F.set == 12) %>%
      dplyr::select(method, RMSE, obj, bandwidth) %>%
      arrange(method)
    tab3 = result %>%
      filter(p == i & n == 40 & F.set == 12) %>%
      dplyr::select(method, RMSE, obj, bandwidth) %>%
      arrange(method)
    tab1[, 1] %>%
      cbind(data.frame(NA)) %>%
      cbind(round(tab1[, 2:3], 2)) %>%
      cbind(round(tab1[, 4], 1)) %>%
      cbind(data.frame(NA)) %>%
      cbind(round(tab2[, 2:3], 2)) %>%
      cbind(round(tab2[, 4], 1)) %>%
      cbind(data.frame(NA)) %>%
      cbind(round(tab3[, 2:3], 2)) %>%
      cbind(round(tab3[, 4], 1))
  } %>%
  xtable::xtable() %>%
  print(include.rownames = F)

p = q = 40
n = 30
f = 12
g = foreach(mt = c('banding', 'tapering', 'GC'), .combine = 'rbind') %do%
  {
    load(paste0('Output/L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq4/', mt, '-ev-p_F', f, '.Rdata'))
    if(mt == 'banding'){
      mt = 'HD-EnKF(BL)'; y = floor(analyse$bandwidth); bw = 10
    }else if(mt == 'tapering'){
      mt = 'HD-EnKF(CZZ)'; y = round(analyse$bandwidth, 1); bw = 10
    }else if(mt == 'GC'){
      mt = 'HD-EnKF(GC)'; y = analyse$bandwidth; bw = 20
    }
    data.frame(
      method = mt,
      x = 4 * 1:500,
      y = zoo::rollmean(analyse$bandwidth, bw, fill = 'extend')
    )
  } %>% ggplot()+
  geom_line(aes(x = x, y = y))+
  facet_wrap(vars(method), ncol = 1, scale = 'free_y')+
  labs(x = 'Time step', y = 'Selected bandwidth')+
  theme_bw()
print(g)
ggsave('Output/L96_final_1125/Result/result_bandwidth_2.png', g, width = 8, height = 4)  

g = foreach(mt = c('standard', 'oracle', 'inflation', 'GC-ev-p'), .combine = 'rbind') %do%
  {
    load(paste0('Output/L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq4/', mt, '_F', f, '.Rdata'))
    if(mt == 'GC-ev-p'){mt = 'HD-EnKF'}
    data.frame(
      method = mt,
      x = 1:2000,
      y = zoo::rollmean(analyse$error, 20, fill = 'extend')
    )
  } %>% 
  mutate(method = factor(method, c('standard', 'oracle', 'inflation', 'HD-EnKF'))) %>%
  ggplot()+
  geom_line(aes(x = x, y = y, col = method))+
  labs(x = 'Time step', y = latex2exp::TeX('RMSE$_i$'))+
  # ggsci::scale_color_npg()+
  scale_color_manual(values = c('#8B0000', '#CC5500', '#006400', '#00008B'))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
print(g)
ggsave('Output/L96_final_1125/Result/result_RMSE_t_2.png', g, width = 5, height = 5) 
