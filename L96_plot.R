rm(list = ls())
setwd("E:/Project/EnKF_GTE")
library(tidyverse)
library(foreach)
library(doParallel)
library(ggpubr)


## fig 2 ####
result1 = read.csv('Output/L96_final_1125/Result/result-ev-p_231220(3).csv')
result2 = read.csv('Output/L96_final_1125/Result/result-all_1129.csv')
result3 = read.csv('Output/L96_final_1125/Result/result-GC-only_231220.csv')
# result = rbind(result1, result2 %>% filter(method %in% c('standard', 'oracle', 'inflation'))) %>% rbind(result3)
result = rbind(result1, result2 %>% filter(method %in% c('standard', 'inflation'))) %>% rbind(result3)

g = result %>%
  # filter(method %in% c('standard', 'oracle', 'inflation', 'GC-only', 'HD-EnKF(GC)')) %>%
  filter(method %in% c('standard', 'oracle', 'inflation', 'GC-only', 'HD-EnKF(GC)')) %>%
  mutate(method = {
    x = method
    x[which(x == 'GC-only')] = 'taper'
    x[which(x == 'HD-EnKF(GC)')] = 'HD-EnKF'
    x = factor(x, c('standard', 'inflation', 'taper', 'HD-EnKF', 'oracle'))
    x
  }) %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, color = method))+
  facet_wrap(vars(p, n), ncol = 3,
             labeller = label_bquote(p==.(p)~n==.(n)))+
  labs(x = 'Force term', y = 'RMSE')+
  # ggsci::scale_color_npg()+
  scale_color_manual(values = c('#8B0000', '#CC5500', '#006400', '#00008B', '#800080'))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
# ggsave('Output/L96_plot/Fig2_all_RMSE.png', g, width = 8, height = 8, dpi = 600)
# ggsave('Output/L96_plot/Fig2_all_RMSE.tiff', g, width = 8, height = 8, dpi = 600)
ggsave('Output/L96_plot/Slides_all_RMSE.pdf', g, width = 5, height = 5, dpi = 600)

## fig 3 ####
dat.plot = foreach(p = c(40, 100, 200), .combine = 'rbind') %do%
  {
    q = p
    n = 30
    f = 12
    S.seq = 4
    # foreach(mt = c('standard', 'oracle', 'inflation', 'GC-only', 'GC-ev-p'), .combine = 'rbind') %do%
    foreach(mt = c('standard', 'inflation', 'GC-only', 'GC-ev-p'), .combine = 'rbind') %do%
      {
        if(mt == 'GC-only'){
          rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
            {
              load(file = paste0('Output/L96_set26_taper-only/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/',
                                 mt, '_F', f, '_B', b, '.Rdata'))
              
              data.frame(
                boot = b,
                RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
              )
            }
          b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
          load(paste0('Output/L96_set26_taper-only/L96_p', p, '_q', q, '_n', n, '_seq4/', mt, '_F', f, '_B', b, '.Rdata'))
        }else if(mt == 'GC-ev-p'){
          rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
            {
              load(file = paste0('Output/L96_set25_taper-ev-p_scale/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/',
                                 mt, '_F', f, '_B', b, '.Rdata'))
              
              data.frame(
                boot = b,
                RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
              )
            }
          b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
          load(paste0('Output/L96_set25_taper-ev-p_scale/L96_p', p, '_q', q, '_n', n, '_seq4/', mt, '_F', f, '_B', b, '.Rdata'))
        }else{
          load(paste0('Output/L96_final_1125/L96_p', p, '_q', q, '_n', n, '_seq4/', mt, '_F', f, '.Rdata'))
        }
        if(mt == 'GC-only'){mt = 'taper'}
        if(mt == 'GC-ev-p'){mt = 'HD-EnKF'}
        data.frame(
          method = mt,
          p = p, n = n, F.set = f,
          x = 1:2000,
          y = zoo::rollmean(analyse$error, 20, fill = 'extend')
        )
      } 
  }
g = dat.plot %>% 
  # mutate(method = factor(method, c('standard', 'inflation', 'taper', 'HD-EnKF', 'oracle'))) %>%
  mutate(method = factor(method, c('standard', 'inflation', 'taper', 'HD-EnKF'))) %>%
  ggplot()+
  geom_line(aes(x = x, y = y, col = method))+
  facet_wrap(vars(p, n, F.set), ncol = 1, labeller = label_bquote(p==.(p)~~n==.(n)~~'F'==.(F.set)))+
  # labs(x = 'Time step', y = latex2exp::TeX('RMSE$_i$'))+
  labs(x = 'Time step', y = 'Assimilation error')+
  # ggsci::scale_color_npg()+
  ylim(c(0, 8.7))+
  scale_color_manual(values = c('#8B0000', '#CC5500', '#006400', '#00008B', '#800080'))+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
# ggsave('Output/L96_plot/Fig3_RMSE_F12.png', g, width = 6, height = 6, dpi = 600) 
# ggsave('Output/L96_plot/Fig3_RMSE_F12.tiff', g, width = 6, height = 6, dpi = 600) 
ggsave('Output/L96_plot/Slides_RMSE_F12.pdf', g, width = 4.5, height = 4.5, dpi = 600)

## taper-const ####
dat.plot = foreach(p = c(40, 100, 200), .combine = 'rbind') %do%
  {
    q = p
    n = 30
    S.seq = 4
    tmp1 = foreach(f = 4:12, .combine = 'rbind') %do%
      {
        mt = 'banding-ev-p'
        dir.open = paste0('Output/L96_set25_taper-ev-p_scale/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
        rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
          {
            load(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
            data.frame(
              boot = b,
              RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
            )
          }
        b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
        load(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
        data.frame(
          p = p, 
          n = n,
          F.set = f,
          boot = b,
          method = 'HD-EnKF',
          RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
        )
      }
    tmp2 = foreach(f = 4:12, .combine = 'rbind') %do%
      {
        mt = 'BL-constant'
        dir.open = paste0('Output/taper-const/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
        rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
          {
            load(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
            data.frame(
              boot = b,
              RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
            )
          }
        b = rmse.all %>% arrange(RMSE) %>% slice(4) %>% pull(boot)
        # b = 1
        load(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
        data.frame(
          p = p, 
          n = n,
          F.set = f,
          boot = b,
          method = 'HD-EnKF(constant)',
          RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
        )
      }
    rbind(tmp1, tmp2)
  } 
g1 = dat.plot %>%
  mutate(method = factor(method, c('HD-EnKF', 'HD-EnKF(constant)'))) %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, color = method, linetype = method))+
  facet_wrap(vars(p, n), ncol = 1,
             labeller = label_bquote(p==.(p)~n==.(n)))+
  # scale_color_manual(values = c( '#008080', '#FF8C00','#00008B', '#8B4500', '#800080', '#8B2300'))+
  scale_color_manual(values = c('#00008B', '#8B4500'))+
  labs(x = 'Force term', y = 'RMSE')+
  theme_bw()+
  theme(legend.position = 'bottom',
        legend.title = element_blank())
print(g1)

dat.plot = foreach(p = c(40, 100, 200), .combine = 'rbind') %do%
  {
    q = p
    n = 30
    S.seq = 4
    f = 12
    b = 1
    
    mt = 'BL'
    dir.open = paste0('Output/taper-const/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
    load(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
    tmp1 = data.frame(
      p = p, 
      n = n,
      method = 'HD-EnKF',
      time = 1:2000,
      RMSE = zoo::rollmean(analyse$error, 20, fill = 'extend')
    )
    
    mt = 'BL-constant'
    dir.open = paste0('Output/taper-const/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
    load(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
    tmp2 = data.frame(
      p = p, 
      n = n,
      method = 'HD-EnKF(constant)',
      time = 1:2000,
      RMSE = zoo::rollmean(analyse$error, 20, fill = 'extend')
    )
    
    rbind(tmp1, tmp2)
  } 
g2 = dat.plot %>%
  mutate(method = factor(method, c('HD-EnKF', 'HD-EnKF(constant)'))) %>%
  ggplot()+
  geom_line(aes(x = time, y = RMSE, color = method, linetype = method))+
  facet_wrap(vars(p, n), ncol = 1,
             labeller = label_bquote(p==.(p)~n==.(n)))+
  # scale_color_manual(values = c( '#008080', '#FF8C00','#00008B', '#8B4500', '#800080', '#8B2300'))+
  scale_color_manual(values = c('#00008B', '#8B4500'))+
  labs(x = 'Time step', y = 'Assimilation error')+
  theme_bw()+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(),
        legend.title = element_blank())
print(g2)

g = ggarrange(g1, g2, common.legend = TRUE, legend = 'bottom')
print(g)
# ggsave('Output/L96_plot/Fig4_taper_const.png', g, width = 8, height = 6, dpi = 600)
# ggsave('Output/L96_plot/Fig4_taper_const.tiff', g, width = 8, height = 6, dpi = 600)
ggsave('Output/L96_plot/Slides_taper_const.pdf', g, width = 7, height = 5.5, dpi = 600)

## sum-eigen ####
read.csv('Output/L96_final_1125/Result/result-all_231225.csv') %>%
  filter(p == 200 & n == 20) %>%
  select(method, F.set, RMSE) %>%
  spread(F.set, RMSE)

dir.open = 'Output/L96_set31_sumeigen/L96_p200_q200_n20_seq4/'
result = data.frame()
for(taper in c('BL', 'CZZ', 'GC')){
  for(seg in c(0.5, 0.9, 0.95, 0.99)){
    for(f in 4:12){
      b = 1
      load(file = paste0(dir.open, taper, '_sumeigen', seg, '_F', f, '_B', b, '.Rdata'))
      tmp = data.frame(
        taper = taper,
        sumeigen = seg,
        F.set = f,
        RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
      )
      result = rbind(result, tmp)
    }
  }
}
result %>%
  group_by(taper, sumeigen) %>%
  spread(F.set, RMSE)
dat1 = read.csv('Output/L96_final_1125/Result/result-all_231225.csv') %>%
  filter(p == 200 & n == 20 & method %in% c('standard', 'oracle', 'HD-EnKF(GC)')) %>%
  select(method, F.set, RMSE)
dat1$method[which(dat1$method == 'HD-EnKF(GC)')] = '100%'
dat2 = result %>%
  filter(taper == 'GC') %>%
  mutate(method = paste0('HD-EnKF(', sumeigen * 100, '%)')) %>%
  select(method, F.set, RMSE)
dat.plot = result %>%
  filter(taper == 'GC') %>%
  mutate(method = paste0(sumeigen * 100, '%')) %>%
  select(method, F.set, RMSE) %>%
  rbind(dat1) %>%
  mutate(method = factor(method, paste0(c(50, 90, 95, 99, 100), '%'))) %>%
  filter(method %in% paste0(c(50, 90, 95, 99, 100), '%'))
g = dat.plot %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, color = method))+
  scale_color_manual(name = 'proportion of the total  \n sum of eigenvalues',
                     values = c('#8B4500',  '#8B6B00', '#004500', '#5A4A8B', '#00008B'))+
  labs(x = 'Force term', y = 'RMSE')+
  theme_bw()+
  theme(legend.position = 'bottom')
print(g)
# ggsave(paste0('Output/L96_plot/FigA1_sumeigen.png'), g, width = 5.2, height = 4, dpi = 600)
# ggsave(paste0('Output/L96_plot/FigA1_sumeigen.tiff'), g, width = 5.2, height = 4, dpi = 600)
g = g + geom_line(aes(x = F.set, y = RMSE, color = method), linewidth = 0.5)
ggsave(paste0('Output/L96_plot/Slides_sumeigen.pdf'), g, width = 5.2, height = 4, dpi = 600)

## taper-block ####
dir.open = 'Output/L96_set32_block/L96_p200_q200_n20_seq4/'
result = data.frame()
for(taper in c('BL', 'CZZ', 'GC')){
  for(number.block in c(5)){
    for(overlap.block in c(1, 3)){
      for(f in 4:12){
        rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
          {
            load(file = paste0(dir.open, taper, '_block_a_', number.block, '_overlap', overlap.block, '_F', f, '_B', b, '.Rdata'))
            data.frame(
              boot = b,
              RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
            )
          }
        b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
        load(file = paste0(dir.open, taper, '_block_a_', number.block, '_overlap', overlap.block, '_F', f, '_B', b, '.Rdata'))
        tmp = data.frame(
          taper = taper,
          number.block = number.block,
          overlap.block = overlap.block,
          F.set = f,
          RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
        )
        result = rbind(result, tmp)
      }
    }
  }
}
result %>%
  group_by(taper, number.block, overlap.block) %>%
  spread(F.set, RMSE)
read.csv('Output/L96_final_1125/Result/result-all_231225.csv') %>%
  filter(p == 200 & n == 20) %>%
  select(method, F.set, RMSE) %>%
  spread(F.set, RMSE)
dat1 = read.csv('Output/L96_final_1125/Result/result-all_231225.csv') %>%
  filter(p == 200 & n == 20 & method %in% c('standard', 'oracle', 'HD-EnKF(GC)')) %>%
  select(method, F.set, RMSE)
dat1$method[which(dat1$method == 'HD-EnKF(GC)')] = 'HD-EnKF'
dat2 = result %>%
  filter(taper == 'GC' & number.block == 5 & overlap.block == 1) %>%
  mutate(method = 'HD-EnKF(block)') %>%
  select(method, F.set, RMSE)
dat.plot = rbind(dat1, dat2)  %>%
  mutate(method = factor(method, c('standard', 'HD-EnKF', 'HD-EnKF(block)', 'oracle'))) %>%
  filter(method %in% c('HD-EnKF', 'HD-EnKF(block)'))
g1 = dat.plot %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, color = method))+
  scale_color_manual(values = c('#00008B','#8B4500'))+
  labs(x = 'Force term', y = 'RMSE')+
  theme_bw()+
  theme(legend.position = 'bottom',
        legend.title = element_blank())
print(g1)

p = q = 200
n = 20
S.seq = 4
f = 12
number.block = 5
overlap.block = 1
dir.open = 'Output/L96_set32_block/L96_p200_q200_n20_seq4/'
rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
  {
    load(file = paste0(dir.open, taper, '_block_a_', number.block, '_overlap', overlap.block, '_F', f, '_B', b, '.Rdata'))
    data.frame(
      boot = b,
      RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
    )
  }
b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
load(file = paste0(dir.open, taper, '_block_a_', number.block, '_overlap', overlap.block, '_F', f, '_B', b, '.Rdata'))
tmp1 = data.frame(
  taper = taper,
  method = 'HD-EnKF(block)',
  time = 1:2000,
  RMSE = analyse$error
)

dir.open = 'Output/L96_set32_block/L96_p200_q200_n20_seq4/'
rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
  {
    load(file = paste0(dir.open, taper, '_block_a_', number.block, '_overlap', overlap.block, '_F', f, '_B', b, '.Rdata'))
    data.frame(
      boot = b,
      RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
    )
  }
b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
load(file = paste0(dir.open, taper, '_block_a_', number.block, '_overlap', overlap.block, '_F', f, '_B', b, '.Rdata'))
tmp1 = data.frame(
  method = 'HD-EnKF(blockwise)',
  time = 1:2000,
  RMSE = zoo::rollmean(analyse$error, 20, fill = 'extend')
)
mt = 'GC-ev-p'
dir.open = paste0('Output/L96_set25_taper-ev-p_scale/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
rmse.all = foreach(b = 1:10, .combine = 'rbind') %do%
  {
    load(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
    data.frame(
      boot = b,
      RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
    )
  }
b = rmse.all %>% arrange(RMSE) %>% slice(1) %>% pull(boot)
load(paste0(dir.open, mt, '_F', f, '_B', b, '.Rdata'))
tmp2 = data.frame(
  method = 'HD-EnKF',
  time = 1:2000,
  RMSE = zoo::rollmean(analyse$error, 20, fill = 'extend')
)
g2 = rbind(tmp1, tmp2) %>%
  mutate(method = factor(method, c('HD-EnKF', 'HD-EnKF(blockwise)'))) %>%
  ggplot()+
  geom_line(aes(x = time, y = RMSE, color = method))+
  scale_color_manual(values = c('#00008B', '#8B4500'))+
  labs(x = 'Time step', y = 'Assimilation error')+
  theme_bw()+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(),
        legend.title = element_blank())
print(g2)

g = ggarrange(g1, g2, common.legend = TRUE, legend = 'bottom')
print(g)
# ggsave('Output/L96_plot/FigA2_blockwise.png', g, width = 10, height = 4, dpi = 600)
# ggsave('Output/L96_plot/FigA2_blockwise.tiff', g, width = 10, height = 4, dpi = 600)
ggsave('Output/L96_plot/Slides_blockwise.pdf', g, width = 8, height = 3.2, dpi = 600)
