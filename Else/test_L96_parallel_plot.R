write.csv(result, paste0(dir.save, 'result_all.csv'), row.names = F)

methods = c('oracle', 'Wang & Bishop', 'inflation', 'inflation & iterative', 'banding', 'banding & inflation', 'banding & iteration', 'banding & inflation & iteration')
result = read.csv(paste0(dir.save, 'result_all.csv'))

g = result %>%
  group_by(method, F.set) %>%
  summarise(RMSE = quantile(RMSE, 0.5, na.rm = T)) %>%
  # summarise(RMSE = min(RMSE, na.rm = T)) %>%
  mutate(method = factor(method, levels = methods)) %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, color = method))+
  labs(x = 'Force', y = 'RMSE')+
  theme_bw()+
  theme(legend.position = 'bottom',
        legend.title = element_blank())
print(g)
ggsave('Output/Simulation_L96_RMSE_1007.png', g, width = 8, height = 6, dpi = 600)


A = result %>%
  filter(F.set %in% c(4, 8, 12)) %>%
  group_by(method, F.set) %>%
  arrange(RMSE, .by_group = T) %>%
  slice(2)
# summarise(boot = which.min(na.omit(RMSE))) 
dat.RMSE = foreach(mt = A$method, f = A$F.set, b = A$boot, .combine = 'rbind') %do% 
  {
    filename = paste0(mt, '_F', f, '_B', b, '.Rdata')
    load(paste0(dir.save, mt, '/', filename))
    
    data.frame(
      method = mt,
      F.set = f,
      boot = b,
      X = 1:2000,
      RMSE = analyse$error
    )
  } 
dat1.RMSE = dat %>%
  filter(X > 1000) %>%
  mutate(method = factor(method, levels = methods)) %>%
  group_by(method, F.set) %>%
  # summarise(RMSE = sqrt(mean(RMSE ^ 2))) %>%
  summarise(RMSE = mean(RMSE)) %>%
  spread(F.set, RMSE)
dat.obj = foreach(mt = A$method, f = A$F.set, b = A$boot, .combine = 'rbind') %do% 
  {
    filename = paste0(mt, '_F', f, '_B', b, '.Rdata')
    load(paste0(dir.save, mt, '/', filename))
    
    data.frame(
      method = mt,
      F.set = f,
      boot = b,
      X = 1:500,
      obj = analyse$likelihood
    )
  } 
dat2 = dat.obj %>%
  filter(X > 250) %>%
  mutate(method = factor(method, levels = methods)) %>%
  group_by(method, F.set) %>%
  summarise(obj = mean(obj)) %>%
  spread(F.set, obj)
dat1[, 1] %>% cbind(data.frame(' ')) %>%
  cbind(dat1[, 2]) %>% cbind(dat2[, 2]) %>% cbind(data.frame(' ')) %>%
  cbind(dat1[, 3]) %>% cbind(dat2[, 3]) %>% cbind(data.frame(' ')) %>%
  cbind(dat1[, 4]) %>% cbind(dat2[, 4]) %>%
  xtable() %>%
  print(include.rownames = F)

dat1 = dat.RMSE %>%
  filter(F.set == 4) %>%
  filter(method %in% c('oracle', 'Wang & Bishop', 'inflation & iterative', 'banding & inflation'))
dat2 = dat.RMSE %>%
  filter(F.set == 8) %>%
  filter(method %in% c('oracle', 'Wang & Bishop','inflation', 'banding', 'banding & iteration'))
dat3 = dat.RMSE %>%
  filter(F.set == 12) %>%
  filter(method %in% c('oracle', 'Wang & Bishop','inflation & iteration', 'banding & inflation & iteration'))

g = dat1 %>% rbind(dat2) %>% rbind(dat3) %>%
  mutate(method = factor(method, levels = methods)) %>%
  mutate(F.set = factor(paste0('F = ', F.set), levels = c('F = 4', 'F = 8', 'F = 12'))) %>%
  group_by(F.set, method) %>%
  mutate(RMSE = zoo::rollmean(RMSE, 4, fill = NA)) %>%
  ggplot()+
  geom_line(aes(x = X, y = RMSE, color = method))+
  facet_wrap(vars(F.set), ncol = 1)+
  labs(x = 'time', y = 'RMSE')+
  theme_bw()+
  theme(legend.position = 'bottom',
        legend.title = element_blank())
print(g)
ggsave('Output/Simulation_L96_RMSE_t_1007.png', g, width = 8, height = 6, dpi = 600)

dat1 = dat.obj %>%
  filter(F.set == 4) %>%
  filter(method %in% c('oracle', 'Wang & Bishop', 'inflation & iterative', 'banding & inflation'))
dat2 = dat.obj %>%
  filter(F.set == 8) %>%
  filter(method %in% c('oracle', 'Wang & Bishop','inflation', 'banding', 'banding & iteration'))
dat3 = dat.obj %>%
  filter(F.set == 12) %>%
  filter(method %in% c('oracle', 'Wang & Bishop','inflation & iteration', 'banding & inflation & iteration'))

g = dat1 %>% rbind(dat2) %>% rbind(dat3) %>%
  mutate(method = factor(method, levels = methods)) %>%
  mutate(F.set = factor(paste0('F = ', F.set), levels = c('F = 4', 'F = 8', 'F = 12'))) %>%
  group_by(F.set, method) %>%
  ggplot()+
  geom_line(aes(x = X, y = obj, color = method))+
  facet_wrap(vars(F.set), ncol = 1, scale = 'free_y')+
  labs(x = 'time', y = 'objective function')+
  theme_bw()+
  theme(legend.position = 'bottom',
        legend.title = element_blank())
print(g)
ggsave('Output/Simulation_L96_obj_t_1007.png', g, width = 8, height = 6, dpi = 600)

result %>% 
  mutate(method = factor(method, levels = methods)) %>% 
  group_by(F.set, method) %>% 
  summarise(converge = mean(cv)) %>% 
  arrange(method) %>% 
  spread(F.set, converge) %>%
  xtable() %>%
  print(include.rownames = F)
