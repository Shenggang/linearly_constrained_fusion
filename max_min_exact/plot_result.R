source('benchmark.R')

# load ts
lpday = 24*10 # number of entries per day
test_days = 7
total_rows = dim(results[[1]]$result)[1]

# load observation constraints
obs_raw = read.csv('./observation_subset.csv', header=T)
obs_length = nrow(obs_raw)
constraints = obs_raw[(obs_length-test_days*24*2+1):obs_length, 2]
baseline = constraints[1:((total_rows-5)*48)]


# extract results
fitted_quantile = matrix(nrow = 2, ncol=48*(total_rows-5))
for (row in 6:total_rows)
{
  for (i in 1:48)
  {
    fitted_quantile[1,i+(row-6)*48] = min(apply(results[[i]]$result[row,,], 1, mean))
    fitted_quantile[2,i+(row-6)*48] = max(apply(results[[i]]$result[row,,], 1, mean))
  }
}


# read true data
gt = read.csv('./target_subset.csv', header=T)
discarded_idx = 1:(nrow(gt) - test_days*2*24)
gt = data.matrix(gt[-discarded_idx,])


library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggpubr)

merged_max_min = cbind(fitted_quantile[1,], fitted_quantile[2,], gt[1:((total_rows-5)*48),2:3], baseline)
colnames(merged_max_min) = c('min_pred', 'max_pred', 'gt_max', 'gt_min', 'baseline')
data = tibble(data.frame(merged_max_min)) %>% mutate(t=1:((total_rows-5)*48))

min_tibble = data[c('min_pred', 'gt_min', 'baseline', 't')] %>% 
  pivot_longer(cols = c('min_pred', 'gt_min', 'baseline'), names_to = 'group', values_to='min_value' )

g_min = ggplot(data=min_tibble, aes(color=group, shape=group)) + 
  geom_line(aes(t, min_value, group=group), size=1) + 
  geom_point(aes(t, min_value, group=group), size=2) + 
  labs(x='', y='Min Value', tag='(a)') +  scale_colour_discrete(labels =c('Baseline', 'True Value', 'Predicted'))+
  scale_shape_discrete(labels = c('Baseline', 'True Value', 'Predicted'))+
  theme(legend.title =  element_blank(), text = element_text(size = 20))
g_min
#ggsave('min_plot.jpg', g_min, width=width, height=height, device='jpg', dpi=300)

max_tibble = data[c('max_pred', 'gt_max', 'baseline', 't')] %>% 
  pivot_longer(cols = c('max_pred', 'gt_max', 'baseline'), names_to = 'group', values_to='max_value' )

g_max = ggplot(data=max_tibble, aes(color=group, shape=group)) + 
  geom_line(aes(t, max_value, group=group), size=1) + 
  geom_point(aes(t, max_value, group=group), size=2) + 
  labs(x='Time Index', y='Max Value', tag='(b)') +  scale_colour_discrete(labels =c('Baseline', 'True Value', 'Predicted'))+
  scale_shape_discrete(labels = c('Baseline', 'True Value', 'Predicted'))+
  theme(legend.title =  element_blank(), text = element_text(size = 20))
g_max
#ggsave('max_plot.jpg', g_max, width=width, height=height, device='jpg', dpi=300)
arranged = ggarrange(g_min, g_max, nrow=2, common.legend = TRUE, legend='bottom')
arranged

width=16
height=12
ggsave(paste('min_max_plot_',test_days,'.jpg',sep = ""), arranged, width=width, height=height, device = 'jpg', dpi=300)
