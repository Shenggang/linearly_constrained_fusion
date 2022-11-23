library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggpubr)

plot_res = function(mean_est, gt, upper, lower, ymax, filename)
{
  temp = c()
  for (i in 1:3)
  {
    seg = cbind(1:14, mean_est[,i], gt[,i], upper[,i], lower[,i],i)
    temp=rbind(temp,seg)
  }
  colnames(temp) = c("Day", "Mean", "GT", "Upper_Qtl", "Lower_Qtl", "interval")
  
  data = tibble(data.frame(temp))
  
  legend_labs = c('Mean', 'Ground Truth', '2.5%-quantile')
  linetypes=  c("solid", "solid", "longdash")
  plotcolors = c("red", "blue", "darkgreen")
  pointshapes = c(5, 0, 1)
  ylab = c("Consumption", "","")
  title_op = c("7.00-15.00", "15.00-23.00", "23.00-7.00")
  
  unconsylims = c(min(lower),max(upper*1.4))
  
  pics = list()
  for (i in 1:3)
  {
    pics[[i]] = data %>% filter(interval==i) %>%
      pivot_longer(cols = c("Mean", "GT", "Upper_Qtl", "Lower_Qtl"), names_to="group", values_to="consumption")%>%
      mutate(aes = factor(case_when(group=="Mean" ~ 1,
                                    group=="GT" ~ 2,
                                    TRUE ~ 3))) %>%
      ggplot(aes(x=Day)) +  ylim(0, ymax) + labs(title= title_op[i], x="Day", y=ylab[i]) +
      geom_line(aes(y=consumption, group=group, color=aes, linetype=aes), size=1)+
      geom_point(aes(y=consumption, group=group, color=aes, shape=aes, size=aes)) +
      scale_linetype_manual("Consumption", values=linetypes, labels=legend_labs) +
      scale_color_manual("Consumption", values=plotcolors, labels=legend_labs)+
      scale_shape_manual("Consumption", values = pointshapes, labels=legend_labs)+
      scale_size_manual("Consumption", values=c(2.5,2.5,0), labels=legend_labs) +
      theme(text = element_text(size = 30))
  }
  arranged = ggarrange(pics[[1]], pics[[2]], pics[[3]], nrow=1, common.legend = TRUE, legend='bottom')
  arranged
  width=18
  height=6
  ggsave(filename, arranged, width=width, height=height, device = 'jpg', dpi=300)
}


# data from       pred_consumption     dim = total_days, 3, total_samples
# true data from  consY                dim = total_days, 3
pred_result = sim_result$pred_consumption[8:21,,]

mean_est = apply(pred_result, c(1,2), mean)
var_est = apply(pred_result, c(1,2), var)


q = 0.025 # 2.5%
upper = apply(pred_result, c(1,2), function(array){ return(quantile(array, prob=1-q))})
lower = apply(pred_result, c(1,2), function(array){ return(quantile(array, prob=q))})

source("unconstrained_part2_glg.R")
pf_mean_est = apply(pf_pred, c(1,2), mean)
pf_upper = apply(pf_pred, c(1,2), function(array){ return(quantile(array, prob=1-q))})
pf_lower = apply(pf_pred, c(1,2), function(array){ return(quantile(array, prob=q))})
ymax = max(pf_upper)

gt = consY[1,,]

plot_res(mean_est = mean_est, gt = gt, upper = upper, lower = lower, 
         ymax = ymax, filename = './plots/constrained_exact.jpg')
plot_res(mean_est = pf_mean_est, gt = gt, upper = pf_upper, lower = pf_lower, 
         ymax = ymax, filename = './plots/unconstrained_pf.jpg')


