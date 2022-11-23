result_set = list()

######## less-skewed ###############
source("base_settings_glg.R")

a = c(3,3,3)
b = c(1.2,2,2)
gam = c(1,2,2)
C = c(0,0,0)

N=1e5

total_time=0.1
source("mcf_exact.R")
mcf_result = results

stepsize = 1
source("rw_mcmc.R")
mh_result = results

source("linear_pf.R")
linear_pf_result = results

result_set[[1]] = list(params=list(a=a,b=b,gamma=gam,C=C),  total_time=total_time,
                       mcf_exact_result = mcf_result,
                       mh_result = mh_result,
                       linear_pf_result=linear_pf_result)

######################################

######## positive-skewed ###############
source("base_settings_glg.R")

a = c(3,3,3)
b = c(0.4,0.4,0.4)
gam = c(2,1,1)
C = c(-5,-2,-3)

N=1e5

total_time=0.1
source("mcf_exact.R")
mcf_result = results

stepsize = 1
source("rw_mcmc.R")
mh_result = results

source("linear_pf.R")
linear_pf_result = results

result_set[[1]] = list(params=list(a=a,b=b,gamma=gam,C=C),  total_time=total_time, 
                       mcf_exact_result = mcf_result,
                       mh_result = mh_result,
                       linear_pf_result=linear_pf_result)


######################################


save(result_set, file="convergence_result_glg.RData")





