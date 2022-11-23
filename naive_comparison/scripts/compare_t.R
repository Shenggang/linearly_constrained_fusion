result_set = list()
total_time=0.1

#######################################
# small uncertainty
source("base_settings_t.R")

deg = c(5,5,5)
mu = c(-2, 3, 5)

N=1e5

source("mcf_exact.R")
mcf_result = results

stepsize = 1
source("rw_mcmc.R")
mh_result = results

source("linear_pf.R")
linear_pf_result = results

result_set[[1]] = list(deg=deg, mu=mu, total_time=total_time, 
                       mcf_exact_result = mcf_result, 
                       linear_pf_result = linear_pf_result,
                       mh_result = mh_result)

######################################

######################################
# normal uncertainty
source("base_settings_t.R")

deg = c(2+1e-2,2+1e-2,2+1e-2)
mu = c(-2, 3, 5)

N=1e5
source("mcf_exact.R")
mcf_result = results

stepsize = 1
source("rw_mcmc.R")
mh_result = results

source("linear_pf.R")
linear_pf_result = results

result_set[[2]] = list(deg=deg, mu=mu, total_time=total_time, 
                       mcf_exact_result = mcf_result, 
                       linear_pf_result = linear_pf_result,
                       mh_result = mh_result)

######################################

save(result_set, file="convergence_result_t.RData")

