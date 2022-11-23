# for 2d-integral
library(rmutil)

print_table = function(vec, collapse= " & ")
{
  cat(paste(round(vec, digits=2), collapse=collapse),"\\\\\n")
}

normalize_weights = function(weights)
{
  nW = weights - max(weights)
  nW = exp(nW)/sum(exp(nW))
  return(nW)
}

idx = 1

cur_result = result_set[[idx]]

deg = cur_result$deg
mu = cur_result$mu
rhs = 10

cons_density = function(x,y)
{
  d = dt(x-mu[1], df=deg[1], log=TRUE) +
    dt(y-mu[2], df=deg[2], log=TRUE) +
    dt(rhs-x-y-mu[3], df=deg[3], log=TRUE)
  return(exp(d))
}
Z = int2(cons_density, a=c(-Inf,-Inf), b=c(Inf,Inf))

true_cons_density = function(x,y)
{
  d = dt(x-mu[1], df=deg[1], log=TRUE) +
    dt(y-mu[2], df=deg[2], log=TRUE) +
    dt(rhs-x-y-mu[3], df=deg[3], log=TRUE)
  return(exp(d)/Z)
}

m1x = int2(function(x,y){return(true_cons_density(x,y)*x)}, a=c(-Inf,-Inf), b=c(Inf,Inf))
m2x = int2(function(x,y){return(true_cons_density(x,y)*(x-m1x)^2)}, a=c(-Inf,-Inf), b=c(Inf,Inf))
m1y = int2(function(x,y){return(true_cons_density(x,y)*y)},  a=c(-Inf,-Inf), b=c(Inf,Inf))
m2y = int2(function(x,y){return(true_cons_density(x,y)*(y-m1y)^2)},  a=c(-Inf,-Inf), b=c(Inf,Inf))
m1z = int2(function(x,y){return(true_cons_density(x,y)*(rhs-x-y))}, a=c(-Inf,-Inf), b=c(Inf,Inf))
m2z = int2(function(x,y){return(true_cons_density(x,y)*(rhs-x-y-m1z)^2)}, a=c(-Inf,-Inf), b=c(Inf,Inf))

benchmark = c(m1x, m1y, m1z, m2x, m2y, m2z)
result_set[[idx]]$benchmark = benchmark

# benchmark = cur_result$benchmark

mcf_exact_particles = cur_result$mcf_exact_result$result
pf_particles = cur_result$linear_pf_result$particles
pf_weights = cur_result$linear_pf_result$weights
mh_particles = cur_result$mh_result[[1]]

mcf_N = dim(mcf_exact_particles)[1]
pf_N = dim(pf_particles)[1]

levels = c(0.01,0.1,1)

mcf_exact_levels = list()
pf_levels = list()
mh_levels = list()

for (k in 1:length(levels))
{
  level = levels[k]
  
  mcf_exact_est = 1:6
  for (i in 1:3)
  {
    mcf_exact_est[i] = mean(mcf_exact_particles[1:(level*mcf_N),i])
    mcf_exact_est[i+3] = mean((mcf_exact_particles[1:(level*mcf_N),i]-mcf_exact_est[i])^2)
  }
  mcf_exact_levels[[k]] = mcf_exact_est
  
  mh_est = 1:6
  for (i in 1:3)
  {
    mh_est[i] = mean(mh_particles[1:(level*mcf_N),i])
    mh_est[i+3] = mean((mh_particles[1:(level*mcf_N),i]-mh_est[i])^2)
  }
  mh_levels[[k]] = mh_est
  
  pf_est = 1:6
  for (i in 1:3)
  {
    pf_est[i] = pf_particles[1:(level*pf_N),i]%*%pf_weights[1:(level*pf_N)]/level
    pf_est[i+3] = (pf_particles[1:(level*pf_N),i]-pf_est[i])^2%*%pf_weights[1:(level*pf_N)]/level
  }
  pf_levels[[k]] = pf_est
}

# for(k in 1:3)
# {
#   print((abs(mcf_levels[[k+1]]-benchmark)-abs(pf_levels[[k]]-benchmark))/abs(benchmark))
# }


for (k in 1:3)
{
  print_table(-(abs(mcf_exact_levels[[k]]-benchmark)-abs(pf_levels[[k]]-benchmark))/abs(benchmark)*100, collapse = "\\% & ")
}

for (k in 1:3)
{
  print_table(-(abs(mcf_exact_levels[[k]]-benchmark)-abs(mh_levels[[k]]-benchmark))/abs(benchmark)*100, collapse = "\\% & ")
}
