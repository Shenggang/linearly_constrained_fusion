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

idx=1

cur_result = result_set[[idx]]

consMat = matrix(c(1,1,1), nrow=1)

a = cur_result$params$a
b = cur_result$params$b
gamma = cur_result$params$gamma
C = cur_result$params$C
rhs = 10

dglg = function(x,a,b,gam,C, log=FALSE)
{
  psi = exp((x-C)/gam)
  ipsi = exp(-(x-C)/gam)
  if(log)
  {
    d = -log(gam*beta(a,b)) - a*log(1+ipsi) - b*log(1+psi)
  } else
  {
    d = 1/((gam*beta(a,b))*(1+ipsi)^a*(1+psi)^b)
  }
  return(d)
}

cons_density = function(x,y)
{
  d = dglg(x,a[1],b[1],gamma[1], C[1]) *
    dglg(y,a[2],b[2],gamma[2], C[2]) *
    dglg(rhs-x-y, a[3],b[3],gamma[3], C[3])
  return(d)
}
Z = int2(cons_density, a=c(-Inf,-Inf), b=c(Inf,Inf))

true_cons_density = function(x,y)
{
  d = dglg(x,a[1],b[1],gamma[1], C[1]) *
    dglg(y,a[2],b[2],gamma[2], C[2]) *
    dglg(rhs-x-y, a[3],b[3],gamma[3], C[3])
  return(d/Z)
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
  
  pf_est = 1:6
  for (i in 1:3)
  {
    pf_est[i] = pf_particles[1:(level*pf_N),i]%*%pf_weights[1:(level*pf_N)]/level
    pf_est[i+3] = (pf_particles[1:(level*pf_N),i]-pf_est[i])^2%*%pf_weights[1:(level*pf_N)]/level
  }
  pf_levels[[k]] = pf_est
  
  mh_est = 1:6
  for (i in 1:3)
  {
    mh_est[i] = mean(mh_particles[1:(level*mcf_N),i])
    mh_est[i+3] = mean((mh_particles[1:(level*mcf_N),i]-mh_est[i])^2)
  }
  mh_levels[[k]] = mh_est
}

for (k in 1:3)
{
  print_table(-(abs(mcf_exact_levels[[k]]-benchmark)-abs(pf_levels[[k]]-benchmark))/abs(benchmark)*100, collapse = "\\% & ")
}


for (k in 1:3)
{
  print_table(-(abs(mcf_exact_levels[[k]]-benchmark)-abs(mh_levels[[k]]-benchmark))/abs(benchmark)*100, collapse = "\\% & ")
}