# initialize
results = list()
particles = matrix(nrow = N, ncol = ncomps)

x0 = rlcnorm(1, consMat, rhs)

prev_W = 0
for (j in 1:ncomps)
{
  if (type=="gamma")
    prev_W = prev_W  + dgamma(x0[j], shape=shape[j], rate=rate[j], log = TRUE)
  if (type=="t")
    prev_W = prev_W + dt(x0[j]-mu[j], df=deg[j],log=TRUE)
  if (type=="glg")
    prev_W = prev_W + dglg(x0[j], a[j], b[j], gam[j], C[j], log=TRUE)
}
cur_X = x0
for (i in 1:(N+10000))
{
  if (i %% 2000 == 0)
    cat(strftime(Sys.time(), "%H:%M:%S")," --- MCMC Sample ",i,"-\n")
  cur_step = rlcnorm(1,consMat,0)*stepsize
  
  next_X = cur_X+cur_step
  cur_W = 0
  for (j in 1:ncomps)
  {
    if (type=="gamma")
      cur_W = cur_W  + dgamma(next_X[j], shape=shape[j], rate=rate[j], log = TRUE)
    if (type=="t")
      cur_W = cur_W + dt(next_X[j]-mu[j], df=deg[j],log=TRUE)
    if (type=="glg")
      cur_W = cur_W + dglg(next_X[j], a[j], b[j], gam[j], C[j], log=TRUE)
  }
  acc = cur_W-prev_W
  if (log(runif(1)) < acc)
  {
    if (i > 10000)
      particles[i-10000,] = next_X
    prev_W = cur_W
    cur_X = next_X
  } else
  {
    if (i > 10000)
      particles[i-10000,] = cur_X
  }
}

results = list(particles)