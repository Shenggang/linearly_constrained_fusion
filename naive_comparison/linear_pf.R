# initialize
particles = matrix(nrow = N, ncol = ncomps)

# components
# approximate by Gaussian
mu_temp=c()
std=c()
if (type=="gamma")
{
  for (i in 1:ncomps)
  {
    mu_temp[i] = shape[i]/rate[i]
    std[i] = sqrt(shape[i]/rate[i]/rate[i])
  }
  positive=TRUE
} else if (type=="t")
{
  for (i in 1:ncomps)
  {
    mu_temp[i] = mu[i]
    std[i] = sqrt(deg[i]/(deg[i]-2))
    if (deg[i] <=2)
      std[i] = abs(mu[i])
  }
  positive=FALSE
} else if (type=="glg")
{
  for (i in 1:ncomps)
  {
    mu_temp[i] = C[i] + gam[i]*(digamma(a[i])- digamma(b[i]))
    std[i] = gam[i]*sqrt(trigamma(a[i]) +trigamma(b[i]))
  }
  positive=FALSE
}


sqrtSigma = diag(std)

c = rhs-consMat%*%mu_temp
B = consMat%*%sqrtSigma


# importance sampling via direct Gaussian approximation
weights = c()
for (i in 1:N)
{
  if (i %% 2000 == 0)
    cat(strftime(Sys.time(), "%H:%M:%S")," --- Sample ",i,"-\n")
  
  trials = 0
  
  # generate new sample
  # reject if negative
  y =  sqrtSigma%*%t(rlcnorm(1, B, c)) + mu_temp
  while (!all(y>0) && (trials < 100) && positive)
  {
    y =  sqrtSigma%*%t(rlcnorm(1, B, c)) + mu_temp
    trials = trials + 1
  }
  particles[i, ] = y
  
  # compute weight
  if (!all(y>0) && positive)
  {
    weights[i]=-Inf
  } else{
    weights[i] = 0
    for (j in 1:ncomps)
    {
      if (type=="gamma")
        weights[i] = weights[i] + dgamma(y[j], shape=shape[j], rate=rate[j], log = TRUE)-dnorm(y[j], mean=mu_temp[j], sd=std[j], log=TRUE)
      if (type=="t")
        weights[i] = weights[i] + dt(y[j]-mu_temp[j], df=deg[j],log=TRUE) - dnorm(y[j], mean=mu_temp[j], sd=std[j], log=TRUE)
      if (type=="glg")
        weights[i] = weights[i] + dglg(y[j], a[j], b[j], gam[j], C[j], log=TRUE) - dnorm(y[j], mean=mu_temp[j], sd=std[j], log=TRUE)
    }
  }
}

nW = weights - max(weights)
nW = exp(nW)/sum(exp(nW))
weights=nW
# Compute ESS
ESS = sum(nW)^2/sum(nW^2)
print(ESS)

results = list(particles=particles, weights=weights, ESS=ESS)