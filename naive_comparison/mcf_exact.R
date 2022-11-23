# initialize
results = list()
particles = matrix(nrow = N, ncol = ncomps)

# components
comps = list()
if (type=="gamma")
{
  for (i in 1:ncomps)
  {
    comps[[i]] =  GammaComp$new(alpha = shape[i], beta = rate[i])
    comps[[i]]$ai=c(1:5)*1e-2
    comps[[i]]$d = 1e-2
  }
  positive=TRUE
} else if (type=="t")
{
  for (i in 1:ncomps)
  {
    comps[[i]] =  StudentTComp$new(mean = mu[i], nu = deg[i])
  }
  positive=FALSE
} else if (type=="glg")
{
  for (i in 1:ncomps)
  {
    comps[[i]] = GenLogisticComp$new(alpha=a[i], beta=b[i], gamma=gam[i], C=C[i])
  }
  positive=FALSE
}

results = mcf_cons(N, total_time, comps, consMat, rhs, positive=positive)