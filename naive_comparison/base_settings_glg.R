setwd("../MCF/")
source("./load.R")

setwd("../naive_comparison/")
source("./categorical_distribution.R")

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

######### default settings #################
type="glg"

a = c(3,3,3)
b = c(0.2,0.3,0.5)
gam = c(1,1,1)
C = c(0,0,0)

ncomps = length(a)

# constraint
consMat = matrix(c(1,1,1), nrow=1)
rhs = 10

# samples
N = 10000

#####################################
