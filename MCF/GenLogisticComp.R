library(R6)

GenLogisticComp = R6Class("GenLogisticComp", list(
  alpha = 3,
  beta = 3,
  gamma = 1,
  C = 0,
  total_time = 1,
  ai = c(1:10)*0.1,
  d = 0.1,
  
  initialize = function(alpha, beta, gamma, C)
  {
    self$alpha = alpha
    self$beta = beta
    self$gamma = gamma
    self$C = C
  },
  
  generateSample = function(n = 1)
  {
    x1 = rgamma(n, shape=self$alpha, rate=1)
    x2 = rgamma(n, shape=self$beta, rate=1)
    y = self$gamma * log(x1/x2) + self$C
    return(y)
  },
  
  density = function(x, log=TRUE)
  {
    a = self$alpha
    b = self$beta
    gam = self$gamma
    psi = exp((x-self$C)/gam)
    ipsi = exp(-(x-self$C)/gam)
    if(log)
    {
      d = -log(gam*beta(a,b)) - a*log(1+ipsi) - b*log(1+psi)
    } else
    {
      d = 1/((gam*beta(a,b))*(1+ipsi)^a*(1+psi)^b)
    }
    return(d)
  },
  
  phi = function(x)
  {
    a = self$alpha
    b = self$beta
    gam = self$gamma
    psi = exp((x-self$C)/gam)
    ipsi = exp(-(x-self$C)/gam)
    return (a^2*(1/(1+psi))^2 + b^2*(1/(1+ipsi))^2 - (2*a*b+a+b)*(2+psi+ipsi)^(-1))/(2*gam^2)
  },
  
  lowerBound = function(interval)
  {
    xrange = seq(interval[1], interval[2], length.out = 100)
    return(min(self$phi(xrange)))
  },
  
  upperBound = function(interval)
  {
    xrange = seq(interval[1], interval[2], length.out = 100)
    return(max(self$phi(xrange)))
  }
)
)