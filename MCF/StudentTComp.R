library(R6)

StudentTComp = R6Class("StudentTComp", list(
  mean = 0,
  nu = 1,
  sigma = 1,
  total_time = 3,
  ai = c(1:10)*0.1,
  d = 0.1,
  
  initialize = function(mean = 0, nu = 4, sigma = 1)
  {
    self$mean = mean
    self$nu = nu
    self$sigma = sigma
    sol = sqrt((nu+4)*nu*sigma^2/(nu+2))
    solutions = c(sol+mean,-sol+mean)
    private$solutions = solutions
    private$peaks = self$phi(sol+mean)
    private$lb = self$phi(mean)
  },
  
  generateSample = function(total_samples=1)
  {
    return(rt(total_samples,self$nu)*self$sigma+self$mean)
  },
  
  density = function(x, log=TRUE)
  {
    nu = self$nu
    mean = self$mean
    s2 = self$sigma^2
    if (log)
    {
      d = log(gamma((nu+1)/2))-0.5*log(nu*pi*s2)-log(gamma((nu)/2))-(nu+1)/2*log(1+(x-mean)^2/nu/s2)
    } else
    {
      d = gamma((nu+1)/2)*(1+(x-mean)^2/nu/s2)^(-(nu+1)/2)/sqrt(nu*pi*s2)/gamma(nu/2)
    }
    return(d)
  },
  
  phi = function(x)
  {
    nu = self$nu
    mean = self$mean
    s2 = self$sigma^2
    return(0.5*(nu+1) * ((x-mean)^2 + nu*s2)^(-2) * ((nu+2)*(x-mean)^2 - nu*s2))
  },
  
  lowerBound = function(interval)
  {
    return(private$lb)
  },
  
  upperBound = function(interval)
  {
    sol = private$solutions[1]
    within = (sol < interval[1]) || (-sol > interval[2])
    if (within==1)
    {
      return(max(self$phi(interval)))
    } else
    {
      return(private$peaks[1])
    }
  }),
  private = list(solutions = vector(),
                 peaks = vector(), lb = 0)
)