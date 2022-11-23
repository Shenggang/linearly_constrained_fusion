library(moments)

# glg definition:
# For X1, X2 both Gamma distributions with shape alpha and beta and rate 1,
# Y is glg(alpha, beta, gamma,C) is defined as gamma*log(X1/X2)+C


glg_moments = function(x)
{
  # input x = [alpha,beta,gamma]
  k2 = x[3]^2*(trigamma(x[1])+trigamma(x[2]))
  k3 = x[3]^3*(psigamma(x[1],2) - psigamma(x[2],2))
  k4 = x[3]^4*(psigamma(x[1],3) + psigamma(x[2],3))
  return(c(k2,k3,k4))
}


glg_lr = function(X, Y, allow_intercept=TRUE)
{
  # X ---- Observation/explainatory variables, d x n
  # Y ---- Response, n x 1 vector
  # n ---- number of individual data entries, no need to input
  stopifnot(dim(X)[2]== dim(Y)[1])
  
  d = dim(X)[1]; n = dim(X)[2]
  
  
  # compute regression parameters
  df = data.frame(cbind(Y, t(X)))
  X_temp = rbind(array(1, dim = c(1,n)), X)
  starting = solve(X_temp%*%t(X_temp), X_temp%*%Y)
  if (allow_intercept)
  {
    beta = starting
    btx = t(X_temp)%*%beta
  } else
  {
    beta = solve(X%*%t(X), X%*%Y)
    btx = t(X)%*%beta
  }
  
  # compute distribution parameters
  residue = Y-btx
  target = c(var(residue), skewness(residue), kurtosis(residue)-3)
  target[3] = max(0, target[3])
  
  sol = optim(c(2,2,1), function(x){return(sum((glg_moments(x)-target)^2)+1e-3*x[1]^2+1e-3*x[2]^2+1e-6*x[3]^2)}, 
              method="L-BFGS-B", lower = c(0.1,0.1,0), upper=c(10,10,1e6))
  sol = sol$par
  # sol = nleqslv(c(2,2,1), function(x){return(glg_moments(x)-target)}, method="Newton", control = list(maxit=1000, allowSingular=TRUE))
  # sol = sol$x
  # if (any(sol<0))
  #   print(target)
  C = -sol[3]*(digamma(sol[1])-digamma(sol[2]))
  glg_param = list(a=sol[1], b=sol[2], gamma=sol[3], C=C)

  return(list(beta = beta, glg_param = glg_param))
}

glg_ar = function(Y,order = 1, xreg= NaN, allow_intercept=TRUE)
{
  # Y     ---- Response time series, n x L
  # n     ---- number of individual data entries
  # L     ---- number of time points
  # order ---- order of the autoregressive model, 0 is equiv to linear regression
  # xreg  ---- additional covariates, n*L x d
  
  stopifnot(order >= 0, round(order)==order)
  
  if (is.null(dim(Y)))
  {
    Y = matrix(Y, nrow=1)
  }
  
  iXreg = !any(is.nan(xreg))
  d = 0
  if (iXreg)
  {
    d = dim(xreg)[2]
  }
  
  n = dim(Y)[1]; L = dim(Y)[2]
  # recreate observation matrix X
  ext = L-order
  Xnew = array(NaN, dim = c(d+order, n*ext))
  Ynew = array(NaN, dim = c(n*ext,1))
  
  for (i in 1:n)
  {
    Ynew[1:ext+(i-1)*ext,] = Y[i, (order+1):L]
    if (order > 0)
    {
      for (j in 1:ext)
      {
        Xnew[1:order, j+(i-1)*ext] = Y[i, j:(j+order-1)]
      }
    }
    if (d > 0)
      Xnew[1:d+order, 1:ext+(i-1)*ext] = t(xreg[((order+1):L)+(i-1)*L, ]) 
  }
  
  result = glg_lr(Xnew, Ynew, allow_intercept)
  glg_param = result$glg_param
  
  alpha = NaN
  beta = NaN
  intercept=NaN
  if (allow_intercept)
  {
    intercept = result$beta[1]
  }
  if (order > 0)
  {
    alpha = result$beta[1:order+allow_intercept]
  }
  if (!any(is.nan(xreg)))
  {
    beta = result$beta[-(1:(allow_intercept+order))]
  }
  return(list(intercept=intercept, alpha = alpha, beta=beta, glg_param=glg_param))
}