rcat = function(n, prob)
{
  stopifnot(all(prob>=0))
  prob = prob/sum(prob)
  steps = cumsum(prob)
  result = array(NaN, dim = n)
  for (i in 1:n)
  {
    result[i] = sum(runif(1) > steps)+1
  }
  return(result)
}