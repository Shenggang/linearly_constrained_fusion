# The submission will be scored based on a skill score 
# which will be the ratio of the RMSE error of your submission 
# and the RMSE of the benchmark model (see below for the benchmark model). 

benchmark = function(gt, estimate)
{
  error = gt - estimate
  return(sqrt(mean(error^2)))
}