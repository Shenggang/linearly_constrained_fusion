total_days = 14
total_samples = 10000
houses = 1
shift = 15

# extract constraints
consY = array(consumption_array[(1:houses+n+shift), 1:total_days+order,], dim = c(houses, total_days, 3))
total_sum = t(apply(consY, c(2,3), sum)) # dim = c(3, days)
indi_sum = apply(consY, c(1,2), sum) # dim = c(houses, days)

pf_pred = array(NaN, dim = c(houses, total_days+order, 3, total_samples))
predictor = t(array(survey_50[(1:houses+n+shift), ], dim = c(houses, dim(survey_50)[2]))) # d x houses

#initialize first 7 days
for (i in 1:total_samples)
{
  pf_pred[, 1:order, , i] = consumption_array[1:houses+n, 1:order, ]
}

# Predict data
for (day in 1:total_days+order)
{
  for (count in 1:total_samples)
  {
    # extract time series
    day_range = (day-order):(day-1)
    cur_Y = array(pf_pred[ , day_range, , count], dim = c(houses, length(day_range), 3))
    
    # compute mean
    means = t(t(beta)%*% predictor) + intercept # houses x 3
    for (k in 1:3)
    {
      means[, k] = means[ ,k] + alpha[, k]%*% t(array(cur_Y[ , , k], dim = c(houses,length(day_range))))
    }
    
    # initialise components
    outcome = array(NaN, dim = c(houses, 3))
    for (idx in 1:houses)
    {
      for (r in 1:3)
      {
        glg_param = glg_params[[r]]
        comp = GenLogisticComp$new(glg_param$a, glg_param$b, glg_param$gamma, glg_param$C+means[idx,r])
        outcome[idx,r] = max(0,comp$generateSample())
      }
    }
    # draw samples
    pf_pred[, day, ,count] = outcome
  }
}

pf_pred = pf_pred[,1:total_days+order, ,]
