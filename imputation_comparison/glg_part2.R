max_trial = 1000

total_days = 14
total_samples = 10000
max_samples = 1e5
houses = 1
shift = 15
time_pool = c(0.01, 0.1, 0.5, 1, 1.5, 2)

# extract constraints
consY = array(consumption_array[(1:houses+n+shift), 1:total_days+order,], dim = c(houses, total_days, 3))
total_sum = t(apply(consY, c(2,3), sum)) # dim = c(3, days)
indi_sum = apply(consY, c(1,2), sum) # dim = c(houses, days)

pred_consumption = array(NaN, dim = c(houses, total_days+order, 3, max_samples))
predictor = t(array(survey_50[(1:houses+n+shift), ], dim = c(houses, dim(survey_50)[2]))) # d x houses

#initialize first 7 days
for (i in 1:max_samples)
{
  pred_consumption[, 1:order, , i] = consumption_array[1:houses+n, 1:order, ]
}

# initialise constraint matrix
consMat = array(0, dim = c(houses, 3*houses))

for (i in 1:houses)
{
  # individual house over the whole day
  consMat[i, 0:2*houses+i] = 1
}

# extract constaints
rhs = array(NaN, dim = c(dim(consMat)[1], total_days))

rhs[1:houses,] = indi_sum

for (i in 1:houses)
{
  # individual house over the whole day
  consMat[i, 0:2*houses+i] = 1
}


# Predict data
count = 0
total = 0
while ((count < total_samples) && (total < max_samples))
{
  cur_idx = total+1
  if ((cur_idx %% 2000 == 0) || (cur_idx==1))
    cat(strftime(Sys.time(), "%H:%M:%S")," --- Sample ",cur_idx,"-\n")
  
  for (day in 1:total_days+order)
  {
    cur_rhs = rhs[, day-order]
    # extract time series
    day_range = (day-order):(day-1)
    cur_Y = array(pred_consumption[ , day_range, , cur_idx], dim = c(houses, length(day_range), 3))
    
    # compute mean
    means = t(t(beta)%*% predictor) + intercept # houses x 3
    for (k in 1:3)
    {
      means[, k] = means[ ,k] + alpha[, k]%*% t(array(cur_Y[ , , k], dim = c(houses,length(day_range))))
    }
    
    # construct components
    comps = list()
    comps[[3*houses]] = NULL
    for (idx in 1:houses)
    {
      for (r in 1:3)
      {
        glg_param = glg_params[[r]]
        comps[[idx + (r-1)*houses]] = GenLogisticComp$new(glg_param$a, glg_param$b, glg_param$gamma, glg_param$C+means[idx,r])
      }
    }
    
    # produce sample
    for (tc in time_pool)
    {
      outcome = mcf_cons_ctrl(1, tc, comps, consMat, cur_rhs, positive = TRUE, suppressed = TRUE)
      if (!outcome$failed)
        break;
    }
    if (outcome$failed)
    {
      total=total+1
      break;
    }
    pred_consumption[, day, ,cur_idx] = outcome$samples
  }
  if (day==total_days+order)
  {
    count=count+1
    total=total+1
  }
}

pred_consumption = pred_consumption[1,,,]
pred_consumption = pred_consumption[,,!is.nan(apply(pred_consumption, 3, sum))]

filename = "Simulation_exact_glg.RData"
sim_result = list(houses = houses, total_days=total_days, pred_consumption=pred_consumption, consY=consY, total=total, count=count)
save(sim_result, file = filename)
