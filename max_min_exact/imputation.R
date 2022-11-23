impute_energy = function(total_samples, total_steps, comps_per_step, initial_data, parameters,
                         positive = TRUE, margin = 10, max_trial=1000, xreg=NaN, idx=1)
{
  order = parameters$order
  
  # initialize result matrix
  max_samples = total_samples*10
  pred_consumption = array(NaN, dim = c(total_steps+order, comps_per_step, max_samples))
  for (i in 1:max_samples)
  {
    pred_consumption[1:order, ,i] = initial_data
  }
  
  #predict data
  count = 0
  total = 0
  while ((count < total_samples) && (total < max_samples))
  {
    cur_idx = total+1
    if (cur_idx %% 1000 == 0)
      cat(strftime(Sys.time(), "%H:%M:%S")," --- Slot ",idx,"--- Sample ",cur_idx,"-- accepted = ", count, "\n")
    
    for (day in 1:total_steps+order)
    {
      cur_rhs = parameters$consY[day-order]
      # extract time series
      day_range = (day-order):(day-1)
      cur_Y = array(pred_consumption[day_range, , cur_idx], dim = c(length(day_range), comps_per_step))
      
      # construct components
      if (!any(is.nan(xreg)))
      {
        x_in = xreg[,day-order]
      }
      comps = create_comps(cur_Y, comps_per_step, parameters, xreg=x_in)
      
      # produce sample
      for (tc in parameters$time_pool)
      {
        outcome = mcf_cons_ctrl(1, tc, comps, parameters$consMat, 
                                cur_rhs, positive = positive, suppressed = TRUE, max_trial=max_trial)
        if (!outcome$failed)
          break;
      }
      if (outcome$failed)
      {
        total=total+1
        break;
      }
      pred_consumption[day, ,cur_idx] = outcome$samples
    }
    if (day==total_steps+order && !outcome$failed)
    {
      count=count+1
      total=total+1
    }
  }
  pred_consumption = pred_consumption[,,!is.nan(apply(pred_consumption, 3, sum))]
  return(list(result = pred_consumption, total, count))
}

MHFC_get_sample = function(comps, consMat, cur_rhs, total_time, positive = TRUE, margin = 10, max_trial=10)
{
  # produce samples
  trials = 0
  rst = MHFC_gen_sample_layered(comps, consMat, cur_rhs, total_time, positive = TRUE, margin = 10)
  while ((is.infinite(rst$est) || is.nan(rst$est)) && (trials < max_trial))
  {
    rst = MHFC_gen_sample_layered(comps, consMat, cur_rhs, total_time, positive = TRUE, margin = 10)
    trials = trials+1
  }
  return(rst)
}

create_comps = function(cur_Y, comps_per_step, parameters, xreg=NaN)
{
  comps = list()
  comps[[comps_per_step]] = NULL
  glg_params = parameters$glg_params
  alpha = parameters$alpha
  intercept = parameters$intercept
  for (r in 1:comps_per_step)
  {
    means = alpha[, r]%*% cur_Y[ , r] + intercept[r]
    if (!any(is.nan(xreg)))
    {
      means = means + parameters$beta[, r] %*% xreg
    }
    glg_param = glg_params[[r]]
    comps[[r]] = GenLogisticComp$new(glg_param$a, glg_param$b, glg_param$gamma, glg_param$C+means)
    comps[[r]]$ai = c(1:10)*0.1
    comps[[r]]$d = 0.1
  }
  return(comps)
}


