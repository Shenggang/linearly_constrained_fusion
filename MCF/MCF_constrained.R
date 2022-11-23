check_ppp_pass = function(ppp, bb, comp)
{
  time_points = ppp[1,]
  values = ppp[2,]
  for (r1 in consc(1, length(time_points)))
  {
    t = time_points[r1]
    for (r2 in consc(r1,length(bb[1,])))
    {
      if (bb[1,r2] == t)
      {
        if (comp$phi(bb[2,r2])  > values[r1])
        {
          return(TRUE)# rejected = TRUE
        }
      }
    }
  }
  return(FALSE)
}

mcf_cons_generate_y = function(total_time, total_comps, A, b, x, positive = FALSE)
{
  y = (sqrt(total_time)*diag(total_comps))%*%t(rlcnorm(1, A, b))+x
  count = 1
  if (positive)
  {
    if (!all(y>0))
      return(NaN)
  }
  return(y)
}
  
# run LCF with maximum trial number, may not fulfill the sample number requirement if the acceptance rate is too low. 
mcf_cons_ctrl = function(total_samples, total_time, comps, A, b, positive = FALSE, suppressed = TRUE, max_trial=1e3)
{
  total_comps = length(comps)
  i = 0; trials = 0; result = matrix(nrow = total_samples, ncol = total_comps)
  failed=FALSE
  ptm = proc.time()
  while (i < total_samples)
  {
    trials = trials + 1
    if (trials > max_trial*total_samples)
    {
      failed=TRUE
      break;
    }
    # simulate x_i, y
    x = vector()
    y = vector()
    for(j in c(1:total_comps))
    {
      x[j] = comps[[j]]$generateSample()
    }
    B = A%*%(sqrt(total_time)*diag(total_comps))
    alpha = b - A%*%x
    y = mcf_cons_generate_y(total_time, total_comps, B, alpha, x, positive)
    if (is.nan(y[1]))
    {
      sprint("Failed to find y", suppressed)
      next
    }
    # check AP1
    if (runif(1) <= norm_cnst(x,total_time*diag(total_comps),A,b))
    {
      rejected = rep(FALSE, total_comps)
      #for each component
      for (j in c(1:total_comps))
      {
        ai = comps[[j]]$ai; d = comps[[j]]$d
        I = layer(total_time, x[j], y[j], ai, d)
        K = fetch_bound(I-1, ai, d); L = fetch_bound(I, ai, d)
        interval = c(min(x[j],y[j])-L, max(x[j],y[j])+L)
        #simulate PPP
        lb = comps[[j]]$lowerBound(interval)
        ub = comps[[j]]$upperBound(interval) - lb
        if (ub*total_time > 1000)
        {
          rejected[j]=TRUE
          break;
        }
        ppp = poisson_point_process(total_time, ub)
        if (length(ppp) == 1)
        {
          rejected[j] = TRUE
        } else{
          #simulate BB
          time_points = ppp[1,]; values = ppp[2,] + lb
          bb = layered_BB(total_time, x[j], y[j], time_points, K, L)
          #check if pass
          rejected[j] = check_ppp_pass(rbind(time_points, values), bb, comps[[j]])
          if (rejected[j])
          {
            break;
          }
        }
      }
      if (sum(rejected) == 0)
      {
        i = i+1; result[i,] = y
        sprint(i, suppressed)
        if (i%%2000==0)
          print(paste(Sys.time(),i,sep = "-----"))
      } else
      {
        sprint("rejected", suppressed)
      }
    }else
    {
      sprint("rejected in AP1", suppressed)
    }
  }
  return(list(samples=result, failed=failed, trials = trials, rate = total_samples/trials, time =  proc.time()-ptm))
}

mcf_cons = function(total_samples, total_time, comps, A, b, positive = FALSE, suppressed = TRUE)
{
  total_comps = length(comps)
  i = 0; trials = 0; result = matrix(nrow = total_samples, ncol = total_comps)
  ptm = proc.time()
  while (i < total_samples)
  {
    trials = trials + 1
    # simulate x_i, y
    x = vector()
    y = vector()
    for(j in c(1:total_comps))
    {
      x[j] = comps[[j]]$generateSample()
    }
    B = A%*%(sqrt(total_time)*diag(total_comps))
    alpha = b - A%*%x
    y = mcf_cons_generate_y(total_time, total_comps, B, alpha, x, positive)
    if (is.nan(y[1]))
    {
      sprint("Failed to find y", suppressed)
      next
    }
    # check AP1
    if (runif(1) <= norm_cnst(x,total_time*diag(total_comps),A,b))
    {
      rejected = rep(FALSE, total_comps)
      #for each component
      for (j in c(1:total_comps))
      {
        ai = comps[[j]]$ai; d = comps[[j]]$d
        I = layer(total_time, x[j], y[j], ai, d)
        K = fetch_bound(I-1, ai, d); L = fetch_bound(I, ai, d)
        interval = c(min(x[j],y[j])-L, max(x[j],y[j])+L)
        #simulate PPP
        lb = comps[[j]]$lowerBound(interval)
        ub = comps[[j]]$upperBound(interval) - lb
        if (ub*total_time > 5000)
        {
          rejected[j]=TRUE
          break;
        }
        ppp = poisson_point_process(total_time, ub)
        if (length(ppp) == 1)
        {
          rejected[j] = TRUE
        } else{
          #simulate BB
          time_points = ppp[1,]; values = ppp[2,] + lb
          bb = layered_BB(total_time, x[j], y[j], time_points, K, L)
          #check if pass
          rejected[j] = check_ppp_pass(rbind(time_points, values), bb, comps[[j]])
          if (rejected[j])
          {
            break;
          }
        }
      }
      if (sum(rejected) == 0)
      {
        i = i+1; result[i,] = y
        sprint(i, suppressed)
        if (i%%2000==0)
          print(paste(Sys.time(),i,sep = "-----"))
      } else
      {
        sprint("rejected", suppressed)
      }
    }else
    {
      sprint("rejected in AP1", suppressed)
    }
  }
  return(list(result=result,trials = trials, rate = total_samples/trials, time =  proc.time()-ptm))
}