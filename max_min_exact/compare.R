source('benchmark.R')

# load ts
lpday = 24*10 # number of entries per day
test_days = 7
total_rows = dim(results[[1]]$result)[1]

# load observation constraints
obs_raw = read.csv('./observation_subset.csv', header=T)
obs_length = nrow(obs_raw)
constraints = obs_raw[(obs_length-test_days*24*2+1):obs_length, 2]
baseline = constraints[1:((total_rows-5)*48)]


# extract results
fitted_quantile = matrix(nrow = 2, ncol=48*(total_rows-5))
for (row in 6:total_rows)
{
  for (i in 1:48)
  {
    fitted_quantile[1,i+(row-6)*48] = min(apply(results[[i]]$result[row,,], 1, mean))
    fitted_quantile[2,i+(row-6)*48] = max(apply(results[[i]]$result[row,,], 1, mean))
  }
}


# read true data
gt = read.csv('./target_subset.csv', header=T)
discarded_idx = 1:(nrow(gt) - test_days*2*24)
gt = data.matrix(gt[-discarded_idx,])
gt = gt[1:((total_rows-5)*48),2:3]

merged_max_min = cbind(fitted_quantile[1,], fitted_quantile[2,], gt[1:((total_rows-5)*48),], baseline)

min_diff = fitted_quantile[1,] - gt[,2]
max_diff = fitted_quantile[2,] - gt[,1]

min_diff_p = min_diff/gt[,2]*100
max_diff_p = max_diff/gt[,1]*100

plot(min_diff_p, max_diff_p)

cat("Average max difference percentage = ", mean(abs(max_diff_p)), '\n')

# benchmark : use observation as both max and min estimation
# compute performance as ratio  between benchmark and current model
base_score = benchmark(gt, cbind(baseline,baseline))
my_score = benchmark(gt, t(fitted_quantile[c(2,1),]))

print(my_score/base_score)
