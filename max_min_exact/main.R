# try higher order models
# look at holt-winter model
# https://www.analyticsvidhya.com/blog/2021/08/holt-winters-method-for-time-series-analysis/
# Forecasting seasonals and trends by exponentially weighted moving averages
#
# Gamma AR basics
# Product Autoregression: A Time-Series Characterization of the Gamma Distribution



# Turn high-freqeuncy real data into separate time series
interval_length = 6
lpday = 24*60/interval_length # number of entries per day

raw_data = read.csv(paste('./',interval_length,'min_subset.csv', sep=''), header=T)
data_mat = matrix(raw_data$value, nrow=lpday)
weather_data = read.csv('./weather_mh_subset.csv', header=T)
weather_data$temperature = scale(weather_data$temperature)
weather_repeat = weather_data[rep(seq_len(nrow(weather_data)), each=60/interval_length),]

scaling = 20 # scaling to make fitting more stable
total_length = ncol(data_mat)
train_index = 1:(total_length-test_days)
training_set = data_mat[,train_index]*scaling
training_length = ncol(training_set)

temperature_mat = matrix(weather_repeat[,'temperature'], nrow=lpday)
humidity_mat = matrix(weather_repeat[, 'spec_humidity'], nrow=lpday)

# fit Gamma AR on every series
order = 5
intercept = c()
phi = matrix(NA, nrow=order, ncol=lpday)
beta = matrix(NA, nrow=2, ncol=lpday)
glg_params = list()

for(i in 1:lpday)
{
  result = glg_ar(training_set[i,], order = order,
                    xreg = t(rbind(temperature_mat[i,train_index], humidity_mat[i,train_index])))
  intercept[i] = result$intercept
  beta[,i] = result$beta
  phi[,i] = result$alpha
  glg_params[[i]] = result$glg_param
}

# load observation constraints
obs_raw = read.csv('./observation_subset.csv', header=T)
obs_length = nrow(obs_raw)
constraints = matrix(obs_raw[(obs_length-test_days*24*2+1):obs_length, 2], nrow = 48)

# imputation
total_samples = 10000
per_steps = 5
time_pool = c(0.01, 0.1, 0.5, 1, 1.5, 2)
#time_pool = c(1,2,4,10,20)
consMat = array(1, dim = c(1,per_steps))/per_steps

total_steps=1
results = list()
for (idx in 1:48)
{
  cat('=============Slot ', idx, '===============\n')
  slots = (1:per_steps)+(idx-1)*per_steps
  parameters = list(alpha=phi[, slots], glg_params=glg_params[slots], intercept=intercept[slots], beta=beta[, slots],
                    order=order, consMat=consMat, consY=constraints[idx, ]*scaling, time_pool=c(1))
  imputed = impute_energy(total_samples, total_steps, per_steps,
                         t(training_set[slots, (training_length-order+1):training_length]), 
                         parameters, xreg = rbind(temperature_mat[i,-train_index], humidity_mat[i,-train_index]),
                         max_trial=1000, positive=TRUE)
  imputed$result =imputed$result/scaling
  results[[idx]] = imputed
}


save(results, file=paste('result_n', total_samples, '_d', total_steps,  
                  '_k',per_steps,'_', as.character(Sys.Date()),'.RData', sep=''))

