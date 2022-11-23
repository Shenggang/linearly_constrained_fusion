# run data_extraction.R to preprocess the dataset
load("data_true.RData")

# define data
n = 30
L = 20
order = 7
Y = consumption_array[1:n, 1:L, ]
X = survey_50[1:n,]

# fit order k linear models

intercept = c()
beta = array(NaN, dim = c(8,3))
alpha = array(NaN, dim = c(order,3))
glg_params = list()

for (i in 1:3)
{
  #print(c(k,i)) #for each period in day
  result = glg_ar(Y[,,i], order = order,
                    xreg = X[rep(1:n, each=L), ])
  intercept[i] = result$intercept
  beta[,i] = result$beta
  alpha[,i] = result$alpha
  glg_params[[i]] = result$glg_param
}