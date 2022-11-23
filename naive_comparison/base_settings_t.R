setwd("../MCF/")
source("./load.R")

setwd("../naive_comparison/")
source("./categorical_distribution.R")


######### student t #################
type="t"

deg = c(2+1e-2,2+1e-1,3)
mu = c(-2, 3, 5)

ncomps = length(deg)

# constraint
consMat = matrix(c(1,1,1), nrow=1)
rhs = 10

# samples
N = 10000

#####################################
