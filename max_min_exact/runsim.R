#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

source('./load.R')

setwd("../imputation_comparison/")
source("./glg_regresssion.R")
source("./categorical_distribution.R")

setwd("../max_min_exact/")
source("./imputation.R")

if (length(args)==0)
{
  test_days=7
} else
{
  test_days=as.integer(args[1])
}
source('./main.R')
