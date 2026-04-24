library(BDgraph)
library(MASS)
library(tidyverse)
library(BPlane) #installed locally from C++ files
library(furrr)
options(future.wait.timeout = 1e9)

setwd("path)

source("R_code/bplane_fcns.R")
source("R_code/sim_fcns.R")
source("R_code/R_fcns-metrics.R")
source("R_code/sim_prior_fcns.R")

out_dir <- "path"
n <- 50
p <- 2000
props <- c(0, .25, .50, .75, 1)
cl_sz <- 500
rep <- 1:500

plan(multisession, workers=5)
gc()


future_walk(
  .x = rep, 
  .f = sim_prior_1,
  n, p, props, cl_sz, out_dir,
  .options = furrr_options(seed=TRUE, scheduling=Inf)
)
