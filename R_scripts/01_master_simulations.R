library(igraph)
library(BDgraph)
library(BH)         #for RJWWA algorithm
library(Rcpp)       #for RJWWA algorithm
library(RcppBlaze)
library(MASS)
library(huge)
library(tidyverse)
library(BPlane) #installed locally from C++ files
library(furrr)
library(tictoc)
options(future.wait.timeout = 1e9)

setwd("path")

source("R_code/RJWWA_fcns.R")
source("R_code/bplane_fcns.R")
source("R_code/sim_fcns.R")
source("R_code/R_fcns-metrics.R")
source("R_code/glasso_fcns.R")
source("R_code/bagus_fcns.R")


rep <- 1:100
n <- c(50)
p <- c(10,25,50,100,200,400)
rjwwa_plim <- 101 #WWA method is only run for p values below this
bagus_plim <- 401 #BAGUS method is only run for p values below this
out_dir <- "path" 
# the output (one row per method, containing metrics) from each simulation 
# is written to a separate file labeled by rep, n, and p

sim_settings <- expand_grid(rep=rep, n = n, p = p)

plan(strategy = multisession, workers=12)
tic()
nothing <- furrr::future_pmap(
  .l = sim_settings, 
  .f = safely(sim_1),
  rjwwa_plim,
  bagus_plim,
  out_dir = out_dir,
  .progress = TRUE, 
  .options =furrr_options(seed=TRUE)
)
toc()
plan(strategy=sequential)
