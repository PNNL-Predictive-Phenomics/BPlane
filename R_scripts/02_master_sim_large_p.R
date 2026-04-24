library(BDgraph)
library(MASS)
library(glassoFast)
library(tidyverse)
library(BPlane) #installed locally from C++ files
library(furrr)
options(future.wait.timeout = 1e9)

setwd("path1")

source("R_code/bplane_fcns.R")
source("R_code/sim_fcns.R")
source("R_code/R_fcns-metrics.R")
source("R_code/glasso_fcns.R")

out_dir <- "path2"

rep <- 1:100
n <- 50
p <- c(2000,4000,8000)

for(i in seq_along(rep)){
  print(paste0("rep=",rep[i], " n=", n, " p=", p))
  nothing <- sim_large_p(rep[i], n, p, out_dir)
}
