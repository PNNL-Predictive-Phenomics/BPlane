library(tidyverse)
setwd("path1")
library(BPlane)
library(tictoc)
library(furrr)
library(igraph)
source("code/R_code/bplane_fcns.R")
load("data/Data/processed.RData")
load("code/results/data_results/opt_hyp.RData")
out_dir <- "path2"

dataname <- c("it", "mck", "wa")
rep <- 1:100

input <- expand_grid(rep, dataname)
gc()

bootstrap_1 <- function(rep, dataname, opt_hyp, out_dir, it, mck, wa)
{
  runif(n=1,min=.01, max=.5) |> Sys.sleep()
  if(
    !file.exists(file.path(out_dir, paste0(dataname, rep, ".RData")))
  ){
    if(dataname=="it") Y <- it
    else if(dataname=="mck") Y <- mck
    else if(dataname=="wa") Y <- wa
    n <- nrow(Y)
    p <- ncol(Y)
    idx <- sample.int(n, size=n, replace=TRUE)
    boot_data <- Y[idx,]
    dn <- dataname
    rm(it,mck,wa,Y,idx)
    gc()
    omega <- bplane_1gp_r(
      boot_data, n, p, 
      v0=filter(opt_hyp, dataname==dn)$v0, 
      v1=filter(opt_hyp, dataname==dn)$v1, 
      tau=.5 / sqrt(n * log(p)), 
      prior=array(.5, dim=c(p,p)), 
      eps=.001, maxiter=100)$omega
    adj_mat <- ifelse(omega > .5, 1, 0)
    graph <- igraph::graph_from_adjacency_matrix(adj_mat, mode="undirected", diag=FALSE)
    output <- list(name=dataname, graph=graph)
    save(output, file=file.path(out_dir, paste0(dataname, rep, ".RData")))
    rm(boot_data,omega,adj_mat,graph,output)
    gc()
  }
}

plan(strategy = multisession, workers=10)
nothing <- future_pmap(
  .l = input,
  .f = bootstrap_1, 
  opt_hyp, out_dir, it, mck, wa,
  .progress = TRUE, 
  .options = furrr_options(seed=TRUE)
)

# add labels
gp <- c("it", "mck", "wa")
for(i in 1:100){
  for(j in seq_along(gp)){
    load(file.path(out_dir, paste0(gp[j], i, ".RData")))
    output$graph <- set_vertex_attr(output$graph, "label", value = colnames(get(gp[j])))
    save(output, file=file.path(out_dir, paste0(gp[j], i, ".RData")))
  }
}

