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

detect_bplane_data_1 <- function(dataname, v0, v1)
{
  load("C:/Users/adriand1/OneDrive - Grand Valley State University/Desktop/Sabbatical/data/Data/processed.RData")
  Y <- get(dataname)
  n <- nrow(Y)
  p <- ncol(Y)
  omega <- bplane_1gp_r(Y, n, p, v0, v1, 
                           tau=.5 / sqrt(n * log(p)), prior=array(.5, dim=c(p,p)), 
                           eps=.001, maxiter=500)$omega
  adj_mat <- ifelse(omega > .5, 1, 0)
  graph <- igraph::graph_from_adjacency_matrix(adj_mat, mode="undirected", diag=FALSE)
  probs <- list(dataname, omega)
  output <- list(dataname, graph)
  save(probs, file=file.path(out_dir, paste0("probs_", dataname, ".RData")))
  output
}

plan(strategy = multisession, workers=3)
graphs <- furrr::future_pmap(
  .l = opt_hyp, 
  .f = detect_bplane_data_1,
  .progress = TRUE, 
  .options = furrr_options(seed=TRUE)
)
plan(sequential)

it_graph <- graphs[[1]][[2]]
mck_graph <- graphs[[2]][[2]]
wa_graph <- graphs[[3]][[2]]

#add protein labels
it_graph <- set_vertex_attr(it_graph, "label", value = colnames(it))
mck_graph <- set_vertex_attr(mck_graph, "label", value = colnames(mck))
wa_graph <- set_vertex_attr(wa_graph, "label", value = colnames(wa))

save(it_graph, mck_graph, wa_graph, file=file.path(out_dir, "covid_graphs_detected.RData"))


