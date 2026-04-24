setwd("path1")
library(ggupset)
library(tidyverse)
library(igraph)
load("code/results/data_results/covid_graphs_detected.RData")
fig_dir <- "path2"

edge_num <- function(x){
  i <- x[1]; j <- x[2]; p <- 6968
  (i-1)*(p-i/2) + (j-i)
}
tot_edg <- 6968*6987/2

wa_edg <- as_edgelist(wa_graph) |> apply(1, edge_num) 
it_edg <- as_edgelist(it_graph) |> apply(1, edge_num) 
mk_edg <- as_edgelist(mck_graph) |> apply(1, edge_num)


edges_tf <- array(FALSE, dim = c(tot_edg, 3))
wa_idx <- it_idx <- mk_idx <- 1
for(i in 1:tot_edg){
  if(i == wa_edg[wa_idx] & wa_idx <=length(wa_edg)){
    edges_tf[i,1] <- TRUE
    wa_idx <- wa_idx + 1
  }
  if(i == it_edg[it_idx] & it_idx <=length(it_edg)){
    edges_tf[i,2] <- TRUE
    it_idx <- it_idx + 1
  }
  if(i == mk_edg[mk_idx] & mk_idx <=length(mk_edg)){
    edges_tf[i,3] <- TRUE
    mk_idx <- mk_idx + 1
  }
  if(i %% 1e6 == 0) print(i)
}
colnames(edges_tf) <- c("Washington", "Italy", "Mock")
graph_dat <- edges_tf |>
  as_tibble(rownames = "Edge") |>
  pivot_longer(cols = c(Washington, Italy, Mock), names_to = "Strain", values_to = "Member") |>
  filter(Member) |>
  select(-Member) |> 
  group_by(Edge) |>
  summarize(Strains = list(Strain))

(upset_plot <- graph_dat |>
  ggplot(aes(x = Strains)) + 
  geom_bar() + 
  scale_x_upset() + 
  ylim(c(NA,22000)) + 
  geom_text(stat = "count", aes(label = after_stat(count)), vjust=-0.5) + 
  ylab("Number of Edges") + 
  xlab("SARS-Cov2 Strain (Combination)") + 
  theme_bw()
)
pdf(file=file.path(fig_dir, "upset.pdf"), width=8, height=3.5)
  upset_plot
dev.off()
