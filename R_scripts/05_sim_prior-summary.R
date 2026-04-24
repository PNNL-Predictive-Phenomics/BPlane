library(tidyverse)
out_dir <- "C:/Users/adriand1/OneDrive - Grand Valley State University/Desktop/Sabbatical/code/results/sim_prior"
all_output <- NULL

all_files <- list.files(out_dir)
all_files <- all_files[all_files!="trial"]
for(file in all_files){
  load(file.path(out_dir, file))
  all_output <- bind_rows(all_output, output)
}

reshape <- all_output |> 
  select(!c(n,p,cl_sz)) |>
  pivot_longer(
    cols = c(auc, pplus, pminus, sens, spec, f1, mcc,inf_pplus,noninf_pplus,inf_sens,noninf_sens),
    names_to = "measure", 
    values_to = "value"
  ) |> 
  mutate(prop = factor(prop, levels = c(0, .25, .5, .75, 1)))

ggplot(
  data = reshape,
  mapping = aes(x = value, color = prop, fill=prop)
) + 
  geom_density(alpha=.1) + 
  facet_wrap(~measure, scales="free") + 
  theme_bw()

reshape |>
  group_by(measure, prop) |>
  summarize(
    mean = mean(value), sd = sd(value), se = sd / sqrt(n())
  ) |> 
  print(n=Inf)
  