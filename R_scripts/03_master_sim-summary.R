library(tidyverse)
library(cowplot)
#collect files
out_dir_small_p <- "path1"
out_dir_large_p <- "path2"
all_output <- NULL
fig_dir <- "path3"

#load small_p files
all_files <- list.files(out_dir_small_p)
for(file in all_files){
  load(file.path(out_dir_small_p, file))
  all_output <- bind_rows(all_output, output)
}
all_output <- all_output |>
  mutate(
    auc = ifelse(p>=1000, NA, auc),
    pplus = ifelse(p>=1000, NA, pplus),
    pminus = ifelse(p>=1000, NA, pminus),
    sens = ifelse(p>=1000, NA, sens),
    spec = ifelse(p>=1000, NA, spec),
    f1 = ifelse(p>=1000, NA, f1),
    mcc = ifelse(p>=1000, NA, mcc),
  )

#load large_p files
all_files <- list.files(out_dir_large_p)
for(file in all_files){
  load(file.path(out_dir_large_p, file))
  all_output <- bind_rows(all_output, output)
}

#summarize data
reshape <- all_output |>
  pivot_longer(
    cols = c(time, auc, pplus, pminus, sens, spec, f1, mcc),
    names_to = "measure", 
    values_to = "value"
  )
summary <- reshape |> 
  group_by(n, p, method, measure) |>
  summarize(reps = sum(!is.na(value)), mean = mean(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE)) |>
  mutate(
    se = sd / sqrt(reps),
    lower = mean - 2*se,
    upper = mean + 2*se, 
  ) |> 
  filter(!(measure == "time" & reps == 0)) |>
  rename(Algorithm = method) |>
  mutate(
    Algorithm = Algorithm |>
      replace_values(
        "glasso" ~ 'GLASSO', 
        "RJWWA" ~ "WWA"
    ), 
    measure = measure |>
      replace_values(
        "auc" ~ "AUC", 
        "pplus" ~ "Pr+", 
        "pminus" ~ "Pr-",
        "f1" ~ "F1", 
        "mcc" ~ "MCC", 
        "sens" ~ "Sensitivity", 
        "spec" ~ "Specificity"
      )
  )

#time plot
options(scipen=999)
time_plot <- ggplot(
  data = summary |> 
    filter(measure == "time"),
  mapping = aes(x = p, y = mean, color = Algorithm)
) + 
  geom_line(linewidth=.1) + 
  geom_point() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Algorithm), alpha=.1, linewidth=.1) + 
  scale_x_log10(breaks = c(10,25,50,100,200,400,1000,2000)) + 
  scale_y_log10(breaks = c(.01, 0.1, 1, 10, 100, 1000, 10000)) + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank()) + 
  xlab("Nodes (p)") + 
  ylab("Compuation time (sec)")

pdf(file=file.path(fig_dir, "comp_time_rev.pdf"), width=8, height=3.5)
  time_plot
dev.off()

#detection measures
(auc_plot <- summary |> 
  filter(measure == "AUC") |> 
  ggplot(
  mapping = aes(x = p, y = mean, color = Algorithm)
) + 
  geom_line() + geom_point() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Algorithm), alpha=.1, linewidth=.25) + 
  scale_x_log10(breaks = c(10,25,50,100,200,400,1000,2000,4000,8000)) + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        plot.margin = margin_part(b=0, l=0),
        legend.position="none"
  ) + 
  ylab("AUC")
)
(prplus_plot <- summary |> 
    filter(measure == "Pr+") |> 
    ggplot(
      mapping = aes(x = p, y = mean, color = Algorithm)
    ) + 
    geom_line() + geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Algorithm), alpha=.1, linewidth=.25) + 
    scale_x_log10(breaks = c(10,25,50,100,200,400,1000,2000,4000,8000)) + 
    theme_bw() + 
    theme(panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin_part(b=0, l=0), 
          legend.position="none"
    ) + 
    ylab(expression("Pr"^"+"))
)
(prminus_plot <- summary |> 
    filter(measure == "Pr-") |> 
    ggplot(
      mapping = aes(x = p, y = mean, color = Algorithm)
    ) + 
    geom_line() + geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Algorithm), alpha=.1, linewidth=.25) + 
    scale_x_log10(breaks = c(10,25,50,100,200,400,1000,2000,4000,8000)) + 
    theme_bw() + 
    theme(panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin_part(b=0, l=0),
          legend.position = "inside", 
          legend.position.inside = c(.75, .75), 
          legend.background = element_rect(color = "grey50")
    ) + 
    ylab(expression("Pr"^" _")) 
)

edge_det <- plot_grid(auc_plot, prplus_plot, prminus_plot, nrow=1)

pdf(file=file.path(fig_dir, "auc.pdf"), width=8, height=3.5)
  edge_det
dev.off()

(sens_plot <- summary |> 
  filter(measure == "Sensitivity") |>
  ggplot(
    mapping = aes(x = p, y = mean, color = Algorithm)
  ) + 
  geom_line() + geom_point() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Algorithm), alpha=.1, linewidth=.25) + 
  scale_x_log10(breaks = c(10,25,50,100,200,400,1000,2000,4000,8000)) + 
  theme_bw() + 
  ylab("Sensitivity") + 
  theme(
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = "none"
  ))
(spec_plot <- summary |> 
  filter(measure == "Specificity") |>
  ggplot(
    mapping = aes(x = p, y = mean, color = Algorithm)
  ) + 
  geom_line() + geom_point() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Algorithm), alpha=.1, linewidth=.25) + 
  scale_x_log10(breaks = c(10,25,50,100,200,400,1000,2000,4000,8000)) + 
  theme_bw() + 
  ylab("Specificity") + 
  theme(
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = "inside",
    legend.position.inside = c(.8, .45), 
    legend.background = element_rect(color = "grey50")
  ))
(mcc_plot <- summary |> 
    filter(measure == "MCC") |>
    ggplot(
      mapping = aes(x = p, y = mean, color = Algorithm)
    ) + 
    geom_line() + geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Algorithm), alpha=.1, linewidth=.25) + 
    scale_x_log10(breaks = c(10,25,50,100,200,400,1000,2000,4000,8000)) + 
    theme_bw() + 
    ylab("Matthews Correlation Coefficient") + 
    theme(
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
      legend.position = "none"
    ))
(f1_plot <- summary |> 
    filter(measure == "F1") |>
    ggplot(
      mapping = aes(x = p, y = mean, color = Algorithm)
    ) + 
    geom_line() + geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Algorithm), alpha=.1, linewidth=.25) + 
    scale_x_log10(breaks = c(10,25,50,100,200,400,1000,2000,4000,8000)) + 
    theme_bw() + 
    ylab("F1 score") + 
    theme(
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
      legend.position = "none"
    ))

(edge_det_supp <- plot_grid(
  sens_plot, spec_plot, mcc_plot, f1_plot,
  nrow = 2, ncol = 2
))

pdf(file=file.path(fig_dir, "edge_det_supp.pdf"), width=8, height=5)
  edge_det_supp
dev.off()




