setwd("path")
library(tidyverse)
out <- readRDS("Data/prot_data_allmissing.RDS") # "bc" means batched corrected

4 * 3 * 5 * 6 #total combinations of 4 viruses, 3 donors, 5 replications, 6 time points

edata <- out$e_data
fdata <- out$f_data
#emeta <- out$e_meta
rm(out)

fdata |> distinct(Strain)
fdata |> distinct(Virus)
fdata |> count(Virus, Donor, Replicate)
fdata |> count(Virus, Donor, Replicate, TimePt)
fdata |> count(Plex)
fdata |> count(Virus, TimePt)

# how many samples for each protein?
edata |> 
  pivot_longer(-1) |> 
  rename(SampleID = name) |> 
  full_join(fdata, by = join_by(SampleID)) |> 
  filter(TimePt %in% c("48Hr", "60Hr", "72Hr")) |> 
  filter(Virus != "NL63") |>
  group_by(Protein) |> 
  summarize(N = sum(!is.na(value))) |> 
  count(N) |> print(n=110)

# susbset to 6,968 proteins with all samples at t=48,60, and 72 hr
full_dat <- edata |> 
  pivot_longer(-1) |> 
  rename(SampleID = name) |> 
  full_join(fdata, by = join_by(SampleID)) |> 
  filter(TimePt %in% c("48Hr", "60Hr", "72Hr")) |> 
  filter(Virus != "NL63") |>
  group_by(Protein) |>
  mutate(N = sum(!is.na(value))) |> 
  filter(N >= 129) |> # N = 129 if no samples are missing
  select(!N)

n_by_p <- full_dat |> 
  group_by(Protein,Virus) |> 
  mutate(mean = mean(value), sd = sd(value)) |> 
  mutate(std_ab = (value - mean) / sd) |> 
  select(SampleID,  Virus, Protein, std_ab) |> 
  pivot_wider(id_cols = c(SampleID, Virus), names_from = Protein, values_from = std_ab)

sep_virus_matrix <- function(n_by_p, virus_name)
{
  n_by_p |> 
    filter(Virus == virus_name) |> 
    ungroup() |>
    select(-c(SampleID,Virus)) |> 
    as.matrix()
}

wa <- sep_virus_matrix(n_by_p, "CoV2-WA")
it <- sep_virus_matrix(n_by_p, "CoV2-IT")
mck <- sep_virus_matrix(n_by_p, "MCK")
#nl63 <- sep_virus_matrix(n_by_p, "NL63")

save(wa, it, mck, file="Data/processed.RData")
