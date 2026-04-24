# BPlane code

For "Computationally efficient Bayesian estimation of graphical networks for omics data" article in Journal of Proteome Research. Contact Dan Adrian (adriand1@gvsu.edu) with questions.

All code was implemented in R, with some C++ backend code.

The code is organized in 3 folders:

* `cpp_code`

* `R_fcns`

* `R_scripts`

## Cpp code

The folder contains C++ code files. I had to put them in an installation of a local R package to get them to work properly with the `{furrr}` package (for parallel computation). The files are:

* `bagus.cpp`: contains modified code from [GemBag Github repo](https://github.com/xinming104/GemBag) used to run BAGUS algorithm.

* `bplane.cpp`: contains self-written C++ code for BPlane algorithm.

* `rjwwa.cpp`: contains C++ code copied from `ggm_new5.cpp` from [Bayesian Structure Learning in GGMs Github repo](https://github.com/lucasvogels33/Review-paper-Bayesian-Structure-Learning-in-GGMs) used to run WWA.

## R_fcns

These files are called through the `source()` function by R scripts. Desriptions of the files:

* `bagus_fcns.R`: fits BAGUS algorithm for a single dataset

* `bplane_fcns.R`: fits BPlane algorithm for a single dataset

* `glasso_fcns.R`: fits GLASSO algorithm (via `{huge}` package) for a single dataset

* `R_fcns_metrics.R`: calculates metrics of edge detection performance, including AUC, Pr+, Pr-, sensitivity, specificity, MCC, and F1 score. Adapted from ["metric_functions.R" at Bayesian Structure Learning in GGMs Github repo](https://github.com/lucasvogels33/Review-paper-Bayesian-Structure-Learning-in-GGMs/blob/main/Section%204%20-%20Empirical%20comparison/metric_functions.R).

* `RJWWA_fcns.R`: fits WWA algorithm for a single dataset. Adapted from [RJWWA_functions.R at Bayesian Structure Learning in GGMs Github repo](https://github.com/lucasvogels33/Review-paper-Bayesian-Structure-Learning-in-GGMs/blob/main/Section%204%20-%20Empirical%20comparison/2.RJWWA_functions.R).

* `sim_fcns.R`: contains functions that simulate graphs and datasets via `{BDgraph}` package. Functions used in two types of simulations: (1) "small p" simulations which are done on a single processor within `future_map()` functions; (2) "large p" simulations that use multiple processors for tuning job for the same data.

* `sim_prior_fcns.R`: functions used in simulations with prior information

## R_scripts

These R script files are ordered numerically by their order of mention in the journal article. These utilize parallel computing through the `{furrr}` package. Because of this, the approach is to take the output from each simulated dataset and write it to a separate file.

* `01_master_simulations.R`: Performs "small p" simulations (p < 1000)

* `02_master_sim_large_p.R`: Performs "large p" simulations (p >= 1000)

* `03_master_sim-summary.R`: Collect multiple files from simulations and produces Figures 1 and 2 and supp. Figure 2.

* `04_sim_prior.R`: Performs simulation study applying different proportions of informative prior information to true edges.

* `05_sim_prior-summary.R`: Collects files from simulation and produces Table 1.

* `06_data_start.R`: Wrangles raw SARS-Cov2 dataset into usable form

* `07_data_detection_parallel.R`: Estimates graph for SARS-Cov2 dataset using parallel computing (1 core per strain)

* `08_data_upset_edges.R`: Produces Figure 3 UpSet plot using `{ggupset}` package

* `09_sim_bootstrap.R`: Estimates graphs for bootstrap samples from original data for stability analysis


