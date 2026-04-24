#run within 04_sim_prior.R
sim_prior_1 <- function(rep, n, p, props, cl_sz, out_dir)
{
  L <- sim_graph_large_p(n=n, p=p, cl_sz=cl_sz)
  G <- L$G
  Y <- L$Y
  v0 <- c(.25,.50,.75,1) / sqrt(n * log(p))
  v1 <- c(2.5,5,7.5,10) / sqrt(n * log(p))
  vs <- expand_grid(v0, v1)
  output <- NULL
  file.path(out_dir, paste0("prog_rep=",rep)) |> dir.create()
  for(i in seq_along(props)){
    prior <- make_prior_mat(p, G, props[i])
    bics <- pmap_dfr(
      .l = vs, 
      .f = tune_1_prior,
      prior, n, p, Y
    )
    opt_tune <- bics |>
      slice_min(bic) |>
      dplyr::select(!bic)
    omega <- bplane_1gp_r(
      Y, n, p, v0=opt_tune$v0, v1=opt_tune$v1, 
      tau=.5 / sqrt(n * log(p)), 
      prior=prior, eps=.001, maxiter=100
    )$omega
    out1 <- tibble( #output
      prop = props[i],
      metrics_prior(omega, G, prior)
    )
    output <- bind_rows(output, out1)
    save("done", 
      file=file.path(out_dir, paste0("prog_rep=",rep), paste0("prop=",props[i], ".RData"))
    )
  }
  output <- tibble(rep=rep, n=n, p=p, cl_sz=cl_sz, output)
  save(output,
    file=file.path(out_dir, paste0("sim_prior_rep=",rep,"_n=",n,"_p=",p,".RData"))
  )
}



tune_1_prior <- function(v0, v1, prior, n, p, Y)
{
  Thetahat <- bplane_1gp_r(
    Y, n, p, v0=v0, v1=v1, tau=.5 / sqrt(n * log(p)), 
    prior = prior, eps=.002, maxiter=100)$Theta
  bic <- bic_bplane_1gp(n, p, Y, Thetahat)
  tibble(v0=v0, v1=v1, bic=bic)
}

make_prior_mat <- function(p, G, prop)
{
  if(prop==0){
    prior <- array(0.5, dim=c(p,p))
    return(prior)
  }
  else if(prop==1){
    edg_idx <- which(G==1 & upper.tri(G))
    prior <- array(0.5, dim=c(p,p))
    prior[edg_idx] <- 0.9
    prior[arrayInd(edg_idx, .dim=c(p,p))[,2:1]] <- 0.9
    return(prior)
  }
  else if(prop > 0 & prop < 1){
    edg_idx <- which(G==1 & upper.tri(G))
    if(prop*length(edg_idx) != floor(prop*length(edg_idx)))
      stop(paste("proportion ",prop, "of edges is not an integer"))
    inf_idx <- sample(edg_idx, size=prop*length(edg_idx))
    prior <- array(0.5, dim=c(p,p))
    prior[inf_idx] <- 0.9
    prior[arrayInd(inf_idx, .dim=c(p,p))[,2:1]] <- 0.9
    return(prior)
  }
  else stop("prop must be [0,1]")
}

metrics_prior <- function(omega, G, prior)
{
  predictor = omega[upper.tri(omega)]
  response = G[upper.tri(G)]
  auc <- calc_AUC_ROC(predictor,response)
  out <- calc_pplus_pmin(predictor,response)
  conf <- calc_confusion_prob(G, omega)
  m <- calc_metrics(conf) #R_fcns-metrics
  edg_idx <- which(G==1 & upper.tri(G))
  inf_idx <- which(prior != 0.5 & upper.tri(G))
  noninf_idx <- setdiff(edg_idx, inf_idx)
  tibble(
    auc = auc, 
    pplus = out$p_plus,
    pminus = out$p_minus,
    sens = m$sens, 
    spec = m$spec,
    f1 = m$f1, 
    mcc = m$mcc,
    inf_pplus = mean(omega[inf_idx]),
    noninf_pplus = mean(omega[noninf_idx]),
    inf_sens = mean(omega[inf_idx] > .5),
    noninf_sens = mean(omega[noninf_idx] > .5)
  )
}
