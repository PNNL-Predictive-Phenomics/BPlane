# run within furrr::future_pmap() on 01_master_simulations.ratio
# for simulations with p < 1000
sim_1 <- function(rep, n, p, rjwwa_plim, bagus_plim, out_dir){
  runif(n=1, min=0, max=5) |> Sys.sleep()
  if(
	!file.exists(file.path(out_dir, paste0("sim_rep=",rep, "_n=",n, "_p=",p,".RData")))
  ){
	print(paste0("rep=",rep, " n=", n, " p=", p))
	n_iter <- 30000
	out_rjwwa <- out_glasso <- out_bplane <- out_bagus <- NULL
	rjwwa <- glasso <- bplane <- bagus <- TRUE
	if(p > rjwwa_plim){
		rjwwa <- FALSE
		if(p > bagus_plim){
			bagus <- FALSE
		}
	}
	L <- bdgraph.sim(p=p, graph="scale-free", n=n, size=p-1)
	Y <- L$data
	G <- L$G
	if(rjwwa){
		out_rjwwa <- tibble(method = "RJWWA", RJWWA_1(Y, G, n_iter=30000))
	}
	if(glasso){
		out_glasso <- tibble(method = "glasso", glasso_1(Y, G))
	}
	if(bplane){
		out_bplane <- tibble(method = "BPlane", bplane_1(Y, G))
	}
	if(bagus){
		out_bagus <- tibble(method = "BAGUS", bagus_1(Y, G))
	}
	res <- bind_rows(out_rjwwa, out_glasso, out_bplane, out_bagus)
	output <- tibble(rep=rep, n=n, p=p, res)
 	save(
 		output,
 		file=file.path(out_dir, paste0("sim_rep=",rep, "_n=",n, "_p=",p,".RData"))
 	)
  }
}

# run within 02_master_sim_large_p.R
sim_large_p <- function(rep, n, p, out_dir)
{
	if(
		!file.exists(file.path(out_dir, paste0("sim_rep=",rep, "_n=",n, "_p=",p,".RData")))
	){
		L <- sim_graph_large_p(n, p, cl_sz=500)
		res <- tuning_large_p(L$Y, L$G)
		output <- tibble(rep=rep, n=n, p=p, res)
		save(
			output,
			file=file.path(out_dir, paste0("sim_rep=",rep, "_n=",n, "_p=",p,".RData"))
		)
	}
	rm(L,res,output)
	gc()
}

sim_graph_large_p <- function(n, p, cl_sz)
# simulates large graphs in multiple clusters 
# - because bdgraph.sim() is very slow for p>=1000
{
  if(p %% cl_sz != 0) stop("p must be a multiple of cl_sz!")
  n_cl <- as.integer(p/cl_sz)
  G <- sigma <- array(0,dim=c(p,p))
  for(i in 1:n_cl){
    idx <- (i-1)*cl_sz + 1:cl_sz
    L <- bdgraph.sim(p=cl_sz, n=n, graph="scale-free", vis=FALSE)
    G[idx,idx] <- L$G
    sigma[idx,idx] <- L$sigma
  }
  rm(L,n_cl,idx)
  gc()
  list(
    G=G, 
    Y = mvrnorm(n=n, mu=rep(0,p), Sigma=sigma)
  )
}

tuning_large_p <- function(Y, G)
# tuning for both GLASSO and BPlane methods using parallel processing
{
	n <- nrow(Y)
	p <- ncol(Y)
	S <- cov(Y)
	lambdas <- calc_lambdas(S, nlambda=16, lambda.min.ratio=.01) #glasso_fcns.R
	v0 <- c(.25,.50,.75,1) / sqrt(n * log(p))
	v1 <- c(2.5,5,7.5,10) / sqrt(n * log(p))
	vs <- expand_grid(v0, v1)
	input_tuning <- bind_rows(
		tibble(method="glasso", tune_par1=lambdas, tune_par2=NA),
		tibble(method="BPlane", tune_par1=vs$v0, tune_par2=vs$v1)
	)
	rm(S,lambdas,v0,v1,vs)
	plan(sequential)
	gc()
	plan(strategy = multisession, workers=6)
	#plan(sequential)
	bics <- future_pmap_dfr(
		.l = input_tuning, 
		.f = tuning_1, 
		n, p, Y, 
		.progress = TRUE,
		.options =furrr_options(
			scheduling=Inf,seed=TRUE,
			packages=c("glassoFast", "BPlane")
		)
	)
	opt_tune <- bics |>
		group_by(method) |>
		slice_min(bic) |>
		dplyr::select(!bic)
	plan(strategy = multisession, workers=2)
	#plan(sequential)
	future_pmap_dfr(
		.l = opt_tune, 
		.f = fit_1, 
		G, Y, n, p,  
		.progress = TRUE,
		.options =furrr_options(scheduling=Inf,seed=TRUE)
	)
}

tuning_1 <- function(method, tune_par1, tune_par2, n, p, Y)
# evaluates for one set of parameters
# used in future_pmap_dfr() function within tuning_large_p() function above
{
	Sys.sleep(runif(1, 0, 2))
	t0 <- proc.time()
	if(method=="glasso"){
		Sys.setenv(OMP_NUM_THREADS = "1")
		Sys.setenv(MKL_NUM_THREADS = "1")
		Thetahat <- glassoFast(S=cov(Y), rho=tune_par1)$wi
		bic <- bic_bplane_1gp(n, p, Y, Thetahat)
		return(
			tibble(method=method, tune_par1=tune_par1, tune_par2=NA, bic=bic)
		)
		t1 <- proc.time()
		comp_time <- (t1-t0)[[3]]
		print(paste0(method, ", par1=", tune_par1, " par2=", tune_par2, " time=", comp_time))
	}
	else if(method=="BPlane"){
		Thetahat <- bplane_1gp_r(
			Y, n, p, v0=tune_par1, v1=tune_par2, tau=.5 / sqrt(n * log(p)), 
			prior = array(.5, dim=c(p,p)), 
			eps=.002, maxiter=100)$Theta
		bic <- bic_bplane_1gp(n, p, Y, Thetahat)
		return(
			tibble(method=method, tune_par1=tune_par1, tune_par2=tune_par2, bic=bic)
		)
		t1 <- proc.time()
		comp_time <- (t1-t0)[[3]]
		print(paste0(method, ", par1=", tune_par1, " par2=", tune_par2, " time=", comp_time))
	}
	else stop("method is not glasso or bplane")
	rm(Thetahat)
	gc()
}

fit_1 <- function(method, tune_par1, tune_par2, G, Y, n, p)
# fits based on one method (used in parallel processing)
{
	if(method=="glasso"){
		Thetahat <- glassoFast(S=cov(Y), rho=tune_par1)$wi
		G_est = ifelse(Thetahat==0, 0, 1)
		out <- calc_pplus_pminus_glasso(G, G_est)
		conf <- calc_confusion_glasso(G, G_est)
		m <- calc_metrics(conf) #R_fcns-metrics.R
		auc <- (m$sens + m$spec) / 2
	}
	else if(method=="BPlane"){
		omega <- bplane_1gp_r(
			Y, n, p, v0=tune_par1, v1=tune_par2, 
			tau=.5 / sqrt(n * log(p)), 
			prior=array(.5, dim=c(p,p)), eps=.001, maxiter=100
			)$omega
		predictor = omega[upper.tri(omega)]
		response = G[upper.tri(G)]
		auc <- calc_AUC_ROC(predictor,response)
		out <- calc_pplus_pmin(predictor,response)
		conf <- calc_confusion_prob(G, omega)
		m <- calc_metrics(conf) #R_fcns-metrics
	}
	else stop("method is not glasso or bplane")
	tibble( #output
		method=method,
		time = NA,
		auc = auc, 
		pplus = out$p_plus,
		pminus = out$p_minus,
		sens = m$sens, 
		spec = m$spec, 
		f1 = m$f1, 
		mcc = m$mcc
	)
}
