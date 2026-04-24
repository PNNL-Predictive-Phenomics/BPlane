bplane_1 <- function(Y, G)
# Estimates graph using BPlane algorithm for one dataset
# Inputs: 
# - Y: n x p dataset
# - G: p x p matrix containing true graph (0s for non-edges, 1s for edges)
# Output: Single row of metrics
{
	n <- nrow(Y)
	p <- ncol(Y)
	if(p<=25) t0 <- Sys.time()
	else t0 <- proc.time()
	vs <- bplane_sel_hyp(Y)
    omega <- bplane_1gp_r(
      Y, n, p, v0=vs[1], v1=vs[2], tau=.5 / sqrt(n * log(p)), 
	  prior=array(.5, dim=c(p,p)), eps=.001, maxiter=100
    )$omega
	if(p<=25){ 
	  t1 <- Sys.time()
	  comp_time <- as.numeric(t1-t0)
	} else {
	  t1 <- proc.time()
	  comp_time <- (t1-t0)[[3]]
	}
	predictor = omega[upper.tri(omega)]
	response = G[upper.tri(G)]
	auc <- calc_AUC_ROC(predictor,response) #R_fcns-metrics
	out <- calc_pplus_pmin(predictor,response) #R_fcns-metrics
	conf <- calc_confusion_prob(G, omega) #R_fcns-metrics
	m <- calc_metrics(conf) #R_fcns-metrics
	tibble( #output
		time = comp_time,
		auc = auc,
		pplus = out$p_plus,
		pminus = out$p_minus,
		sens = m$sens, 
		spec = m$spec, 
		f1 = m$f1, 
		mcc = m$mcc
	)
}

bplane_sel_hyp <- function(Y)
{
	eps <- .002
	maxiter <- 100
	n <- nrow(Y)
	p <- ncol(Y)
	tau <- .5 / sqrt(n * log(p))
	v0 <- c(.25,.50,.75,1) / sqrt(n * log(p))
	v1 <- c(2.5,5,7.5,10) / sqrt(n * log(p))
	bic <- array(dim = c(length(v0), length(v1)))
	prior <- array(.5, dim=c(p,p))
	for(i in seq_along(v0)){
        for(j in seq_along(v1)){
			Thetahat <- bplane_1gp_r(Y, n, p, v0[i], v1[j], tau, prior, eps, maxiter)$Theta
			bic[i,j] <- bic_bplane_1gp(n, p, Y, Thetahat)
		}
	}
	min_idx <- which(bic == min(bic, na.rm=T), arr.ind=TRUE)
	c(v0[min_idx[1]], v1[min_idx[2]])
}

# interface with bplane.cpp for fitting BPlane algorithm
bplane_1gp_r <- function(Y, n, p, v0, v1, tau, prior, eps, maxiter){
  # Y is a matrix, prior is a matrix
  out <- bplane_r_wrapper( #Rcpp function
    n_vec_r = n, p, v0, v1, tau, init_iter = 0,
    max_iter = maxiter, omega_eps = eps, Y_list = list(Y),
    p1_mat = array(1, dim=c(p,p)), p2_list = list(prior)
  )
  list(
    Theta = out[[1]][[1]], 
    omega = out[[2]][[1]]
  )
}

# interface with bplane.cpp for calculating BIC-type measure
bic_bplane_1gp <- function(n, p, Y, Theta){
# Y and Theta are matrices
  bplane_bic(
    K=1, n_vec_r = n, p, Y_list = list(Y), Theta_list = list(Theta)
  )[[1]]
}

