# Adapts code from RJWWA_functions.R at 
# https://github.com/lucasvogels33/Review-paper-Bayesian-Structure-Learning-in-GGMs/blob/main/Section%204%20-%20Empirical%20comparison/2.RJWWA_functions.R

RJWWA_1 <- function(Y, G, n_iter)
# Estimates graph using BPlane algorithm for one dataset
# Inputs: 
# - Y: n x p dataset
# - G: p x p matrix containing true graph (0s for non-edges, 1s for edges)
# Output: Single row of metrics
{
	n = nrow(Y)
	p = ncol(Y)
	adj <- matrix(0L, nrow = p, ncol = p) #initial graph of the Markov chain
	adj_cum <- matrix(0L, nrow = p, ncol = p)
	#intialize paramters
	edge_prob = 0.2 #prior edge inclusion prob
	df_0 = 3.0 #degree for G-wishart distribution
	U <- t(Y) %*% Y
	t0 <- proc.time()
	for (s in 1:n_iter) {
		adj = update_G(adj = adj, edge_prob = edge_prob, df_0 = df_0, U = U, n = n)
		adj_cum = adj_cum + adj
	}
	t1 <- proc.time()
	comp_time <- (t1-t0)[[3]]
	plinks <- adj_cum / n_iter #output edge probabilities
	predictor = plinks[upper.tri(plinks)]
	response = G[upper.tri(G)]
	auc <- calc_AUC_ROC(predictor,response)
	out <- calc_pplus_pmin(predictor,response)
	conf <- calc_confusion_prob(G, plinks)
	m <- calc_metrics(conf)
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


############################
###supporting  functions#######
############################

#this function produces an update of the graph using the WWA algorithm
update_G <- function(adj, edge_prob, df_0, U, n) {
  p <- nrow(adj)
  
  return(update_G_Rcpp(
    adj = adj,
    edge_prob_mat = matrix(edge_prob, nrow = p, ncol = p),
    df = df_0 + n,
    df_0 = df_0,
    rate = diag(p) + U,
    n_edge = p,  # Number of single edge updates attempted in this MCMC update
    seed = sample.int(n = .Machine$integer.max, size = 1),
    loc_bal = FALSE
  ))
}

# this function samples from the G-Wishart distribution
rgwish <- function(adj, df = 3, rate = NULL) {
  p <- nrow(adj)
  if (is.null(rate)) rate <- diag(p)
  
  return(rgwish_Rcpp(
    adj = adj,
    df = df,
    rate = rate,
    seed = sample.int(n = .Machine$integer.max, size = 1)
  ))
}