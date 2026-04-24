#modified from "1.Glasso_functions.R" at 
#https://github.com/lucasvogels33/Review-paper-Bayesian-Structure-Learning-in-GGMs/tree/main/Section%204%20-%20Empirical%20comparison

glasso_1 <- function(Y, G)
# Estimates graph using GLASSO algorithm for one dataset
# Inputs: 
# - Y: n x p dataset
# - G: p x p matrix containing true graph (0s for non-edges, 1s for edges)
# Output: Single row of metrics
{
	n <- nrow(Y)
	p <- ncol(Y)
	nlambda <- 16
	lambda.min.ratio <- .01
	if(p<=25) t0 <- Sys.time()
	else t0 <- proc.time()
	huge_output <- huge(
		x=Y,
		lambda = NULL,
		nlambda = nlambda,
		lambda.min.ratio = lambda.min.ratio,
		method = "glasso",
		verbose = FALSE
	)
	ric_out <- huge.select(huge_output,criterion="ric", verbose=FALSE)
	huge_output_ric <- huge(x=Y,lambda = ric_out$opt.lambda,method="glasso", verbose=FALSE)
	if(p<=25){ 
	  t1 <- Sys.time()
	  comp_time <- as.numeric(t1-t0)
	} else {
	  t1 <- proc.time()
	  comp_time <- (t1-t0)[[3]]
	}
	G_est = huge_output_ric$path[[1]]
	out <- calc_pplus_pminus_glasso(G, G_est)
	conf <- calc_confusion_glasso(G, G_est)
	m <- calc_metrics(conf) #R_fcns-metrics.R
	auc <- (m$sens + m$spec) / 2
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

calc_pplus_pminus_glasso <- function(G, G_est)
{
	predictor <- G_est[upper.tri(G_est)]
	response <- G[upper.tri(G)]
	ones = sum(response)
	zeroes = length(response)-ones
	pplus <- sum(predictor[which(response==1)])/ones
	pminus <- sum(predictor[which(response==0)])/zeroes
	list(p_plus = pplus, p_minus = pminus)
}

calc_confusion_glasso <- function(G, G_est)
{
	conf <- rep(0, 4)
	conf[1] <- (G_est == 1 & G == 1 & upper.tri(G)) |> sum() #tp
	conf[2] <- (G_est == 1 & G == 0 & upper.tri(G)) |> sum() #fp
	conf[3] <- (G_est == 0 & G == 0 & upper.tri(G)) |> sum() #tn
	conf[4] <- (G_est == 0 & G == 1 & upper.tri(G)) |> sum() #fn
	conf
}

calc_lambdas <- function(S, nlambda, lambda.min.ratio)
{
	lambda_max <- S[upper.tri(S)] |> abs() |> max()
	lambda_min <- lambda.min.ratio * lambda_max
	log10_lambda <- seq(from=log10(lambda_min), to=log10(lambda_max), length=nlambda)
	10^log10_lambda #returns lambdas
}