#include "my_project_headers.h"
#include <math.h>
#include <cstring>
#include <R.h>
#include <Rinternals.h>
using namespace std;
using namespace arma;
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

// Soft thresholding function
int sgn(double v) {
  return (v < 0) ? -1 : 1;
}

double SoftThreshold(double a,double b,double lambda, double c = 0){
  return -c+sgn(c-b/a)*fmax(abs(c-b/a)-lambda/a,0.0);
}

arma::Mat<double> Cor_desc(
    arma::Mat<double> Theta_11_inv,
    double W_22,
    arma::Mat<double> Theta_12,
    arma::Mat<double> S_12,
    arma::Mat<double> v_0_12,
    arma::Mat<double> v_1_12,
    arma::Mat<double> P_12,
    int n
){
  int p = Theta_11_inv.n_rows;
  Mat<double> Theta_12_2 = Theta_12;
  for( int iter = 0; iter < 1000; iter++){
    Mat<double> Theta_12_1 = Theta_12;
    for(int i = 0; i<p ; i++){
      double b = n*(as_scalar(Theta_11_inv.row(i)*Theta_12)*W_22- \
                    (Theta_11_inv(i,i)*Theta_12(i))*W_22+S_12(i));
      double a = n*Theta_11_inv(i,i)*W_22;
      double lambda = P_12(i)/v_1_12(i) + (1-P_12(i))/v_0_12(i);
      Theta_12(i) = SoftThreshold(a,b,lambda);
    }

    if((abs(Theta_12_1 - Theta_12)).max() < 0.001) return(Theta_12);
    if(abs(Theta_12).max() > 10) return(Theta_12_2);
  }
  return Theta_12_2;
}

double posterior_fcn( //used for monitoring convergence of the EM algorithm
    int k,
    int p,
    NumericVector n_l,
    List Theta_l,
    List S_l,
    double tau,
    Mat<double> p_1,
    List p_2,
    arma::mat v_1, // TODO: matricize
    arma::mat v_0 // TODO: matricize
){
  double LL = 0, pen_L1 = 0, pen_GB = 0;
  double p2 = 0, theta = 0;
  double v_0_scalar = 0, v_1_scalar = 0;
  double temp1 = 1.0, temp2 = 1.0; //first and second terms of Pen_GB products
  mat temp_mat = ones(p, p);
  for(int i = 0; i < k; i++){
    int n = n_l[i];
    mat Theta = as<arma::mat>(Theta_l[i]);
    mat S = as<arma::mat>(S_l[i]);
    LL = LL + ((double(n) / 2) * (real(arma::log_det(Theta)) - (arma::trace(S * Theta))));  //n_k * (-logdet(Theta_k) + tr(S_k %*% Theta_k))
    for(int j=0; j<p; ++j)
      pen_L1 = pen_L1 + Theta(j,j);
  }
  for(int i = 0; i < p-1; i++){
    for(int j = i+1; j < p; j++){
      temp1 = temp2 = 1.0;
      for(int l = 0; l < k; l++){
        p2 = as<arma::mat>(p_2[l])(i,j);
        theta = as<arma::mat>(Theta_l[l])(i,j);
        v_0_scalar = v_0(i,j);
        v_1_scalar = v_1(i,j);
        temp1 = temp1 * (p2 / (2.0*v_1_scalar) * exp(-abs(theta)/v_1_scalar) + (1 - p2) / (2.0*v_0_scalar) * exp(-abs(theta)/v_0_scalar)); // no problem with matricizing v_0 and v_1 here
        temp2 = temp2 * 1.0 / (2.0*v_0_scalar) * exp(-abs(theta)/v_0_scalar); // no problem with matricizing v_0 and v_1 here
      }
      temp_mat(i,j) = -1.0*log(p_1(i,j) * temp1 + (1 - p_1(i,j)) * temp2);
      pen_GB = pen_GB + temp_mat(i,j);
    }
  }

  return -1.0 * LL + tau * pen_L1 + pen_GB;
}



List calculate_p_l(
    List Theta_l,
    arma::Mat<double> P_1_mat,
    List P_2_mat,
    arma::Mat<double> v_0,
    arma::Mat<double> v_1,
    double p,
    int k,
    bool should_use_p1
) {
  mat p2k = ones(p, p); //matrix of prior probs p2 for each group
  mat shrk = ones(p, p);
  // Calculate P_l
  if(should_use_p1) {
    mat tmp = ones(p, p);
    for(int i = 0; i < k; i++) {
      p2k = as<arma::mat>(P_2_mat[i]);
      tmp = tmp % ((p2k % v_0) / v_1) % exp(abs(as<arma::mat>(Theta_l[i])) % ((1.0 / v_0) - (1.0 / v_1))) + (1 - p2k);
    }
    mat shrk = 1.0 / ( 1.0 + (1.0 - P_1_mat) / P_1_mat / tmp);
  }


  List P_l(k);
  for(int i = 0; i < k; i++) {
    p2k = as<arma::mat>(P_2_mat[i]);
    mat term1 = (v_1/v_0) % exp(-abs(as<arma::mat>(Theta_l[i])) / v_0 + abs(as<arma::mat>(Theta_l[i])) / v_1);
    mat term2 = ((1-p2k)/p2k);

    P_l[i] = shrk % (1.0 / (1.0 + term1 % term2));
  }

  return P_l;
}

double max_diff(int p, int K, List P_new, List P_old)
{
	int i, j, k;
	double max = 0.0;
	for(k=0; k<K; ++k){
		mat oldP = as<arma::mat>(P_old[k]);
		mat newP = as<arma::mat>(P_new[k]);
		for(i=0; i<p-1; ++i)
			for(j=i+1; j<p; ++j)
				max = fmax(max,
					abs(newP(i,j) - oldP(i,j))
					);
	}
	return max;
}

// [[Rcpp::export]]
List gembag_cpp(
    List S_l,
    NumericVector n_l,
    arma::mat v_0,
    arma::mat v_1,
    int maxiter,
    arma::mat P_1_mat,
    List P_2_mat,
    double tau,
    double eps
){
  // Input:
  //   S_l: sample covariance matrices
  //   n_l: sample sizes of the classes
  //   v_0: spike prior parameter
  //   v_1: slab prior parameter
  //   maxiter: maximum of number of iterations
  //   p_1: parameter of Bernoulli prior on binary latent indicator $\gamma_{ij}$
  //   tau: parameter of exponential prior on diagonal entries of precision matrix
  //
  //  Return:
  //    Theta_l: estimated precision matrices
  //    P_l: estimated posterior inclusion probabilities
  //    W_l: estimated covariance matrices

  int curr_iter = 0;
  double p = as<arma::mat>(S_l[0]).n_rows; //p is the number of nodes (biomolecules) in the graph
  int k = S_l.size(), conv = 0;
  //double obj_old, obj_new = 0; //objective function: the posterior
  double max_pdiff;
  bool should_use_p1 = arma::approx_equal(arma::mean(P_1_mat), arma::mat(1, 1, arma::fill::ones), "absdiff", 1e-6);
  mat p2k = ones(p, p); //matrix of prior probs p2 for each group
  List W_l(k);
  List Theta_l(k);
  List P_l(k);
  for(int i = 0; i < k; i++) {
    W_l[i] = eye<mat>(p, p);
    Theta_l[i] = eye<mat>(p, p);
  }

  // loop until max iteration is reached
  // where is the convergence criteria? maybe this is saved for later once we get it working for a small number of iterations?

  //while loop here

  //obj_new = posterior_fcn(k, p, n_l, Theta_l, S_l, tau, P_1_mat, P_2_mat, v_1, v_0);
  //obj_old = obj_new + (eps + 100.0);

  while(conv == 0 && (curr_iter < maxiter)){

    //update previous posterior
    //obj_old = obj_new;

    // E-step
    // Calculate P_l
    P_l = calculate_p_l(Theta_l, P_1_mat, P_2_mat,  v_0, v_1, p, k, should_use_p1);

    for(int i = 0; i < k; i++) {
      mat Theta = as<arma::mat>(Theta_l[i]);
      mat W = as<arma::mat>(W_l[i]);
      mat P = as<arma::mat>(P_l[i]);
      mat S = as<arma::mat>(S_l[i]);
      int n = n_l[i];

      //for(int o=0; o<5; o++){ //not sure why this is here, not present in the algorithm in Gan 2019
      // M-step
      for(int i = 0; i<p; i++){
		//if (i % 10 == 0) Rprintf("column=%d\n", i); // to keep track of progress

        uvec idx(p); // unsigned integer vector of length p
        std::iota(std::begin(idx),std::end(idx),0); // fill vector with increasing integers starting from 0
        idx.shed_row(i); // remove the ith row for idx, so its a length p-1 vector

        uvec indice(1); // devind a length 1 unsigned integer vector
        indice << i; //initialize it with a value of i
        Mat<double> W_11 = W.submat(idx, idx); //get the submatrix of W which excludes row and column i

        Mat<double> W_12 = W.submat(idx, indice); // get the submatrix of W which excludes row i, and only includes column i, so column i excluding row i
        double W_22 = W(i,i); // get the ith diagonal entry

        Mat<double> Theta_12 = Theta.submat(idx, indice); // analogous to W_12
        Mat<double> Theta_11 = Theta(idx,idx); // analogous to W_11
        double Theta_22 = Theta(i,i); // analogous to W_22

        Mat<double> P_12 = P.submat(idx, indice); // anologous to W_12
        Mat<double> S_12 = S.submat(idx, indice); // analogous to W_12
        Mat<double> v_0_12 = v_0.submat(idx, indice); // analogous to W_12
        Mat<double> v_1_12 = v_1.submat(idx, indice); // analogous to W_12

        Mat<double> Theta_11_inv = W_11 - W_12*W_12.t()/W_22;
        W_22 = S(i,i) + 2.0/n*tau;
        Theta_12 = Cor_desc(Theta_11_inv, W_22, Theta_12, S_12, v_0_12, v_1_12, P_12, n);
        Theta_22 = 1.0/W_22 + as_scalar(Theta_12.t()*Theta_11_inv*Theta_12);

        Theta(indice, idx) = Theta_12.t();
        Theta(idx, indice) = Theta_12;
        Theta(i,i) = Theta_22;

        double tmp = as_scalar(Theta_22 - Theta_12.t()*Theta_11_inv*Theta_12);
        Mat<double> tmp1 = Theta_11_inv *  Theta_12;
        W_11 = Theta_11_inv + tmp1 * tmp1.t()/tmp;
        W_12 = -tmp1/tmp;

        W(idx,idx) = W_11;
        W(indice, idx) = W_12.t();
        W(idx, indice) = W_12;
        W(i,i)  = W_22;
      }
      //}
      Theta_l[i] = Theta;
      W_l[i] = W;
      P_l[i] = P;


      //calculate posterior
      //obj_new = posterior_fcn(k, p, n_l, Theta_l, S_l, tau, P_1_mat, P_2_mat, v_1, v_0);

    }

	List P_l_new = calculate_p_l(Theta_l, P_1_mat, P_2_mat,  v_0, v_1, p, k, should_use_p1);
	max_pdiff = max_diff(p, k, P_l_new, P_l);
	//Rprintf("iter = %d max_diff = %f\n", curr_iter, max_pdiff);
	if(max_pdiff < eps)
		conv = 1;
	P_l = P_l_new;

    curr_iter++;
  } //closes EM iteration loop
  // Calculate P_l


  List returnlist(4);
  returnlist(0) = Theta_l;
  returnlist(1) = P_l;
  returnlist(2) = W_l;
  returnlist(3) = curr_iter;

  return returnlist;

}

