#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include "my_project_headers.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// needed for interface with R and C++
void copy_matrix_list_to_3d(const List& mat_list,
                            std::vector<std::vector<std::vector<double>>>& dest,
                            int K, std::vector<int>& n_row, int ncol)
	// nrow is a K-vector
{
  int n;
  for (int k = 0; k < K; ++k) {
    arma::mat curr_mat = as<arma::mat>(mat_list[k]);
	n = n_row[k];
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < ncol; ++j) {
        dest[k][i][j] = curr_mat(i, j);
      }
    }
  }
}

int max_int(const std::vector<int>& vec, int len) {
  int i;
  int max = vec[0];
  for(i=0; i<len; ++i){
    if(vec[i] > max)
      max = vec[i];
  }
  return max;
}

double res_sumsq(
  int i, int k, int p, int n,
  const std::vector<std::vector<std::vector<double>>>& Y,
  const std::vector<std::vector<std::vector<double>>>& Theta
) {
  int idx, j; //idx is the sum over 1, ..., n
  double del_sum = 0.0, rss = 0.0;
  for(idx=0; idx<n; ++idx){
    for(j=0; j<p; ++j){
      if(j != i)
        del_sum += Theta[k][i][j] / Theta[k][i][i] * Y[k][idx][j];
    }
    rss += pow(Y[k][idx][i] + del_sum, 2.0);
    del_sum = 0.0;
  }
  return rss;
}

int nonzero_thetas(
  int i, int p, int k,
  const std::vector<std::vector<std::vector<double>>>& Theta
)
{
  int j;
  int count_nonzero = 0;
  for(j=i+1; j<p; ++j){
    if(fabs(Theta[k][i][j]) > .000001)
      ++count_nonzero;
  }
  return count_nonzero;
}

double bic_fcn(
  int K,
  const std::vector<int>& n_vec,
  int p,
  const std::vector<std::vector<std::vector<double>>>& Y,
  const std::vector<std::vector<std::vector<double>>>& Theta
){
  int i, k, n;
  double bic = 0.0;
  for(k=0; k<K; ++k){
    n = n_vec[k];
    for(i=0; i<p; ++i){
      bic += res_sumsq(i, k, p, n, Y, Theta) +
        log( (double) n) * nonzero_thetas(i, p, k, Theta);
    }
  }
  return bic;
}

// [[Rcpp::export]]
List bplane_bic(
  int K,
  IntegerVector n_vec_r,
  int p,
  List Y_list,
  List Theta_list
) {
  double bic;
  std::vector<int> n_vec(n_vec_r.begin(), n_vec_r.end());
  std::vector<int> p_vec(K, p); // K-elements of p
  if (Y_list.size() != K || Theta_list.size() != K) {
    stop("Number of matrices in Y_list and Theta_list must match length of n_vec");
  }
  int n_max = max_int(n_vec, K);
  std::vector<std::vector<std::vector<double>>> Y(K,
                                                  std::vector<std::vector<double>>(n_max, std::vector<double>(p, 0.0)));
  std::vector<std::vector<std::vector<double>>> Theta(K,
                                                      std::vector<std::vector<double>>(p, std::vector<double>(p, 0.0)));
  copy_matrix_list_to_3d(Y_list, Y, K, n_vec, p);
  copy_matrix_list_to_3d(Theta_list, Theta, K, p_vec, p);

  bic = bic_fcn(K, n_vec, p, Y, Theta);

  List result;
  result["bic"] = bic;
  return result;
}


void fill_in_omega(int K, int p, std::vector<std::vector<std::vector<double>>>& omega) {
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < p; ++i) {
      for (int j = i + 1; j < p; ++j) {
        omega[k][j][i] = omega[k][i][j];
      }
      omega[k][i][i] = 1.0;
    }
  }
}

void print_3array(const std::vector<std::vector<std::vector<double>>>& array,
	int d1, int d2, int d3, std::string name)
{
	int i, j, k;
	for(k=0; k<d1; ++k){
		Rprintf("%s[%d]\n", name.c_str(), k);
		for(i=0; i<d2; ++i){
			for(j=0; j<d3; ++j)
				Rprintf("%f\t", array[k][i][j]);
			Rprintf("\n");
		}
		Rprintf("\n");
	}
}

double soft_threshold(double x, double lambda) {
  return std::copysign(1.0, x) * std::max(std::fabs(x) - lambda, 0.0);
}

double quad_formula(double a, double b, double c) {
  return (-1.0 * b + std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
}

void copy_3array_uppertri(const std::vector<std::vector<std::vector<double>>>& copy_from,
                          int d1, int d2, int d3,
                          std::vector<std::vector<std::vector<double>>>& copy_to) {
  for (int i = 0; i < d1; ++i)
    for (int j = 0; j < d2 - 1; ++j)
      for (int k = j + 1; k < d3; ++k)
        copy_to[i][j][k] = copy_from[i][j][k];
}

double max_diff_uppertri(const std::vector<std::vector<std::vector<double>>>& a,
                         const std::vector<std::vector<std::vector<double>>>& b,
                         int d1, int d2, int d3) {
  double max = 0.0;
  for (int i = 0; i < d1; ++i)
    for (int j = 0; j < d2 - 1; ++j)
      for (int k = j + 1; k < d3; ++k)
        max = std::max(max, std::fabs(a[i][j][k] - b[i][j][k]));

  return max;
}

void print_Theta_mats(int iter, int k, int p,
                      const std::vector<std::vector<std::vector<double>>>& Theta) {
  Rcout << "Iter= " << iter << std::endl;
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++)
      Rcout << Theta[k][i][j] << "\t";
    Rcout << std::endl;
  }
  Rcout << std::endl;
}

void less_np(const std::vector<int>& n_vec, int p, int K,
	std::vector<char>& np_less)
{
	int k;
	for(k=0; k<K; ++k){
		if(n_vec[k] < p)
			np_less[k] = 'n';
		else
			np_less[k] = 'p';
	}
}

// from RC_interface.c

void copy_1d_to_3d(const std::vector<double>& src,
                   std::vector<std::vector<std::vector<double>>>& dest,
                   int K, int p1, int p2) {
  int idx = 0;
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < p1; ++i) {
      for (int j = 0; j < p2; ++j) {
        dest[k][i][j] = src[idx++];
      }
    }
  }
}

void copy_numericmatrix_to_2d(const NumericMatrix src,
	std::vector<std::vector<double>>& dest,
	int d1, int d2)
{
	int i, j;
	for(i=0; i<d1; ++i)
		for(j=0; j<d2; ++j)
			dest[i][j] = src(i,j);
}

List convert_3d_to_matrix_list(const std::vector<std::vector<std::vector<double>>>& src,
                               int K, int d1, int d2) {
  List result(K);
  for (int k = 0; k < K; ++k) {
    arma::mat curr_mat(d1, d2);
    for (int i = 0; i < d1; ++i) {
      for (int j = 0; j < d2; ++j) {
        curr_mat(i, j) = src[k][i][j];
      }
    }
    result[k] = curr_mat;
  }
  return result;
}

// headers

void calc_S(const std::vector<int>& n_vec, int p, int K,
	const std::vector<std::vector<std::vector<double>>>& Y,
	std::vector<std::vector<std::vector<double>>>& S);

void CM_step(int K, const std::vector<int>& n_vec, int p, double tau, double v0, double v1,
             const std::vector<std::vector<std::vector<double>>>& S,
             std::vector<std::vector<std::vector<double>>>& Theta,
             const std::vector<std::vector<std::vector<double>>>& omega,
			 std::vector<std::vector<std::vector<double>>>& R,
			 const std::vector<std::vector<std::vector<double>>>& Y,
			 const std::vector<char>& np_less);

void update_theta_ij(int i, int j, int k, int p, int n,
                     const std::vector<std::vector<std::vector<double>>>& S,
                     std::vector<std::vector<std::vector<double>>>& Theta,
                     const std::vector<std::vector<std::vector<double>>>& lambda,
					 std::vector<std::vector<std::vector<double>>>& R,
					 const std::vector<std::vector<std::vector<double>>>& Y,
					 const std::vector<char> np_less);

void update_theta_ii(int i, int k, int p, int n, double tau,
                     const std::vector<std::vector<std::vector<double>>>& S,
                     std::vector<std::vector<std::vector<double>>>& Theta,
					 std::vector<std::vector<std::vector<double>>>& R,
					 const std::vector<std::vector<std::vector<double>>>& Y,
					 const std::vector<char> np_less);

void calc_omega(int K, int p, double v0, double v1,
                const std::vector<std::vector<double>>& p1,
                const std::vector<std::vector<std::vector<double>>>& p2,
                const std::vector<std::vector<std::vector<double>>>& Theta,
                std::vector<std::vector<std::vector<double>>>& omega);

void starting_values(int K, const std::vector<int>& n_vec, int p, double tau, double v0, double v1,
                     const std::vector<std::vector<std::vector<double>>>& S,
                     std::vector<std::vector<std::vector<double>>>& Theta,
					 const std::vector<std::vector<double>>& p1,
                     const std::vector<std::vector<std::vector<double>>>& p2,
					 std::vector<std::vector<std::vector<double>>>& omega,
					 std::vector<std::vector<std::vector<double>>>& R,
					 const std::vector<std::vector<std::vector<double>>>& Y,
					 const std::vector<char>& np_less);

void bplane(const std::vector<int>& n_vec, int p, int K, double v0, double v1, double tau,
              int init_iter, int max_iter, double omega_eps,
			  const std::vector<std::vector<double>>& p1,
              const std::vector<std::vector<std::vector<double>>>& p2,
			  const std::vector<std::vector<std::vector<double>>>& Y,
              std::vector<std::vector<std::vector<double>>>& Theta,
              std::vector<std::vector<std::vector<double>>>& omega,
              int& n_iter);


// [[Rcpp::export]]
List bplane_r_wrapper(
    IntegerVector n_vec_r,
    int p,
    double v0,
    double v1,
    double tau,
    int init_iter,
    int max_iter,
    double omega_eps,
    List Y_list,
	NumericMatrix p1_mat,
    List p2_list
) {
  std::vector<int> n_vec(n_vec_r.begin(), n_vec_r.end());
  int K = n_vec.size();
  std::vector<int> p_vec(K, p); // K-elements of p

  if (Y_list.size() != K || p2_list.size() != K) {
    stop("Number of matrices in Y_list and prior_list must match length of n_vec");
  }

  int n_max = max_int(n_vec, K);

  std::vector<std::vector<std::vector<double>>> Y(K,
                                                  std::vector<std::vector<double>>(n_max, std::vector<double>(p, 0.0)));
  std::vector<std::vector<std::vector<double>>> p2(K,
                                                      std::vector<std::vector<double>>(p, std::vector<double>(p, 0.0)));
  std::vector<std::vector<std::vector<double>>> Theta(K,
                                                      std::vector<std::vector<double>>(p, std::vector<double>(p, 0.0)));
  std::vector<std::vector<std::vector<double>>> omega(K,
                                                      std::vector<std::vector<double>>(p, std::vector<double>(p, 0.0)));
  std::vector<std::vector<double>> p1(p, std::vector<double>(p, 0.0));

  copy_matrix_list_to_3d(Y_list, Y, K, n_vec, p);
  copy_matrix_list_to_3d(p2_list, p2, K, p_vec, p);
  copy_numericmatrix_to_2d(p1_mat, p1, p, p);

  int n_iter;
  bplane(n_vec, p, K, v0, v1, tau, init_iter, max_iter, omega_eps,
           p1, p2, Y, Theta, omega, n_iter);

  List Theta_list = convert_3d_to_matrix_list(Theta, K, p, p);
  List omega_list = convert_3d_to_matrix_list(omega, K, p, p);

  List result;
  result["Theta"] = Theta_list;
  result["omega"] = omega_list;
  result["n_iter"] = n_iter;

  return result;
}

void bplane(const std::vector<int>& n_vec, int p, int K, double v0, double v1, double tau,
              int init_iter, int max_iter, double omega_eps,
			  const std::vector<std::vector<double>>& p1,
              const std::vector<std::vector<std::vector<double>>>& p2,
			  const std::vector<std::vector<std::vector<double>>>& Y,
              std::vector<std::vector<std::vector<double>>>& Theta,
              std::vector<std::vector<std::vector<double>>>& omega,
              int& n_iter) {
  int iter = 0, conv = 0;

  double max = 0.0;
  int n_max = max_int(n_vec, K);
  std::vector<std::vector<std::vector<double>>> R(K,
                                                          std::vector<std::vector<double>>(n_max, std::vector<double>(p, 0.0)));
  std::vector<std::vector<std::vector<double>>> omega_old(K,
                                                          std::vector<std::vector<double>>(p, std::vector<double>(p, 0.0)));
  std::vector<std::vector<std::vector<double>>> S(K,
                                                          std::vector<std::vector<double>>(p, std::vector<double>(p, 0.0)));
  std::vector<char> np_less(K);
  less_np(n_vec, p, K, np_less);
  calc_S(n_vec, p, K, Y, S);
  starting_values(K, n_vec, p, tau, v0, v1, S, Theta, p1, p2, omega, R, Y, np_less);

  calc_omega(K, p, v0, v1, p1, p2, Theta, omega);

  while (conv == 0 && iter < max_iter) {
    if (iter >= init_iter)
      copy_3array_uppertri(omega, K, p, p, omega_old);

    CM_step(K, n_vec, p, tau, v0, v1, S, Theta, omega, R, Y, np_less);
    calc_omega(K, p, v0, v1, p1, p2, Theta, omega);

    if (iter >= init_iter) {
      max = max_diff_uppertri(omega, omega_old, K, p, p);
	  //Rprintf("iter=%d , max_diff=%f\n", iter, max);
      if (max < omega_eps)
        conv = 1;
    }
    ++iter;
  }

  n_iter = iter;
  fill_in_omega(K, p, omega);
}

void calc_S(const std::vector<int>& n_vec, int p, int K,
	const std::vector<std::vector<std::vector<double>>>& Y,
	std::vector<std::vector<std::vector<double>>>& S)
{
	int i, j, k, n, idx;
	double temp = 0.0;
	for(k=0; k<K; ++k){
		n = n_vec[k];
		for(i=0; i<p; ++i){
			for(j=i; j<p; ++j){
				for(idx=0; idx<n; ++idx)
					temp += Y[k][idx][i] * Y[k][idx][j];
				S[k][i][j] = S[k][j][i] = temp / n;
				temp = 0.0;
			}
		}
	}
}

void starting_values(int K, const std::vector<int>& n_vec, int p, double tau, double v0, double v1,
                     const std::vector<std::vector<std::vector<double>>>& S,
                     std::vector<std::vector<std::vector<double>>>& Theta,
					 const std::vector<std::vector<double>>& p1,
                     const std::vector<std::vector<std::vector<double>>>& p2,
					 std::vector<std::vector<std::vector<double>>>& omega,
					 std::vector<std::vector<std::vector<double>>>& R,
					 const std::vector<std::vector<std::vector<double>>>& Y,
					 const std::vector<char>& np_less) {

  int k, n, i, idx, j;
  for (k = 0; k < K; ++k) { // R = Y and Theta = I
	n = n_vec[k];
    for (i = 0; i < p; ++i) {
	  for(idx=0; idx<n; ++idx)
		R[k][idx][i] = Y[k][idx][i];
      for (j = i; j < p; ++j) {
		omega[k][i][j] = omega[k][j][i] = p1[i][j] * p2[k][i][j];
        if (i == j)
          Theta[k][i][j] = 1.0;
        else
          Theta[k][i][j] = Theta[k][j][i] = 0.0;
      }
    }
  }

  CM_step(K, n_vec, p, tau, v0, v1, S, Theta, omega, R, Y, np_less);
}

double calc_1_omega(int K, int k, int i, int j, double v1, double v0,
	const std::vector<std::vector<double>>& p1,
    const std::vector<std::vector<std::vector<double>>>& p2,
    const std::vector<std::vector<std::vector<double>>>& Theta)
{
  double eta1, eta2, S1, S2, temp, sum=0.0;
  //calculate eta1
  if(fabs(p1[i][j] - 1.0) <  1e-8){
	  eta1 = 1;
  }
  else{
	  // first S1
    for(int l=0; l<K; ++l){
	  temp = p2[l][i][j] / (2.0*v1) * exp(-1.0 * std::fabs(Theta[l][i][j]) / v1) +
		(1.0 - p2[l][i][j]) / (2.0*v0) * exp(-1.0 * std::fabs(Theta[l][i][j]) / v0);
	  sum += log(temp);
    }
    S1 = exp(sum);
	//next S2
    sum = 0.0;
    for(int l=0; l<K; ++l){
	  temp = 0.5/v0* exp(-1.0 * std::fabs(Theta[l][i][j]) / v0);
	  sum += log(temp);
    }
    S2 = exp(sum);
    eta1 = p1[i][j] * S1 / (p1[i][j] * S1 + (1-p1[i][j]) * S2);
  }
  //calculate eta2
  double num = p2[k][i][j] / (2.0*v1) * std::exp(-1.0 * std::fabs(Theta[k][i][j]) / v1);
  double den = num + (1 - p2[k][i][j]) / (2.0*v0) * std::exp(-1.0 * std::fabs(Theta[k][i][j]) / v0);
  eta2 = num / den;
  return eta1 * eta2;
}

void calc_omega(int K, int p, double v0, double v1,
                const std::vector<std::vector<double>>& p1,
                const std::vector<std::vector<std::vector<double>>>& p2,
                const std::vector<std::vector<std::vector<double>>>& Theta,
                std::vector<std::vector<std::vector<double>>>& omega) {
  for (int k = 0; k < K; ++k)
    for (int i = 0; i < p-1; ++i)
      for (int j = i + 1; j < p; ++j)
        omega[k][i][j] = calc_1_omega(K, k, i, j, v1, v0, p1, p2, Theta);
}

void calc_lambda(int K, int p, const std::vector<std::vector<std::vector<double>>>& omega,
                 double v0, double v1,
                 std::vector<std::vector<std::vector<double>>>& lambda) {
  for (int k = 0; k < K; ++k)
    for (int i = 0; i < p - 1; ++i)
      for (int j = i + 1; j < p; ++j)
        lambda[k][i][j] = omega[k][i][j] / v1 + (1 - omega[k][i][j]) / v0;
}



void CM_step(int K, const std::vector<int>& n_vec, int p, double tau, double v0, double v1,
             const std::vector<std::vector<std::vector<double>>>& S,
             std::vector<std::vector<std::vector<double>>>& Theta,
             const std::vector<std::vector<std::vector<double>>>& omega,
			 std::vector<std::vector<std::vector<double>>>& R,
			 const std::vector<std::vector<std::vector<double>>>& Y,
			 const std::vector<char>& np_less) {
  std::vector<std::vector<std::vector<double>>> lambda(K,
                                                       std::vector<std::vector<double>>(p, std::vector<double>(p, 0.0)));

  calc_lambda(K, p, omega, v0, v1, lambda);


  for (int k = 0; k < K; ++k) {
    int n = n_vec[k];
	for (int i = 0; i < p - 1; ++i) {
      for (int j = i + 1; j < p; ++j) {
        update_theta_ij(i, j, k, p, n, S, Theta, lambda, R, Y, np_less);
		//if((i==0)) print_Theta_mats(0, 0, p, Theta);
      }
    }
	for(int i = 0; i < p; ++i)
		update_theta_ii(i, k, p, n, tau, S, Theta, R, Y, np_less);
  }
}

void update_theta_ij(int i, int j, int k, int p, int n,
                     const std::vector<std::vector<std::vector<double>>>& S,
                     std::vector<std::vector<std::vector<double>>>& Theta,
                     const std::vector<std::vector<std::vector<double>>>& lambda,
					 std::vector<std::vector<std::vector<double>>>& R,
					 const std::vector<std::vector<std::vector<double>>>& Y,
					 const std::vector<char> np_less) {

  double sum1 = 0.0, sum2 = 0.0, temp1 = 0.0, temp2 = 0.0;
  double theta_ij_old, mult_i, mult_j;
  theta_ij_old = Theta[k][i][j];
  // two ways of calculating sum(theta[i,]*S[j,]) depending on
  //	whether n or p is smaller
  //Rprintf("n=%d, p=%d\n", n, p);
  if(np_less[k] == 'p'){ // sum of p calculation
    for (int idx = 0; idx < p; ++idx) {
      if (idx != j)
        sum1 += Theta[k][i][idx] * S[k][j][idx];
      if (idx != i)
        sum2 += Theta[k][j][idx] * S[k][i][idx];
    }
  }
  if(np_less[k] == 'n'){ // calculation with t(Y) * r
	for(int idx=0; idx<n; ++idx){
	  temp1 += Y[k][idx][j] * R[k][idx][i];
	  temp2 += Y[k][idx][i] * R[k][idx][j];
	}
	sum1 = -1.0 * Theta[k][i][j] * S[k][j][j] + Theta[k][i][i] / n * temp1;
	sum2 = -1.0 * Theta[k][i][j] * S[k][i][i] + Theta[k][j][j] / n * temp2;
  }
  /*if((i==0)){
	  printf("j=%d\n", j);
	  printf("sum1=%f\tsum2=%f\n", sum1, sum2);
	  printf("lambda=%f, n=%d\n", lambda[0][0][1], n);
  }*/
  Theta[k][i][j] = Theta[k][j][i] = soft_threshold(
    -1.0 * (sum1 + sum2), lambda[k][i][j] / n
  ) / (S[k][i][i] + S[k][j][j]);
  if(np_less[k] == 'n'){ // update r_i and r_j
	mult_i = (Theta[k][i][j] - theta_ij_old) / Theta[k][i][i];
	mult_j = (Theta[k][i][j] - theta_ij_old) / Theta[k][j][j];
	for(int idx=0; idx<n; ++idx){
	  R[k][idx][i] = R[k][idx][i] + mult_i * Y[k][idx][j];
	  R[k][idx][j] = R[k][idx][j] + mult_j * Y[k][idx][i];
	}
  }
}

void update_theta_ii(int i, int k, int p, int n, double tau,
                     const std::vector<std::vector<std::vector<double>>>& S,
                     std::vector<std::vector<std::vector<double>>>& Theta,
					 std::vector<std::vector<std::vector<double>>>& R,
					 const std::vector<std::vector<std::vector<double>>>& Y,
					 const std::vector<char> np_less) {
  double sum = 0.0, temp = 0.0, mult;
  double theta_ii_old = Theta[k][i][i];

  if(np_less[k] == 'p') // sum of p calculation
    for (int idx = 0; idx < p; ++idx)
      if (idx != i)
        sum += Theta[k][i][idx] * S[k][i][idx];

  if(np_less[k] == 'n'){ // calculation with t(Y) * r
	for(int idx=0; idx<n; ++idx)
	  temp += Y[k][idx][i] * R[k][idx][i];
	sum = -1.0 * Theta[k][i][i] * S[k][i][i] + Theta[k][i][i] / n * temp;
  }

  double b = sum - tau;
  Theta[k][i][i] = quad_formula(S[k][i][i], b, -1.0);

  if(np_less[k] == 'n'){ // update r_i;
	mult = theta_ii_old / Theta[k][i][i];
	for(int idx=0; idx<n; ++idx)
	  R[k][idx][i] = (R[k][idx][i] - Y[k][idx][i]) * mult + Y[k][idx][i];
  }
}

double posterior(int K, int p, const std::vector<int>& n_vec, double tau, double v0, double v1,
                 const std::vector<std::vector<std::vector<double>>>& Theta,
                 const std::vector<std::vector<std::vector<double>>>& S,
                 const std::vector<std::vector<std::vector<double>>>& prior) {
  double temp1 = 0.0, temp3 = 0.0;
  double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;

  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < p; ++i) {
      temp1 += std::log(Theta[k][i][i]);
      term2 += Theta[k][i][i];
    }
    term1 += n_vec[k] * temp1;
    temp1 = 0.0;
  }

  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int l = 0; l < p; ++l) {
          temp3 += Theta[k][i][j] * Theta[k][i][l] * S[k][j][l];
        }
      }
    }
    term3 += n_vec[k] * temp3;
    temp3 = 0.0;
  }

  for (int k = 0; k < K; ++k)
    for (int i = 0; i < p; ++i)
      for (int j = i + 1; j < p; ++j)
        term4 += std::log(
          prior[k][i][j] / (2.0 * v1) * std::exp(-1.0 * std::fabs(Theta[k][i][j]) / v1) +
            (1 - prior[k][i][j]) / (2.0 * v0) * std::exp(-1.0 * std::fabs(Theta[k][i][j]) / v0)
        );

  return -1.0 * term1 + tau * term2 + 0.5 * term3 - term4;
}
