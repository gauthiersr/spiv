#include <RcppArmadillo.h>

using namespace arma;

// Anderson-Rubin statistic for SP‑IV
// [[Rcpp::export]]
double compute_AR_stat(const arma::vec& b,
                       const arma::mat& y_res,
                       const arma::mat& Y_res,
                       const arma::mat& Pz,
                       const arma::mat& Mz,
                       int H,
                       int T_ess,
                       int d_AR) {

  // u_res_b = y_res - (b' %x% I_H) * Y_res
  mat I_H = eye(H, H);
  mat kron_b_I = kron(b.t(), I_H);
  mat u_res_b = y_res - kron_b_I * Y_res;

  // Quadratic form: u_res_b * Pz * u_res_b' * inv(u_res_b * Mz * u_res_b')
  mat uPu = u_res_b * Pz * u_res_b.t();
  mat uMu = u_res_b * Mz * u_res_b.t();
  mat inner = uPu * solve(uMu, eye(H, H));

  double trace_val = trace(inner);
  double AR_stat = (T_ess - d_AR) * trace_val;

  return AR_stat;
}

// Kleibergen LM statistic for SP‑IV
// [[Rcpp::export]]
double compute_KLM_stat(const arma::vec& b,
                        const arma::mat& y_res,
                        const arma::mat& Y_res,
                        const arma::mat& v_res,
                        const arma::mat& Mz,
                        const arma::mat& Pz,
                        const arma::mat& R,
                        int H,
                        int T_ess,
                        int d_K) {

  // 1. u_res_b
  mat I_H = eye(H, H);
  mat kron_b_I = kron(b.t(), I_H);
  mat u_res_b = y_res - kron_b_I * Y_res;

  // 2. Projected residuals
  mat u_check_b = u_res_b * Mz;
  mat v_check   = v_res   * Mz;

  // 3. Xi = u_res_b * Mz * u_res_b'
  mat Xi = u_res_b * Mz * u_res_b.t();

  // 4. Y_check_b
  mat u_check_cross = u_check_b * u_check_b.t();
  mat A = solve(u_check_cross, u_res_b * Pz);
  mat Y_check_b = Y_res * Pz - v_check * u_check_b.t() * A;

  // 5. Score matrix and vector S
  mat Score_mat = solve(Xi, u_res_b * Y_check_b.t());
  vec S = R.t() * vectorise(Score_mat);

  // 6. Variance V
  mat inner = solve(Xi, u_res_b * u_res_b.t() * solve(Xi, eye(H, H)));
  mat V_inner = kron(Y_check_b * Y_check_b.t(), inner);
  mat V = R.t() * V_inner * R;

  // 7. KLM statistic
  double stat = (T_ess - d_K) * as_scalar(S.t() * solve(V, S));
  return stat;
}
