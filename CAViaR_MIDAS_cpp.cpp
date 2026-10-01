// CAViaR_MIDAS_cpp.cpp
// ============================================================================
// Rcpp functions for joint VaR-ES CAViaR-RMIDAS estimation
//
// Beta weighting: controlled by the fix_w1 flag
//   fix_w1 = true  (default): w1 = 1 fixed, theta contains w2 only
//   fix_w1 = false:           theta contains (w1, w2), both free
//
// Parameter vectors (fix_w1 = true):
//   Linear  SAV: (phi, psi_1..psi_K, w2, alpha, beta_ar, gamma0)
//   Linear  AS:  (phi, psi_1..psi_K, w2, alpha_pos, alpha_neg, beta_ar, gamma0)
// Parameter vectors (fix_w1 = false):
//   Linear  SAV: (phi, psi_1..psi_K, w1, w2, alpha, beta_ar, gamma0)
//   Linear  AS:  (phi, psi_1..psi_K, w1, w2, alpha_pos, alpha_neg, beta_ar, gamma0)
// (Exponential variants follow the same pattern.)
//
// Convention: VaR and ES are negative for tau < 0.5; mu = 0
// ============================================================================

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// ============================================================================
// Helper: normalised Beta lag weights
// ============================================================================
static std::vector<double> beta_weights_vec(int S, double w1, double w2) {
  std::vector<double> wt(S);
  double max_log = -1.0e30;
  for (int s = 0; s < S; s++) {
    double u = (double)(s + 1) / (double)S;
    wt[s] = (w1 - 1.0) * std::log(u + 1.0e-15) +
            (w2 - 1.0) * std::log(1.0 - u + 1.0e-15);
    if (wt[s] > max_log) max_log = wt[s];
  }
  double sum_wt = 0.0;
  for (int s = 0; s < S; s++) {
    wt[s] = std::exp(wt[s] - max_log);
    sum_wt += wt[s];
  }
  for (int s = 0; s < S; s++) wt[s] /= sum_wt;
  return wt;
}

// ============================================================================
// Helper: monthly long-run component
// factors is pre-lagged: row t = month t-1 value.
// s=0 reads factors[t], s=1 reads factors[t-1], etc.
// Returns vector of length (T_month - trim), trimming the first trim months.
// ============================================================================
static std::vector<double> compute_LR_monthly(
    const double* factors, int T_month, int K, int S,
    const double* psi, const std::vector<double>& wt,
    double phi, bool is_exponential, int trim)
{
  std::vector<double> LR(T_month - trim);
  for (int t = trim; t < T_month; t++) {
    double inner = 0.0;
    for (int s = 0; s < S; s++) {
      int lag_idx = t - s;
      if (lag_idx < 0) lag_idx = 0;
      double factor_sum = 0.0;
      for (int k = 0; k < K; k++) {
        double x_val = factors[lag_idx + k * T_month];
        factor_sum += psi[k] * (is_exponential ? x_val : std::fabs(x_val));
      }
      inner += wt[s] * factor_sum;
    }
    LR[t - trim] = is_exponential ? phi * std::exp(inner) : phi + inner;
  }
  return LR;
}

// ============================================================================
// Helper: unpack (w1, w2) from theta at position idx
// ============================================================================
static bool unpack_w(const NumericVector& theta, int& idx,
                     bool fix_w1, double& w1, double& w2) {
  if (fix_w1) { w1 = 1.0; w2 = theta[idx++]; }
  else        { w1 = theta[idx++]; w2 = theta[idx++]; }
  return (w1 > 0.01 && w1 <= 200.0 && w2 >= 1.0 && w2 <= 200.0);
}


// ============================================================================
// 1. Fused neg AL log-likelihood — Linear SAV
// ============================================================================
// [[Rcpp::export]]
double neg_al_loglik_midas_linear_sav_cpp(
    NumericVector theta, NumericVector y, NumericMatrix factors,
    IntegerVector month_id, double tau, double q_init,
    int S, int n_burn, bool fix_w1, int trim)
{
  const int TT = y.size(), T_month = factors.nrow(), K = factors.ncol();
  int idx = 0;
  double phi = theta[idx++];
  std::vector<double> psi(K);
  for (int k = 0; k < K; k++) psi[k] = theta[idx++];
  double w1, w2;
  if (!unpack_w(theta, idx, fix_w1, w1, w2)) return 1.0e12;
  double alpha = theta[idx++], beta_ar = theta[idx++];
  if (beta_ar < 0.0 || beta_ar >= 1.0) return 1.0e12;
  double gamma0 = theta[idx++], es_mult = 1.0 + std::exp(gamma0);

  auto wt   = beta_weights_vec(S, w1, w2);
  auto LR_m = compute_LR_monthly(factors.begin(), T_month, K, S, psi.data(), wt, phi, false, trim);
  int LR_len = (int)LR_m.size();
  for (int t = 0; t < LR_len; t++) if (!std::isfinite(LR_m[t])) return 1.0e12;

  const int start = n_burn + 1;
  // month_id is 1-based into LR_m (which is already trimmed)
  double SR_prev = q_init - LR_m[month_id[0] - 1], neg_ll = 0.0;
  for (int t = 1; t < TT; t++) {
    double LR_t = LR_m[month_id[t] - 1];
    double SR_t = alpha * std::fabs(y[t-1]) + beta_ar * SR_prev;
    double Q_t  = LR_t + SR_t; SR_prev = SR_t;
    if (t < start) continue;
    double ES_t = es_mult * Q_t;
    if (ES_t >= 0.0) return 1.0e12;
    double resid = y[t] - Q_t, hit = (y[t] <= Q_t) ? 1.0 : 0.0;
    neg_ll -= std::log((tau - 1.0) / ES_t) + resid * (tau - hit) / (tau * ES_t);
    if (!std::isfinite(neg_ll)) return 1.0e12;
  }
  return neg_ll;
}

// ============================================================================
// 2. Fused neg AL log-likelihood — Linear AS
// ============================================================================
// [[Rcpp::export]]
double neg_al_loglik_midas_linear_as_cpp(
    NumericVector theta, NumericVector y, NumericMatrix factors,
    IntegerVector month_id, double tau, double q_init,
    int S, int n_burn, bool fix_w1, int trim)
{
  const int TT = y.size(), T_month = factors.nrow(), K = factors.ncol();
  int idx = 0;
  double phi = theta[idx++];
  std::vector<double> psi(K);
  for (int k = 0; k < K; k++) psi[k] = theta[idx++];
  double w1, w2;
  if (!unpack_w(theta, idx, fix_w1, w1, w2)) return 1.0e12;
  double alpha_pos = theta[idx++], alpha_neg = theta[idx++], beta_ar = theta[idx++];
  if (beta_ar < 0.0 || beta_ar >= 1.0) return 1.0e12;
  double gamma0 = theta[idx++], es_mult = 1.0 + std::exp(gamma0);

  auto wt   = beta_weights_vec(S, w1, w2);
  auto LR_m = compute_LR_monthly(factors.begin(), T_month, K, S, psi.data(), wt, phi, false, trim);
  int LR_len = (int)LR_m.size();
  for (int t = 0; t < LR_len; t++) if (!std::isfinite(LR_m[t])) return 1.0e12;

  const int start = n_burn + 1;
  double SR_prev = q_init - LR_m[month_id[0] - 1], neg_ll = 0.0;
  for (int t = 1; t < TT; t++) {
    double LR_t = LR_m[month_id[t] - 1], ylag = y[t-1];
    double ypos = (ylag > 0.0) ? ylag : 0.0, yneg = (ylag < 0.0) ? -ylag : 0.0;
    double SR_t = alpha_pos * ypos + alpha_neg * yneg + beta_ar * SR_prev;
    double Q_t  = LR_t + SR_t; SR_prev = SR_t;
    if (t < start) continue;
    double ES_t = es_mult * Q_t;
    if (ES_t >= 0.0) return 1.0e12;
    double resid = y[t] - Q_t, hit = (y[t] <= Q_t) ? 1.0 : 0.0;
    neg_ll -= std::log((tau - 1.0) / ES_t) + resid * (tau - hit) / (tau * ES_t);
    if (!std::isfinite(neg_ll)) return 1.0e12;
  }
  return neg_ll;
}

// ============================================================================
// 3. Fused neg AL log-likelihood — Exponential SAV
// ============================================================================
// [[Rcpp::export]]
double neg_al_loglik_midas_exp_sav_cpp(
    NumericVector theta, NumericVector y, NumericMatrix factors,
    IntegerVector month_id, double tau, double q_init,
    int S, int n_burn, bool fix_w1, int trim)
{
  const int TT = y.size(), T_month = factors.nrow(), K = factors.ncol();
  int idx = 0;
  double phi = theta[idx++];
  std::vector<double> psi(K);
  for (int k = 0; k < K; k++) psi[k] = theta[idx++];
  double w1, w2;
  if (!unpack_w(theta, idx, fix_w1, w1, w2)) return 1.0e12;
  double alpha = theta[idx++], beta_ar = theta[idx++];
  if (beta_ar < 0.0 || beta_ar >= 1.0) return 1.0e12;
  double gamma0 = theta[idx++], es_mult = 1.0 + std::exp(gamma0);

  auto wt   = beta_weights_vec(S, w1, w2);
  auto LR_m = compute_LR_monthly(factors.begin(), T_month, K, S, psi.data(), wt, phi, true, trim);
  int LR_len = (int)LR_m.size();
  for (int t = 0; t < LR_len; t++) if (!std::isfinite(LR_m[t])) return 1.0e12;

  const int start = n_burn + 1;
  double SR_prev = q_init - LR_m[month_id[0] - 1], neg_ll = 0.0;
  for (int t = 1; t < TT; t++) {
    double LR_t = LR_m[month_id[t] - 1];
    double SR_t = alpha * std::fabs(y[t-1]) + beta_ar * SR_prev;
    double Q_t  = LR_t + SR_t; SR_prev = SR_t;
    if (t < start) continue;
    double ES_t = es_mult * Q_t;
    if (ES_t >= 0.0) return 1.0e12;
    double resid = y[t] - Q_t, hit = (y[t] <= Q_t) ? 1.0 : 0.0;
    neg_ll -= std::log((tau - 1.0) / ES_t) + resid * (tau - hit) / (tau * ES_t);
    if (!std::isfinite(neg_ll)) return 1.0e12;
  }
  return neg_ll;
}

// ============================================================================
// 4. Fused neg AL log-likelihood — Exponential AS
// ============================================================================
// [[Rcpp::export]]
double neg_al_loglik_midas_exp_as_cpp(
    NumericVector theta, NumericVector y, NumericMatrix factors,
    IntegerVector month_id, double tau, double q_init,
    int S, int n_burn, bool fix_w1, int trim)
{
  const int TT = y.size(), T_month = factors.nrow(), K = factors.ncol();
  int idx = 0;
  double phi = theta[idx++];
  std::vector<double> psi(K);
  for (int k = 0; k < K; k++) psi[k] = theta[idx++];
  double w1, w2;
  if (!unpack_w(theta, idx, fix_w1, w1, w2)) return 1.0e12;
  double alpha_pos = theta[idx++], alpha_neg = theta[idx++], beta_ar = theta[idx++];
  if (beta_ar < 0.0 || beta_ar >= 1.0) return 1.0e12;
  double gamma0 = theta[idx++], es_mult = 1.0 + std::exp(gamma0);

  auto wt   = beta_weights_vec(S, w1, w2);
  auto LR_m = compute_LR_monthly(factors.begin(), T_month, K, S, psi.data(), wt, phi, true, trim);
  int LR_len = (int)LR_m.size();
  for (int t = 0; t < LR_len; t++) if (!std::isfinite(LR_m[t])) return 1.0e12;

  const int start = n_burn + 1;
  double SR_prev = q_init - LR_m[month_id[0] - 1], neg_ll = 0.0;
  for (int t = 1; t < TT; t++) {
    double LR_t = LR_m[month_id[t] - 1], ylag = y[t-1];
    double ypos = (ylag > 0.0) ? ylag : 0.0, yneg = (ylag < 0.0) ? -ylag : 0.0;
    double SR_t = alpha_pos * ypos + alpha_neg * yneg + beta_ar * SR_prev;
    double Q_t  = LR_t + SR_t; SR_prev = SR_t;
    if (t < start) continue;
    double ES_t = es_mult * Q_t;
    if (ES_t >= 0.0) return 1.0e12;
    double resid = y[t] - Q_t, hit = (y[t] <= Q_t) ? 1.0 : 0.0;
    neg_ll -= std::log((tau - 1.0) / ES_t) + resid * (tau - hit) / (tau * ES_t);
    if (!std::isfinite(neg_ll)) return 1.0e12;
  }
  return neg_ll;
}

// ============================================================================
// 5–6. Standalone CAViaR-MIDAS recursions (unchanged)
// ============================================================================
// [[Rcpp::export]]
List caviar_midas_sav_cpp(double alpha, double beta_ar, NumericVector y,
                          NumericVector LR_daily, double q_init) {
  const int TT = y.size(); NumericVector Q(TT), SR(TT);
  SR[0] = q_init - LR_daily[0]; Q[0] = q_init;
  for (int t = 1; t < TT; t++) {
    SR[t] = alpha * std::fabs(y[t-1]) + beta_ar * SR[t-1];
    Q[t]  = LR_daily[t] + SR[t];
  }
  return List::create(Named("Q")=Q, Named("SR")=SR);
}

// [[Rcpp::export]]
List caviar_midas_as_cpp(double alpha_pos, double alpha_neg, double beta_ar,
                         NumericVector y, NumericVector LR_daily, double q_init) {
  const int TT = y.size(); NumericVector Q(TT), SR(TT);
  SR[0] = q_init - LR_daily[0]; Q[0] = q_init;
  for (int t = 1; t < TT; t++) {
    double ylag = y[t-1];
    double ypos = (ylag > 0.0) ? ylag : 0.0, yneg = (ylag < 0.0) ? -ylag : 0.0;
    SR[t] = alpha_pos * ypos + alpha_neg * yneg + beta_ar * SR[t-1];
    Q[t]  = LR_daily[t] + SR[t];
  }
  return List::create(Named("Q")=Q, Named("SR")=SR);
}

// ============================================================================
// 7. Standalone Beta weights
// ============================================================================
// [[Rcpp::export]]
NumericVector beta_weights_cpp(int S, double w1, double w2) {
  auto wt = beta_weights_vec(S, w1, w2);
  NumericVector out(S);
  for (int s = 0; s < S; s++) out[s] = wt[s];
  return out;
}

// ============================================================================
// 8. Monthly LR component (post-estimation)
// ============================================================================
// [[Rcpp::export]]
NumericVector compute_LR_monthly_cpp(
    NumericMatrix factors, NumericVector psi_vec,
    double w1, double w2, int S, double phi, bool is_exponential, int trim)
{
  const int T_month = factors.nrow(), K = factors.ncol();
  std::vector<double> psi(K);
  for (int k = 0; k < K; k++) psi[k] = psi_vec[k];
  auto wt = beta_weights_vec(S, w1, w2);
  auto LR = compute_LR_monthly(factors.begin(), T_month, K, S, psi.data(), wt, phi, is_exponential, trim);
  int LR_len = (int)LR.size();
  NumericVector out(LR_len);
  for (int t = 0; t < LR_len; t++) out[t] = LR[t];
  return out;
}

// ============================================================================
// 9. Fused check-loss — Linear SAV
// ============================================================================
// [[Rcpp::export]]
double rq_loss_midas_linear_sav_cpp(
    NumericVector theta_var, NumericVector y, NumericMatrix factors,
    IntegerVector month_id, double tau, double q_init,
    int S, int n_burn, bool fix_w1, int trim)
{
  const int TT = y.size(), T_month = factors.nrow(), K = factors.ncol();
  int idx = 0; double phi = theta_var[idx++];
  std::vector<double> psi(K);
  for (int k = 0; k < K; k++) psi[k] = theta_var[idx++];
  double w1, w2;
  if (!unpack_w(theta_var, idx, fix_w1, w1, w2)) return 1.0e12;
  double alpha = theta_var[idx++], beta_ar = theta_var[idx++];
  if (beta_ar < 0.0 || beta_ar >= 1.0) return 1.0e12;
  auto wt   = beta_weights_vec(S, w1, w2);
  auto LR_m = compute_LR_monthly(factors.begin(), T_month, K, S, psi.data(), wt, phi, false, trim);
  const int start = n_burn + 1;
  double SR_prev = q_init - LR_m[month_id[0]-1], loss = 0.0;
  for (int t = 1; t < TT; t++) {
    double SR_t = alpha * std::fabs(y[t-1]) + beta_ar * SR_prev;
    double Q_t  = LR_m[month_id[t]-1] + SR_t; SR_prev = SR_t;
    if (t < start) continue;
    double u = y[t] - Q_t;
    loss += u * (tau - ((u < 0.0) ? 1.0 : 0.0));
  }
  return loss;
}

// ============================================================================
// 10. Fused check-loss — Linear AS
// ============================================================================
// [[Rcpp::export]]
double rq_loss_midas_linear_as_cpp(
    NumericVector theta_var, NumericVector y, NumericMatrix factors,
    IntegerVector month_id, double tau, double q_init,
    int S, int n_burn, bool fix_w1, int trim)
{
  const int TT = y.size(), T_month = factors.nrow(), K = factors.ncol();
  int idx = 0; double phi = theta_var[idx++];
  std::vector<double> psi(K);
  for (int k = 0; k < K; k++) psi[k] = theta_var[idx++];
  double w1, w2;
  if (!unpack_w(theta_var, idx, fix_w1, w1, w2)) return 1.0e12;
  double alpha_pos = theta_var[idx++], alpha_neg = theta_var[idx++], beta_ar = theta_var[idx++];
  if (beta_ar < 0.0 || beta_ar >= 1.0) return 1.0e12;
  auto wt   = beta_weights_vec(S, w1, w2);
  auto LR_m = compute_LR_monthly(factors.begin(), T_month, K, S, psi.data(), wt, phi, false, trim);
  const int start = n_burn + 1;
  double SR_prev = q_init - LR_m[month_id[0]-1], loss = 0.0;
  for (int t = 1; t < TT; t++) {
    double ylag = y[t-1];
    double ypos = (ylag>0.0)?ylag:0.0, yneg = (ylag<0.0)?-ylag:0.0;
    double SR_t = alpha_pos*ypos + alpha_neg*yneg + beta_ar*SR_prev;
    double Q_t  = LR_m[month_id[t]-1] + SR_t; SR_prev = SR_t;
    if (t < start) continue;
    double u = y[t] - Q_t;
    loss += u * (tau - ((u < 0.0) ? 1.0 : 0.0));
  }
  return loss;
}
