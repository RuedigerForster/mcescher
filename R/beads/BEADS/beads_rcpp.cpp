// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// beadsbaseline.cpp is compiled as a separate translation unit by sourceCpp;
// include only the header here.
#include "beadsbaseline.h"

//' Baseline Estimation and Denoising with Sparsity (BEADS)
//'
//' Wraps the C++ BEADS algorithm (Ning, Selesnick, Duval 2014) for use in R.
//'
//' @param signal  Numeric vector — raw chromatogram.
//' @param d       Filter order (1 or 2).
//' @param fc      Filter cut-off frequency, cycles/sample (0 < fc < 0.5).
//' @param r       Asymmetry / positivity bias ratio (default 6).
//' @param lam0    Regularisation parameter lambda_0 (sparsity of x).
//' @param lam1    Regularisation parameter lambda_1 (sparsity of Dx).
//' @param lam2    Regularisation parameter lambda_2 (sparsity of D^2 x).
//'
//' @return Numeric vector — estimated baseline, same length as \code{signal}.
//'
// [[Rcpp::export]]
Rcpp::NumericVector beads_baseline(
    Rcpp::NumericVector signal,
    int    d  = 1,
    double fc = 0.1,
    double r  = 6.0,
    double lam0 = 0.4,
    double lam1 = 4.0,
    double lam2 = 3.2
) {
    std::vector<double> y(signal.begin(), signal.end());
    beadsBaseline bl;
    // beads() returns {x (corrected signal), f (baseline)}
    std::vector<std::vector<double>> res = bl.beads(y, d, fc, r, lam0, lam1, lam2);
    return Rcpp::wrap(res[1]);   // res[1] = f, the estimated baseline
}
