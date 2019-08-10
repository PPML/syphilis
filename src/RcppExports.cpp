// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mv_mult
arma::vec mv_mult(arma::mat& lhs, arma::vec& rhs);
RcppExport SEXP _syphilis_mv_mult(SEXP lhsSEXP, SEXP rhsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type lhs(lhsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type rhs(rhsSEXP);
    rcpp_result_gen = Rcpp::wrap(mv_mult(lhs, rhs));
    return rcpp_result_gen;
END_RCPP
}
// mat_mult
arma::mat mat_mult(arma::mat& a, arma::mat& b);
RcppExport SEXP _syphilis_mat_mult(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_mult(a, b));
    return rcpp_result_gen;
END_RCPP
}
// syphSim
List syphSim(List x, double dt, List cm, NumericVector pabx, NumericVector rep_count, int nYrs, NumericVector initPop, NumericVector n_sa, arma::mat births, arma::mat births_sa, arma::mat births_nsa, arma::mat aging, arma::mat aging_nsa);
RcppExport SEXP _syphilis_syphSim(SEXP xSEXP, SEXP dtSEXP, SEXP cmSEXP, SEXP pabxSEXP, SEXP rep_countSEXP, SEXP nYrsSEXP, SEXP initPopSEXP, SEXP n_saSEXP, SEXP birthsSEXP, SEXP births_saSEXP, SEXP births_nsaSEXP, SEXP agingSEXP, SEXP aging_nsaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< List >::type cm(cmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pabx(pabxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rep_count(rep_countSEXP);
    Rcpp::traits::input_parameter< int >::type nYrs(nYrsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initPop(initPopSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_sa(n_saSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type births(birthsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type births_sa(births_saSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type births_nsa(births_nsaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aging(agingSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aging_nsa(aging_nsaSEXP);
    rcpp_result_gen = Rcpp::wrap(syphSim(x, dt, cm, pabx, rep_count, nYrs, initPop, n_sa, births, births_sa, births_nsa, aging, aging_nsa));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_syphilis_mv_mult", (DL_FUNC) &_syphilis_mv_mult, 2},
    {"_syphilis_mat_mult", (DL_FUNC) &_syphilis_mat_mult, 2},
    {"_syphilis_syphSim", (DL_FUNC) &_syphilis_syphSim, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_syphilis(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
