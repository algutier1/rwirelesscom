// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sinc
NumericVector sinc(NumericVector x);
RcppExport SEXP rwirelesscom_sinc(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    __result = Rcpp::wrap(sinc(x));
    return __result;
END_RCPP
}
// rcosine
NumericVector rcosine(NumericVector x, NumericVector B, NumericVector Ns);
RcppExport SEXP rwirelesscom_rcosine(SEXP xSEXP, SEXP BSEXP, SEXP NsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ns(NsSEXP);
    __result = Rcpp::wrap(rcosine(x, B, Ns));
    return __result;
END_RCPP
}
// sqrtrcosine
NumericVector sqrtrcosine(NumericVector x, NumericVector B, NumericVector Ns);
RcppExport SEXP rwirelesscom_sqrtrcosine(SEXP xSEXP, SEXP BSEXP, SEXP NsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ns(NsSEXP);
    __result = Rcpp::wrap(sqrtrcosine(x, B, Ns));
    return __result;
END_RCPP
}
