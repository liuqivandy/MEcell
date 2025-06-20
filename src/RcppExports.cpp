// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// find_knn_rcpp
List find_knn_rcpp(NumericMatrix mat, IntegerMatrix candidates_mat);
RcppExport SEXP _MEcell_find_knn_rcpp(SEXP matSEXP, SEXP candidates_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type candidates_mat(candidates_matSEXP);
    rcpp_result_gen = Rcpp::wrap(find_knn_rcpp(mat, candidates_mat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MEcell_find_knn_rcpp", (DL_FUNC) &_MEcell_find_knn_rcpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MEcell(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
