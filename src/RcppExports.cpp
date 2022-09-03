// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// readsas
Rcpp::List readsas(const char * filePath, const bool debug, Nullable<IntegerVector> selectrows_, Nullable<CharacterVector> selectcols_, const bool empty_to_na);
RcppExport SEXP _readsas_readsas(SEXP filePathSEXP, SEXP debugSEXP, SEXP selectrows_SEXP, SEXP selectcols_SEXP, SEXP empty_to_naSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char * >::type filePath(filePathSEXP);
    Rcpp::traits::input_parameter< const bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type selectrows_(selectrows_SEXP);
    Rcpp::traits::input_parameter< Nullable<CharacterVector> >::type selectcols_(selectcols_SEXP);
    Rcpp::traits::input_parameter< const bool >::type empty_to_na(empty_to_naSEXP);
    rcpp_result_gen = Rcpp::wrap(readsas(filePath, debug, selectrows_, selectcols_, empty_to_na));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_readsas_readsas", (DL_FUNC) &_readsas_readsas, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_readsas(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
