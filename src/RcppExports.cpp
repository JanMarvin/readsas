// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// readsas
Rcpp::List readsas(const char * filePath, const bool debug, const IntegerVector selectrows, const CharacterVector selectcols);
RcppExport SEXP _readsas_readsas(SEXP filePathSEXP, SEXP debugSEXP, SEXP selectrowsSEXP, SEXP selectcolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char * >::type filePath(filePathSEXP);
    Rcpp::traits::input_parameter< const bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type selectrows(selectrowsSEXP);
    Rcpp::traits::input_parameter< const CharacterVector >::type selectcols(selectcolsSEXP);
    rcpp_result_gen = Rcpp::wrap(readsas(filePath, debug, selectrows, selectcols));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_readsas_readsas", (DL_FUNC) &_readsas_readsas, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_readsas(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
