// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// apply_transform
Rcpp::List read_obj_str(std::vector< std::string > string);
RcppExport SEXP _svgViewR_read_obj_str(SEXP stringSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector< std::string > >::type string(stringSEXP);
    rcpp_result_gen = Rcpp::wrap(read_obj_str(string));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_svgViewR_read_obj_str", (DL_FUNC) &_svgViewR_read_obj_str, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_svgViewR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
