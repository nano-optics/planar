// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// integrand_collection
double integrand_collection(const arma::colvec& rt, const arma::colvec& r2, const double k0, const double psi, const arma::cx_vec& epsilon, const arma::vec& thickness);
RcppExport SEXP _planar_integrand_collection(SEXP rtSEXP, SEXP r2SEXP, SEXP k0SEXP, SEXP psiSEXP, SEXP epsilonSEXP, SEXP thicknessSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< const double >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type thickness(thicknessSEXP);
    rcpp_result_gen = Rcpp::wrap(integrand_collection(rt, r2, k0, psi, epsilon, thickness));
    return rcpp_result_gen;
END_RCPP
}
// cpp_field_collection
arma::vec cpp_field_collection(const arma::mat& r2, const double k0, const double psi, const arma::vec& omega, const arma::cx_vec& epsilon, const arma::vec& thickness, const int maxEval, const double reqAbsError, const double tol, bool progress);
RcppExport SEXP _planar_cpp_field_collection(SEXP r2SEXP, SEXP k0SEXP, SEXP psiSEXP, SEXP omegaSEXP, SEXP epsilonSEXP, SEXP thicknessSEXP, SEXP maxEvalSEXP, SEXP reqAbsErrorSEXP, SEXP tolSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< const double >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type thickness(thicknessSEXP);
    Rcpp::traits::input_parameter< const int >::type maxEval(maxEvalSEXP);
    Rcpp::traits::input_parameter< const double >::type reqAbsError(reqAbsErrorSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_field_collection(r2, k0, psi, omega, epsilon, thickness, maxEval, reqAbsError, tol, progress));
    return rcpp_result_gen;
END_RCPP
}
// cpp_layer_fresnel
Rcpp::List cpp_layer_fresnel(const arma::colvec& k0, const arma::cx_mat& kx, const arma::cx_mat& epsilon, const double& thickness);
RcppExport SEXP _planar_cpp_layer_fresnel(SEXP k0SEXP, SEXP kxSEXP, SEXP epsilonSEXP, SEXP thicknessSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type kx(kxSEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const double& >::type thickness(thicknessSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_layer_fresnel(k0, kx, epsilon, thickness));
    return rcpp_result_gen;
END_RCPP
}
// cpp_recursive_fresnel
Rcpp::List cpp_recursive_fresnel(const arma::colvec& k0, const arma::cx_mat& kx, const arma::cx_mat& epsilon, const arma::colvec& thickness, const int& polarisation);
RcppExport SEXP _planar_cpp_recursive_fresnel(SEXP k0SEXP, SEXP kxSEXP, SEXP epsilonSEXP, SEXP thicknessSEXP, SEXP polarisationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type kx(kxSEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type thickness(thicknessSEXP);
    Rcpp::traits::input_parameter< const int& >::type polarisation(polarisationSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_recursive_fresnel(k0, kx, epsilon, thickness, polarisation));
    return rcpp_result_gen;
END_RCPP
}
// cpp_integrand_gb_ml
arma::colvec cpp_integrand_gb_ml(const arma::colvec& rt, const arma::colvec& r2, const double k0, const double psi, const double alpha, const double w0, const arma::cx_vec& epsilon, const arma::vec& thickness);
RcppExport SEXP _planar_cpp_integrand_gb_ml(SEXP rtSEXP, SEXP r2SEXP, SEXP k0SEXP, SEXP psiSEXP, SEXP alphaSEXP, SEXP w0SEXP, SEXP epsilonSEXP, SEXP thicknessSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< const double >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type thickness(thicknessSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_integrand_gb_ml(rt, r2, k0, psi, alpha, w0, epsilon, thickness));
    return rcpp_result_gen;
END_RCPP
}
// cpp_integrand_gb_layer
arma::colvec cpp_integrand_gb_layer(const arma::colvec& rt, const arma::colvec& r2, const double ki, const double psi, const double alpha, const double w0, const double ni, const double no, const arma::cx_double nl, const double d);
RcppExport SEXP _planar_cpp_integrand_gb_layer(SEXP rtSEXP, SEXP r2SEXP, SEXP kiSEXP, SEXP psiSEXP, SEXP alphaSEXP, SEXP w0SEXP, SEXP niSEXP, SEXP noSEXP, SEXP nlSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< const double >::type ki(kiSEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< const double >::type ni(niSEXP);
    Rcpp::traits::input_parameter< const double >::type no(noSEXP);
    Rcpp::traits::input_parameter< const arma::cx_double >::type nl(nlSEXP);
    Rcpp::traits::input_parameter< const double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_integrand_gb_layer(rt, r2, ki, psi, alpha, w0, ni, no, nl, d));
    return rcpp_result_gen;
END_RCPP
}
// cpp_field_gb_layer
arma::cx_mat cpp_field_gb_layer(const arma::mat& r2, const double k0, const double psi, const double alpha, const double w0, const arma::cx_vec& epsilon, const arma::vec& thickness, const int maxEval, const double reqAbsError, const double tol, bool progress);
RcppExport SEXP _planar_cpp_field_gb_layer(SEXP r2SEXP, SEXP k0SEXP, SEXP psiSEXP, SEXP alphaSEXP, SEXP w0SEXP, SEXP epsilonSEXP, SEXP thicknessSEXP, SEXP maxEvalSEXP, SEXP reqAbsErrorSEXP, SEXP tolSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< const double >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type thickness(thicknessSEXP);
    Rcpp::traits::input_parameter< const int >::type maxEval(maxEvalSEXP);
    Rcpp::traits::input_parameter< const double >::type reqAbsError(reqAbsErrorSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_field_gb_layer(r2, k0, psi, alpha, w0, epsilon, thickness, maxEval, reqAbsError, tol, progress));
    return rcpp_result_gen;
END_RCPP
}
// cpp_field_gb_ml
arma::cx_mat cpp_field_gb_ml(const arma::mat& r2, const double k0, const double psi, const double alpha, const double w0, const arma::cx_vec& epsilon, const arma::vec& thickness, const int maxEval, const double reqAbsError, const double tol, bool progress);
RcppExport SEXP _planar_cpp_field_gb_ml(SEXP r2SEXP, SEXP k0SEXP, SEXP psiSEXP, SEXP alphaSEXP, SEXP w0SEXP, SEXP epsilonSEXP, SEXP thicknessSEXP, SEXP maxEvalSEXP, SEXP reqAbsErrorSEXP, SEXP tolSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< const double >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type thickness(thicknessSEXP);
    Rcpp::traits::input_parameter< const int >::type maxEval(maxEvalSEXP);
    Rcpp::traits::input_parameter< const double >::type reqAbsError(reqAbsErrorSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_field_gb_ml(r2, k0, psi, alpha, w0, epsilon, thickness, maxEval, reqAbsError, tol, progress));
    return rcpp_result_gen;
END_RCPP
}
// cpp_multilayer
Rcpp::List cpp_multilayer(const arma::colvec& k0, const arma::cx_mat& kx, const arma::cx_mat& epsilon, const arma::colvec& thickness, const arma::colvec& z, const double psi, const bool intensity);
RcppExport SEXP _planar_cpp_multilayer(SEXP k0SEXP, SEXP kxSEXP, SEXP epsilonSEXP, SEXP thicknessSEXP, SEXP zSEXP, SEXP psiSEXP, SEXP intensitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type kx(kxSEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type thickness(thicknessSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const bool >::type intensity(intensitySEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_multilayer(k0, kx, epsilon, thickness, z, psi, intensity));
    return rcpp_result_gen;
END_RCPP
}
// cpp_multilayer_field
Rcpp::List cpp_multilayer_field(const double k0, const double kx, const arma::cx_vec& epsilon, const arma::colvec& thickness, const arma::colvec& z, const double psi);
RcppExport SEXP _planar_cpp_multilayer_field(SEXP k0SEXP, SEXP kxSEXP, SEXP epsilonSEXP, SEXP thicknessSEXP, SEXP zSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type k0(k0SEXP);
    Rcpp::traits::input_parameter< const double >::type kx(kxSEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type thickness(thicknessSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_multilayer_field(k0, kx, epsilon, thickness, z, psi));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_planar_integrand_collection", (DL_FUNC) &_planar_integrand_collection, 6},
    {"_planar_cpp_field_collection", (DL_FUNC) &_planar_cpp_field_collection, 10},
    {"_planar_cpp_layer_fresnel", (DL_FUNC) &_planar_cpp_layer_fresnel, 4},
    {"_planar_cpp_recursive_fresnel", (DL_FUNC) &_planar_cpp_recursive_fresnel, 5},
    {"_planar_cpp_integrand_gb_ml", (DL_FUNC) &_planar_cpp_integrand_gb_ml, 8},
    {"_planar_cpp_integrand_gb_layer", (DL_FUNC) &_planar_cpp_integrand_gb_layer, 10},
    {"_planar_cpp_field_gb_layer", (DL_FUNC) &_planar_cpp_field_gb_layer, 11},
    {"_planar_cpp_field_gb_ml", (DL_FUNC) &_planar_cpp_field_gb_ml, 11},
    {"_planar_cpp_multilayer", (DL_FUNC) &_planar_cpp_multilayer, 7},
    {"_planar_cpp_multilayer_field", (DL_FUNC) &_planar_cpp_multilayer_field, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_planar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
