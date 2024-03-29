// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/targeted.h"
#include "../inst/include/targeted_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// NB
Rcpp::List NB(arma::vec y, arma::mat x, arma::uvec xlev, arma::vec ylev, arma::vec weights, double laplacesmooth);
static SEXP _targeted_NB_try(SEXP ySEXP, SEXP xSEXP, SEXP xlevSEXP, SEXP ylevSEXP, SEXP weightsSEXP, SEXP laplacesmoothSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xlev(xlevSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ylev(ylevSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type laplacesmooth(laplacesmoothSEXP);
    rcpp_result_gen = Rcpp::wrap(NB(y, x, xlev, ylev, weights, laplacesmooth));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_NB(SEXP ySEXP, SEXP xSEXP, SEXP xlevSEXP, SEXP ylevSEXP, SEXP weightsSEXP, SEXP laplacesmoothSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_NB_try(ySEXP, xSEXP, xlevSEXP, ylevSEXP, weightsSEXP, laplacesmoothSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// predNB
arma::mat predNB(arma::mat const& X, Rcpp::List const& condprob, Rcpp::List const& xord, arma::uvec multinomial, arma::vec prior, double threshold);
static SEXP _targeted_predNB_try(SEXP XSEXP, SEXP condprobSEXP, SEXP xordSEXP, SEXP multinomialSEXP, SEXP priorSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::List const& >::type condprob(condprobSEXP);
    Rcpp::traits::input_parameter< Rcpp::List const& >::type xord(xordSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type multinomial(multinomialSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(predNB(X, condprob, xord, multinomial, prior, threshold));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_predNB(SEXP XSEXP, SEXP condprobSEXP, SEXP xordSEXP, SEXP multinomialSEXP, SEXP priorSEXP, SEXP thresholdSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_predNB_try(XSEXP, condprobSEXP, xordSEXP, multinomialSEXP, priorSEXP, thresholdSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ode_solve
arma::mat ode_solve(odeptr_t f, arma::mat input, arma::mat init, arma::mat par);
static SEXP _targeted_ode_solve_try(SEXP fSEXP, SEXP inputSEXP, SEXP initSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< odeptr_t >::type f(fSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type input(inputSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type init(initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(ode_solve(f, input, init, par));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_ode_solve(SEXP fSEXP, SEXP inputSEXP, SEXP initSEXP, SEXP parSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_ode_solve_try(fSEXP, inputSEXP, initSEXP, parSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ode_solve2
arma::mat ode_solve2(Rcpp::Function f, arma::mat input, arma::mat init, arma::mat par);
static SEXP _targeted_ode_solve2_try(SEXP fSEXP, SEXP inputSEXP, SEXP initSEXP, SEXP parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type input(inputSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type init(initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type par(parSEXP);
    rcpp_result_gen = Rcpp::wrap(ode_solve2(f, input, init, par));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_ode_solve2(SEXP fSEXP, SEXP inputSEXP, SEXP initSEXP, SEXP parSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_ode_solve2_try(fSEXP, inputSEXP, initSEXP, parSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_logl
arma::vec bin_logl(const arma::vec& y, const arma::vec& a, const arma::mat& x1, const arma::mat& x2, const arma::vec par, const arma::vec& weights, std::string type, bool indiv);
static SEXP _targeted_bin_logl_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type indiv(indivSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_logl(y, a, x1, x2, par, weights, type, indiv));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_bin_logl(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_bin_logl_try(ySEXP, aSEXP, x1SEXP, x2SEXP, parSEXP, weightsSEXP, typeSEXP, indivSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_dlogl
arma::mat bin_dlogl(const arma::vec& y, const arma::vec& a, const arma::mat& x1, const arma::mat& x2, const arma::vec par, const arma::vec& weights, std::string type, bool indiv);
static SEXP _targeted_bin_dlogl_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type indiv(indivSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_dlogl(y, a, x1, x2, par, weights, type, indiv));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_bin_dlogl(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_bin_dlogl_try(ySEXP, aSEXP, x1SEXP, x2SEXP, parSEXP, weightsSEXP, typeSEXP, indivSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_pa
arma::mat bin_pa(const arma::vec& y, const arma::vec& a, const arma::mat& x1, const arma::mat& x2, const arma::vec par, std::string type);
static SEXP _targeted_bin_pa_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_pa(y, a, x1, x2, par, type));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_bin_pa(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP typeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_bin_pa_try(ySEXP, aSEXP, x1SEXP, x2SEXP, parSEXP, typeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_dlogl_c
arma::cx_mat bin_dlogl_c(const arma::cx_vec& y, const arma::cx_vec& a, const arma::cx_mat& x1, const arma::cx_mat& x2, const arma::cx_vec par, const arma::cx_vec& weights, std::string type, bool indiv);
static SEXP _targeted_bin_dlogl_c_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type indiv(indivSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_dlogl_c(y, a, x1, x2, par, weights, type, indiv));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_bin_dlogl_c(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP, SEXP indivSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_bin_dlogl_c_try(ySEXP, aSEXP, x1SEXP, x2SEXP, parSEXP, weightsSEXP, typeSEXP, indivSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_esteq
arma::mat bin_esteq(const arma::vec& y, const arma::vec& a, const arma::mat& x1, const arma::mat& x2, const arma::vec& pr, arma::vec alpha, arma::vec par, const arma::vec& weights, std::string type);
static SEXP _targeted_bin_esteq_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP prSEXP, SEXP alphaSEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pr(prSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_esteq(y, a, x1, x2, pr, alpha, par, weights, type));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_bin_esteq(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP prSEXP, SEXP alphaSEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_bin_esteq_try(ySEXP, aSEXP, x1SEXP, x2SEXP, prSEXP, alphaSEXP, parSEXP, weightsSEXP, typeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bin_esteq_c
arma::cx_mat bin_esteq_c(const arma::cx_vec& y, const arma::cx_vec& a, const arma::cx_mat& x1, const arma::cx_mat& x2, const arma::cx_mat& x3, arma::cx_vec alpha, arma::cx_vec par, const arma::cx_vec& weights, std::string type);
static SEXP _targeted_bin_esteq_c_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP x3SEXP, SEXP alphaSEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type x3(x3SEXP);
    Rcpp::traits::input_parameter< arma::cx_vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::cx_vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_esteq_c(y, a, x1, x2, x3, alpha, par, weights, type));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_bin_esteq_c(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP x3SEXP, SEXP alphaSEXP, SEXP parSEXP, SEXP weightsSEXP, SEXP typeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_bin_esteq_c_try(ySEXP, aSEXP, x1SEXP, x2SEXP, x3SEXP, alphaSEXP, parSEXP, weightsSEXP, typeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// ace_est
Rcpp::List ace_est(const arma::vec& y, const arma::mat& a, const arma::mat& x1, const arma::mat& x2, const arma::vec& theta, const arma::vec& weights, const arma::vec& offset, std::string link);
static SEXP _targeted_ace_est_try(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP thetaSEXP, SEXP weightsSEXP, SEXP offsetSEXP, SEXP linkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< std::string >::type link(linkSEXP);
    rcpp_result_gen = Rcpp::wrap(ace_est(y, a, x1, x2, theta, weights, offset, link));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_ace_est(SEXP ySEXP, SEXP aSEXP, SEXP x1SEXP, SEXP x2SEXP, SEXP thetaSEXP, SEXP weightsSEXP, SEXP offsetSEXP, SEXP linkSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_ace_est_try(ySEXP, aSEXP, x1SEXP, x2SEXP, thetaSEXP, weightsSEXP, offsetSEXP, linkSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// fast_iid
arma::mat fast_iid(const arma::vec& y, const arma::vec& p, const arma::mat& x1, const arma::vec& weights, bool logistic);
static SEXP _targeted_fast_iid_try(SEXP ySEXP, SEXP pSEXP, SEXP x1SEXP, SEXP weightsSEXP, SEXP logisticSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type logistic(logisticSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_iid(y, p, x1, weights, logistic));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_fast_iid(SEXP ySEXP, SEXP pSEXP, SEXP x1SEXP, SEXP weightsSEXP, SEXP logisticSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_fast_iid_try(ySEXP, pSEXP, x1SEXP, weightsSEXP, logisticSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// clusterid
List clusterid(const arma::uvec& id);
static SEXP _targeted_clusterid_try(SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type id(idSEXP);
    rcpp_result_gen = Rcpp::wrap(clusterid(id));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_clusterid(SEXP idSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_clusterid_try(idSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// groupsum
arma::mat groupsum(const arma::mat& x, const arma::uvec& cluster, bool reduce);
static SEXP _targeted_groupsum_try(SEXP xSEXP, SEXP clusterSEXP, SEXP reduceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type cluster(clusterSEXP);
    Rcpp::traits::input_parameter< bool >::type reduce(reduceSEXP);
    rcpp_result_gen = Rcpp::wrap(groupsum(x, cluster, reduce));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_groupsum(SEXP xSEXP, SEXP clusterSEXP, SEXP reduceSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_groupsum_try(xSEXP, clusterSEXP, reduceSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// softmax
arma::mat softmax(arma::mat& lp, bool ref, bool log);
static SEXP _targeted_softmax_try(SEXP lpSEXP, SEXP refSEXP, SEXP logSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type lp(lpSEXP);
    Rcpp::traits::input_parameter< bool >::type ref(refSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    rcpp_result_gen = Rcpp::wrap(softmax(lp, ref, log));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_softmax(SEXP lpSEXP, SEXP refSEXP, SEXP logSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_softmax_try(lpSEXP, refSEXP, logSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// pava
List pava(const arma::vec& y, const NumericVector& x, const NumericVector& weights);
static SEXP _targeted_pava_try(SEXP ySEXP, SEXP xSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(pava(y, x, weights));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_pava(SEXP ySEXP, SEXP xSEXP, SEXP weightsSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_pava_try(ySEXP, xSEXP, weightsSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// nondom
arma::mat nondom(const arma::mat& x);
static SEXP _targeted_nondom_try(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(nondom(x));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _targeted_nondom(SEXP xSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_targeted_nondom_try(xSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _targeted_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Rcpp::List(*.NB)(arma::vec,arma::mat,arma::uvec,arma::vec,arma::vec,double)");
        signatures.insert("arma::mat(*.predNB)(arma::mat const&,Rcpp::List const&,Rcpp::List const&,arma::uvec,arma::vec,double)");
        signatures.insert("arma::mat(*.ode_solve)(odeptr_t,arma::mat,arma::mat,arma::mat)");
        signatures.insert("arma::mat(*.ode_solve2)(Rcpp::Function,arma::mat,arma::mat,arma::mat)");
        signatures.insert("arma::vec(*bin_logl)(const arma::vec&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::vec,const arma::vec&,std::string,bool)");
        signatures.insert("arma::mat(*bin_dlogl)(const arma::vec&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::vec,const arma::vec&,std::string,bool)");
        signatures.insert("arma::mat(*bin_pa)(const arma::vec&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::vec,std::string)");
        signatures.insert("arma::cx_mat(*bin_dlogl_c)(const arma::cx_vec&,const arma::cx_vec&,const arma::cx_mat&,const arma::cx_mat&,const arma::cx_vec,const arma::cx_vec&,std::string,bool)");
        signatures.insert("arma::mat(*bin_esteq)(const arma::vec&,const arma::vec&,const arma::mat&,const arma::mat&,const arma::vec&,arma::vec,arma::vec,const arma::vec&,std::string)");
        signatures.insert("arma::cx_mat(*bin_esteq_c)(const arma::cx_vec&,const arma::cx_vec&,const arma::cx_mat&,const arma::cx_mat&,const arma::cx_mat&,arma::cx_vec,arma::cx_vec,const arma::cx_vec&,std::string)");
        signatures.insert("Rcpp::List(*ace_est)(const arma::vec&,const arma::mat&,const arma::mat&,const arma::mat&,const arma::vec&,const arma::vec&,const arma::vec&,std::string)");
        signatures.insert("arma::mat(*fast_iid)(const arma::vec&,const arma::vec&,const arma::mat&,const arma::vec&,bool)");
        signatures.insert("List(*.clusterid)(const arma::uvec&)");
        signatures.insert("arma::mat(*.groupsum)(const arma::mat&,const arma::uvec&,bool)");
        signatures.insert("arma::mat(*.softmax)(arma::mat&,bool,bool)");
        signatures.insert("List(*.pava)(const arma::vec&,const NumericVector&,const NumericVector&)");
        signatures.insert("arma::mat(*.nondom)(const arma::mat&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _targeted_RcppExport_registerCCallable() { 
    R_RegisterCCallable("targeted", "_targeted_.NB", (DL_FUNC)_targeted_NB_try);
    R_RegisterCCallable("targeted", "_targeted_.predNB", (DL_FUNC)_targeted_predNB_try);
    R_RegisterCCallable("targeted", "_targeted_.ode_solve", (DL_FUNC)_targeted_ode_solve_try);
    R_RegisterCCallable("targeted", "_targeted_.ode_solve2", (DL_FUNC)_targeted_ode_solve2_try);
    R_RegisterCCallable("targeted", "_targeted_bin_logl", (DL_FUNC)_targeted_bin_logl_try);
    R_RegisterCCallable("targeted", "_targeted_bin_dlogl", (DL_FUNC)_targeted_bin_dlogl_try);
    R_RegisterCCallable("targeted", "_targeted_bin_pa", (DL_FUNC)_targeted_bin_pa_try);
    R_RegisterCCallable("targeted", "_targeted_bin_dlogl_c", (DL_FUNC)_targeted_bin_dlogl_c_try);
    R_RegisterCCallable("targeted", "_targeted_bin_esteq", (DL_FUNC)_targeted_bin_esteq_try);
    R_RegisterCCallable("targeted", "_targeted_bin_esteq_c", (DL_FUNC)_targeted_bin_esteq_c_try);
    R_RegisterCCallable("targeted", "_targeted_ace_est", (DL_FUNC)_targeted_ace_est_try);
    R_RegisterCCallable("targeted", "_targeted_fast_iid", (DL_FUNC)_targeted_fast_iid_try);
    R_RegisterCCallable("targeted", "_targeted_.clusterid", (DL_FUNC)_targeted_clusterid_try);
    R_RegisterCCallable("targeted", "_targeted_.groupsum", (DL_FUNC)_targeted_groupsum_try);
    R_RegisterCCallable("targeted", "_targeted_.softmax", (DL_FUNC)_targeted_softmax_try);
    R_RegisterCCallable("targeted", "_targeted_.pava", (DL_FUNC)_targeted_pava_try);
    R_RegisterCCallable("targeted", "_targeted_.nondom", (DL_FUNC)_targeted_nondom_try);
    R_RegisterCCallable("targeted", "_targeted_RcppExport_validate", (DL_FUNC)_targeted_RcppExport_validate);
    return R_NilValue;
}

RcppExport SEXP _rcpp_module_boot_riskregmodel();

static const R_CallMethodDef CallEntries[] = {
    {"_targeted_NB", (DL_FUNC) &_targeted_NB, 6},
    {"_targeted_predNB", (DL_FUNC) &_targeted_predNB, 6},
    {"_targeted_ode_solve", (DL_FUNC) &_targeted_ode_solve, 4},
    {"_targeted_ode_solve2", (DL_FUNC) &_targeted_ode_solve2, 4},
    {"_targeted_bin_logl", (DL_FUNC) &_targeted_bin_logl, 8},
    {"_targeted_bin_dlogl", (DL_FUNC) &_targeted_bin_dlogl, 8},
    {"_targeted_bin_pa", (DL_FUNC) &_targeted_bin_pa, 6},
    {"_targeted_bin_dlogl_c", (DL_FUNC) &_targeted_bin_dlogl_c, 8},
    {"_targeted_bin_esteq", (DL_FUNC) &_targeted_bin_esteq, 9},
    {"_targeted_bin_esteq_c", (DL_FUNC) &_targeted_bin_esteq_c, 9},
    {"_targeted_ace_est", (DL_FUNC) &_targeted_ace_est, 8},
    {"_targeted_fast_iid", (DL_FUNC) &_targeted_fast_iid, 5},
    {"_targeted_clusterid", (DL_FUNC) &_targeted_clusterid, 1},
    {"_targeted_groupsum", (DL_FUNC) &_targeted_groupsum, 3},
    {"_targeted_softmax", (DL_FUNC) &_targeted_softmax, 3},
    {"_targeted_pava", (DL_FUNC) &_targeted_pava, 3},
    {"_targeted_nondom", (DL_FUNC) &_targeted_nondom, 1},
    {"_rcpp_module_boot_riskregmodel", (DL_FUNC) &_rcpp_module_boot_riskregmodel, 0},
    {"_targeted_RcppExport_registerCCallable", (DL_FUNC) &_targeted_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_targeted(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
