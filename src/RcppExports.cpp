// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// agree_C
arma::mat agree_C(arma::mat mat_A, arma::mat mat_B);
RcppExport SEXP ludic_agree_C(SEXP mat_ASEXP, SEXP mat_BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mat_A(mat_ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mat_B(mat_BSEXP);
    rcpp_result_gen = Rcpp::wrap(agree_C(mat_A, mat_B));
    return rcpp_result_gen;
END_RCPP
}
// agree_C_sparse
arma::sp_mat agree_C_sparse(arma::mat mat_A, arma::mat mat_B);
RcppExport SEXP ludic_agree_C_sparse(SEXP mat_ASEXP, SEXP mat_BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mat_A(mat_ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mat_B(mat_BSEXP);
    rcpp_result_gen = Rcpp::wrap(agree_C_sparse(mat_A, mat_B));
    return rcpp_result_gen;
END_RCPP
}
// estep_C_vect
arma::mat estep_C_vect(arma::mat agreemat, double p, arma::colvec m, arma::colvec u);
RcppExport SEXP ludic_estep_C_vect(SEXP agreematSEXP, SEXP pSEXP, SEXP mSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type agreemat(agreematSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(estep_C_vect(agreemat, p, m, u));
    return rcpp_result_gen;
END_RCPP
}
// EMstep_C_sparse_big
List EMstep_C_sparse_big(arma::mat mat_A, arma::mat mat_B, double p, arma::rowvec m, arma::rowvec u);
RcppExport SEXP ludic_EMstep_C_sparse_big(SEXP mat_ASEXP, SEXP mat_BSEXP, SEXP pSEXP, SEXP mSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mat_A(mat_ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mat_B(mat_BSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(EMstep_C_sparse_big(mat_A, mat_B, p, m, u));
    return rcpp_result_gen;
END_RCPP
}
// loglikC_bin
NumericMatrix loglikC_bin(arma::mat Bmat, arma::mat Amat, NumericVector eps_p, NumericVector eps_n, NumericVector piA, NumericVector piB);
RcppExport SEXP ludic_loglikC_bin(SEXP BmatSEXP, SEXP AmatSEXP, SEXP eps_pSEXP, SEXP eps_nSEXP, SEXP piASEXP, SEXP piBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Bmat(BmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eps_p(eps_pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eps_n(eps_nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type piA(piASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type piB(piBSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikC_bin(Bmat, Amat, eps_p, eps_n, piA, piB));
    return rcpp_result_gen;
END_RCPP
}
// strsplitC
std::vector<std::string> strsplitC(std::string s, char sep);
RcppExport SEXP ludic_strsplitC(SEXP sSEXP, SEXP sepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s(sSEXP);
    Rcpp::traits::input_parameter< char >::type sep(sepSEXP);
    rcpp_result_gen = Rcpp::wrap(strsplitC(s, sep));
    return rcpp_result_gen;
END_RCPP
}
// loglikC_bin_wDates
NumericMatrix loglikC_bin_wDates(arma::mat Bmat, arma::mat Amat, StringMatrix Bdates, StringMatrix Adates, NumericVector eps_p, NumericVector eps_n, NumericVector piA, NumericVector piB);
RcppExport SEXP ludic_loglikC_bin_wDates(SEXP BmatSEXP, SEXP AmatSEXP, SEXP BdatesSEXP, SEXP AdatesSEXP, SEXP eps_pSEXP, SEXP eps_nSEXP, SEXP piASEXP, SEXP piBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Bmat(BmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< StringMatrix >::type Bdates(BdatesSEXP);
    Rcpp::traits::input_parameter< StringMatrix >::type Adates(AdatesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eps_p(eps_pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eps_n(eps_nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type piA(piASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type piB(piBSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikC_bin_wDates(Bmat, Amat, Bdates, Adates, eps_p, eps_n, piA, piB));
    return rcpp_result_gen;
END_RCPP
}
// loglikratioC_diff_arbitrary
NumericMatrix loglikratioC_diff_arbitrary(arma::mat Bmat, arma::mat Amat, NumericVector d_max, NumericVector cost);
RcppExport SEXP ludic_loglikratioC_diff_arbitrary(SEXP BmatSEXP, SEXP AmatSEXP, SEXP d_maxSEXP, SEXP costSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Bmat(BmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d_max(d_maxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cost(costSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikratioC_diff_arbitrary(Bmat, Amat, d_max, cost));
    return rcpp_result_gen;
END_RCPP
}
// matchingScore_C
arma::mat matchingScore_C(arma::mat agreemat, arma::vec m, arma::vec u, int nA, int nB);
RcppExport SEXP ludic_matchingScore_C(SEXP agreematSEXP, SEXP mSEXP, SEXP uSEXP, SEXP nASEXP, SEXP nBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type agreemat(agreematSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< int >::type nA(nASEXP);
    Rcpp::traits::input_parameter< int >::type nB(nBSEXP);
    rcpp_result_gen = Rcpp::wrap(matchingScore_C(agreemat, m, u, nA, nB));
    return rcpp_result_gen;
END_RCPP
}
// matchingScore_C_sparse_big
arma::mat matchingScore_C_sparse_big(arma::mat mat_A, arma::mat mat_B, arma::vec m, arma::vec u);
RcppExport SEXP ludic_matchingScore_C_sparse_big(SEXP mat_ASEXP, SEXP mat_BSEXP, SEXP mSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mat_A(mat_ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mat_B(mat_BSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(matchingScore_C_sparse_big(mat_A, mat_B, m, u));
    return rcpp_result_gen;
END_RCPP
}
// matchProbs_rank_full_C
NumericMatrix matchProbs_rank_full_C(NumericMatrix computed_dist, double prop_match);
RcppExport SEXP ludic_matchProbs_rank_full_C(SEXP computed_distSEXP, SEXP prop_matchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type computed_dist(computed_distSEXP);
    Rcpp::traits::input_parameter< double >::type prop_match(prop_matchSEXP);
    rcpp_result_gen = Rcpp::wrap(matchProbs_rank_full_C(computed_dist, prop_match));
    return rcpp_result_gen;
END_RCPP
}
