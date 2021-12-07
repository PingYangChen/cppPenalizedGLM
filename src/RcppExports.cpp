#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

