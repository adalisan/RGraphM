
#include <Rcpp.h>
#include <RcppGSL.h>

#include <iostream>
#include <stdlib.h>

#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

#include "graphm/graph.h"
#include "graphm/experiment.h"
#include "graphm/algorithm.h"

using namespace std;
//using namespace Rcpp;

// declare a dependency on the RcppGSL package; also activates plugin
// (but not needed when ’LinkingTo: RcppGSL’ is used with a package)
//
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::export]]
Rcpp::List run_graph_match(const RcppGSL::Matrix& A, const RcppGSL::Matrix& B, const Rcpp::List& algorithm_params){
  graph graphm_obj_A(A);
  graph graph_test("");
  //graphm_obj_A.set_adjmatrix(A);
  graph graphm_obj_B(B);
  //graphm_obj_B.set_adjmatrix(B);
  experiment exp;
  //exp.read_config(conc_params_string);
  exp.run_experiment();
  return Rcpp::List();

}

