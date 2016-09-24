#include <Rcpp.h>
#include <RcppGSL.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include "graph.h"
#include "experiment.h"
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>


using namespace std;
using namespace Rcpp;

// declare a dependency on the RcppGSL package; also activates plugin
// (but not needed when ’LinkingTo: RcppGSL’ is used with a package)
//
// [[Rcpp::depends(RcppGSL)]]

//
//
// [[Rcpp::export]]
Rcpp::List run_graph_match(const RcppGSL::Matrix & A, const RcppGSL::Matrix & B, const Rcpp::List &algorithm_params ){
  graph::graph graphm_obj_A = new graph();
  graphm_obj_A.set_adjmatrix (&A);
  graph::graph graphm_obj_B = new graph();
  graphm_obj_B.set_adjmatrix (&B);
  experiment exp;
  exp.read_config(conc_params_string);
  exp.run_experiment();


}

