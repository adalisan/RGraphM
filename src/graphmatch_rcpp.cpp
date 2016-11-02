
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
using namespace Rcpp;

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
  CharacterVector param_names = algorithm_params.names();
  int param_count  = param_names.size();

  for (int i=0; i < param_count; ++i) {
  	Rcpp::String key = param_names[i];
  	if (!(Rf_isNull(algorithm_params[key]))){
  	  if ( Rf_isNumeric(algorithm_params[key]) && !Rf_isInteger(algorithm_params[key]) ) {
  	 double value = as<double>( algorithm_params[key]);
     exp.set_param(key,value);
  	  } else if ( Rf_isInteger(algorithm_params[key]) ) {
  	  	int value = as<int> ( algorithm_params[key]);
  	  }  	  	else if ( Rf_isString( (algorithm_params[key]) )) {
  	  	string value = as<string>( algorithm_params[key]);
  	  	exp.set_param(key,value);
  	  }
  	}
  }
  //exp.read_config(conc_params_string);

	try {
  exp.run_experiment(graphm_obj_A,graphm_obj_B);

	} catch(std::exception &ex) {
		forward_exception_to_r(ex);
	} catch(...) {
  return Rcpp::List();

}

