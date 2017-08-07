// # Copyright (C) Sancar Adali 2015-2017

#include <Rcpp.h>
#include <RcppGSL.h>

#include <iostream>
#include <stdlib.h>

#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <stdexcept>
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
// [[Rcpp::depends(RcppGSL,Rcpp)]]
// [[Rcpp::export]]
Rcpp::List run_graph_match(const RcppGSL::Matrix& A, const RcppGSL::Matrix& B, const Rcpp::List& algorithm_params){
	//check if adj matrices are square
	if (A.ncol()!= A.nrow()){
		Rf_error("A is not a square matrix (eq. number of rows and columns)!");
		return Rcpp::List();
	}
	if (B.ncol()!= B.nrow()){
		Rf_error("B is not a square matrix (eq. number of rows and columns)!");
		return Rcpp::List();
	}
  graph graphm_obj_A(A);
	//print("Read in matrix for graph A")
  //graphm_obj_A.printout("graphA.txt");
  //graphm_obj_A.set_adjmatrix(A);
  graph graphm_obj_B(B);
  //print("Read in matrix for graph B")
  //graphm_obj_B.set_adjmatrix(B);
  //graphm_obj_B.printout("graphB.txt");
  // ofstream debug_out("gm_rcpp_debug.txt", ofstream::out);

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
  	  	  exp.set_param(key,value);
  	  } else if ( Rf_isString( (algorithm_params[key]) )) {

  	  	string tmp_str(as<string> (algorithm_params[key]));



  //       debug_out << " \n ";
  // 	  	for( unsigned int i=0; i < tmp_str.length(); i++ ){
  // 	  		debug_out << "i is: " << i << ", the element value is: " << tmp_str[i];
  // 	  		debug_out << "\n";
  // 	  	}

  	  	if ( tmp_str.size() == 1){
  	  	  char value = tmp_str[0];
  	  	  exp.set_param(key,value);
  	  	}	else {
  	  		//string value = as<string>( algorithm_params[key]);
  	  		exp.set_param(key,tmp_str);
  	  	}
  	  }
  	}
  }
  int expected_result_count=0;
  vector<string> algo_list;
  if (!Rf_isNull(algorithm_params["algo"]) ) {
  	string str(as<string>(algorithm_params["algo"]));
  	string buf; // Have a buffer string
  	stringstream ss(str); // Insert the string into a stream
  	 // Create vector to hold our words

  	while (ss >> buf)
  		algo_list.push_back(buf);

  	expected_result_count = algo_list.size();

  } else{
    //algorithm_params["algo"] = "PATH";
  	//algorithm_params["algo_init_sol"] = "unif";

  	exp.set_param("algo",std::string("I"));
  	exp.set_param("algo_init_sol",std::string("unif"));
  	expected_result_count = 1;
  	algo_list.push_back("I");
  }
  Rcpp::Rcout << "Running Graph matching algorithms with following algorithms\n";
  for (vector<string>::iterator algo_it = algo_list.begin() ; algo_it != algo_list.end(); algo_it++)
  	Rcpp::Rcout << *algo_it <<"\n";

  //exp.read_config(conc_params_string);
  //exp.printout("before_Exp.txt","Parameters");
  int P_nr = max(A.nrow(),B.ncol());

  std::vector<Rcpp::NumericMatrix> P_matrix_list ;
  std::vector<Rcpp::NumericVector> P_vector_list ;
  try {

  	exp.run_experiment(graphm_obj_A, graphm_obj_B);

  	int exp_count = exp.get_algo_len();

  	//exp.printout("experiment ran.txt");

    //Rf_warning("exp count is %d",exp_count);

   	for (unsigned int exp_i = 0;  exp_i < exp_count; exp_i++) {
   		Rcpp::NumericMatrix P(P_nr);
   		Rcpp::NumericVector Pv_tmp(P_nr);
   		try{
   		  exp.get_P_result(exp_i);
   		}
   		catch (...){
   			Rf_warning("Unable to get perm mat results from graphm experiment object ");
   		}
   		//RcppGSL::Matrix P_tmp (A.nrow(),B.ncol())	;

   		RcppGSL::Matrix P_tmp ( exp.get_P_result(exp_i) )	;

   		if ( P_tmp != NULL ) {
   			int mat_size = P_tmp.nrow();
   			mat_size *= P_tmp.ncol();
   			if (P_tmp.nrow()==0 || P_tmp.ncol()==0 ) {
   				Rf_error("0 sized P_tmp: Permutation matrix");
   				return Rcpp::List();
   			}

   			try{
   				for (unsigned int j = 0 ; j < P_tmp.nrow(); j++){
   					for (unsigned int k = 0 ; k < P_tmp.ncol(); k++){
   						P(j,k) = P_tmp(j,k);
   						if (P(j,k)==1) {
   							Pv_tmp[j]=k+1;
   						}
   						//debug_out<<P_tmp(j,k)<<std::endl;
   					}
   				}
   				//Rf_PrintValue(P);
   			}
   			catch (...){
   				Rf_error("could not copy graphm result to Rcpp matrix");
   				return Rcpp::List();
  //
   			}
   		}	else {
   			Rf_error("graphm returned NULL as matrix pointer");
   			return Rcpp::List();
   		}
  //
  //
   		P_matrix_list.push_back(P);
   		P_vector_list.push_back(Pv_tmp);

   		P_tmp.free();

   		}
   	Rcpp::List  res = Rcpp::List::create(
   			Rcpp::Named("debugprint_file") = algorithm_params["debugprint_file"],
                                           Rcpp::Named("algo_names") = algo_list,
                                           Rcpp::Named("Pmat") = P_matrix_list,
                                           Rcpp::Named("exp_count") = exp_count,
                                           Rcpp::Named("Pvec") = P_vector_list
   		);
   	// if (debug_out.is_open())
   	// 	debug_out.close();
   	return res;

  } catch( std::exception &ex) {
  	// if (debug_out.is_open())
  	// 	debug_out.close();
  	Rf_error(ex.what());
  } catch (...) {
  	// if (debug_out.is_open())
  	// 	debug_out.close();
  	return Rcpp::List();
  }

}

