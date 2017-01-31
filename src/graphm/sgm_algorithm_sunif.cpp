#include "sgm_algorithm_sunif.h"

sgm_algorithm_sunif::sgm_algorithm_sunif(int num_seeds) :sgm_algorithm(num_seeds)
 {

}
match_result sgm_algorithm_sunif::match(graph &g,graph &h,gsl_matrix* gm_P_i, gsl_matrix* gm_ldh,double dalpha_ldh){
 return(this->match_with_seeds(g, h ,gm_P_i, gm_ldh, dalpha_ldh, 0))	;
}

match_result sgm_algorithm_sunif::match_with_seeds(graph& g, graph& h,gsl_matrix* gm_P_i, gsl_matrix* gm_ldh,double dalpha_ldh, unsigned int m_seeds)
{
	//some duplicate variables
	match_result mres;
	mres.gm_P=gsl_matrix_alloc(N,N);
    gsl_matrix_set_all(mres.gm_P,1.0/(N- m_seeds));
    nonseededPtoseededP(mres.gm_P, m_seeds);

	
	mres.gm_P_exact=NULL;
        //initial score
	mres.vd_trace.push_back(graph_dist(g,h,cscore_matrix));
	//final score
	mres.vd_trace.push_back(graph_dist(g,h,mres.gm_P,cscore_matrix));
	mres.dres=mres.vd_trace.at(1);
	mres.inum_iteration=2;
	return mres;
}