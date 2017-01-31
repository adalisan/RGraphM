#include "sgm_algorithm.h"
#include "algorithm.h"
#include "rpc.h"

sgm_algorithm::sgm_algorithm(std::string fconfig)
: algorithm(fconfig)
{
	gm_ldh=NULL;dalpha_ldh=0;bnosymm=false; this->m =0; 
}
sgm_algorithm::sgm_algorithm()
: algorithm()
{
	gm_ldh=NULL;dalpha_ldh=0;bnosymm=false;df_norm=0;N=0;cdesc_matrix='A';cscore_matrix='A'; this->m=0;
}
sgm_algorithm::sgm_algorithm(int m)
	: algorithm()
{
	gm_ldh=NULL;dalpha_ldh=0;bnosymm=false;df_norm=0;N=0;cdesc_matrix='A';cscore_matrix='A';this->m=m; 
}


void sgm_algorithm::nonseededPtoseededP(gsl_matrix *P, unsigned int _m){
	//*gout<<"Using the first "<<m<<" vertex pairs as seeds"<<std::endl;
	if (_m>0) {

	gsl_matrix_view gmv_P_seed = gsl_matrix_submatrix(P,0,0,_m,_m);
	gsl_matrix_set_identity(&gmv_P_seed.matrix);
	gsl_matrix_view gmv_P_seed_nonseed = gsl_matrix_submatrix(P,0,_m,_m,N-_m);
	gsl_matrix_set_zero(&gmv_P_seed_nonseed.matrix);
	gsl_matrix_view gmv_P_nonseed_seed = gsl_matrix_submatrix(P,_m,0,N-_m,_m);
	gsl_matrix_set_zero(&gmv_P_nonseed_seed.matrix);
	}

}

//common framework for graph matching algorithm
match_result sgm_algorithm::gmatch(graph& g, graph& h,gsl_matrix *gm_P_i,gsl_matrix* _gm_ldh,double _dalpha_ldh)
{
	dalpha_ldh=_dalpha_ldh;
	if (_gm_ldh!=NULL) set_ldhmatrix(_gm_ldh);
	pdebug = get_param("debugprint");
	pdebug_f = get_param("debugprint_file");
	pdebug.strvalue=pdebug_f.strvalue;
	bverbose=(get_param_i("verbose_mode")==1);
	sverbfile=get_param_s("verbose_file");
	cdesc_matrix=get_param_c("cdesc_matrix");
	cscore_matrix=get_param_c("cscore_matrix");
	bnosymm=(get_param_i("nosymm_matrix")==1);
	//if (sverbfile.compare("cout")==0)
	//	gout=&std::cout;
	//else {
	fverbose.open(sverbfile.c_str(),std::ios_base::app);
	this->gout=&fverbose;
		
	//};
	N=g.getN();
	unsigned  int m_seeds = 0;
	int do_seeds = get_param_i("do_sgm");
	 

	gsl_matrix* gm_Ag_d=g.get_descmatrix(cscore_matrix);
	gsl_matrix* gm_Ah_d=h.get_descmatrix(cscore_matrix);
	gsl_matrix* gm_temp=gsl_matrix_alloc(N,N);
	df_norm=pow(gsl_matrix_norm(gm_Ag_d,2),2)+pow(gsl_matrix_norm(gm_Ah_d,2),2);
	df_norm=(df_norm>EPSILON)?df_norm:EPSILON;

	time_t t1=time(NULL);
	match_result mres;
	if (do_seeds==1) {
		
		mres = match_with_seeds(g,h,gm_P_i,_gm_ldh,_dalpha_ldh,m);
	    
	} else {

		mres = match(g,h,gm_P_i,_gm_ldh,_dalpha_ldh);
	}

	mres.dfvalue=f_qcv(gm_Ag_d,gm_Ah_d,mres.gm_P,gm_temp,true);
	if (mres.gm_P_exact!=NULL)
		mres.dfvalue_exact=f_qcv(gm_Ag_d,gm_Ah_d,mres.gm_P_exact,gm_temp,true);
	else
		mres.dfvalue_exact=mres.dfvalue;
	gsl_matrix_free(gm_Ag_d);
	gsl_matrix_free(gm_Ah_d);
	gsl_matrix_free(gm_temp);

	mres.dtime=difftime(time(NULL),t1);
	if (fverbose.is_open()){
		fverbose.close();
	}
	return mres;
}
