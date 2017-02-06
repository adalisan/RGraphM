/***************************************************************************
 *   Copyright (C) 2008 by Mikhail Zaslavskiy   *
 *   mikhail.zaslavskiy@ensmp.fr   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
//#include "algorithm_path.h"
#include "sgm_algorithm_sfaq.h"

sgm_algorithm_sfaq::sgm_algorithm_sfaq(int num_seeds) :sgm_algorithm(num_seeds)
 {

}
match_result sgm_algorithm_sfaq::match(graph &g,graph &h,gsl_matrix* gm_P_i, gsl_matrix* gm_ldh,double dalpha_ldh){
 return(this->match_with_seeds(g, h ,gm_P_i, gm_ldh, dalpha_ldh, 0))	;
}

match_result sgm_algorithm_sfaq::match_with_seeds(graph& g, graph& h, gsl_matrix* gm_P_i, gsl_matrix* gm_ldh,double dalpha_ldh, unsigned int m_seeds)
{
	//int m  = get_param_i("num_seeds");
        bool bblast_match=(get_param_i("blast_match")==1);
        bool bblast_match_end=(get_param_i("blast_match_proj")==1);
	bool bbest_path_proj=(get_param_i("best_path_proj_sol")==1);
	bool bbest_path_blast_proj=(get_param_i("best_path_blast_proj_sol")==1);
	bool bbest_path_greedy=(get_param_i("best_path_greedy_sol")==1);
	bool bbest_path_blast_greedy=(get_param_i("best_path_blast_greedy_sol")==1);
	bool bbest_path=bbest_path_proj or bbest_path_blast_proj or bbest_path_greedy or bbest_path_blast_greedy;
	std::ofstream fverbose;
	double dfw_xeps=get_param_d("algo_fw_xeps");
	double dfw_feps=get_param_d("algo_fw_feps");
	double dlambda_M=get_param_d("qcvqcc_lambda_M");
	double dlambda_min=get_param_d("qcvqcc_lambda_min");
	bool bverbose=(get_param_i("verbose_mode")==1);
	double dhung_max=get_param_d("hungarian_max");
        double bgreedy=(get_param_i("hungarian_greedy")==1);

	if (bverbose)
		*gout<<"Path matching"<<std::endl;
	//some duplicate variables
	gsl_matrix* gm_Ag_d=g.get_descmatrix(cdesc_matrix);
	gsl_matrix* gm_Ah_d=h.get_descmatrix(cdesc_matrix);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ag_d,"Ag",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ah_d,"Ah",pdebug.strvalue);
	//laplacian construction
	//memory allocation
	bool bstop_algo=false;
	gsl_vector* gv_C=gsl_vector_alloc(N*N);
	gsl_vector* gv_temp=gsl_vector_alloc(N*N);
	gsl_vector* gv_temp2=gsl_vector_alloc(N*N);
	gsl_matrix_view gmv_temp2=gsl_matrix_view_array(gv_temp2->data,N,N);
	gsl_matrix* gm_temp2=&gmv_temp2.matrix;
	gsl_matrix_view gmv_temp=gsl_matrix_view_array(gv_temp->data,N,N);
	gsl_matrix* gm_temp=&gmv_temp.matrix;
	gsl_matrix* C;
	gsl_matrix_view gmv_C;
	gsl_vector_view gvv_P,gvv_P_prev,gvv_dP,gvv_P_lambda;
	gsl_matrix* gm_P=gsl_matrix_alloc(N,N);
	gsl_matrix* gm_P_prev=gsl_matrix_alloc(N,N);
	gsl_matrix* gm_P_lambda=gsl_matrix_alloc(N,N);
	gsl_matrix* gm_dP=gsl_matrix_alloc(N,N);
	gsl_matrix_set_zero(gm_P_prev);

	if (gm_P_i==NULL)
	  nonseededPtoseededP(gm_P,m);
		gsl_matrix_set_all(gm_P,1.0/(N-m);
	else
		gsl_matrix_memcpy(gm_P,gm_P_i);

        //perm matrix transformation into vector
	gvv_P=gsl_vector_view_array(gm_P->data,N*N);
	gvv_P_prev=gsl_vector_view_array(gm_P_prev->data,N*N);
	gvv_dP=gsl_vector_view_array(gm_dP->data,N*N);
	//and in opposite direction for gradient
	gmv_C=gsl_matrix_view_vector(gv_C,N,N);
	C=&gmv_C.matrix;
	gsl_vector*gv_debug_trace=gsl_vector_alloc(3);//debug trace information
	gsl_vector_set_zero(gv_debug_trace);
	//extern cycle over dlambda_cvcc
	//some temp variables

	//temp variables to trace the best solution along the algorithm path
	gsl_matrix* gm_P_bp_temp=NULL;
	gsl_matrix* gm_P_bp=NULL;
	double fbest_path=1e+300;

        //*************ALGORITHM***********************
	bool bpath_continue=true;
	double dlambda_add=dlambda_min;
	double dlambda_fix;
	double df_value_old,df_value,dP_norm,ddP_norm,dtemp;
	double dlambda=1;
	double dlambda_M_c=1;
	bool binc_lambda=false;
	bool bmax_lambda;
	gsl_matrix_memcpy(gm_P_prev,gm_P);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Lg_d,"Lg",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Lh_d,"Lh",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ag_d,"Ag",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ah_d,"Ah",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P",pdebug.strvalue);
	double dt0,dt1,dt2,dt3;
	while (bpath_continue){
	//LINEAR COMBINATION
	//new step initialisation

	//main optimization cycle
	int icounter=0;
	ddP_norm=1;
	df_value_old=f_FAQ(gm_Ag_d,gm_Ah_d,gm_P_prev,gm_temp,gm_temp2);
	while(!bstop_algo)
	{
	dt0=clock();
	//the default gsl representation is made by rows
	//gradient estimation: Agh*P, Here instead of P we have to use gvv_P_prev.
	if (pdebug.ivalue) gsl_matrix_printout(&gvv_P_prev.vector,"gvv_P_prev",pdebug.strvalue);
	faq_gradient(gm_Ag_d,gm_Ah_d, gvv_P_prev, gv_C);
	if (pdebug.ivalue) gsl_matrix_printout(gv_C,"gv_C",pdebug.strvalue);

	
	//debug trace
	if (pdebug.ivalue) gsl_vector_set(gv_debug_trace,0,gsl_vector_get(gv_debug_trace,0)+1);
	if (pdebug.ivalue) gsl_vector_set(gv_debug_trace,1,gsl_matrix_norm(C,1));
	if (pdebug.ivalue) gsl_vector_set(gv_debug_trace,2,graph_dist(g,h,gm_P,cscore_matrix));
	if (pdebug.ivalue) gsl_matrix_printout(gv_debug_trace,"debug_trace",pdebug.strvalue);

	//result save
	if (pdebug.ivalue) gsl_matrix_printout(C,"C=gradient",pdebug.strvalue);
        //update_C_hungarian(C);
	double dscale_factor =gsl_matrix_max_abs(C);
	dscale_factor=(dscale_factor>EPSILON)?dscale_factor:EPSILON;
	dscale_factor=dhung_max/dscale_factor;
	gsl_matrix_scale(C,dscale_factor);

	if (pdebug.ivalue) gsl_matrix_printout(C,"scale(C)",pdebug.strvalue);
	//hungarian, before  C matrix must be transposed
	
	dt1=clock();
	gsl_matrix_hungarian(C,gm_P,NULL,NULL,false,(bblast_match?gm_ldh:NULL),bgreedy);
	dt2=clock();
	
	gsl_matrix_scale(C,1/dscale_factor);

	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_P_prev,"gm_P_prev",pdebug.strvalue);
	//line search
	gsl_matrix_memcpy(gm_dP,gm_P);
	gsl_matrix_sub(gm_dP,gm_P_prev);

	if (pdebug.ivalue) gsl_matrix_printout(gm_dP,"gm_dP",pdebug.strvalue);

	double a,b1,b2,b3;
	//transpose all matrices but not Delta
	df_value=f_faq(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,gm_temp2);

	if (a>0)//projection is convex
	{
		double alpha = -b1/(2*a);

		if ((alpha<1) and (alpha>0))
			{
				gsl_matrix_scale(gm_dP,(1-alpha));
				gsl_matrix_sub(gm_P,gm_dP);
			};
		if (!(alpha>0))
			gsl_matrix_memcpy(gm_P,gm_P_prev);
	}
	//if (abs(a)<1e-50)
	//	if (!(b1>0))
	//		gsl_matrix_memcpy(gm_P,gm_P_prev);
	if ((a<0) or (abs(a)<1e-50)) //projection is concave
	{
		if (df_value>df_value_old)
			{gsl_matrix_memcpy(gm_P,gm_P_prev);
			 df_value=df_value_old;};
	};
	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P_step_finish",pdebug.strvalue);


	//stop criterion
	dP_norm=gsl_matrix_norm(gm_P_prev,1);
	df_value=f_qcvqcc(gm_Ag_d,gm_Ah_d,gm_Lg_d,gm_Lh_d,gm_Delta,gm_P,dlambda,gm_temp,gm_temp2);
	gsl_matrix_memcpy(gm_temp,gm_P_prev);
	gsl_matrix_sub(gm_P_prev,gm_P);
	if (df_value>df_value_old)
		int dbg=1;
	ddP_norm=gsl_matrix_norm(gm_P_prev,1);
	bstop_algo=((ddP_norm<dfw_xeps*N*dlambda_M_c) and ((abs(df_value-df_value_old)<dfw_feps*abs(df_value_old)*dlambda_M_c) or (ddP_norm==0)));
	bstop_algo = (bstop_algo or (icounter>pow(N,4)));
	bstop_algo=(bstop_algo or (abs(df_value-df_value_old)<1e-30));
	bstop_algo=(bstop_algo or ((icounter>0) and binc_lambda));//if we are on the increment step do not repreat many times
	icounter++;
	dt3=clock();
	long lnum_constraints=0;
	for (long li=0;li<N*N;li++)
		 if (gm_P->data[li]<1e-30)
			{
			lnum_constraints++;
			};
	if (bverbose)
		*gout<<"iter="<<icounter<<", x="<<dP_norm<<", dx="<<ddP_norm<<", f="<<df_value<<", df="<<df_value_old-df_value<<", "<<"grad="<<(df_value_old-df_value)/ddP_norm<<". #Act.Constr="<<lnum_constraints<<". Timing="<<(dt1-dt0)/CLOCKS_PER_SEC<<" "<<(dt2-dt1)/CLOCKS_PER_SEC<<" "<<(dt3-dt2)/CLOCKS_PER_SEC<<std::endl;
	df_value_old=df_value;

	//now we test different projection to estimate the best permutation
	if (bbest_path_proj)
	{
		//permuation projection
		gsl_matrix_transpose_memcpy(gm_P_bp_temp,gm_P);
		gsl_matrix_scale(gm_P_bp_temp,-10000);
		gsl_matrix_hungarian(gm_P_bp_temp,gm_P_prev,NULL,NULL,false,NULL,bgreedy);
		double df_bp_new=f_qcv(gm_Ag_d,gm_Ah_d,gm_P_prev,gm_temp2,true);
		if (df_bp_new<fbest_path)
		{
			fbest_path=df_bp_new;
			gsl_matrix_memcpy(gm_P_bp,gm_P_prev);
		};
	};

	if (bbest_path_blast_proj)
	{
		//permuation projection
		gsl_matrix_transpose_memcpy(gm_P_bp_temp,gm_P);
		gsl_matrix_scale(gm_P_bp_temp,-10000);
		gsl_matrix_hungarian(gm_P_bp_temp,gm_P_prev,NULL,NULL,false,gm_ldh,bgreedy);
		double df_bp_new=f_qcv(gm_Ag_d,gm_Ah_d,gm_P_prev,gm_temp2,true);
		if (df_bp_new<fbest_path)
		{
			fbest_path=df_bp_new;
			gsl_matrix_memcpy(gm_P_bp,gm_P_prev);
		};
	};

	if (bbest_path_greedy)
	{
	};

	if (bbest_path_blast_greedy)
	{
	};
	gsl_matrix_memcpy(gm_P_prev,gm_P);
	};//end Frank-Wolfe cycle
	bpath_continue=(dlambda>0);//stop if it is a concave function
	//we decrement lambda and repeat the Frank-Wolfe step
	//we are on the lambda-increment step
	if (!(dlambda>0) and (icounter>1) and (dlambda_add>2*dlambda_min))//if dlambda_add is very small then we have to stop
		{
			binc_lambda=true;
			bpath_continue=true;
			dlambda_M_c=dlambda_M;//for the next times
		};
	if (binc_lambda)
	{
		if ((bmax_lambda) and (icounter>1)) //there was two and more steps
		{
			dlambda_add/=2;//return to the previos step
			bmax_lambda=false;
		}
		else{
			if (icounter>1)	dlambda_add/=2;
			if (bmax_lambda) dlambda_add*=2;//continue lambda step
			//we change the current point if it is a limit step or if step is too small

		    };
		if (((!bmax_lambda) and (icounter==1)) or (dlambda_add<dlambda_min)){ binc_lambda=false;dlambda_M_c=1;//the first time we use the original stop criterion
			}; //just  change the lambda_fix
 		dlambda=dlambda_fix-dlambda_add;
		dlambda=(dlambda<0)?0:dlambda;
		if (dlambda==0){ //concave function
			binc_lambda=false;dlambda_M_c=1;
				};
		gsl_matrix_memcpy(gm_P_prev,gm_P_lambda);
	}
	else
	{
		bmax_lambda=true;
		dlambda_fix=dlambda;
		gsl_matrix_memcpy(gm_P_lambda,gm_temp);
		binc_lambda=true;
		dlambda_M_c=dlambda_M;//for the next times
		bool beigen_trace=false;double dmin_eval=0;
		if (beigen_trace)
		{
		gsl_matrix* gm_Hessian=gsl_matrix_alloc(N*N,N*N);
		gsl_matrix* gm_Hessian_2=gsl_matrix_alloc(N*N,N*N);
		gsl_matrix* gm_Hessian_3=gsl_matrix_alloc(N*N,N*N);

		for (int i1=0;i1<N;i1++)
			for (int i2=0;i2<N;i2++)
			   for (int j1=0;j1<N;j1++)
				for (int j2=0;j2<N;j2++)
			  {
				double dvalue=0;
				if (i2==j2)
					dvalue+=gm_Ag_d->data[i1*N+j1];
				if (i1==j1)
					dvalue-=gm_Ah_d->data[i2*N+j2];
				gsl_matrix_set(gm_Hessian,i1+i2*N,j1+j2*N,dvalue);
			  };
		gsl_matrix_transpose_memcpy(gm_Hessian_2,gm_Hessian);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Hessian_2,gm_Hessian,0,gm_Hessian_3);
		gsl_matrix_transpose_memcpy(gm_Hessian,gm_Hessian_3);
		gsl_matrix_scale(gm_Hessian,dlambda);
		for (int i1=0;i1<N;i1++)
			for (int i2=0;i2<N;i2++)
			   for (int j1=0;j1<N;j1++)
				for (int j2=0;j2<N;j2++)
			  {
				double dvalue=0;
				dvalue+=-(1-dlambda)*gm_Lg_d->data[i1*N+j1]*gm_Lh_d->data[i2*N+j2];
				gsl_matrix_set(gm_Hessian,i1+i2*N,j1+j2*N,gsl_matrix_get(gm_Hessian,i1+i2*N,j1+j2*N)+dvalue);
			   };
		gsl_eigen_symm_workspace * gesw= gsl_eigen_symm_alloc (N*N);
		gsl_vector* eval_h=gsl_vector_alloc(N*N);
		gsl_eigen_symm (gm_Hessian, eval_h,gesw);
		gsl_matrix_free(gm_Hessian);
		gsl_matrix_free(gm_Hessian_2);
		gsl_matrix_free(gm_Hessian_3);
		dmin_eval=gsl_vector_min(eval_h);
		gsl_eigen_symm_free (gesw);
		gsl_vector_free(eval_h);
		};
		if (bverbose) *gout<<"lambda_fix="<<dlambda_fix<<", dl_add="<<dlambda_add<<", min eval="<<dmin_eval<<std::endl;
	};
	if (dlambda_add<dlambda_min)
		dlambda_add=dlambda_min;
	//std::cout<<"lambda="<<dlambda<<", dl_add="<<dlambda_add<<std::endl;
	bstop_algo=false;
	};//end pathway cycle
	match_result mres;
	mres.gm_P_exact=gsl_matrix_alloc(N,N);
	gsl_matrix_memcpy(mres.gm_P_exact,gm_P);

	//permuation projection
	gsl_matrix_transpose_memcpy(gm_P_prev,gm_P);
	gsl_matrix_scale(gm_P_prev,-10000);
	gsl_matrix_hungarian(gm_P_prev,gm_P,NULL,NULL,false,(bblast_match_end?gm_ldh:NULL),bgreedy);
	if (bbest_path){
		double df_bp_new=f_qcv(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,true);
		if (df_bp_new>fbest_path)
		{
			*gout<<"Best path solution is used"<<std::endl;
			gsl_matrix_memcpy(gm_P,gm_P_bp);
		};
	};

	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P_projected",pdebug.strvalue);
	//memory free
	gsl_matrix_free(gm_P_lambda);
	gsl_matrix_free(gm_P_prev);
	gsl_matrix_free(gm_dP);
	gsl_matrix_free(gm_Lg_d);
	gsl_matrix_free(gm_Lh_d);
	gsl_matrix_free(gm_Ag_d);
	gsl_matrix_free(gm_Ah_d);
	gsl_matrix_free(gm_Delta_tmp);
	gsl_matrix_free(gm_Delta);
	gsl_vector_free(gv_C);
	gsl_vector_free(gv_temp);
	gsl_vector_free(gv_temp2);
	gsl_vector_free(gv_debug_trace);
	if (bbest_path)
	{
		gsl_matrix_free(gm_P_bp_temp);
		gsl_matrix_free(gm_P_bp);
	};

	mres.gm_P=gm_P;

	//initial score
	mres.vd_trace.push_back(graph_dist(g,h,cscore_matrix));
	//final score
	mres.vd_trace.push_back(graph_dist(g,h,gm_P,cscore_matrix));
	//other output parameters
	mres.dres=mres.vd_trace.at(1);
	//transpose matrix save
	mres.gm_P=gm_P;
	return mres;
}

double f_faq(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,gm_temp2);
void nonseededPtoseededP(gsl_matrix * gmv_P, unsigned int m);
//concave function value
double grad_faq(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,gsl_matrix * gm_temp,gsl_matrix *gm_temp2, unsigned int m = 0);

//gradient for qcvqcc function
void grad_faq(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,gsl_matrix * gm_temp,gsl_matrix *gm_temp2, unsigned int m = 0);

void sgm_algorithm_sfaq::void grad_faq(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_P,gsl_matrix * gm_temp,gsl_matrix *gm_temp2);
f_faq(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,gm_temp2);
	gsl_matrix_view gmv_temp=gsl_matrix_view_array(gv_temp->data,N,N);
	qcv_gradient_opt(gm_Ag_d,gm_Ah_d,gv_P,gv_grad,&gmv_temp.matrix,m);
	gsl_vector_scale(gv_grad,dlambda);
	qcc_gradient_opt(gm_Lg_d, gm_Lh_d, gv_P, gv_grad,-(1-dlambda),&gmv_temp.matrix,m);
}

void sgm_algorithm_sfaq::nonseededPtoseededP(gsl_matrix *P, unsigned int m){
	//*gout<<"Using the first "<<m<<" vertex pairs as seeds"<<std::endl;
	if (m>0) {

	gsl_matrix_view gmv_P_seed = gsl_matrix_submatrix(P,0,0,m,m);
	gsl_matrix_set_identity(&gmv_P_seed.matrix);
	gsl_matrix_view gmv_P_seed_nonseed = gsl_matrix_submatrix(P,0,m,m,N-m);
	gsl_matrix_set_zero(&gmv_P_seed_nonseed.matrix);
	gsl_matrix_view gmv_P_nonseed_seed = gsl_matrix_submatrix(P,m,0,N-m,m);
	gsl_matrix_set_zero(&gmv_P_nonseed_seed.matrix);
	}

}

//optimized version of gradient calculations: tensor product tricks
// Note: returns  the transpose of gradient
void sgm_algorithm_sfaq::qcv_gradient_opt(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_vector* gv_P,  gsl_vector* gv_grad, gsl_matrix * gm_temp, unsigned int m)
{
	
}



//convex-concave function value
double sgm_algorithm_sfaq::f_faq(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,double dlambda,gsl_matrix * gm_temp,gsl_matrix *gm_temp2)
{
	return dres;
}
