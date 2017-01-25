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
#ifndef SGMALGORITHM_H
#define SGMALGORITHM_H
#define EPSILON 1e-100
#include "rpc.h"
#include "algorithm"
#include "graph.h"
#include <math.h>
#include "hungarian.h"
#include "algorithm.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include <vector>
#include <iostream>
#include <ctime>

/*


Parent class for all graph matching algorithms

@author Mikhail Zaslavskiy <mikhail.zaslavskiy@ensmp.fr>
*/

class sgm_algorithm : public algorithm
{
public:
	sgm_algorithm(std::string );
	sgm_algorithm();
	sgm_algorithm(int );
	match_result gmatch(graph& g, graph& h,gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1);//common stuff,
	virtual match_result match(graph& g, graph& h, gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1)=0;//particular method implementation
	virtual match_result match_with_seeds(graph& g, graph& h, gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1, unsigned int m_seeds=0)=0;//particular method implementation
		~sgm_algorithm();


protected:

	//double f_qcv(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_P,gsl_matrix * gm_temp,bool bqcv=false);

	long long m;


};

#endif //SGM_ALGORITHM
