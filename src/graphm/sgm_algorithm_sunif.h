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
#ifndef SGM_ALGORITHM_sUNIF_H
#define SGM_ALGORITHM_sUNIF_H
#include "sgm_algorithm.h"


/**
A path following graph matching algorithm

	@author Sancar Adali <sadali@gmail.com>;Mikhail Zaslavskiy <mikhail.zaslavskiy@ensmp.fr>
*/
class sgm_algorithm_sunif : public sgm_algorithm
{

public:
	  sgm_algorithm_sunif(int);
    match_result match(graph &g,graph &h,gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1);

		match_result match_with_seeds(graph& g, graph& h, gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1, unsigned int m_seeds=0);

    

};

#endif
