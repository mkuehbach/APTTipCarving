/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#ifndef __BICARVING_MATH_H__
#define __BICARVING_MATH_H__

#include "BiCarving_Datatypes.h"




class mathHdl
{
public:
	mathHdl(){}
	~mathHdl(){}

	apt_real closestPointOnTriangle( const tri3d face, const p3d src );
};



#endif
