/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#ifndef __BICARVING_MATH_H__
#define __BICARVING_MATH_H__

#include "BiCarving_Datatypes.h"


struct t3x3
{
	apt_real a11;				//a second order rank tensor with row-column indexing
	apt_real a12;
	apt_real a13;
	apt_real a21;
	apt_real a22;
	apt_real a23;
	apt_real a31;
	apt_real a32;
	apt_real a33;
	t3x3() :	a11(1.0), a12(0.0), a13(0.0),
				a21(0.0), a22(1.0), a23(0.0),
				a31(0.0), a32(0.0), a33(1.0) {}	//initialize to identity tensor
	t3x3(	const apt_real _a11, const apt_real _a12, const apt_real _a13,
			const apt_real _a21, const apt_real _a22, const apt_real _a23,
			const apt_real _a31, const apt_real _a32, const apt_real _a33 ) :
				a11(_a11), a12(_a12), a13(_a13),
				a21(_a21), a22(_a22), a23(_a23),
				a31(_a31), a32(_a32), a33(_a33) {}
	t3x3 premultiplyR1( t3x3 const & R1 );
	t3x3 inverse();
};

std::ostream& operator << (std::ostream& in, t3x3 const & val);


class mathHdl
{
public:
	mathHdl(){}
	~mathHdl(){}

	apt_real closestPointOnTriangle( const tri3d face, const p3d src );
};



#endif
