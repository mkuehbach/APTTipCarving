/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#ifndef __BICARVING_ORIMATH_H__
#define __BICARVING_ORIMATH_H__

#include "BiCarving_Math.h"

struct axisangle
{
	apt_real v1;
	apt_real v2;
	apt_real v3;
	apt_real theta;
	axisangle() : v1(1.f), v2(0.f), v3(0.f), theta(0.f) {}
	axisangle(const apt_real _v1, const apt_real _v2, const apt_real _v3, const apt_real _th) :
		v1(_v1), v2(_v2), v3(_v3), theta(_th) {}
	t3x3 ax2om();
};


#endif
