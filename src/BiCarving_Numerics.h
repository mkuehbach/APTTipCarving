/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_NUMERICS_H__
#define __BICARVING_NUMERICS_H__

#include "BiCarving_Information.h"


//precision
//MK::we utiĺize double precision by default
#define EMPLOY_DOUBLE_PRECISION

#define EPSILON							(1.0e-12)
typedef double apt_real;
typedef size_t apt_int;

//type range
#define UINT64MX						(numeric_limits<size_t>::max())
#define UINT64MI						(numeric_limits<size_t>::lowest())
#define UINT32MX						(numeric_limits<unsigned int>::max())
#define UINT32MI						(numeric_limits<unsigned int>::lowest())
#define F32MX							(numeric_limits<apt_real>::max())
#define F32MI							(numeric_limits<apt_real>::lowest())
#define SIZETMX							(numeric_limits<size_t>::max())

//handling of iontypes and flagging them for inclusion/exclusion in analyses
#define UNKNOWNTYPE						0

//MK::coordinate system is right-handed x,y,z
#define PARAPROBE_XAXIS					0
#define PARAPROBE_YAXIS					1
#define PARAPROBE_ZAXIS					2

//MK::Mersenne initialization
#define MT19937SEED						(-1)
#define MT19937WARMUP					(700000)

#endif
