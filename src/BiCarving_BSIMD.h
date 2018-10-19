/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_BSIMD_H__
#define __BICARVING_BSIMD_H__

#include "BiCarving_Parallelization.h"

//MK::connecting to the NumScale boostSIMD library to access portable vector intrinsics for recon
//#define USE_BOOST

#ifdef USE_BOOST

	#include <boost/simd/pack.hpp>
	namespace bs = boost::simd;

	//include Boost functionality
	#include <boost/dynamic_bitset.hpp>

	//include BoostSIMD
	//##MK

#endif

#endif
