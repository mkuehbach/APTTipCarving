/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_ABINITIO_H__
#define __BICARVING_ABINITIO_H__

#include "BiCarving_Numerics.h"

//natural constants
#define PI								(3.141592653589793238462643383279502884197169399375105820974)
#define kboltzmann						(1.3806488e-23)		//J/Kmol
#define echarge							(1.602176565e-19)	//Coulomb
#define Navogadro						(6.022140857e+23)	//1/mol
#define RGAS							(8.31446154) 		//(Navogadro)*(kboltzmann)) //J/K/mol

//natural beauty
#define SQR(a)							((a)*(a))
#define CUBE(a)							((a)*(a)*(a))
#define MIN(X,Y)						(((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y)						(((X) < (Y)) ? (Y) : (X))
#define CLAMP(x,lo,hi)					(MIN(MAX((x), (lo)), (hi)))


//unit conversions
#define CELSIUS2KELVIN(T)				((273.15) + (T))
#define KELVIN2CELSIUS(T)				((T) - (273.15))

#define DEGREE2RADIANT(theta)			((theta)/(180.0)*(PI))
#define RADIANT2DEGREE(rad)				((rad)/(PI)*(180.0))

//scaling conversions
#define NANOMETER2ANGSTROEM(nm)			((10.0)*(nm))
#define ANGSTROEM2NANOMETER(ang)		((ang)/(10.0))
#define METER2NANOMETER(m)				((1.0e9)*(m))
#define NANOMETER2METER(nm)				((1.0e-9)*(nm))


//crystal lattice identifiers
#define FCC								1


#endif
