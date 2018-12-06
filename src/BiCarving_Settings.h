/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_SETTINGS_H__
#define __BICARVING_SETTINGS_H__

#include "BiCarving_Profiling.h"


enum WHAT_TO_DO {
	E_NOTHING,
	E_SINGLECRYSTAL,
	E_BICRYSTAL
};


enum GB_INPUT_HOW {
	E_GB_NOTHING,
	E_GB_EXPLICIT_TRIANGLEMESH,
	E_GB_PLANENORMAL_AND_POINT
};

enum ORI_INPUT_HOW {
	E_ORI_NOTHING,
	E_ORI_AXISANGLE,
	E_ORI_BUNGEEULER,
	E_ORI_MATRIX
};

//add thirdparty XML reader header library functionality by M. Kalicinski
#include "thirdparty/RapidXML/rapidxml.hpp"

using namespace rapidxml;


class Settings {

public:
	static WHAT_TO_DO WhatToDo;
	static GB_INPUT_HOW GBModel;
	static ORI_INPUT_HOW OriParameterization;

	static string GBTriSurface;
	static string GBPlaneNormal;
	static string GBPointOnPlane;
	static string AxisAngleSX;
	static string BungeEulerGrain1;
	static string BungeEulerGrain2;
	static string Orimatrix1;
	static string Orimatrix2;

	static apt_real GBBoundaryThickness;
	static apt_real CylHeight;
	static apt_real CylRadius;
	static apt_real LatticeConstant;

	static apt_real SXFrustHeight;
	static apt_real SXFrustRadiusB;
	static apt_real SXFrustRadiusT;
	static apt_real SXBlowup;
	static apt_real SXBaseEps;

	static apt_real PRNGThinning;
	static size_t PRNGSeed;
	static size_t PRNGDiscard;

	static unsigned int SimID;
	static unsigned int TAPSimVacIdx;
	static unsigned int TAPSimBseIdx;
	static unsigned int TAPSimTipIdx;


	static bool IOBoundary;
	static bool IOCrystals;

//prototypes
	static void readXML(string filename = "");
	static bool checkUserInput( void );
};

#endif
