/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#include "BiCarving_Settings.h"


WHAT_TO_DO Settings::WhatToDo = E_NOTHING;
string Settings::GBTriSurface = "";
string Settings::AxisAngleSX = "";
string Settings::BungeEulerGrain1 = "";
string Settings::BungeEulerGrain2 = "";
apt_real Settings::CylHeight = 0.f;
apt_real Settings::CylRadius = 0.f;
apt_real Settings::LatticeConstant = 0.f;
apt_real Settings::SXFrustHeight = 0.f;
apt_real Settings::SXFrustRadiusB = 0.f;
apt_real Settings::SXFrustRadiusT = 0.f;
apt_real Settings::SXBlowup = 0.f;
apt_real Settings::SXBaseEps = 0.05; //nm
apt_real Settings::PRNGThinning = 1.f; //no thinning by default
size_t Settings::PRNGSeed = 1234;
size_t Settings::PRNGDiscard = 700000; //see Matsumoto and Saito as well as L'Ecuyer
unsigned int Settings::SimID = 0;
unsigned int Settings::TAPSimVacIdx = 1; //see C. Oberdorfer et al. HowTo TAPSim manual
unsigned int Settings::TAPSimBseIdx = 2;
unsigned int Settings::TAPSimTipIdx = 10;
bool Settings::IOBoundary = false;
bool Settings::IOCrystals = true;


inline apt_real str2real( const string str )
{
	return stod( str );
}

inline size_t str2long( const string str )
{
	return stoul( str );
}


void Settings::readXML(string filename) {
	//find the desired .xml file
	if ( 0 == filename.compare("") )
		filename = string("PARAPROBE.Input.Debug.xml");

	ifstream file( filename );
	if ( file.fail() ) 
		throw runtime_error(string("Unable to locate input file ") + filename);

	stringstream contents;
	contents << file.rdbuf();
	string xmlDocument(contents.str());

	xml_document<> tree;
	tree.parse<0>(&xmlDocument[0]);

	xml_node<>* rootNode = tree.first_node();
	if (0 != strcmp(rootNode->name(), "Parameters")) {
		throw runtime_error("Undefined parameters file!");
	}

	unsigned int mode = 0;
	if (0 != rootNode->first_node("WhatToDo"))
		mode = str2long( rootNode->first_node("WhatToDo")->value() );
	switch (mode)
	{
		case E_SINGLECRYSTAL:
			WhatToDo = E_SINGLECRYSTAL; break;
		case E_BICRYSTAL:
			WhatToDo = E_BICRYSTAL; break;
		default:
			WhatToDo = E_NOTHING;
	}

	if (0 != rootNode->first_node("GBTriangleSurface"))
		GBTriSurface = rootNode->first_node("GBTriangleSurface")->value();
	if (0 != rootNode->first_node("CylinderHeight"))
		CylHeight = str2real( rootNode->first_node("CylinderHeight")->value() );
	if (0 != rootNode->first_node("CylinderRadius"))
		CylRadius = str2real( rootNode->first_node("CylinderRadius")->value() );
	if (0 != rootNode->first_node("LatticeConstant"))
		LatticeConstant = str2real( rootNode->first_node("LatticeConstant")->value() );
	if (0 != rootNode->first_node("AxisAngleSX"))
		AxisAngleSX = rootNode->first_node("AxisAngleSX")->value();
	if (0 != rootNode->first_node("BungeEulerGrain1"))
		BungeEulerGrain1 = rootNode->first_node("BungeEulerGrain1")->value();
	if (0 != rootNode->first_node("BungeEulerGrain2"))
		BungeEulerGrain2 = rootNode->first_node("BungeEulerGrain2")->value();

	if (0 != rootNode->first_node("SXFrustHeight"))
		SXFrustHeight = str2real( rootNode->first_node("SXFrustHeight")->value() );
	if (0 != rootNode->first_node("SXFrustRadiusB"))
		SXFrustRadiusB = str2real( rootNode->first_node("SXFrustRadiusB")->value() );
	if (0 != rootNode->first_node("SXFrustRadiusT"))
		SXFrustRadiusT = str2real( rootNode->first_node("SXFrustRadiusT")->value() );
	if (0 != rootNode->first_node("SXBlowup"))
		SXBlowup = str2real( rootNode->first_node("SXBlowup")->value() );
	if (0 != rootNode->first_node("SXBaseEpsilon"))
		SXBaseEps = str2real( rootNode->first_node("SXBaseEpsilon")->value() );

	if (0 != rootNode->first_node("PRNGSeed"))
		PRNGSeed = str2long( rootNode->first_node("PRNGSeed")->value() );
	if (0 != rootNode->first_node("PRNGDiscard"))
		PRNGDiscard = str2long( rootNode->first_node("PRNGDiscard")->value() );
	if (0 != rootNode->first_node("PRNGThinning"))
		PRNGThinning = str2real( rootNode->first_node("PRNGThinning")->value() );
	if (0 != rootNode->first_node("TAPSIMVacIndex"))
		TAPSimVacIdx = str2long( rootNode->first_node("TAPSIMVacIndex")->value() );
	if (0 != rootNode->first_node("TAPSIMBaseIndex"))
		TAPSimBseIdx = str2long( rootNode->first_node("TAPSIMBaseIndex")->value() );
	if (0 != rootNode->first_node("TAPSIMTipIndex"))
		TAPSimTipIdx = str2long( rootNode->first_node("TAPSIMTipIndex")->value() );

	//##MK::IOBoundary
	//##MK::convert units to SI
}


bool Settings::checkUserInput()
{
	//##MK::check user input for validity and good sense
	switch(Settings::WhatToDo)
	{
		case E_SINGLECRYSTAL:
			cout << "Generating single crystalline pillar" << "\n";
			break;
		case E_BICRYSTAL:
			cout << "Generating bicrystalline pillar with defined boundary geometry" << "\n";
			break;
		default:
			cout << "Nothing to do, will quit now"; return false;
	}

	if ( Settings::CylHeight < EPSILON ) {
		cout << "Cylinder height is too small" << "\n"; return false;
	}
	if ( Settings::CylRadius < EPSILON ) {
		cout << "Cylinder radius is too small" << "\n"; return false;
	}
	if ( Settings::CylRadius > Settings::CylHeight ) {
		cout << "Cylinder should be not larger than height" << "\n"; return false;
	}
	if ( Settings::LatticeConstant < EPSILON ) {
		cout << "Lattice constant unrealistically too small" << "\n"; return false;
	}

	cout << "BicrystalCarving utilizes the following settings..." << "\n";
	cout << "Triangles in\t\t\t" << Settings::GBTriSurface << "\n";
	cout << "CylinderHeight\t\t\t" << Settings::CylHeight << " nm" << "\n";
	cout << "CylinderRadius\t\t\t" << Settings::CylRadius << " nm" << "\n";
	cout << "LatticeConstant\t\t\t" << Settings::LatticeConstant << " nm" << "\n";
	cout << "AxisAngleSX\t\t" << Settings::AxisAngleSX << "\n";
	cout << "BungeEulerGrain1\t\t" << Settings::BungeEulerGrain1 << "\n";
	cout << "BungeEulerGrain2\t\t" << Settings::BungeEulerGrain2 << "\n";

	cout << "SXFrustHeight\t\t" << Settings::SXFrustHeight << " nm" << "\n";
	cout << "SXFrustRadiusB\t\t" << Settings::SXFrustRadiusB << " nm" << "\n";
	cout << "SXFrustRadiusT\t\t" << Settings::SXFrustRadiusT << " nm" << "\n";
	cout << "SXBlowup\t\t\t" << Settings::SXBlowup << " nm" << "\n";

	cout << "PRNGSeed\t\t\t\t" << Settings::PRNGSeed << "\n";
	cout << "PRNGDiscard\t\t" << Settings::PRNGDiscard << "\n";
	cout << "PRNGThinning\t\t\t" << Settings::PRNGThinning << "\n";
	cout << "TAPSimVacuumIndex\t" << Settings::TAPSimVacIdx << "\n";
	cout << "TAPSimBaseIndex\t" << Settings::TAPSimBseIdx << "\n";
	cout << "TAPSimTipatomIndex\t" << Settings::TAPSimTipIdx << "\n";

	if ( Settings::IOBoundary == true )
		cout << "--->Visualizing boundary in VTK" << "\n";
	if ( Settings::IOCrystals == true )
		cout << "--->Visualizing crystals in VTK" << "\n";
	cout << endl;


	return true;
}


/*apt_real Settings::read_real( xml_node<>* const src, const string keyword, const apt_real defaultvalue )
{
	if (0 != src->first_node(keyword))
		return str2real( src->first_node(keyword)->value() );
	else
		return defaultvalue;
}*/

