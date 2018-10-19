/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#include "BiCarving_BiXXHdl.h"

bicrystal::bicrystal()
{
	//calls type default constructs by default behavior
	owner = NULL;
}


bicrystal::bicrystal( solverHdl * own )
{
	owner = own;
}


bicrystal::~bicrystal()
{
	//default constructo destruct solverHdl class object internals
	//do not delete owner, only backreference
}


void bicrystal::define_domain()
{
	//cylinder centroid at coordinate system origin 0,0,0
	world.H = Settings::CylHeight;
	world.R = Settings::CylRadius;
	world.mybox.xmi = -1.f*Settings::CylRadius;
	world.mybox.xmx = +1.f*Settings::CylRadius;
	world.mybox.ymi = -1.f*Settings::CylRadius;
	world.mybox.ymx = +1.f*Settings::CylRadius;
	world.mybox.zmi = -0.5*Settings::CylHeight;
	world.mybox.zmx = +0.5*Settings::CylHeight;
	world.mybox.scale();
	world.center = p3d( 0.f, 0.f, 0.f);
}


void bicrystal::load_grainboundary()
{
	gAB.read_triangulation();
}


void bicrystal::generate_bicrystal()
{
	//##MK::build lattice of crystal gA with atoms in front of the boundary gAB

	/*void crystal::build_lattice( cylinder const & box, boundary const & constraint )
	{
		//##MK::define crystal structure
		//sample over positions in cylinder
		//
	}*/

	//build lattice of crystal gB with atoms behind the boundary gAB
}


void bicrystal::carve_tip()
{
}


void bicrystal::process()
{
}
