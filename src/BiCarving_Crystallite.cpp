/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#include "BiCarving_Crystallite.h"

unitcellaggr::unitcellaggr()
{
	this->a = 0.f;

	this->umin = 0;
	this->umax = 0;
	this->vmin = 0;
	this->vmax = 0;
	this->wmin = 0;
	this->wmax = 0;

	this->a1 = v3d();
	this->a2 = v3d();
	this->a3 = v3d();

	this->base.clear();
}


unitcellaggr::unitcellaggr(const apt_real _a, const aabb3d unitbox, const unsigned int model )
{
	//initialize fcc base atoms
	a = _a;

	//initialize cubic base vectors
	a1 = v3d( _a*1.f, 0.f, 0.f);
	a2 = v3d( 0.f, _a*1.f, 0.f);
	a3 = v3d( 0.f, 0.f, _a*1.f);

	base.push_back( p3d(0.f, 0.f, 0.f) ); //Al 8x 1/8 = 1
	base.push_back( p3d(0.5, 0.5, 0.f) ); //Al
	base.push_back( p3d(0.f, 0.5, 0.5) ); //Al
	base.push_back( p3d(0.5, 0.f, 0.5) ); //Al 6x 1/2 = 3 Al --> 4 units per EZ okay

	//unitbox gives min/max dimensions in nanometer that we have to fill construct on positive sectors of \mathcal{R}^3
	umin = static_cast<int>(floor(unitbox.xmi / _a));
	umax = static_cast<int>(ceil(unitbox.xmx / _a));
	vmin = static_cast<int>(floor(unitbox.ymi / _a));
	vmax = static_cast<int>(ceil(unitbox.ymx / _a));
	wmin = static_cast<int>(floor(unitbox.zmi / _a));
	wmax = static_cast<int>(ceil(unitbox.zmx / _a));

	//unitbox is axis-aligned to standard orientation 0.0, 0.0, 0.0 Bunge Euler fcc crystal lattice
}

unitcellaggr::~unitcellaggr(){}


p3d unitcellaggr::get_atom(const size_t b, const int u, const int v, const int w)
{
	apt_real uu = static_cast<apt_real>(u);
	apt_real vv = static_cast<apt_real>(v);
	apt_real ww = static_cast<apt_real>(w);

	//##MK::implicit origin at 0,0,0
	p3d res = p3d(
			(base[b].x + uu)*a1.u + (base[b].y + vv)*a2.u + (base[b].z + ww)*a3.u,
			(base[b].x + uu)*a1.v + (base[b].y + vv)*a2.v + (base[b].z + ww)*a3.v,
			(base[b].x + uu)*a1.w + (base[b].y + vv)*a2.w + (base[b].z + ww)*a3.w  );

	return res;
}


string unitcellaggr::report_unitcell()
{
	string str = "imi/imx = " + to_string(umin) + ";" + to_string(umax) + ";" + to_string(vmin) + ";" + to_string(vmax) + ";" + to_string(wmin) + ";" + to_string(wmax) + "\n";
	return str;
}


crystal::crystal()
{
	ori = t3x3();
}

crystal::~crystal()
{
}



boundary::boundary()
{
}


boundary::~boundary()
{
}


void boundary::read_triangulation()
{
	triangulation.clear();

	//##MK::open Settings::GBTriSurface leave triangulation empty on error

	simplex.clear();

	//##MK::build simplex
}


/*
inline bool boundary::in_front_of_me( p3d const & p ) const
{
	//##MK::test against if in front of every triangle in triangulation assuming a consistent normal
	for(auto it = simplex.begin(); it != simplex.end(); ++it) {
		if ( it->in_front_of_me( p ) == true )
			continue;
		else
			return false;
	}
	return true;
}
*/

