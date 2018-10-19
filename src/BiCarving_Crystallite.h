/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#ifndef __BICARVING_CRYSTALLITE_H__
#define __BICARVING_CRYSTALLITE_H__


#include "BiCarving_Orimath.h"

class unitcellaggr
{
public:
	unitcellaggr();
	unitcellaggr(const apt_real _a, aabb3d unitbox, const unsigned int model );
	~unitcellaggr();

	p3d get_atom(const size_t b, const int u, const int v, const int w);
	string report_unitcell();

	apt_real a;
	int umin;
	int umax;
	int vmin;
	int vmax;
	int wmin;
	int wmax;

	//orthogonal base vectors
	v3d a1;
	v3d a2;
	v3d a3;

	//base atoms
	vector<p3d> base;
};



class crystal
{
public:
	crystal();
	~crystal();

	//void build_lattice();

	//struct crystalstructure;
	t3x3 ori;
	//vector<p3dm1> atoms;
};


class boundary
{
public:
	boundary();
	~boundary();

	void read_triangulation();

	//inline bool in_front_of_me( p3d const & p ) const;

	vector<tri3d> triangulation;
	vector<plane3d> simplex;
};


#endif
