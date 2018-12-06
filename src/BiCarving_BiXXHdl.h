/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_BIXXHDL_H__
#define __BICARVING_BIXXHDL_H__

#include "BiCarving_SiXHdl.h"

/*
class bicrystal
{
	//top-level construct implementing the worker instance at process rank within the MPI process level parallelism

public:
	bicrystal();
	bicrystal( solverHdl * own );
	~bicrystal();


	void define_domain();
	void load_grainboundary();
	void generate_bicrystal();
	void carve_tip();

	void process();

	solverHdl* owner;

	cylinder world;					//a container of atoms
	crystal gA;						//grain A of atoms
	crystal gB;						//grain B of atoms
	boundary gAB;					//the grain boundary
};
*/


class polycrystal
{
	//top-level construct implementing the worker instance at process rank within the MPI process level parallelism

public:
	polycrystal();
	polycrystal( solverHdl * own );
	~polycrystal();

	unsigned int process_on_which_side( p3d const & in );
	void process_px( t3x3 & o1, t3x3 & o2 );

	solverHdl* owner;

	t3x3 gA;
	t3x3 gB;
	boundary gb;

	//TAPSim employs a cascading mesh scheme to bridge from the nm to the m scale
	tip zoneII;						//cylindrical container with reduced number of meshpoints
	tip zoneI;						//same shape a generate tip full density of meshpoint to get sufficient accuracy of launch angle and early trajectory
	tip pillar;						//zone0 aka the synthetic pillar

	vector<p3dm1> pp3f;				//all mesh points of the tip, zoneI, and zoneII to pass to MeshGen
};

#endif
