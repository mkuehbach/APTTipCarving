/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_SIXHDL_H__
#define __BICARVING_SIXHDL_H__

#include "BiCarving_VTKIO.h"

class solverHdl;

class singlecrystal
{
	//top-level construct implementing the worker instance at process rank within the MPI process level parallelism

public:
	singlecrystal();
	singlecrystal( solverHdl * own );
	~singlecrystal();

	void process_sx( axisangle const & targetori_ax, t3x3 const & targetori_om );

/*
	void define_domain();
	void load_grainboundary();
	void generate_bicrystal();
	void carve_tip();
*/
	solverHdl* owner;

	crystal gA;						//grain A of atoms

	//TAPSim employs a cascading mesh scheme to bridge from the nm to the m scale
	tip zoneII;						//cylindrical container with reduced number of meshpoints
	tip zoneI;						//same shape a generate tip full density of meshpoint to get sufficient accuracy of launch angle and early trajectory
	tip pillar;						//zone0 aka the synthetic pillar

	vector<p3dm1> pp3f;				//all mesh points of the tip, zoneI, and zoneII to pass to MeshGen
};


#endif
