/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_BIXXHDL_H__
#define __BICARVING_BIXXHDL_H__

#include "BiCarving_SiXHdl.h"


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

#endif
