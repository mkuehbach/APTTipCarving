/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_SOLVERHDL_H__
#define __BICARVING_SOLVERHDL_H__

#include "BiCarving_BiXXHdl.h"


class solverHdl
{
	//top-level construct implementing the respective task-specific worker instances at process rank within the MPI process level parallelism

public:
	solverHdl();
	~solverHdl();

	singlecrystal* sx;
	bicrystal* bx;

	profiler tictoc;

	int myRank;						//my MPI ID in the MPI_COMM_WORLD
	int nRanks;						//total number of MPI processes that work in the world
};

#endif
