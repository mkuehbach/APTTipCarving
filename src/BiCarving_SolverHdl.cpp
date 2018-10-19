/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#include "BiCarving_SolverHdl.h"

solverHdl::solverHdl()
{
	//calls type default constructs by default behavior

	myRank = MASTER;
	nRanks = SINGLEPROCESS;
}


solverHdl::~solverHdl()
{
	//default constructor destruct solverHdl class object internals
}
