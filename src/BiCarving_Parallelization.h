/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_PARALLELIZATION_H__
#define __BICARVING_PARALLELIZATION_H__

#include "BiCarving_Abinitio.h"

#include <mpi.h>
#include <omp.h>

#define MASTER									0
#define	SINGLEPROCESS							1

#define	MPI_COMM_WORLD_OMP_GET_NUM_THREADS		1 //12

//implicitly performance affecting choices
#define SEQIO_READ_CACHE						((10)*(1024)*(1024)) //bytes
#define MPIIO_READ_CACHE						((10)*(1024)*(1024)) //bytes

#ifdef EMPLOY_SINGLEPRECISION
	#define SIMDREGISTER_WIDTH					(8) //elements assuming eight 32 bit floats to fit in 256bit wide SIMD register of contemporary processor
#else
	#define SIMDREGISTER_WIDTH					(4) //256bit can take four 64bit double at a time
#endif

#endif
