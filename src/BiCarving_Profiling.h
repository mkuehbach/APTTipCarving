/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#ifndef __BICARVING_PROFILING_H__
#define __BICARVING_PROFILING_H__

#include "BiCarving_CGALInterface.h"


//program profiling should use double precision in general as
//MPI_Wtime() and omp_get_wtime() fires in double precision

//type of computational operations
#define BICRV_XX				0		//default, unspecified
#define BICRV_IO				1		//I/O
#define BICRV_GEO				2		//computational geometry
#define BICRV_CRY				3		//crystallite
#define BICRV_UTL				4		//utilities

#define BICRV_IS_PARALLEL		0
#define BICRV_IS_SEQUENTIAL		1

class plog
{
public:
	plog() : dt(0.0), tstart(0.0), tend(0.0), what(""), typ(BICRV_XX), pll(BICRV_IS_SEQUENTIAL) {}
	plog(const double _dt, const string s, const unsigned short t, const unsigned short p) :
		dt(_dt), tstart(0.0), tend(0.0), what(s), typ(t), pll(p) {}
	plog(const double _ts, const double _te, const string s, const unsigned short t, const unsigned short p)
		: tstart(_ts), tend(_te), what(s), typ(t), pll(p) {
		dt = _te - _ts;
	}
	~plog(){}

	double get_dt(){
		return dt;
	}
	double get_tstart(){
		return tstart;
	}
	double get_tend() {
		return tend;
	}
	string get_what() {
		return what;
	}
	unsigned short get_typ() {
		return typ;
	}
	unsigned short get_pll() {
		return pll;
	}

private:
	double dt;
	double tstart;
	double tend;
	string what;
	unsigned short typ;		//task identifier
	unsigned short pll;		//parallelism identifier
};



class profiler
{
public:
	profiler() {};
	~profiler() {};

	void prof(const string whichenv, const unsigned short category,
			const unsigned short parallelism, const double st, const double en);
	size_t get_nentries( void );
	void spit_profiling( const unsigned int simid, const int rank );

private:
	vector <plog> evn;
};


#endif
