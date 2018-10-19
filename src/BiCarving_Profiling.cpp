/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#include "BiCarving_Profiling.h"


void profiler::prof(const string whichenv,
		const unsigned short category, const unsigned short parallelism,
		const double st, const double en )
{
	evn.push_back( plog(st, en, whichenv, category, parallelism) );
}

size_t profiler::get_nentries( void ){
	return evn.size();
}

bool SortProfLogAscWallClock( plog & first, plog & second )
{
	bool comp = first.get_dt() < second.get_dt();
	return comp;
}

void profiler::spit_profiling( const unsigned int simid, const int rank )
{
	//##MK::further optimization aand convenience tasks: bundle all in one file, incr ID and so forth
	//##MK::suboptimal... one file per rank
	string fn = "PARAPROBE.SimID." + to_string(simid) + ".Rank." + to_string(rank) + ".MyProfiling.csv";

	ofstream csvlog;
	csvlog.open(fn.c_str(), ofstream::out | ofstream::trunc);
	if (csvlog.is_open() == true) {
		//header
		csvlog << "What;Category;ParallelismInfo;WallClock;CumulatedWallClock;CDF;WallClockFraction\n";
		csvlog<< ";;;s;s;1;1\n";
		csvlog << "What;Category;ParallelismInfo;MPI_Wtime;CumulatedWallClock;CDF;WallClockFraction\n";

		//build map of categories
		map<unsigned short, string> categories;
		categories[BICRV_XX] = "BICRV_XX";
		categories[BICRV_IO] = "BICRV_IO";
		categories[BICRV_GEO] = "BICRV_GEO";
		categories[BICRV_CRY] = "BICRY_CRY";
		categories[BICRV_UTL] = "BICRY_UTL";
		map<unsigned short, string> parallelism;
		parallelism[BICRV_IS_PARALLEL] = "PARALLEL";
		parallelism[BICRV_IS_SEQUENTIAL] = "SEQUENTIAL";

		//sort events increasing wallclock time
		sort( evn.begin(), evn.end(), SortProfLogAscWallClock);

		//compute total time
		double dt_total = 0.f;
		for(auto it = evn.begin(); it != evn.end(); ++it) { dt_total += it->get_dt(); }

		//report
		double dt_cumsum = 0.f;
		for (auto it = evn.begin(); it != evn.end(); ++it) {
			dt_cumsum += it->get_dt();
			auto cat = categories.find(it->get_typ());
			auto par = parallelism.find(it->get_pll());
			csvlog << it->get_what() << ";" << cat->second << ";" << par->second << ";" << it->get_dt();
			csvlog << ";" << dt_cumsum << ";" << (dt_cumsum / dt_total) << ";" << (it->get_dt() / dt_total) << endl;
		}

		csvlog.flush();
		csvlog.close();
	}
	else {
		stopping("Unable to write local processing files");
	}
}

//program profiling should use double precision in general as MPI_Wtime() and omp_get_wtime() fires in double precision
