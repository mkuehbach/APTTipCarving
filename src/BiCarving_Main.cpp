/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#include "BiCarving_SolverHdl.h"

//parameter handshake
#define SIMID									1
#define CONTROLFILE								2

#define AXIS1									3
#define AXIS2									4
#define AXIS3									5
#define AXIS4									6

#define OM11	3
#define OM12	4
#define OM13	5
#define OM21	6
#define OM22	7
#define OM23	8
#define OM31	9
#define OM32	10
#define OM33	11

void helloworld ( int pargc, char** pargv )
{
	cout << "Starting up BicrystalCarving v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION;
	if ( VERSION_BETASTAGE == 1 )
		cout << " beta stage" << endl;
	else
		cout << endl;

	cout << "ERRORs are terminating in every case..." << endl;
	if ( pargc < 2 ) {
		std::cout << "\t\tERROR::We need at least a simulation id <unsigned int> and an XML control file <*.xml> before doing something useful!" << endl;
		return;
	}
}


bool init(  int pargc, char** pargv )
{
	Settings::SimID = stoul( pargv[SIMID] );
	try {
		Settings::readXML(pargv[CONTROLFILE]);
	}
	catch (exception& e) { 
		cout << endl << "\t\tERROR::Unable to parse control file! Details:\n" << e.what() << endl; return false;
	}
	if ( Settings::checkUserInput() == false ) {
		cout << endl << "\t\tERROR::Control file settings failed the validity check!" << endl; return false;
	}
	else {
		cout << endl << "Input is valid under SimulationID = " << "SimID." <<  Settings::SimID << endl;
	}
	cout << "All console prompts which follow are intended for debugging purposes only..." << endl;
	cout << endl << endl;

	return true;
}


//genering global functions to report state, warnings and erros
void reporting( const string what ) {
	cout << "VERBOSE::" << what << endl;
}
void reporting( const int rank, const string what ) {
	cout << "VERBOSE::" << rank << what << endl;
}

void complaining( const string what ) {
	cout << "WARNING::" << what << endl;
}

void stopping( const string what ) {
	cout << "ERROR::" << what << endl;
}


int main(int argc, char** argv)
{    
//SETUP PROGRAM AND PARAMETER
	helloworld( argc, argv );
	if ( init( argc, argv ) == false ) {
		return 0;
	}

//INIT MPI and go MPI process parallel with hybrid OpenMP threading capability, funneled means only main thread will make MPI calls
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	int supportlevel_provided;
	MPI_Init_thread( &argc, &argv, supportlevel_desired, &supportlevel_provided);
	//now we are parallel...

	double gtic = MPI_Wtime();
	int nr = 1;
	int r = MASTER;
	if ( supportlevel_provided < supportlevel_desired ) {
		stopping("Insufficient threading capabilities of the MPI library!");
		MPI_Finalize(); //required because threading insufficiency does not imply process is incapable to work at all
		return 0;
	}
	else { 
		MPI_Comm_size(MPI_COMM_WORLD, &nr);
		MPI_Comm_rank(MPI_COMM_WORLD, &r);
	}
	reporting( r, "-th MPI process initialized, we are now MPI_COMM_WORLD parallel using MPI_THREAD_FUNNELED");

//INITIALIZE MPI SOLVER INSTANCE AND LOAD DATASET
	//generate MPI rank solverHdl instance managing which of the parameter configuration rank 0 does
	int localhealth = 1;
	solverHdl* hdl = NULL;
	try {
		hdl = new solverHdl;
		hdl->sx = NULL;
		hdl->bx = NULL;
		if ( Settings::WhatToDo == E_SINGLECRYSTAL ) {
			hdl->sx = new singlecrystal( hdl );
		}
		else { // ( Settings::WhatToDo == E_BICRYSTAL ) {
			hdl->bx = new bicrystal( hdl );
		}

		hdl->myRank = r;
		hdl->nRanks = nr;
	}
	catch (bad_alloc &mpiexc) {
		localhealth = 0;
		stopping("Unable to allocate a solverInstance");
	}
	//were all processes, who should have, able to generate a process-local solver class object instance?
	int globalhealth = 0;
	MPI_Allreduce(&localhealth, &globalhealth, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if ( globalhealth != nr ) {
		stopping("Not all processes were able to allocate a solverInstance");
		MPI_Finalize();
		return 0;
	}

//all processes now have a solver instance, well then set it up

	if ( Settings::WhatToDo == E_SINGLECRYSTAL ) {
/*
		//CRYSTAL LATTICE FCC SYMMETRY ELEMENTS
		axisangle ax_dummy = axisangle( stod(argv[AXIS1]), stod(argv[AXIS2]), stod(argv[AXIS3]), stod(argv[AXIS4]) );
		t3x3 om_dummy = t3x3();
*/

/*
		//SYSTEMATIC ROTATIONS ABOUT X
		//axisangle ax_dummy = axisangle( stod(argv[AXIS1]), stod(argv[AXIS2]), stod(argv[AXIS3]), stod(argv[AXIS4]) );
		//t3x3 om_dummy = t3x3();

		//COMPARISON WITH EXPERIMENTS
		axisangle ax_dummy = axisangle();
		t3x3 om_dummy = t3x3( 	stod(argv[OM11]), stod(argv[OM12]), stod(argv[OM13]),
								stod(argv[OM21]), stod(argv[OM22]), stod(argv[OM23]),
								stod(argv[OM31]), stod(argv[OM32]), stod(argv[OM33]) );
*/

		//COMPARISON WITH EXPERIMENTS USING INVERSE OF INPUT MATRIX FOR ROTATING VECTORS
		axisangle ax_dummy = axisangle();
		t3x3 om_dummy = t3x3( 	stod(argv[OM11]), stod(argv[OM12]), stod(argv[OM13]),
								stod(argv[OM21]), stod(argv[OM22]), stod(argv[OM23]),
								stod(argv[OM31]), stod(argv[OM32]), stod(argv[OM33])   );

		hdl->sx->process( ax_dummy, om_dummy );
	}
	else {
		hdl->bx->process();
	}

//delete solverHdl, deconstruct parallel processes, and exit
	double gtoc = MPI_Wtime();
	string mess = "Elapsed time on process " + to_string(hdl->myRank) + " " + to_string((gtoc-gtic)) + " seconds";
	reporting( mess );
	delete hdl->sx; hdl->sx = NULL;
	delete hdl->bx; hdl->bx = NULL;
	delete hdl; hdl = NULL;

//better have all processes joining to exit program cooperatively
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}


/*
//test analytical inverse of 3x3
t3x3 AA = t3x3( -1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f );
cout << "AA = " << AA << endl;
t3x3 invAA = AA.inverse();
cout << "invAA = " << invAA << endl;
return 0;
*/
