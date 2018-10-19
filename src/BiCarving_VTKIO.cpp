/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#include "BiCarving_VTKIO.h"

void write_meshgen_conformant_result( vector<p3dm1> const & p, const string msh_io_fn )
{
	//MK::write VTK file showing the positions of all ions in reconstructed space
	if ( p.size() == 0 ) {
		reporting("MeshGen reporting finds no p3dm1 objects in p!");
		return;
	}

	ofstream dat;
	dat.open( msh_io_fn.c_str() );
	if ( dat.is_open() == true ) {
		dat << "ASCII " << p.size() << " 0 0\n";
		for(size_t i = 0; i < p.size(); ++i) {
			dat << p[i].x << "\t" <<  p[i].y << "\t" << p[i].z << "\t" << p[i].m << "\n";
		}
		dat << endl;
		dat.flush();
		dat.close();
	}
}

void reconstruction_vtk( vector<p3dm1> const & p, const string vtk_io_fn )
{
	//MK::write VTK file showing the positions of all ions in reconstructed space
	if ( p.size() == 0 ) {
		reporting("VTK reporting finds no p3dm1 objects in p!");
		return;
	}

	ofstream vtk;
	vtk.open( vtk_io_fn.c_str() );
	if ( vtk.is_open() == true ) {
		vtk << "# vtk DataFile Version 2.0\n"; //header
		vtk << "PARAPROBE TAPSIMInputGeometry " << Settings::SimID << "\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		size_t nevents = p.size();
		vtk << "POINTS " << nevents << " double\n";
		for(size_t i = 0; i < p.size(); ++i) {
			vtk << p.at(i).x << " " << p.at(i).y << " " << p.at(i).z << "\n";
		}
		vtk << "\n";
		vtk << "VERTICES " << nevents << " " << 2*nevents << "\n";
		for ( size_t e = 0; e < nevents; ++e ) {
			vtk << 1 << " " << e << "\n";
		}
		//MK::ranged ion type as field data for coloring in Paraview
		vtk << "POINT_DATA " << nevents << "\n";
		vtk << "FIELD FieldData 1\n";
		vtk << "IonType 1 " << nevents << " float\n"; //do not make int because numeric_limits<unsigned int>::max() is not readable then
		for(size_t i = 0; i < p.size(); ++i) {
			vtk << p.at(i).m << "\n";
		}
		vtk << endl;
		vtk.flush();
		vtk.close();
	}
	//cout << "VTK Reconstruction ion locations written in " << (toc - tic) << " seconds" << endl;
}
