/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_VTKIO_H__
#define __BICARVING_VTKIO_H__

#include "BiCarving_Crystallite.h"

void write_meshgen_conformant_result( vector<p3dm1> const & p, const string msh_io_fn );
void reconstruction_vtk( vector<p3dm1> const & p, const string vtk_io_fn );

#endif
