/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_CGALINTERFACE_H__
#define __BICARVING_CGALINTERFACE_H__

#include "BiCarving_BSIMD.h"

//comment this line out to remove all compilation dependencies on and capability of utilizing the CGAL library
//#define UTILIZE_CGAL

#ifdef UTILIZE_CGAL

	//required for Convex_hull_3 functionality
	#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	#include <CGAL/Polyhedron_items_with_id_3.h>
	#include <CGAL/Surface_mesh.h>
	#include <CGAL/convex_hull_3.h>

	typedef CGAL::Exact_predicates_inexact_constructions_kernel			K;
	typedef CGAL::Polyhedron_3<K,CGAL::Polyhedron_items_with_id_3>		Polyhedron_3;
	typedef K::Point_3													Point_3;
	typedef CGAL::Surface_mesh<Point_3>									Surface_mesh;


	//required for Alpha_shape_3 functionality
	//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	#include <CGAL/Delaunay_triangulation_3.h>
	#include <CGAL/Alpha_shape_3.h>

	typedef CGAL::Alpha_shape_vertex_base_3<K>							Vb;
	typedef CGAL::Alpha_shape_cell_base_3<K>							Fb;
	typedef CGAL::Triangulation_data_structure_3<Vb,Fb>					Tds;
	typedef CGAL::Delaunay_triangulation_3<K,Tds,CGAL::Fast_location>	Delaunay;
	typedef CGAL::Alpha_shape_3<Delaunay>								Alpha_shape_3;

	typedef K::Point_3													Point;
	typedef Alpha_shape_3::Alpha_iterator								Alpha_iterator;
	typedef Alpha_shape_3::NT											NT;

#endif


#endif
