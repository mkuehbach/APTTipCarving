/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#include "BiCarving_Crystallite.h"

unitcellaggr::unitcellaggr()
{
	this->a = 0.f;

	this->umin = 0;
	this->umax = 0;
	this->vmin = 0;
	this->vmax = 0;
	this->wmin = 0;
	this->wmax = 0;

	this->a1 = v3d();
	this->a2 = v3d();
	this->a3 = v3d();

	this->base.clear();
}


unitcellaggr::unitcellaggr(const apt_real _a, const aabb3d unitbox, const unsigned int model )
{
	//initialize fcc base atoms
	a = _a;

	//initialize cubic base vectors
	a1 = v3d( _a*1.f, 0.f, 0.f);
	a2 = v3d( 0.f, _a*1.f, 0.f);
	a3 = v3d( 0.f, 0.f, _a*1.f);

	base.push_back( p3d(0.f, 0.f, 0.f) ); //Al 8x 1/8 = 1
	base.push_back( p3d(0.5, 0.5, 0.f) ); //Al
	base.push_back( p3d(0.f, 0.5, 0.5) ); //Al
	base.push_back( p3d(0.5, 0.f, 0.5) ); //Al 6x 1/2 = 3 Al --> 4 units per EZ okay

	//unitbox gives min/max dimensions in nanometer that we have to fill construct on positive sectors of \mathcal{R}^3
	umin = static_cast<int>(floor(unitbox.xmi / _a));
	umax = static_cast<int>(ceil(unitbox.xmx / _a));
	vmin = static_cast<int>(floor(unitbox.ymi / _a));
	vmax = static_cast<int>(ceil(unitbox.ymx / _a));
	wmin = static_cast<int>(floor(unitbox.zmi / _a));
	wmax = static_cast<int>(ceil(unitbox.zmx / _a));

	//unitbox is axis-aligned to standard orientation 0.0, 0.0, 0.0 Bunge Euler fcc crystal lattice
}

unitcellaggr::~unitcellaggr(){}


p3d unitcellaggr::get_atom(const size_t b, const int u, const int v, const int w)
{
	apt_real uu = static_cast<apt_real>(u);
	apt_real vv = static_cast<apt_real>(v);
	apt_real ww = static_cast<apt_real>(w);

	//##MK::implicit origin at 0,0,0
	p3d res = p3d(
			(base[b].x + uu)*a1.u + (base[b].y + vv)*a2.u + (base[b].z + ww)*a3.u,
			(base[b].x + uu)*a1.v + (base[b].y + vv)*a2.v + (base[b].z + ww)*a3.v,
			(base[b].x + uu)*a1.w + (base[b].y + vv)*a2.w + (base[b].z + ww)*a3.w  );

	return res;
}


string unitcellaggr::report_unitcell()
{
	string str = "imi/imx = " + to_string(umin) + ";" + to_string(umax) + ";" + to_string(vmin) + ";" + to_string(vmax) + ";" + to_string(wmin) + ";" + to_string(wmax) + "\n";
	return str;
}


crystal::crystal()
{
	ori = t3x3();
}

crystal::~crystal()
{
}



boundary::boundary()
{
	container = aabb3d();
}


boundary::~boundary()
{
}


void boundary::read_plane( v3d const & normal, p3d const & ptest )
{
	simple.ounormal = normal;
	simple.ptest = ptest;
	//container is irrelevant boundary has infinite extent
}


void boundary::read_triangulation( const string vtk_io_fn )
{
	trimesh.clear();

	//MK::read 3D hull of a single grain from vtkfname VTK file into local class container object LeanGrain, i.e. pass geometry to actual grain
	ifstream vtkfile;
	string vtkline;
	string datapiece;
	vector<p3d> trianglepatch_vertices;
	aabb3d meshbox = aabb3d();

	vtkfile.open( vtk_io_fn.c_str() );
	if ( vtkfile.is_open() == true ) { //read in file
		unsigned int nheaderlines = 5;
		unsigned int j = 0; //jump over header
		while ( vtkfile.good() == true && j < nheaderlines ) {
			getline( vtkfile, vtkline);
			++j;
		}

		if ( vtkfile.good() == true ) { //read number of disjoint vertices
			getline( vtkfile, vtkline );
			size_t nvertices = 0;
			if ( vtkline.size() > 0 ) { //not empty
				istringstream line( vtkline );
				getline(line, datapiece, ' ');
				getline(line, datapiece, ' ');	nvertices = stol( datapiece.c_str() );
				getline(line, datapiece, ' ');

cout << "GB triangle patch has " << nvertices << " meshpoints" << endl;
				trianglepatch_vertices.reserve(nvertices);

				j = 0;
				double xx = F32MX;
				double yy = F32MX;
				double zz = F32MX; //##MK::absolute coordinates in micrometer!
				while ( vtkfile.good() == true && j < nvertices ) {
					getline( vtkfile, vtkline);
					stringstream line(vtkline);
					//###MK::skip checking format of each line
					line >> xx >> yy >> zz;
					//getline(line, datapiece, ' ');	xx = stod( datapiece.c_str() );
					//getline(line, datapiece, ' ');	yy = stod( datapiece.c_str() );
					//getline(line, datapiece, ' ');	zz = stod( datapiece.c_str() );

					trianglepatch_vertices.push_back( p3d( xx, yy, zz) );

					meshbox.xmi = min(meshbox.xmi, xx);
					meshbox.xmx = max(meshbox.xmx, xx);
					meshbox.ymi = min(meshbox.ymi, yy);
					meshbox.ymx = max(meshbox.ymx, yy);
					meshbox.zmi = min(meshbox.zmi, zz);
					meshbox.zmx = max(meshbox.zmx, zz);
	//cout << xx << "\t\t" << yy << "\t\t" << zz << "\n";

					++j;
				}

				if ( vtkfile.good() == true ) { //read number of disjoint triangles of grain hull
					getline ( vtkfile, vtkline );
					size_t nfaces = 0;
					if ( vtkline.size() > 0 ) { //not an empty line
						istringstream line( vtkline );
						getline(line, datapiece, ' ');
						getline(line, datapiece, ' '); nfaces = stol( datapiece.c_str() );
						getline(line, datapiece, ' ');

cout << "GB triangle patch has " << nfaces << " triangles" << endl;
						j = 0;
						unsigned int skip = 0;
						unsigned int uu = UINT32MX;
						unsigned int vv = UINT32MX;
						unsigned int ww = UINT32MX;
						while( vtkfile.good() == true && j < nfaces ) {
							getline( vtkfile, vtkline);
							stringstream line(vtkline);
							//###MK::skip VTK relevant info how many vertices per type
							line >> skip >> uu >> vv >> ww;

							/*
							getline( vtkfile, vtkline );
							istringstream line ( vtkline );
							getline(line, datapiece, ' '); //skip VTK piece of information how many vertices for thew polygon
							//MK::GraGLeS utilizes clock-wise ordering of triangles in a right-handed coordinate system
							getline(line, datapiece, ' '); uu = stol( datapiece.c_str() );
							getline(line, datapiece, ' '); vv = stol( datapiece.c_str() );
							getline(line, datapiece, ' '); ww = stol( datapiece.c_str() );
							*/

							trimesh.push_back( tri3d(
									trianglepatch_vertices[uu].x,
									trianglepatch_vertices[uu].y,
									trianglepatch_vertices[uu].z,
									trianglepatch_vertices[vv].x,
									trianglepatch_vertices[vv].y,
									trianglepatch_vertices[vv].z,
									trianglepatch_vertices[ww].x,
									trianglepatch_vertices[ww].y,
									trianglepatch_vertices[ww].z ) );
							++j;
						} //done reading id triplets
						//##MK::do not read other potential content from VTK file
					}
				} //triangles read
			} //non-empty vtk file interpreted
		} //finish vertices
		vtkfile.close();
cout << "VTK file " << vtk_io_fn << " loaded successfully vertices/triangles " << trianglepatch_vertices.size() << "/" << trimesh.size() << endl;

		container = meshbox;
		container.scale();
cout << "Mesh bounding box determined " << container << endl;
	}
	else {
		trimesh.clear();
cout << "ERROR::Unable to load file " << vtk_io_fn << endl;
	}
}

void boundary::compute_consistent_ounormals()
{
	cutplanes.reserve( trimesh.size() );

	for ( auto it = trimesh.begin(); it != trimesh.end(); ++it ) {
		v3d v31 = v3d( it->x3 - it->x1, it->y3 - it->y1, it->z3 - it->z1);
		v3d v21 = v3d( it->x2 - it->x1, it->y2 - it->y1, it->z2 - it->z1);
		v3d v31x21 = v31.cross( v21 );
		v31x21.normalize();
		p3d itc = it->center();
		cutplanes.push_back( plane3d( itc, v31x21 ) );
cout << v31x21.u << "\t\t" << v31x21.v << "\t\t" << v31x21.w << "\n";
	}

cout << "Cutplanes with consistent normals defined in total there are " << cutplanes.size() << endl;
}


unsigned int boundary::robust_relative_position_plane( p3d const & p )
{
	//units in nanometer
	double distance = 	simple.ounormal.u * (p.x - simple.ptest.x) +
						simple.ounormal.v * (p.y - simple.ptest.y) +
						simple.ounormal.w * (p.z - simple.ptest.z);
	if ( distance >= (Settings::GBBoundaryThickness + 0.f) ) //##MK::was BOUNDARYEPSDISTANCE + 0.f
		return BOUNDARYFRONT;
	else {
		if ( distance < (0.f - Settings::GBBoundaryThickness) )
			return BOUNDARYBEHIND;
	}
	//]BOUNDARYEPSDISTANCE 0.f BOUNDARYEPSDISTANCE[
	return BOUNDARYINSIDE;
}


unsigned int boundary::robust_relative_position_trianglepatch( p3d const & p )
{
	//for arbitrarily shaped surface patches the rel. pos (behind on or in front) of a point to a surface is nontrivial
	//namely even if local normals are outer unit normals only and consistent but body has convex + concave sections
	//point can lay in front of some triangle cutplanes while behind most others
	//on top of that most surface triangle meshers strictly speaking yield numerically a triangle soup only
	//that is the triangles do not form a necessarily numerically accurately connected network with collapsing edges
	//even though during rendering "the mesh seems qualitatively apparently appears connected" as order of magnitude machine precision
	//this is expected as close to machine precision differences cannot be detected in the visualization
	//therefore we have to not only test if point is behind or in front of all triangles of the patch but
	//instead where the point gets projected onto the triangle and whether this projected position
	//falls in the triangle or not
	//in case it is inside we can check for front or back, in case the normals have to be fully consistent
	//further we need to account for triangle soup mentioned afore
	//given that input triangle mesh not necessarily has full double precision or even if we do not get
	//sufficient precision to solve the challenge using multiprecision (beyond double precision) floating
	//point arithmetics we solve greedily by enlarging the triangles by a guard zone eps
	//specifically we compute for each triangle edge the 1) its outer unit normal
	//next 2) we displace the edge by an amount eps in the direction of the ounormal
	//finally we compute the new triangle vertices as intersections of the eps displaced line segments 01, 12, 20

	//####MK::furthermore we can easily validate that for wavy boundary shapes below scheme requires modifications
	//as projection may fall in triangle though despite consistent normals point is outside...
	//####MK::in a nutshell::so far only convex shapes
	//##MK::here so far account only for projection into triangle
	//####MK::it remains to enlarge the triangle do so during consistent_ounormal check

	//units in nanometer
	//##MK::accelerate procedure by pruning candidate triangles...
	for ( auto it = trimesh.begin(); it != trimesh.end(); ++it ) {
		//https://gamedev.stackexchange.com/questions/28781/easy-way-to-project-point-onto-triangle-or-plane
		v3d uu = v3d( it->x2 - it->x1, it->y2 - it->y1, it->z2 - it->z1 );
		v3d vv = v3d( it->x3 - it->x1, it->y3 - it->y1, it->z3 - it->z1 );
		v3d nn = uu.cross( vv );
		v3d ww = v3d( p.x - it->x1, p.y - it->y1, p.z - it->z1 );

		// Barycentric coordinates of the projection P′of P onto T:
		// gamma=[(u cross w)⋅n]/n^2
		// beta=[(w cross v)⋅n]/n^2
		// alpha = 1-beta-gamma
		// The point P star lies inside T if:
			    /*return ((0 <= alpha) && (alpha <= 1) &&
				            (0 <= beta)  && (beta  <= 1) &&
				            (0 <= gamma) && (gamma <= 1));*/
		//MK::given that point is only inside if [0,1] forall three angles we can exit early which is in particular
		//using no triangle candidate pruning is beneficial as most points are not in the triangle
		//therefore high chance of early reject i.e. triangle does not participate in the voting where
		//the point is relative to surface
		double gamma = uu.cross(ww).dot(nn) / nn.dot(nn);
		if ( gamma > 0.f || gamma > 1.f )
			continue;

		double beta = ww.cross(vv).dot(nn) / nn.dot(nn);
		if ( beta > 0.f || beta > 1.f )
			continue;

		double alpha = 1.f - gamma - beta;
		if ( alpha > 0.f || alpha > 1.f )
			continue;

		//not continued so projected point is in volume of infinitely high triangular prism it participates in voting
		//so on which side of the boundary?


		//##MK::assumption the point is behind the boundary if it is not in front of every cutplane of the triangle patch
		//this requires consistent outer unit normals
		nn.normalize();
		if ( nn.dot(ww) >= (Settings::GBBoundaryThickness + 0.f) )
			continue; //this guy voted in favor that it is still possible that the point is in front of all other triangles
		else
			return BOUNDARYBEHIND;
	}
	//not returned BEHIND so all who vote did so in favor of INFRONT
	return BOUNDARYFRONT;

	//##MK::there is another peculiarity to watch out for namely the case of point exactly on the boundary...
}


/*
inline bool boundary::in_front_of_me( p3d const & p ) const
{
	//##MK::test against if in front of every triangle in triangulation assuming a consistent normal
	for(auto it = simplex.begin(); it != simplex.end(); ++it) {
		if ( it->in_front_of_me( p ) == true )
			continue;
		else
			return false;
	}
	return true;
}
*/

