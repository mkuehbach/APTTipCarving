/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#include "BiCarving_BiXXHdl.h"

polycrystal::polycrystal()
{
	//calls type default constructs by default behavior
	owner = NULL;
}


polycrystal::polycrystal( solverHdl * own )
{
	owner = own;
	zoneII = tip();
	zoneI = tip();
	pillar = tip();
}


polycrystal::~polycrystal()
{
	//default constructo destruct solverHdl class object internals
	//do not delete owner, only backreference
}


unsigned int polycrystal::process_on_which_side( p3d const & in )
{
	if ( Settings::GBModel == E_GB_EXPLICIT_TRIANGLEMESH ) {
		return gb.robust_relative_position_trianglepatch( in );
	}
	else if ( Settings::GBModel == E_GB_PLANENORMAL_AND_POINT ) {
		return gb.robust_relative_position_plane( in );
	}
	return BOUNDARYINSIDE;
}


void polycrystal::process_px( t3x3 & o1, t3x3 & o2 )
{
	//all measures in nanometer
	//user has to take care that boundary is inside the tip volume because reposition can be arbitrary
	//even with a simple planar boundary now already sufficient setup to sample all inclinations of planar boundaries
	//define zone geometry: zoneII (vacuum capsule, zoneI (guard zone about tip), and the tip itself

	apt_real eps = 0.01;
	apt_real defaultr = 30.0;
	//##MK::15.0 for crystal symmetry operators and systematic rot about x
	//##MK::30.0 for comparison with experiment
	zoneII = tip( 64.35 - defaultr - eps, defaultr - eps, defaultr - eps );
	zoneI = tip( Settings::SXFrustHeight, Settings::SXFrustRadiusB + Settings::SXBlowup, Settings::SXFrustRadiusT + Settings::SXBlowup );
	pillar = tip( Settings::SXFrustHeight, Settings::SXFrustRadiusB, Settings::SXFrustRadiusT );

	cout << zoneII.report_tip() << endl;
	cout << zoneI.report_tip() << endl;
	cout << pillar.report_tip() << endl;

	//are such defined zoneI and pillar in zone II with sufficient guard space remaining to the zoneII bounds?
	//##MK::reconsider
	apt_real debug_guard = 2.0; //nm
	if ( zoneII.mybox.is_inside_box_xy( zoneI.mybox, debug_guard ) == false ) {
		reporting("Insufficient guard between zoneI and zoneII vacuum capsule!"); return;
	}
	if ( zoneII.mybox.is_inside_box_xy( pillar.mybox, debug_guard ) == false ) {
		reporting("Insufficient guard between tip and zoneII vacuum capsule!"); return;
	}

	cout << zoneII.report_tip() << endl;
	cout << zoneI.report_tip() << endl;
	cout << pillar.report_tip() << endl;

	//relocate the three primitives and the boundary to the center of the coordinate system
	p3d zoneIImid0 = zoneII.mybox.mid();
	cout << zoneIImid0 << endl;
	apt_real offsetz = -1.f * zoneIImid0.z;
	zoneII.relocate_z( offsetz );
	zoneI.relocate_z( offsetz );
	pillar.relocate_z( offsetz );

	//relocate the boundary
	gb.container.zmi += offsetz;
	gb.container.zmx += offsetz;
	gb.simple.ptest.z += offsetz; //the normal of a plane does not change if the entire plane is translated through Euclidian space

	//##MK::no relocating of the boundary triangle mesh as user has enforced that it is
	//within a bounding box that is z-axis centered and in positive i.e. the +z halfspace at/above xy plane

	//zoneII, zoneI, and pillar center are now at coordinate system origin (0,0,0)
	cout << "ZoneII, zoneI, and pillar relocated to coordinate center" << endl;

	//define a circumsphere about zoneII centered, as zoneII at right-handed xyz coordinate system origin 0,0,0
	apt_real sphR = sqrt(2.f) / 2.f * max(max(zoneII.mybox.xsz, zoneII.mybox.ysz), zoneII.mybox.zsz);
	aabb3d sphbox = aabb3d( -sphR, +sphR, -sphR, +sphR, -sphR, +sphR );
	sphbox.scale();

	cout << "sphR/sphbox " << sphR << "\t\t" << sphbox << endl;

	//define an aggregate of crystal unit cells in world coordinate system that is large enough to be inside sphbox
	unitcellaggr rve = unitcellaggr( Settings::LatticeConstant, sphbox, FCC );

	cout << rve.report_unitcell() << endl;

	//find all atoms of this aggregate with positions inside this circumsphere
	//##MK::in what follows use this circumsphere-inscribed lattice and rotate to get any crystal orientation we want
	//##MK::such we are no longer dependent on VESTA and create the heavy data, i.e. atomic positions in-place
	//##MK::this avoids I/O, thus improves efficiency, given that by circumsphere is centered at the origin we safe some arithmetics in the in-sphere test
	apt_real SQRsphR = SQR(sphR);
cout << "SQRsphR = " << SQRsphR << endl;

cout << "R1 = " << o1 << endl;
	t3x3 om1f = o1.inverse();
cout << "inv(R1) = " << om1f << endl;

cout << "R2 = " << o2 << endl;
	t3x3 om2f = o2.inverse();
cout << "inv(R2) = " << om2f << endl;

	//identify in which zone the atoms of this rotated atom ball are, if any, and assign corresponding point marks, following the concept of TAPSim to handle marked points
	synstats2 vbs;
	vector<p3dm1> pp3f;
	//now find in which zone and add marks, relocate back just that base plate is at x,y,z=0
	//because we rotated the lattice and not the geometric primitives and zoneII, zoneI, and pillar are centered at origin
	//the offset to relocate all positions is
	//p3d offset_z = p3d( 0.f, 0.f, 0.5*(pillar.FH+pillar.SH) ); //the same as for the pillar

	//random thinning of integration point in zoneII is achieved by random thinning using a PRNG with long period
	mt19937 dice;
	dice.seed( Settings::PRNGSeed );
	//specifically the MersenneTwister algorithm, which should be warmed up
	uniform_real_distribution<apt_real> unifrnd(0.f, 1.f);
	for (size_t i = 0; i < Settings::PRNGDiscard; ++i)
		apt_real discard = unifrnd(dice);
	reporting("MersenneTwister warmed up");

	for ( size_t b = 0; b < rve.base.size(); ++b) {
		for ( int w = rve.wmin; w <= rve.wmax; ++w ) {
			for ( int v = rve.vmin; v <= rve.vmax; ++v ) {
				for ( int u = rve.umin; u <= rve.umax; ++u ) {
					//MK::both crystal halfs are generated from the initial two lattices centered at the same origin!
					p3d ap = rve.get_atom(b, u, v, w); //##MK::for fixed w, v, and base y and z are the same, hence further optimization potential

					//##MK::implicitly assuming that the sphere is large enough about zoneII !
					if ( (SQR(ap.x)+SQR(ap.y)+SQR(ap.z)) <= SQRsphR ) { //inside the circumsphere, so consider at all
						//##MK::first approach
						p3d ap1 = ap;
						p3d ap2 = ap;
						//we rotate only the points in zone 0 and 1
						//accepting either rotated point r1 or r2 based on decision on which site of the boundary

						//MK::the boundary cuts two ideal crystal lattices
						//points rotated by om1f are only accepted if in front of the boundary
						//points rotated by om2f are only accepted if behind the boundary
						//if numerically doubt and both points on the boundary random decision which one is chosen

						//MK::the boundary is either a spatially finite sized triangle mesh or a specific infinite plane
						//in either case the boundary is defined beyond the triangle representation by planes
						//in the first case by checking for positioning in front or behind all triangles, i.e.
						//	also the triangle patch bounding triangles
						//in the second case the situation is trivial

						//##MK::improvements possible here...
						p3d r1 = ap1.active_rotation_relocate( om1f );
						unsigned int decision1 = process_on_which_side( r1 );

						p3d r2 = ap2.active_rotation_relocate( om2f );
						unsigned int decision2 = process_on_which_side( r2 );

if ( pp3f.size() % 100000 == 0 )
	cout << "x/y/z dec1/dec2\t\t" << ap.x << ";" << ap.y << ";" << ap.z << "\t\t" << decision1 << ";" << decision2 << "\n";

						if ( decision1 == BOUNDARYFRONT ) { //##MK ### && ### && decision2 == BOUNDARYBEHIND ) { //CSL points ?
							//take rotated lattice r1 for everything in front of the boundary
							if ( r1.z > (zoneII.mybox.zmi + Settings::SXBaseEps) ) { //"normal" points, not forming the base plate
								if ( zoneII.is_inside_tip( r1 ) == true ) { //most likely case inside zoneII
									if ( pillar.is_inside_tip( r1 ) == true ) {
										vbs.tipatomsA++;
										pp3f.push_back( p3dm1(r1.x, r1.y, r1.z + zoneIImid0.z, Settings::TAPSimTipIdx ) ); //+100 in-place relocating in z-direction to correct for offset
									}
									else { //potentially in zoneI
										if ( zoneI.is_inside_tip( r1 ) == true ) {
											vbs.zoneIatomsA++;
											pp3f.push_back( p3dm1(r1.x, r1.y, r1.z + zoneIImid0.z, Settings::TAPSimVacIdx) );
										}
										else { //really in zoneII
											if ( unifrnd(dice) > Settings::PRNGThinning ) { //every eta-th gets through so every 1-eta not, if eta = 0.0, almost every dice result > 0.0 discard, if eta --> 1 almost no dice result is > 1-small number
												vbs.thinnedoutA++;
											}
											else {
												vbs.zoneIIatomsA++;
												pp3f.push_back( p3dm1(r1.x, r1.y, r1.z + zoneIImid0.z, Settings::TAPSimVacIdx) );
											}
										}
									}
								}
								//else outside zoneII not considered, MeshGen employs own meshing of zones III and IV
							}
							else { //potentially close to base plate points require special treatment
								if ( zoneII.is_inside_tip( r1 ) == true ) { //most likely case inside zoneII
									if ( pillar.is_inside_tip( r1 ) == true ) {
										vbs.tipatomsA++;
										pp3f.push_back( p3dm1(r1.x, r1.y, r1.z + zoneIImid0.z, Settings::TAPSimTipIdx ) ); //+100
										//no marking with DUMMYBASE because no rebasing of z coordinates to 0.f for points inside tip
									}
									else { //potentially in zoneI
										if ( zoneI.is_inside_tip( r1 ) == true ) {
											vbs.baseatomsA++;
											pp3f.push_back( p3dm1(r1.x, r1.y, r1.z + zoneIImid0.z, UINT32MX) );
										}
										else { //really in zoneII but not thinning, because base plate!
											vbs.baseatomsA++;
											pp3f.push_back( p3dm1(r1.x, r1.y, r1.z + zoneIImid0.z, UINT32MX) );
										}
									}
								}
								//else outside zoneII not considered, MeshGen employs own meshing of zones III and IV
							}
						} //not in clearly in front of boundary


						if ( decision2 == BOUNDARYBEHIND ) { //##MK####decision1 == BOUNDARYBEHIND && decision2 == BOUNDARYFRONT ) { //most likely now clearly behind
							//take rotated lattice r2 if behind
							if ( r2.z > (zoneII.mybox.zmi + Settings::SXBaseEps) ) { //"normal" points, not forming the base plate
								if ( zoneII.is_inside_tip( r2 ) == true ) { //most likely case inside zoneII
									if ( pillar.is_inside_tip( r2 ) == true ) {
										vbs.tipatomsB++;
										pp3f.push_back( p3dm1(r2.x, r2.y, r2.z + zoneIImid0.z, Settings::TAPSimTipIdx ) ); //+200 in-place relocating in z-direction to correct for offset
									}
									else { //potentially in zoneI
										if ( zoneI.is_inside_tip( r2 ) == true ) {
											vbs.zoneIatomsB++;
											pp3f.push_back( p3dm1(r2.x, r2.y, r2.z + zoneIImid0.z, Settings::TAPSimVacIdx) );
										}
										else { //really in zoneII
											if ( unifrnd(dice) > Settings::PRNGThinning ) { //every eta-th gets through so every 1-eta not, if eta = 0.0, almost every dice result > 0.0 discard, if eta --> 1 almost no dice result is > 1-small number
												vbs.thinnedoutB++;
											}
											else {
												vbs.zoneIIatomsB++;
												pp3f.push_back( p3dm1(r2.x, r2.y, r2.z + zoneIImid0.z, Settings::TAPSimVacIdx) );
											}
										}
									}
								}
								//else outside zoneII not considered, MeshGen employs own meshing of zones III and IV
							}
							else { //potentially close to base plate points require special treatment
								if ( zoneII.is_inside_tip( r2 ) == true ) { //most likely case inside zoneII
									if ( pillar.is_inside_tip( r2 ) == true ) {
										vbs.tipatomsB++;
										pp3f.push_back( p3dm1(r2.x, r2.y, r2.z + zoneIImid0.z, Settings::TAPSimTipIdx ) ); //+200
										//no marking with DUMMYBASE because no rebasing of z coordinates to 0.f for points inside tip
									}
									else { //potentially in zoneI
										if ( zoneI.is_inside_tip( r2 ) == true ) {
											vbs.baseatomsB++;
											pp3f.push_back( p3dm1(r2.x, r2.y, r2.z + zoneIImid0.z, UINT32MX) );
										}
										else { //really in zoneII but not thinning, because base plate!
											vbs.baseatomsB++;
											pp3f.push_back( p3dm1(r2.x, r2.y, r2.z + zoneIImid0.z, UINT32MX) );
										}
									}
								}
								//else outside zoneII not considered, MeshGen employs own meshing of zones III and IV
							}
						} //neither clearly in front nor behind

						//###########MK::special handling required THINK ABOUT IT
						/*else {}*/
					} //done with processing lattice aggregate point inside sphere
				} //x
			} //y
		} //z for every base atom
	}

	//pp3i defines a collection of intermediate points, some of which we interpret as atoms, some of which as mesh supporting points within zoneII for TAPSim
	cout << "pp3f.size() " << pp3f.size() << endl;
	cout << vbs << endl;

	//define base plate geometry and get all dimensions into nanometer
	for( auto it = pp3f.begin(); it != pp3f.end(); ++it ) {
		//is not of base plate? only SI scaling and continue
		//is of base plate define, SI scaling and continue
		apt_real xx = NANOMETER2METER(it->x);
		apt_real yy = NANOMETER2METER(it->y);
		apt_real zz = NANOMETER2METER(it->z);
		if ( it->m != UINT32MX ) { //no rebasing
			it->x = xx;
			it->y = yy;
			it->z = zz;
		}
		else {
			it->x = xx;
			it->y = yy;
			it->z = 0.f;
			it->m = Settings::TAPSimBseIdx;
		}
	}

	 string mshgenfn = "TAPSIMInputGeometry.SimID." + to_string(Settings::SimID) + ".dat";
	 write_meshgen_conformant_result( pp3f, mshgenfn );

	cout << "Rebase, SI unit scaling, and base plate fixed!" << endl;

	if ( Settings::IOCrystals == true ) {
		string vtkfn = "TAPSIMInputGeometry.SimID." + to_string(Settings::SimID) + ".vtk";
		reconstruction_vtk( pp3f, vtkfn );
	}

	//##MK:: ##### ??? ##### zoneII, zoneI, and pillar aabb3ds not relocated!
}

/*
void bicrystal::define_domain()
{
	//cylinder centroid at coordinate system origin 0,0,0
	world.H = Settings::CylHeight;
	world.R = Settings::CylRadius;
	world.mybox.xmi = -1.f*Settings::CylRadius;
	world.mybox.xmx = +1.f*Settings::CylRadius;
	world.mybox.ymi = -1.f*Settings::CylRadius;
	world.mybox.ymx = +1.f*Settings::CylRadius;
	world.mybox.zmi = -0.5*Settings::CylHeight;
	world.mybox.zmx = +0.5*Settings::CylHeight;
	world.mybox.scale();
	world.center = p3d( 0.f, 0.f, 0.f);
}


void bicrystal::load_grainboundary()
{
	gAB.read_triangulation();
}


void bicrystal::generate_bicrystal()
{
	//##MK::build lattice of crystal gA with atoms in front of the boundary gAB

	//void crystal::build_lattice( cylinder const & box, boundary const & constraint )
	{
		//##MK::define crystal structure
		//sample over positions in cylinder
		//
	//}

	//build lattice of crystal gB with atoms behind the boundary gAB
}


void bicrystal::carve_tip()
{
}


void bicrystal::process()
{
}
*/
