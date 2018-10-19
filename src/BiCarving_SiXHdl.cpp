/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#include "BiCarving_SiXHdl.h"

singlecrystal::singlecrystal()
{
	//calls type default constructs by default behavior
	owner = NULL;
	gA = crystal();
	zoneII = tip();
	zoneI = tip();
	pillar = tip();
}

singlecrystal::singlecrystal( solverHdl * own )
{
	owner = own;
	gA = crystal();
	zoneII = tip();
	zoneI = tip();
	pillar = tip();
}


singlecrystal::~singlecrystal()
{
	//do not delete owner, only backreference
}


/*
void singlecrystal::process( axisangle const & targetori )
{
	//define crystal orientation and geometry of zoneII, zoneI, and the tip
	gA.ori = t3x3( 	1.f, 0.f, 0.f,
					0.f, 1.f, 0.f,
					0.f, 0.f, 1.f  ); //##MK::import from xml
	zoneII = cylinder( Settings::CylHeight, Settings::CylRadius );

	//all measures in nanometer first
	zoneI = tip( Settings::SXFrustHeight, Settings::SXFrustRadiusB+Settings::SXBlowup, Settings::SXFrustRadiusT+Settings::SXBlowup );
	pillar = tip( Settings::SXFrustHeight, Settings::SXFrustRadiusB, Settings::SXFrustRadiusT );

	//are such defined zoneI and pillar in zone II with sufficient guard space remaining to the zoneII bounds?
	apt_real debug_guard = 5.0; //nm
	if ( zoneII.mybox.is_inside_box( zoneI.mybox, debug_guard ) == false ) {
		reporting("Insufficient guard between zoneI and zoneII vacuum capsule!"); return;
	}
	if ( zoneII.mybox.is_inside_box( pillar.mybox, debug_guard ) == false ) {
		reporting("Insufficient guard between tip and zoneII vacuum capsule!"); return;
	}

	cout << gA.ori << endl;
	cout << zoneII.report_cylinder() << endl;
	cout << zoneI.report_tip() << endl;
	cout << pillar.report_tip() << endl;

	//zoneII, zoneI, and pillar center at coordinate system origin (0,0,0)

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
	vector<p3d> pp3i;
	for ( size_t b = 0; b < rve.base.size(); ++b) {
		for ( int w = rve.wmin; w <= rve.wmax; ++w ) {
			for ( int v = rve.vmin; v <= rve.vmax; ++v ) {
				for ( int u = rve.umin; u <= rve.umax; ++u ) {
					p3d ap = rve.get_atom(b, u, v, w); //##MK::for fixed w, v, and base y and z are the same, hence further optimization potential
					if ( (SQR(ap.x)+SQR(ap.y)+SQR(ap.z)) <= SQRsphR ) {
						pp3i.push_back( p3d( ap.x, ap.y, ap.z ) );
//						if ( pp3i.size() % 100000 != 0 )
//							continue;
//						else
//							cout << "b/w/v/u/sz = " << b << "\t\t" << w << "\t\t" << v << "\t\t" << u << "\t\t" << pp3i.size() << endl;
					}
				}
			}
		} //for every base atom
	}
	//pp3i defines a collection of intermediate points, some of which we interpret as atoms, some of which as mesh supporting points within zoneII for TAPSim
	cout << "pp3i.size() " << pp3i.size() << endl;

	//next we rotate all points into the desired final orientation of the pillar we seek to mill out

//	istringstream line( Settings::AxisAngleSX );
//	string datapiece = "";
//	getline( line, datapiece, ';'); apt_real i1 = stof(datapiece.c_str());
//	getline( line, datapiece, ';'); apt_real i2 = stof(datapiece.c_str());
//	getline( line, datapiece, ';'); apt_real i3 = stof(datapiece.c_str());
//	getline( line, datapiece, ';'); apt_real i4 = stof(datapiece.c_str());
//	//##MK::assuming i1,i2,i3 making a normalized
//cout << "i1234 = " << i1 << "\t\t" << i2 << "\t\t" << i3 << "\t\t" << i4 << endl;
//	axisangle Rtarget = axisangle( i1, i2, i3, DEGREE2RADIANT(i4) );

	axisangle Rtarget = targetori;
	t3x3 R = Rtarget.ax2om();
	//t3x3 R = gA.ori;
cout << "R = " << R << endl;
	apt_real alpha = 0.5 / 180.0 * PI;
	t3x3 R_active_about_x = t3x3(   1.f, 0.f, 0.f,
								0.f, cos(alpha), -sin(alpha),
								0.f, sin(alpha), cos(alpha) ); //##MK::more efficiency with quaternions
	t3x3 R_active_about_y = t3x3(   cos(alpha), 0.f, sin(alpha),
			  	  	  	  	  	0.f, 1.f, 0.f,
								-sin(alpha), 0.f, cos(alpha) );
	t3x3 R_active_about_z = t3x3(  cos(alpha), -sin(alpha), 0.f,
								sin(alpha), cos(alpha), 0.f,
								0.f, 0.f, 1.f );
	t3x3 RxR = R.premultiplyR1( R_active_about_x );
	t3x3 RyxR = RxR.premultiplyR1( R_active_about_y );
	t3x3 RzyxR = RyxR.premultiplyR1( R_active_about_z );
cout << "Rfinal = " << RzyxR << endl;

	for( auto it = pp3i.begin(); it != pp3i.end(); ++it ) {
		apt_real xx = it->x;
		apt_real yy = it->y;
		apt_real zz = it->z;
		//"active" repositioning of position in fixed coordinate system by premultiply with orientation rotation matrix
		it->x = RzyxR.a11 * xx + RzyxR.a12 * yy + RzyxR.a13 * zz;
		it->y = RzyxR.a21 * xx + RzyxR.a22 * yy + RzyxR.a23 * zz;
		it->z = RzyxR.a31 * xx + RzyxR.a32 * yy + RzyxR.a33 * zz;
	}

	//identify in which zone the atoms of this rotated atom ball are, if any, and assign corresponding point marks, following the concept of TAPSim to handle marked points
	synstats vbs;

	//random thinning of integration point in zoneII is achieved by random thinning using a PRNG with long period
	mt19937 dice;
	dice.seed( Settings::PRNGSeed );
	//specifically the MersenneTwister algorithm, which should be warmed up
	uniform_real_distribution<apt_real> unifrnd(0.f, 1.f);
	for (size_t i = 0; i < Settings::PRNGDiscard; ++i)
		apt_real discard = unifrnd(dice);

	reporting("MersenneTwister warmed up");

	//now find in which zone and add marks, relocate back just that base plate is at x,y,z=0
	//because we rotated the lattice and not the geometric primitives and zoneII, zoneI, and pillar are centered at origin
	//the offset to relocate all positions is
	p3d offset_z = p3d( 0.f, 0.f, 0.5*(pillar.FH+pillar.SH) ); //the same as for the pillar

	for (auto it = pp3i.begin(); it != pp3i.end(); ++it) {
		if ( zoneII.is_inside_cylinder( *it ) == true ) { //most likely case inside zoneII
			if ( pillar.is_inside_tip( *it ) == true ) { //relocate
				if ( (it->z + offset_z.z) > Settings::SXBaseEps ) {
					pp3f.push_back( p3dm1(
						NANOMETER2METER(it->x+offset_z.x),
						NANOMETER2METER(it->y+offset_z.y),
						NANOMETER2METER(it->z+offset_z.z), //no z clamping for tip atoms
						Settings::TAPSimTipIdx ) );
				}
				else {
					pp3f.push_back( p3dm1(
						NANOMETER2METER(it->x+offset_z.x),
						NANOMETER2METER(it->y+offset_z.y),
						0.f,
						Settings::TAPSimTipIdx ) );
				}
				vbs.tipatoms++;
			}
			else { //potentially in zoneI
				if ( zoneI.is_inside_tip( *it ) == true ) {
					if ( (it->z + offset_z.z) > Settings::SXBaseEps ) {
						pp3f.push_back( p3dm1(
								NANOMETER2METER(it->x + offset_z.x),
								NANOMETER2METER(it->y + offset_z.y),
								NANOMETER2METER(it->z + offset_z.z),
								Settings::TAPSimVacIdx ) );
					}
					else { //rebase z coordinate in case z close to zero
						pp3f.push_back( p3dm1(
							NANOMETER2METER(it->x + offset_z.x),
							NANOMETER2METER(it->y + offset_z.y),
							0.f,
							( ((SQR(it->x+offset_z.x)+SQR(it->y+offset_z.y)) > SQR(pillar.RB)) ?
									Settings::TAPSimBseIdx : Settings::TAPSimTipIdx)  )); //##MK::else case here set to tip!
					}
					vbs.zoneIatoms++;
				}
				else { //is in zone II
					if ( unifrnd(dice) > Settings::PRNGThinning ) { //every eta-th gets through so every 1-eta not, if eta = 0.0, almost every dice result > 0.0 discard, if eta --> 1 almost no dice result is > 1-small number
						vbs.thinnedout++;
						continue;
					}
					else {
						apt_real lift_zoneII = it->z + offset_z.z;
						if ( lift_zoneII >= 0.f ) { //cut cylinder part in the negative
							if ( lift_zoneII > Settings::SXBaseEps ) {
								pp3f.push_back( p3dm1(
									NANOMETER2METER(it->x + offset_z.x),
									NANOMETER2METER(it->y + offset_z.y),
									NANOMETER2METER(it->z + offset_z.z),
									Settings::TAPSimVacIdx ) );
							}
							else {
								pp3f.push_back( p3dm1(
									NANOMETER2METER(it->x + offset_z.x),
									NANOMETER2METER(it->y + offset_z.y),
									0.f,
									( ((SQR(it->x+offset_z.x)+SQR(it->y+offset_z.y)) > SQR(pillar.RB)) ?
										Settings::TAPSimBseIdx : Settings::TAPSimTipIdx) ));
							}
						}
						vbs.zoneIIatoms++;
					}
				}
			}
		}
		//else { continue; } //not considered, MeshGen employs own meshing of zones III and IV
	}

	cout << "pp3f.size() " << pp3f.size() << endl;
	cout << vbs << endl;


	string mshgenfn = "TAPSIMInputGeometry.SimID." + to_string(Settings::SimID) + ".dat";
	write_meshgen_conformant_result( pp3f, mshgenfn );

	cout << "Rebase, SI unit scaling, and base plate fixed!" << endl;

	if ( Settings::IOCrystals == true ) {
		string vtkfn = "TAPSIMInputGeometry.SimID." + to_string(Settings::SimID) + ".vtk";
		reconstruction_vtk( pp3f, vtkfn );
	}
}
*/

void singlecrystal::process(  axisangle const & targetori_ax, t3x3 const & targetori_om  )
{
	//define geometry: zoneII (vacuum capsule, zoneI (guard zone about tip), and the tip itself
	//all measures in nanometer
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

	//relocate the three primitives to the center of the coordinate system
	p3d zoneIImid0 = zoneII.mybox.mid();
	cout << zoneIImid0 << endl;
	zoneII.relocate_z( -1.f * zoneIImid0.z );
	zoneI.relocate_z( -1.f * zoneIImid0.z );
	pillar.relocate_z( -1.f * zoneIImid0.z );

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
	vector<p3d> pp3i;
	for ( size_t b = 0; b < rve.base.size(); ++b) {
		for ( int w = rve.wmin; w <= rve.wmax; ++w ) {
			for ( int v = rve.vmin; v <= rve.vmax; ++v ) {
				for ( int u = rve.umin; u <= rve.umax; ++u ) {
					p3d ap = rve.get_atom(b, u, v, w); //##MK::for fixed w, v, and base y and z are the same, hence further optimization potential
					if ( (SQR(ap.x)+SQR(ap.y)+SQR(ap.z)) <= SQRsphR ) {
						pp3i.push_back( p3d( ap.x, ap.y, ap.z ) );
//						if ( pp3i.size() % 100000 != 0 )
//							continue;
//						else
//							cout << "b/w/v/u/sz = " << b << "\t\t" << w << "\t\t" << v << "\t\t" << u << "\t\t" << pp3i.size() << endl;
					}
				}
			}
		} //for every base atom
	}
	//pp3i defines a collection of intermediate points, some of which we interpret as atoms, some of which as mesh supporting points within zoneII for TAPSim
	cout << "pp3i.size() " << pp3i.size() << endl;

	//next we rotate all points into the desired final orientation of the pillar we seek to mill out
	cout << "Next now rotate all point into the desired final ori" << endl;

/*
 	//CRYSTAL SYMMETRY OPERATORS
	axisangle Rtarget = targetori_ax;
	t3x3 R = Rtarget.ax2om();
cout << "Rinput = " << R << endl;
	apt_real th = 1.0 / 180.0 * PI;
	t3x3 R_active_about_x = t3x3(   1.f, 0.f, 0.f,				0.f, cos(th), -sin(th),		0.f, sin(th), cos(th) ); //##MK::more efficiency with quaternions
	t3x3 R_active_about_y = t3x3(   cos(th), 0.f, sin(th),		0.f, 1.f, 0.f,				-sin(th), 0.f, cos(th) );
	t3x3 R_active_about_z = t3x3(  cos(th), -sin(th), 0.f,		sin(th), cos(th), 0.f,		0.f, 0.f, 1.f );
	t3x3 RxR = R.premultiplyR1( R_active_about_x );
	t3x3 RyxR = RxR.premultiplyR1( R_active_about_y );
	t3x3 RzyxR = RyxR.premultiplyR1( R_active_about_z );
	R = RzyxR;
cout << "Rfinal = " << R << endl;
*/

/*	//consistency of systematic about x rotated
	axisangle Rtarget = targetori_ax;
cout << "uvwtheta = " << targetori_ax.v1 << ";" << targetori_ax.v2 << ";" << targetori_ax.v3 << ";" << targetori_ax.theta << endl;
	t3x3 R = Rtarget.ax2om();
*/

/*
	//SYSTEMATIC ROTATION ABOUT X
	axisangle RR = targetori_ax;
	//for systematic rotation about x
	t3x3 R = RR.ax2om();
*/
/*
	//COMPARISON AGAINST EXP
	//for comparison with experiment
	t3x3 R = targetori_om;
cout << "Rtarget = " << R << endl;
*/

	//COMPARISON AGAINST EXP WITH INVERSE
	t3x3 Rin = targetori_om;
cout << "Rinput = " << Rin << endl;
	t3x3 R = Rin.inverse();
cout << "Rinput-inverse = " << R << endl;

	//size_t dummy = 0;
	for( auto it = pp3i.begin(); it != pp3i.end(); ++it ) {
		apt_real xx = it->x;
		apt_real yy = it->y;
		apt_real zz = it->z;
		//"active" repositioning of position in fixed coordinate system by premultiply with orientation rotation matrix
/*
		it->x = RzyxR.a11 * xx + RzyxR.a12 * yy + RzyxR.a13 * zz;
		it->y = RzyxR.a21 * xx + RzyxR.a22 * yy + RzyxR.a23 * zz;
		it->z = RzyxR.a31 * xx + RzyxR.a32 * yy + RzyxR.a33 * zz;
*/
//cout << R.a11 << ";" << R.a12 << ";" << R.a13 << "__" << R.a21 << ";" << R.a22 << ";" << R.a23 << "__" << R.a31 << ";" << R.a32 << ";" << R.a33 << "\n";

		it->x = R.a11 * xx + R.a12 * yy + R.a13 * zz;
		it->y = R.a21 * xx + R.a22 * yy + R.a23 * zz;
		it->z = R.a31 * xx + R.a32 * yy + R.a33 * zz;

/*
		dummy++;
		if ( dummy % 1000 ) {
		//if ( it == pp3i.begin() || it == pp3i.begin()+1 || it == pp3i.begin()+2  ) {
cout << "xyz before\t\t" << xx << ";" << yy << ";" << zz << endl;
cout << "xyz rotated\t\t" << it->x << ";" << it->y << ";" << it->z << endl;
		//}
		}
*/
	}

	//identify in which zone the atoms of this rotated atom ball are, if any, and assign corresponding point marks, following the concept of TAPSim to handle marked points
	synstats vbs;

	//random thinning of integration point in zoneII is achieved by random thinning using a PRNG with long period
	mt19937 dice;
	dice.seed( Settings::PRNGSeed );
	//specifically the MersenneTwister algorithm, which should be warmed up
	uniform_real_distribution<apt_real> unifrnd(0.f, 1.f);
	for (size_t i = 0; i < Settings::PRNGDiscard; ++i)
		apt_real discard = unifrnd(dice);

	reporting("MersenneTwister warmed up");

	//now find in which zone and add marks, relocate back just that base plate is at x,y,z=0
	//because we rotated the lattice and not the geometric primitives and zoneII, zoneI, and pillar are centered at origin
	//the offset to relocate all positions is
	p3d offset_z = p3d( 0.f, 0.f, 0.5*(pillar.FH+pillar.SH) ); //the same as for the pillar

	unsigned int dummymark = numeric_limits<unsigned int>::max();
	for (auto it = pp3i.begin(); it != pp3i.end(); ++it) {
		if ( it->z > (zoneII.mybox.zmi + Settings::SXBaseEps) ) { //"normal" points, not forming the base plate
			if ( zoneII.is_inside_tip( *it ) == true ) { //most likely case inside zoneII
				if ( pillar.is_inside_tip( *it ) == true ) {
					vbs.tipatoms++;
					pp3f.push_back( p3dm1(it->x, it->y, it->z + zoneIImid0.z, Settings::TAPSimTipIdx) ); //in-place relocating in z-direction to correct for offset
				}
				else { //potentially in zoneI
					if ( zoneI.is_inside_tip( *it ) == true ) {
						vbs.zoneIatoms++;
						pp3f.push_back( p3dm1(it->x, it->y, it->z + zoneIImid0.z, Settings::TAPSimVacIdx) );
					}
					else { //really in zoneII
						if ( unifrnd(dice) > Settings::PRNGThinning ) { //every eta-th gets through so every 1-eta not, if eta = 0.0, almost every dice result > 0.0 discard, if eta --> 1 almost no dice result is > 1-small number
							vbs.thinnedout++;
							continue;
						}
						else {
							vbs.zoneIIatoms++;
							pp3f.push_back( p3dm1(it->x, it->y, it->z + zoneIImid0.z, Settings::TAPSimVacIdx) );
						}
					}
				}
			}
			//else outside zoneII not considered, MeshGen employs own meshing of zones III and IV
		}
		else { //potentially close to base plate points require special treatment
			if ( zoneII.is_inside_tip( *it ) == true ) { //most likely case inside zoneII
				if ( pillar.is_inside_tip( *it ) == true ) {
					vbs.tipatoms++;
					pp3f.push_back( p3dm1(it->x, it->y, it->z + zoneIImid0.z, Settings::TAPSimTipIdx) );
					//no marking with DUMMYBASE because no rebasing of z coordinates to 0.f for points inside tip
				}
				else { //potentially in zoneI
					if ( zoneI.is_inside_tip( *it ) == true ) {
						vbs.baseatoms++;
						pp3f.push_back( p3dm1(it->x, it->y, it->z + zoneIImid0.z, dummymark) );
					}
					else { //really in zoneII but not thinning, because base plate!
						vbs.baseatoms++;
						pp3f.push_back( p3dm1(it->x, it->y, it->z + zoneIImid0.z, dummymark) );
					}
				}
			}
			//else outside zoneII not considered, MeshGen employs own meshing of zones III and IV
		}
	}

	cout << "pp3f.size() " << pp3f.size() << endl;
	cout << vbs << endl;

	//define base plate geometry and get all dimensions into nanometer
	for( auto it = pp3f.begin(); it != pp3f.end(); ++it ) {
		//is not of base plate? only SI scaling and continue
		//is of base plate define, SI scaling and continue
		apt_real xx = NANOMETER2METER(it->x);
		apt_real yy = NANOMETER2METER(it->y);
		apt_real zz = NANOMETER2METER(it->z);
		if ( it->m != dummymark ) { //no rebasing
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

	//##MK::zoneII, zoneI, and pillar aabb3ds not relocated!
}



/*
void singlecrystal::define_domain()
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


void singlecrystal::load_grainboundary()
{
	gAB.read_triangulation();
}


void singlecrystal::generate_bicrystal()
{
	//##MK::build lattice of crystal gA with atoms in front of the boundary gAB


	//build lattice of crystal gB with atoms behind the boundary gAB
}


void singlecrystal::carve_tip()
{
}

*/
