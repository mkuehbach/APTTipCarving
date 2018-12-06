/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/



#include "BiCarving_Datatypes.h"



ostream& operator<<(ostream& in, p3d const & val)
{
	in << val.x << ";" << val.y << ";" << val.z << endl;
	return in;
}


ostream& operator<<(ostream& in, p3dm1 const & val)
{
	in << val.x << ";" << val.y << ";" << val.z << "--" << val.m << endl;
	return in;
}

/*
inline apt_real v3d::len() const
{
	return (this->SQR_len > EPSILON) ? sqrt(this->SQR_len) : 0.f;
}
*/


t3x3 t3x3::premultiplyR1( t3x3 const & R1)
{
	return t3x3(	R1.a11*this->a11 + R1.a12*this->a21 + R1.a13*this->a31,
					R1.a11*this->a12 + R1.a12*this->a22 + R1.a13*this->a32,
					R1.a11*this->a13 + R1.a12*this->a23 + R1.a13*this->a33,

					R1.a21*this->a11 + R1.a22*this->a21 + R1.a23*this->a31,
					R1.a21*this->a12 + R1.a22*this->a22 + R1.a23*this->a32,
					R1.a21*this->a13 + R1.a22*this->a23 + R1.a23*this->a33,

					R1.a31*this->a11 + R1.a32*this->a21 + R1.a33*this->a31,
					R1.a31*this->a12 + R1.a32*this->a22 + R1.a33*this->a32,
					R1.a31*this->a13 + R1.a32*this->a23 + R1.a33*this->a33  );
}


t3x3 t3x3::inverse()
{
	apt_real det = 	1.f /
					( + this->a11*(this->a22*this->a33 - this->a32*this->a23)
					- this->a21*(this->a12*this->a33 - this->a32*this->a13)
					+ this->a31*(this->a12*this->a23 - this->a22*this->a13)    );

	return t3x3(	det*(this->a22*this->a33 - this->a23*this->a32),
					det*(this->a13*this->a32 - this->a12*this->a33),
					det*(this->a12*this->a23 - this->a13*this->a22),

					det*(this->a23*this->a31 - this->a21*this->a33),
					det*(this->a11*this->a33 - this->a13*this->a31),
					det*(this->a13*this->a21 - this->a11*this->a23),

					det*(this->a21*this->a32 - this->a22*this->a31),
					det*(this->a12*this->a31 - this->a11*this->a32),
					det*(this->a11*this->a22 - this->a12*this->a21)     );
}


ostream& operator << (ostream& in, t3x3 const & val) {
	in << val.a11 << ";" << val.a12 << ";" << val.a13 << "\n";
	in << val.a21 << ";" << val.a22 << ";" << val.a23 << "\n";
	in << val.a31 << ";" << val.a32 << ";" << val.a33 << endl;
	return in;
}


p3d p3d::active_rotation_relocate( t3x3 const & R )
{
	return p3d( 	R.a11 * this->x + R.a12 * this->y + R.a13 * this->z,
					R.a21 * this->x + R.a22 * this->y + R.a23 * this->z,
					R.a31 * this->x + R.a32 * this->y + R.a33 * this->z   );
}

/*
apt_real v3d::SQRLen()
{
	return (SQR(this->u)+SQR(this->v)+SQR(this->w));
}
*/


apt_real v3d::dot( v3d const & b)
{
	return (this->u * b.u + this->v * b.v + this->w * b.w);
}


v3d v3d::cross( v3d const & b )
{
	return v3d( this->v * b.w - this->w * b.v,
				this->w * b.u - this->u * b.w,
				this->u * b.v - this->v * b.u );
}

void v3d::normalize()
{
	//##MK::safety
	apt_real SQRLen = sqrt(SQR(this->u) + SQR(this->v) + SQR(this->w));
	this->u /= SQRLen;
	this->v /= SQRLen;
	this->w /= SQRLen;
}


void v3d::orientnormal( v3d const & reference )
{
	apt_real dot = (this->u * reference.u) + (this->v * reference.v) + (this->w * reference.w);
	if ( dot < 0.f ) { //flip normal
		this->u *= -1.f;
		this->v *= -1.f;
		this->w *= -1.f;
	}
}


ostream& operator<<(ostream& in, v3d const & val)
{
	in << val.u << ";" << val.v << ";" << val.w << endl; // "----" << val.SQR_len << endl;
	return in;
}


p3d tri3d::center()
{
	return p3d(
			(this->x1 + this->x2 + this->x3) / 3.f,
			(this->y1 + this->y2 + this->y3) / 3.f,
			(this->z1 + this->z2 + this->z3) / 3.f );
}


inline bool tri3d::in_front_of_me( p3d const & p ) const
{
	//##MK::dot product with [111] if larger epsilon on or in front else behind
/*
 	//##MK::dot product with consistent outer unit normal, if negative behind, if positive in front or on plane
	//p3d discriminator = p3d( 	p.x - this->center.x,
	//							p.y - this->center.y,
	//							p.z - this->center.z );
	apt_real dot = ((p.x - this->center.x) * this->ounormal.u) + ((p.y - this->center.y) * this->ounormal.v) + ((p.z - this->center.z) * this->ounormal.w);
	if ( dot < 0.f )
		return false;
	else
		return true;
 */
	return true;
}


ostream& operator<<(ostream& in, tri3d const & val)
{
	in << val.x1 << ";" << val.y1 << ";" << val.z1 << "\n";
	in << val.x2 << ";" << val.y2 << ";" << val.z2 << "\n";
	in << val.x3 << ";" << val.y3 << ";" << val.z3 << endl;
	return in;
}

ostream& operator<<(ostream& in, triref3d const & val)
{
	in << val.v1 << "\t\t" << val.v2 << "\t\t" << val.v3 << endl;
	return in;
}


void aabb3d::scale()
{
	this->xsz = this->xmx - this->xmi;
	this->ysz = this->ymx - this->ymi;
	this->zsz = this->zmx - this->zmi;
}


void aabb3d::blowup( const apt_real f )
{
	this->xmi -= f;
	this->xmx += f;
	this->ymi -= f;
	this->ymx += f;
	this->zmi -= f;
	this->zmx += f;
	this->scale();
}


apt_real aabb3d::diag()
{
	return sqrt(SQR(this->xmx-this->xmi)+SQR(this->ymx-this->ymi)+SQR(this->zmx-this->zmi));
}


bool aabb3d::is_inside_box_xy( aabb3d const & reference, apt_real guard )
{
	if ( 	(reference.xmi - guard) > this->xmi &&
			(reference.xmx + guard) < this->xmx &&
			(reference.ymi - guard) > this->ymi &&
			(reference.ymx + guard) < this->ymx )
		return true;
	else
		return false;
	//&& (reference.zmi - guard) > this->zmi && (reference.zmx + guard) < this->zmx )
}


p3d aabb3d::mid()
{
	return p3d( 	0.5 * (this->xmi + this->xmx),
					0.5 * (this->ymi + this->ymx),
					0.5 * (this->zmi + this->zmx) );
}


ostream& operator<<(ostream& in, aabb3d const & val)
{
	in << val.xmi << ";" << val.xmx << "--" << val.ymi << ";" << val.ymx << "--" << val.zmi << ";" << val.zmx << "----" << val.xsz << ";" << val.ysz << ";" << val.zsz << endl;
	return in;
}



size_t sqb::where( const p3dm1 p )
{
	//3d implicit x+y*nx+z*nx*ny
	apt_real ix = floor((p.x - this->box.xmi) / (this->box.xmx - this->box.xmi) * static_cast<apt_real>(this->nx));
	apt_real iy = floor((p.y - this->box.ymi) / (this->box.ymx - this->box.ymi) * static_cast<apt_real>(this->ny));
	apt_real iz = floor((p.z - this->box.zmi) / (this->box.zmx - this->box.zmi) * static_cast<apt_real>(this->nz));
	size_t res = static_cast<size_t>(ix) + static_cast<size_t>(iy) * this->nx + static_cast<size_t>(iz) * this->nxy;
	return res;
}


ostream& operator<<(ostream& in, sqb const & val)
{
	in << "BinningNXYZ = " << val.nx << ";" << val.ny << ";" << val.nz << "\t\t" << val.nxy << ";" << val.nxyz << "\t\t" << endl;
	in << "BinningBoundingBox = " << val.box << endl;
	return in;
}


ostream& operator<<(ostream& in, synstats const & val)
{
	in << "Zone0/ZoneI/ZoneII/Base/Thinnedout = " << val.tipatoms << ";" << val.zoneIatoms << ";" << val.zoneIIatoms << ";" << val.baseatoms << ";" << val.thinnedout << endl;
	return in;
}


ostream& operator<<(ostream& in, synstats2 const & val)
{
	in << "Zone0(A/B) = " << val.tipatomsA << ";" << val.tipatomsB << "\n";
	in << "ZoneI(A/B) = " << val.zoneIatomsA << ";" << val.zoneIatomsB << "\n";
	in << "ZoneII(A/B) = " << val.zoneIIatomsA << ";" << val.zoneIIatomsB << "\n";
	in << "Base(A/B) = " << val.baseatomsA << ";" << val.baseatomsB << "\n";
	in << "ThinnedOut(A/B) = "<< val.thinnedoutA << ";" << val.thinnedoutB << endl;
	return in;
}


cylinder::cylinder()
{
	this->H = 0.f;
	this->R = 0.f;

	this->mybox.xmi = 0.f;
	this->mybox.xmx = 0.f;
	this->mybox.ymi = 0.f;
	this->mybox.ymx = 0.f;
	this->mybox.zmi = 0.f;
	this->mybox.zmx = 0.f;
	this->mybox.scale();

	this->center.x = 0.f;
	this->center.y = 0.f;
	this->center.z = 0.f;
}


cylinder::cylinder(const apt_real _h, const apt_real _r )
{
	this->H = _h;
	this->R = _r;

	this->mybox.xmi = -1.f * this->R; //center cylinder at origin by default
	this->mybox.xmx = +1.f * this->R;
	this->mybox.ymi = -1.f * this->R;
	this->mybox.ymx = +1.f * this->R;
	this->mybox.zmi = -0.5 * this->H;
	this->mybox.zmx = +0.5 * this->H;
	this->mybox.scale();

	this->center.x = 0.f;
	this->center.y = 0.f;
	this->center.z = 0.f;

	cout << "Cylinder height/radius = " << this->H << ";" << this->R << endl;
}


bool cylinder::is_inside_cylinder( p3d const & p)
{
	//MK::check if position p is inside cylinder volume defined by mybox

	//outside bounding box
	if ( 	p.x < this->mybox.xmi ||
			p.x > this->mybox.xmx ||
			p.y < this->mybox.ymi ||
			p.y > this->mybox.ymx ||
			p.z < this->mybox.zmi ||
			p.z > this->mybox.zmx )
		return false;

	//cylinder axis is parallel to z, so outside 2D circle?
	if ( SQR(p.x-this->center.x)+SQR(p.y-this->center.y) <= SQR(this->R) )
		return true;
	else
		return false;
}

string cylinder::report_cylinder()
{
	string str = "Cylinder H/R = " + to_string(this->H) + ";" + to_string(this->R) + "\n";
	str += "imi/imx = " + to_string(this->mybox.xmi) + ";" + to_string(this->mybox.xmx) + ";" + to_string(this->mybox.ymi) + ";" + to_string(this->mybox.ymx) + ";" + to_string(this->mybox.zmi) + ";" + to_string(this->mybox.zmx) + "\n";
	return str;
}


tip::tip()
{
	FH = 0.f;
	RB = 0.f;
	RT = 0.f;
	SH = 0.f;
	mybox = aabb3d();
	center = p3d();
}


tip::tip(const apt_real _fh, const apt_real _rb, const apt_real _rt )
{
	FH = _fh;
	RB = _rb;
	RT = _rt;
	SH = _rt;
	mybox = aabb3d( -RB, +RB, -RB, +RB, 0.f, FH + SH );
	mybox.scale();
	center = mybox.mid();
}


bool tip::is_inside_tip( p3d const & p )
{
	//MK::check if position p is inside tip geometry (conical frustum with spherical cap on top)

	//outside bounding box
	if ( 	p.x < this->mybox.xmi ||
			p.x > this->mybox.xmx ||
			p.y < this->mybox.ymi ||
			p.y > this->mybox.ymx ||
			p.z < this->mybox.zmi ||
			p.z > this->mybox.zmx )
		return false;

	if ( (p.z - mybox.zmi) <= FH ) { //in conical frustum section ?
		apt_real Rz = RB - (p.z - mybox.zmi)/FH * (RB - RT);
		if ( (SQR(p.x-this->center.x)+SQR(p.y-this->center.y)) <= SQR(Rz) ) {
			return true; //inside frustum part so inside, for TAPSim simplified no concave section at the bottom for TAPSIM
		}
		else
			return false; //in mybox but not in frustum means outside
	}
	else { //high enough in z to be in halfsphere potentially
		apt_real Hz = (p.z - mybox.zmi - FH);
		if ( Hz <= (RT - EPSILON) ) { //##MK::works because we assume a halfsphere on top, otherwise test against SHhalfsphere on top
			//apt_real Raz = RT * sin(acos(Hz/RT));
			//apt_real Raz = RT * sqrt(1.f - SQR(Hz/RT));
			apt_real SQRRaz = SQR(RT) * (1.f - SQR(Hz/RT));
			if ( (SQR(p.x-this->center.x)+SQR(p.y-this->center.y)) <= SQRRaz )
				return true; //inside halfsphere
			else
				return false;
		}
		else { //too much above
			return false;
		}
	}
}


void tip::relocate_z( const apt_real offsetz )
{
	this->mybox.zmi += offsetz;
	this->mybox.zmx += offsetz;
	this->mybox.scale();
	this->center = this->mybox.mid();
}


string tip::report_tip()
{
	string str = "FH/RB/RT/SH = " + to_string(FH) + ";" + to_string(RB) + ";" + to_string(RT) + ";" + to_string(SH) + "\n";
	str += "imi/imx = " + to_string(mybox.xmi) + ";" + to_string(mybox.xmx) + ";" + to_string(mybox.ymi) + ";" + to_string(mybox.ymx) + ";" + to_string(mybox.zmi) + ";" + to_string(mybox.zmx) + "\n";
	return str;
}


/*
void secondphasemodel::reportParticlesVTK( const unsigned int simid, const int rank )
{
	//MK::write VTK file showing barycenter of all synthesized cluster in reconstructed space
	//includes particles outside actual tip
	string vtk_io_fn = "PARAPROBE.SimID." + to_string(simid) + ".Rank." + to_string(rank) + ".SyntheticTipCluster.vtk";

	ofstream vtk;
	vtk.open(vtk_io_fn.c_str(), ofstream::out | ofstream::trunc);
	if (vtk.is_open() == true) {
		//header
		vtk << "# vtk DataFile Version 2.0\n";
		vtk << "PARAPROBE Synthetic tip cluster placed\n";
		vtk << "ASCII\n";
		vtk << "DATASET POLYDATA\n";
		vtk << "POINTS " << particles.size() << " double\n";
		for( auto it = particles.begin(); it != particles.end(); ++it ) {
			vtk << it->center.x << " " << it->center.y << " " << it->center.z << "\n";
		}
		vtk << "\n";
		vtk << "VERTICES " << particles.size() << " " << 2*particles.size() << "\n";
		for( size_t i = 0; i < particles.size(); ++i ) {
			vtk << 1 << " " << i << "\n";
		}
		vtk << "\n";
		vtk << "POINT_DATA " << particles.size() << "\n";
		vtk << "FIELD FieldData 1\n";
		vtk << "Radius 1 " << particles.size() << " float\n";
		for( auto it = particles.begin(); it != particles.end(); ++it ) {
			vtk << it->radius << "\n";
		}
		vtk.flush();
		vtk.close();
		cout << "VTK file of synthesized particles positions written to file" << endl;
	}
	else {
		cout << "VTK file of synthesized particles positions was not written" << endl;
	}
}
*/
