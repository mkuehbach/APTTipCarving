/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/


#ifndef __BICARVING_DATATYPES_H__
#define __BICARVING_DATATYPES_H__

#include "BiCarving_Settings.h"


struct p3d
{
	apt_real x;
	apt_real y;
	apt_real z;

	p3d() : x(0.f), y(0.f), z(0.f) {}
	p3d(const apt_real _x, const apt_real _y, const apt_real _z) :
		x(_x), y(_y), z(_z) {}
};

ostream& operator<<(ostream& in, p3d const & val);


struct p3dm1
{
	apt_real x;
	apt_real y;
	apt_real z;
	unsigned int m;

	p3dm1() : x(0.f), y(0.f), z(0.f), m(UNKNOWNTYPE) {}
	p3dm1(const apt_real _x, const apt_real _y, const apt_real _z, const unsigned int _m ) :
		x(_x), y(_y), z(_z), m(_m) {}
};

ostream& operator<<(ostream& in, p3dm1 const & val);


struct v3d
{
	apt_real u;
	apt_real v;
	apt_real w;
	apt_real SQR_len;
	v3d() : u(0.f), v(0.f), w(0.f) {} //, SQR_len(0.f) {}
	v3d( const apt_real _u, const apt_real _v, const apt_real _w ) :
		u(_u), v(_v), w(_w) {} //, SQR_len( SQR(_u)+SQR(_v)+SQR(_w) ) {}

	//inline apt_real len() const;
	void normalize();
	void orientnormal( v3d const & reference );
};

ostream& operator<<(ostream& in, v3d const & val);


inline bool SortSQRLenAscending( const v3d &aa1, const v3d &aa2)
{
	return aa1.SQR_len < aa2.SQR_len;
}


struct tri3d
{
	apt_real x1;
	apt_real y1;
	apt_real z1;

	apt_real x2;
	apt_real y2;
	apt_real z2;

	apt_real x3;
	apt_real y3;
	apt_real z3;
	tri3d() : x1(0.f), y1(0.f), z1(0.f), x2(0.f), y2(0.f), z2(0.f), x3(0.f), y3(0.f), z3(0.f) {}
	tri3d( const apt_real _x1, const apt_real _y1, const apt_real _z1,
		const apt_real _x2, const apt_real _y2, const apt_real _z2,
		const apt_real _x3, const apt_real _y3, const apt_real _z3 ) :
		x1(_x1), y1(_y1), z1(_z1),
		x2(_x2), y2(_y2), z2(_z2),
		x3(_x3), y3(_y3), z3(_z3) {}
	p3d center();
	inline bool in_front_of_me( p3d const & p ) const;
};

ostream& operator<<(ostream& in, tri3d const & val);


struct plane3d
{
	p3d center;
	v3d ounormal;
	plane3d() : center(p3d()), ounormal(v3d()) {}
	plane3d(const p3d _c, const v3d _oun) :
		center(_c), ounormal(_oun) {}
};


struct triref3d
{
	size_t v1;
	size_t v2;
	size_t v3;
	triref3d() : v1(0), v2(0), v3(0) {}
	triref3d( const size_t _v1, const size_t _v2, const size_t _v3 ) :
		v1(_v1), v2(_v2), v3(_v3) {}
};

ostream& operator<<(ostream& in, triref3d const & val);


struct aabb3d
{
	apt_real xmi;
	apt_real xmx;
	apt_real ymi;
	apt_real ymx;
	apt_real zmi;
	apt_real zmx;
	apt_real xsz;
	apt_real ysz;
	apt_real zsz;
	aabb3d() :
		xmi(F32MX), xmx(F32MI), ymi(F32MX), ymx(F32MI), zmi(F32MX), zmx(F32MI), xsz(0.f), ysz(0.f), zsz(0.f)  {}
	aabb3d(const apt_real _xmi, const apt_real _xmx, const apt_real _ymi, const apt_real _ymx, const apt_real _zmi, const apt_real _zmx) :
		xmi(_xmi), xmx(_xmx), ymi(_ymi), ymx(_ymx), zmi(_zmi), zmx(_zmx), xsz(_xmx-_xmi), ysz(_ymx-_ymi), zsz(_zmx-_zmi) {}

	void scale();
	void blowup( const apt_real f );
	apt_real diag();
	bool is_inside_box_xy( aabb3d const & reference, apt_real guard );
	p3d mid();
};

ostream& operator<<(ostream& in, aabb3d const & val);



struct sqb
{
	//MK::add reject if binning is too fine thereby exceeding UINT32
	size_t nx;
	size_t ny;
	size_t nz;

	size_t nxy;
	size_t nxyz;
	
	apt_real width;
	aabb3d box;

	sqb() : nx(1), ny(1), nz(1), nxy(1), nxyz(1), width(F32MX), box(aabb3d()) {}
	sqb(const size_t _nx, const size_t _ny, const size_t _nz, const apt_real _w, const aabb3d _bx) :
		nx(_nx), ny(_ny), nz(_nz), nxy(_nx*_ny), nxyz(_nx*_ny*_nz), width(_w), box(_bx) {}

	size_t where( const p3dm1 p );
};

ostream& operator<<(ostream& in, sqb const & val);


struct synstats
{
	apt_int tipatoms;
	apt_int zoneIatoms;
	apt_int zoneIIatoms;
	apt_int baseatoms;
	apt_int thinnedout;
	synstats() : tipatoms(0), zoneIatoms(0), zoneIIatoms(0),
			baseatoms(0), thinnedout(0) {}
	synstats( const apt_int _t, const apt_int _zI, const apt_int _zII, const apt_int _ba,
			const apt_int _to) : tipatoms(_t), zoneIatoms(_zI), zoneIIatoms(_zII),
					baseatoms(_ba), thinnedout(_to) {}
};

ostream& operator<<(ostream& in, synstats const & val);


class cylinder
{
public:
	cylinder();
	cylinder(const apt_real _h, const apt_real _r );
	bool is_inside_cylinder( p3d const & p );
	string report_cylinder();

	apt_real H;			//height
	apt_real R;			//radius

	aabb3d mybox;
	p3d center;
};


class tip
{
public:
	tip();
	tip(const apt_real _fh, const apt_real _rb, const apt_real _rt );
	bool is_inside_tip( p3d const & p );
	void relocate_z( const apt_real offsetz );
	string report_tip();

	apt_real FH;		//conical frustum height
	apt_real RB; 		//conical frustum bottom radius
	apt_real RT;		//conical frustum top radius
	apt_real SH;		//spherical cap on top radius

	aabb3d mybox;
	p3d center;			//center of geometrical primitive spherical cap above FH, axial symmetry of setup, aligned primitives to z-axis
};


#endif
