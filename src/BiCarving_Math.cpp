/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#include "BiCarving_Math.h"


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




/* 
	Distance Between Point and Triangle in 3D
	David Eberly
	Geometric Tools, LLC
	http://www.geometrictools.com/
	Copyright 
	c 1998-2016. All Rights Reserved.
	Created: September 28, 1999
	Last Modified: March 1, 2008
	vector3 closesPointOnTriangle( const vector3 *triangle, const vector3 &sourcePosition )
*/


apt_real mathHdl::closestPointOnTriangle( const tri3d face, const p3d src )
{
	v3d edge0( face.x2-face.x1, face.y2-face.y1, face.z2-face.z1 );
	v3d edge1( face.x3-face.x1, face.y3-face.y1, face.z3-face.z1 );
	v3d v0( face.x1-src.x, face.y1-src.y, face.z1-src.z );

	apt_real a = SQR(edge0.u) + SQR(edge0.v) + SQR(edge0.w); //edge0.dot( edge0 );
	apt_real b = edge0.u*edge1.u + edge0.v*edge1.v + edge0.w*edge1.w; //edge0.dot( edge1 );
	apt_real c = SQR(edge1.u) + SQR(edge1.v) + SQR(edge1.w); //edge1.dot( edge1 );
	apt_real d = edge0.u*v0.u + edge0.v*v0.v + edge0.w*v0.w; //edge0.dot( v0 );
	apt_real e = edge1.u*v0.u + edge1.v*v0.v + edge1.w*v0.w; //edge1.dot( v0 );

	apt_real det = a*c - b*b;
	apt_real s = b*e - c*d;
	apt_real t = b*d - a*e;

	if ( s + t < det ) {
		if ( s < 0.f ) {
			if ( t < 0.f ) { //region 4
				if ( d < 0.f ) {
					s = CLAMP( -d/a, 0.f, 1.f );
					t = 0.f;
				}
				else {
					s = 0.f;
					t = CLAMP( -e/c, 0.f, 1.f );
				}
			}
			else {//region 3
				s = 0.f;
				t = CLAMP( -e/c, 0.f, 1.f );
			}
		}
		else if ( t < 0.f ) { //region 5
			s = CLAMP( -d/a, 0.f, 1.f );
			t = 0.f;
		}
		else { //region 0
			apt_real invDet = 1.f / det;
			s *= invDet;
			t *= invDet;
		}
	}
	else {
		if ( s < 0.f ) { //region 2
			apt_real tmp0 = b+d;
			apt_real tmp1 = c+e;
			if ( tmp1 > tmp0 ) {
				apt_real numer = tmp1 - tmp0;
				apt_real denom = a-2*b+c;
				s = CLAMP( numer/denom, 0.f, 1.f );
				t = 1-s;
			}
			else {
				t = CLAMP( -e/c, 0.f, 1.f );
				s = 0.f;
			}
		}
		else if ( t < 0.f ) { //region 6
			if ( a+d > b+e ) {
				apt_real numer = c+e-b-d;
				apt_real denom = a-2*b+c;
				s = CLAMP( numer/denom, 0.f, 1.f );
				t = 1-s;
			}
			else {
				s = CLAMP( -e/c, 0.f, 1.f );
				t = 0.f;
			}
		}
		else { //region 1
			apt_real numer = c+e-b-d;
			apt_real denom = a-2*b+c;
			s = CLAMP( numer/denom, 0.f, 1.f );
			t = 1.f - s;
		}
	}

	//closest point is cp
	p3d cp( face.x1 + s*edge0.u + t*edge1.u, face.y1 + s*edge0.v + t*edge1.v, face.z1 + s*edge0.w + t*edge1.w );

	//so return SQR of difference to avoid sqrt in the code
	return ( SQR(cp.x-src.x) + SQR(cp.y-src.y) + SQR(cp.z-src.z) );
}
