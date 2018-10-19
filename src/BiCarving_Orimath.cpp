/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#include "BiCarving_Orimath.h"

//all conventions following D. J. Rowenhorst et al., MSMSE, Tutorial on Consistent Rotation conventions


t3x3 axisangle::ax2om()
{
	apt_real norm = sqrt(SQR(this->v1)+SQR(this->v2)+SQR(this->v3));
	apt_real n1 = this->v1 / norm;
	apt_real n2 = this->v2 / norm;
	apt_real n3 = this->v3 / norm;

	apt_real c = cos(this->theta / 180.0 * PI);
    apt_real s = sin(this->theta / 180.0 * PI);
    apt_real P = -1.f; //sign convention discriminator
    apt_real omc = 1.f - c;

    return t3x3( 	c + omc*SQR(n1),		omc*n1*n2 + s*n3,		omc*n1*n3 - s*n2,
    				omc*n1*n2 - s*n3,		c + omc*SQR(n2),		omc*n2*n3 + s*n1,
					omc*n1*n3 + s*n2,		omc*n2*n3 - s*n1, 		c + omc*SQR(n3) 	);
}



