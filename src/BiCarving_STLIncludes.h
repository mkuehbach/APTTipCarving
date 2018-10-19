/*
	Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"sseldorf
	Markus K\"uhbach, 08/2018
*/

#ifndef __BICARVING_STLINCLUDES_H__
#define __BICARVING_STLINCLUDES_H__

//C++ STL
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <stdexcept>
#include <string>
#include <iomanip>
#include <limits>
#include <vector>
#include <algorithm>
#include <cmath>
#include <list>
#include <cassert>
#include <map>
#include <iterator>
#include <utility>
#include <random>
#include <set>
#include <assert.h>
#include <stdint.h>

//forward declaration for global scope
using namespace std;

//generic global functions to report state, warnings, and erros
void reporting( const string what );
void reporting( const int rank, const string what );
void complaining( const string what );
void stopping( const string what );


#endif
