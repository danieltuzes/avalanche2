// corr_utils.h : stores user defined utilities for corr.cpp
//

#define CORR_UTILS_VERSION 4.2
#define CORR_VERSION 2.3

// What's new?

// CORR_VERSION 2.3: fillUpperHalf applied to corr_type::F output map
// CORR_VERSION 2.2: corr_type D calculates the ratio as well, prints autoheader, treats nans as 1
// CORR_VERSION 2.1: corr_type D is implemented
// CORR_VERSION 2:

// CORR_UTILS_VERSION 4.2: fillUpperHalf is introduced
// CORR_UTILS_VERSION 4.1: corr_type::D is implemented
// CORR_UTILS_VERSION 4:

#pragma once
#include "stdafx.h"
#include "../gen_utils/gen_siomanip.h"
#include "../gen_utils/gen_fmanip.h"
#include "../gen_utils/gen_version.h"
#include "../gen_utils/gen_param.h"
#include "../gen_utils/gen_fftw3.h"
#include "../gen_utils/FT_block_data.h"

#pragma region corr_type IO
enum corr_type
{
	X, //cross-correlation
	A, //auto-correlation
	F,  //frequencies, fourier space (abs square)
	D // frequency differences, files must be given in pairs
	
};

istream& operator >> (istream& i, corr_type& c)
{
	string s;
	i >> s;
	if (s == "X")
		c = corr_type::X;
	else if (s == "A")
		c = corr_type::A;
	else if (s == "F")
		c = corr_type::F;
	else if (s == "D")
		c = corr_type::D;
	else
	{
		cerr << "Unknown corr_type " << s << " Program terminates." << endl;
	}
	return i;
}

ostream& operator << (ostream& o, corr_type c)
{
	if (c == corr_type::X)
		o << "X";
	else if (c == corr_type::A)
		o << "A";
	else if (c == corr_type::F)
		o << "F";
	else if (c == corr_type::D)
		o << "D";
	else
		cerr << "Unsupported corr_type" << endl;

	return o;
}

template <> string param<corr_type>::getValType() { return "corr_type"; }

#pragma endregion

void intro() {
	stringstream version;
	version << "VERSION:\n"
		<< "    " << "corr:          " << CORR_VERSION << "\n"
		<< "    " << "corr_utils:    " << XSTR(CORR_UTILS_VERSION) << "\n"
		<< "    " << "FT_block_data: " << XSTR(FT_BLOCK_DATA_VERSION) << "\n"
		<< "    " << "gen_siomanip:  " << XSTR(GEN_SIOMANIP_VERSION) << "\n"
		<< "    " << "gen_fmanip:    " << XSTR(GEN_FMANIP_VERSION) << "\n"
		<< "    " << "gen_math:      " << XSTR(GEN_MATH_VERSION) << "\n"
		<< "    " << "gen_version:   " << XSTR(GEN_VERSION_VERSION) << "\n";

	intro("This is a correlation-calculation program for avalanche. ", version.str());
}

void fillUpperHalf(block_data<double>& matrix)
{
	int size = matrix.getSize();
	for (int x = 0; x < size; ++x)
	{
		for (int y = size / 2 + 1; y < size; ++y)
		{
			matrix(x, y) = matrix(x, size - y);
		}
	}
}
