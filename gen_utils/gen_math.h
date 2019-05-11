// gen_fmanip.h : stores user defined general function relatad to file and data_block manipulation
//
#define GEN_MATH_VERSION 0.01

// GEN_MATH_VERSION 0.01: linvec has member printAsCol

#pragma once
#include "stdafx.h"
using namespace std;

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286L
#define SQRT2 1.414213562373095048801688724209698078569671875376948073176679737990732478L
#define DINFTY numeric_limits<double>::infinity()
#define SMALLEST_NUMBER nextafter(0.,DINFTY)

int getPower(int number); //number must be a power of two
int lb(int lbof);
int pow2(int exponent);

int getLinPos(int x, int y, int power);
int getX(int linPos, int power);
int getY(int linPos, int mask);

int nofdigit(int number);

template <typename T = double> class linvec //vector with defined operators
{
protected:
	bool mprintAsCol = false;
public:
	vector<T> element;
	linvec();
	linvec(vector<T> element);

	bool consist(const T& compare) const;

	linvec<T>& operator+=(const linvec& toadd);

	linvec<T> operator+(const linvec& toadd) const;

	linvec<T>& operator/=(double divisor);

	linvec<T> operator/(double divisor) const;

	template <typename U> friend ostream& operator<<(ostream& o, const linvec<U>& myvec);

	linvec<T>& printAsCol(bool col);

	bool printAsCol() const;
};

int fint(double value); //gives the biggest integer that is smaller than value

double dist(const vector<double>& a, const vector<double>& b, int sysSize); //the Euclidean distance with periodic boundary conditions, http://en.wikipedia.org/wiki/Euclidean_distance

int dist(int x1, int x2, int sysSize);
double dist(double x1, double x2, int sysSize);

double dist(double x1, double y1, double x2, double y2, int sysSize);

int round_int(double d);