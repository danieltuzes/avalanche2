
#include "stdafx.h"
#include "gen_math.h"

using namespace std;


int getPower(int number) //number must be a power of two
{
	return lb(number);
}

int lb(int number)
{
	int power = 0;
	while (number >>= 1)
		++power;
	return power;
}

int pow2(int exponent)
{
	return (1 << exponent);
}

int getLinPos(int x, int y, int power) { return (x << power) + y; }

int getX(int linPos, int power) { return (linPos >> power); }

int getY(int linPos, int mask) { return (linPos & mask); }

int nofdigit(int number)
{
	int digit = 0;
	while (number != 0) { number /= 10; ++digit; }
	return digit;
}

int fint(double value)
{
	if (value < 0)
		return (static_cast<int>(value)-1);
	else
		return static_cast<int>(value);
}

template <typename T> linvec<T>::linvec() {}
template <typename T> linvec<T>::linvec(vector<T> element) : element(element) {}

template <typename T> bool linvec<T>::consist(const T& compare) const
{
	for (auto x : element)
		if (x == compare)
			return true;

	return false;
}

template <typename T> linvec<T>& linvec<T>::operator+=(const linvec<T>& toadd)
{
	transform(element.begin(),element.end(),toadd.element.begin(),element.begin(),[](T a, T b){return a+b;});
	return *this;
}

template <typename T> linvec<T> linvec<T>::operator+(const linvec<T>& toadd) const
{
	return linvec<T>(*this)+=toadd;
}

template <typename T> linvec<T>& linvec<T>::operator/=(double divisor)
{
	transform(element.begin(),element.end(),element.begin(),[divisor](T elemval){return elemval/static_cast<T>(divisor);});
	return *this;
}

template <typename T> linvec<T> linvec<T>::operator/(double divisor) const
{
	return linvec<T>(*this)/=divisor;
}

template <typename T> ostream& operator<<(ostream& o, const linvec<T>& myvec)
{
	for (auto it = myvec.element.cbegin(); it != prev(myvec.element.cend()); ++it)
	{
		if (!myvec.printAsCol())
			o << *it << "\t";
		else
			o << *it << "\n";

	}
	o << myvec.element.back();
	return o;
}

template <typename T> linvec<T>& linvec<T>::printAsCol(bool col)
{
	mprintAsCol = col;
	return *this;
}

template <typename T> bool linvec<T>::printAsCol() const
{
	return mprintAsCol;
}

template class linvec<int>;
template ostream& operator<<(ostream& o, const linvec<int>& myvec);
template class linvec<double>;
template ostream& operator<<(ostream& o, const linvec<double>& myvec);

double dist(const vector<double>& a, const vector<double>& b, int sysSize) //the Euclidean distance with periodic boundary conditions, http://en.wikipedia.org/wiki/Euclidean_distance
{
	double ret = sqrt(inner_product(
		a.begin(),
		a.end(),
		b.begin(),
		0.,
		std::plus<double>(),
		[sysSize](double ai, double bi){return dist(ai, bi, sysSize)*dist(ai, bi, sysSize); }));
	return ret;
}

int dist(int x1, int x2, int sysSize) // periodic distance of x1 and x2
{
	return min(abs(x1 - x2), sysSize - abs(x1 - x2));
}

double dist(double x1, double x2, int sysSize) // periodic distance of x1 and x2
{
	return min(abs(x1 - x2), sysSize - abs(x1 - x2));
}

double dist(double x1, double y1, double x2, double y2, int sysSize) // periodic 2D distance of (x1,y1) and (x2,y2)
{
	double distx = dist(x1, x2, sysSize);
	double disty = dist(y1, y2, sysSize);

	return sqrt(distx*distx + disty*disty);
}

int round_int(double d)
{
	return static_cast<int>(d + 0.5);
}