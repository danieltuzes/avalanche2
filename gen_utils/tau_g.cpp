#include "stdafx.h"
#include "tau_g.h"

using namespace std;


tau_g_line::tau_g_line(double tau_ext, double sG, double sGn, double av_size, double av_size_n, int random_times)
	: tau_ext(tau_ext), sG(sG), sGn(sGn), av_size(av_size), av_size_n(av_size_n), random_times(random_times) {}

tau_g_line::tau_g_line() : tau_g_line(0, 0, 0, 0, 0, 0) {}

tau_g_line::tau_g_line(istream& s)
{
	s >> *this;
}

tau_g_line::tau_g_line(string line)
{
	stringstream(line) >> *this;
}


tau_g_line& tau_g_line::extend(const tau_g_line& extender)
{
	tau_ext = max(tau_ext, extender.tau_ext);
	av_size = max(av_size, extender.av_size);
	av_size_n = max(av_size_n, extender.av_size_n);
	sG = max(sG, extender.sG);
	sGn = max(sGn, extender.sGn);
	random_times = max(random_times, extender.random_times);

	return *this;
}

double tau_g_line::DsG() const
{
	return sG - sGn;
}

double tau_g_line::tDsG(double springConstant) const // adding external stress / Young modulus * sys size **2
{
	return (sG - sGn) + get_tau_ext() / springConstant;
}

double tau_g_line::DG() const
{
	return av_size - av_size_n;
}

double tau_g_line::get_tau_ext() const
{
	return tau_ext;
}


tau_g_line& tau_g_line::operator+=(const tau_g_line& addum)
{
	tau_ext += addum.tau_ext;
	sG += addum.sG;
	sGn += addum.sGn;
	av_size += addum.av_size;
	av_size_n += addum.av_size_n;
	random_times += addum.random_times;
	return *this;
}

tau_g_line tau_g_line::operator+(const tau_g_line& addum) const
{
	tau_g_line ret(*this);
	return ret += addum;
}

tau_g_line& tau_g_line::operator*=(double factor)
{
	sG *= factor;
	sGn *= factor;
	av_size *= factor;
	av_size_n *= factor;
	return *this;
}

tau_g_line tau_g_line::operator*(double factor) const
{
	tau_g_line ret(*this);
	return ret *= factor;
}

tau_g_line& tau_g_line::operator/=(double divisor)
{
	return *this *= (1 / divisor);
}

tau_g_line tau_g_line::operator/(double divisor) const
{
	tau_g_line ret(*this);
	return ret /= divisor;
}




ostream& operator<< (ostream& o, const tau_g_line& line)
{
	o << line.tau_ext << "\t"
		<< line.sG << "\t"
		<< line.sGn << "\t"
		<< line.av_size << "\t"
		<< line.av_size_n << "\t"
		<< line.random_times;
	return o;
}

istream& operator>> (istream& i, tau_g_line& line)
{
	i >> line.tau_ext
		>> line.sG
		>> line.sGn
		>> line.av_size
		>> line.av_size_n
		>> line.random_times;
	return i;
}

tau_g_line_pos::tau_g_line_pos(int col) : col(col)
{
	if (col < 6)
	{
		cout << "The position of the column containing the position of the avalanche supposed to be at the least the 6th column.";
		exit(1);
	}

}

tau_g_line_pos::tau_g_line_pos(double tau_ext, double sG, double sGn, double av_size, double av_size_n, int random_times, int startPosx, int startPosy) :
tau_g_line(tau_ext, sG, sGn, av_size, av_size_n, random_times), posx(startPosx), posy(startPosy) {}

istream& operator>> (istream& i, tau_g_line_pos& line)
{
	i >> static_cast<tau_g_line &>(line);
	for (int colindex = 6; colindex != line.col; ++colindex)
		i.ignore(numeric_limits<streamsize>::max(), '\t');
	i >> line.posx
		>> line.posy;
	return i;
}

ostream& operator<< (ostream& o, const tau_g_line_pos& line)
{
	o << static_cast<tau_g_line &>(const_cast<tau_g_line_pos &>(line)) << "\t"
		<< line.posx << "\t"
		<< line.posy;
	return o;
}