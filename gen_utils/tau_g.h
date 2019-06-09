#include "stdafx.h"

using namespace std;

class tau_g_line
{
public:
	double tau_ext;
	double sG;
	double sGn;
	double av_size;
	double av_size_n;
	int random_times;

	tau_g_line(double tau_ext, double sG, double sGn, double av_size, double av_size_n, int random_times);
	tau_g_line();
	tau_g_line(string line);
	tau_g_line(istream& i);

	tau_g_line& extend(const tau_g_line& extender);
	double DsG() const;
	double tDsG(double springConstant) const; // adding external stress / Young modulus * sys size **2
	double DG() const;
	double get_tau_ext() const;

	tau_g_line& operator+=(const tau_g_line& addum);
	tau_g_line operator+(const tau_g_line& addum) const;
	tau_g_line& operator/=(double divisor);
	tau_g_line operator/(double divisor) const;
	tau_g_line& operator*=(double factor);
	tau_g_line operator*(double factor) const;
};

ostream& operator<<(ostream& o, const tau_g_line& line);
istream& operator>>(istream& i, tau_g_line& line);

class tau_g_line_pos : public tau_g_line
{
public:
	tau_g_line_pos(int col); // : tau_g_line() {};
	tau_g_line_pos(double tau_ext, double sG, double sGn, double av_size, double av_size_n, int startPosx, int startPosy, int random_times);
	int col;
	double posx;
	double posy;
};

istream& operator>>(istream& i, tau_g_line_pos& line);
ostream& operator<<(ostream& o, const tau_g_line_pos& line);