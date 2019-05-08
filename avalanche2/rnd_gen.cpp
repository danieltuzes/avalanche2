


#include "stdafx.h"
#include "rnd_gen.h"

using namespace std;

ostream& operator<< (ostream& o, distribution distr)
{
	if (distr==distribution::exponential)
		o << "exponential";
	else if (distr == distribution::uniform)
		o << "uniform";
	else if (distr == distribution::uniform_int)
		o << "uniform_int";
	else if (distr==distribution::normal)
		o << "normal";
	else if (distr==distribution::kloster)
		o << "kloster";
	else if (distr==distribution::poisson)
		o << "poisson";
	else if (distr == distribution::weibull)
		o << "weibull";
	else if (distr == distribution::weibullAB)
		o << "weibullAB";
	else if (distr == distribution::lognormal)
		o << "lognormal";
	else
		cerr << "Error: ostream& operator<< not defined on the distribution provided." << endl;

	return o;
}

istream& operator>> (istream& i, distribution& distr)
{
	string name;
	i >> name;
	if (name == "exponential")
		distr = distribution::exponential;
	else if (name == "uniform")
		distr = distribution::uniform;
	else if (name == "uniform_int")
		distr = distribution::uniform_int;
	else if (name == "normal")
		distr = distribution::normal;
	else if (name == "kloster")
		distr = distribution::kloster;
	else if (name == "poisson")
		distr = distribution::poisson;
	else if (name == "weibull")
		distr = distribution::weibull;
	else if (name == "weibullAB")
		distr = distribution::weibullAB;
	else if (name == "lognormal")
		distr = distribution::lognormal;
	else
	{
		cerr << "Error: distribution " << name << " is not supported" << endl;
		exit(1);
	}
	return i;
}