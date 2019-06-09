// rnd_gen.h : random number generation
//
#define RND_GEN_VERSION 0.04



#include "stdafx.h"
#pragma once

using namespace std;

enum distribution
{
	exponential,
	uniform,
	uniform_int,
	normal,
	kloster,
	poisson,
	weibull,
	weibullAB,
	lognormal
};

ostream& operator<< (ostream& o, distribution distr);

istream& operator>> (istream& i, distribution& distr);


class universal_rnd
{
protected:
	distribution distr;
	string params;
	double multiA; // the random number is multiplicated with this number, it effectively multiplies the average with multiA
	double multiB; // the random number is multiplicated with this number, it effectively multiplies the average with multiB, with multiA can be used to create "double Weibull" distribution
	double weight; // the weight for multiA so that multiA and multiB doesn't modify the average
	std::mt19937 generator;
	int seed;
	int times_advanced;
	bool good;

	normal_distribution<double> norm_distr;
	exponential_distribution<double> exp_distr;
	weibull_distribution<double> weib_distr;
	uniform_real_distribution<double> unif_distr;
	lognormal_distribution<double> logn_distr;

	double rand_gen() //todo: ezt a részt talán ifek nélkül is meg lehetne oldani valahogy virtuális függvényekkel és osztályszármaztatással
	{
		if (distr==distribution::normal)
			return norm_distr(generator);
		else if (distr==distribution::exponential)
			return exp_distr(generator);
		else if (distr==distribution::weibull)
			return weib_distr(generator);
		else if (distr == distribution::weibullAB)
		{
			double x = weib_distr(generator);
			double ret = weib_distr(generator);
			if (x < pow(-log(1 - weight), 1 / weib_distr.a()) * weib_distr.b())
				return ret * multiA;
			else
				return ret* multiB;
		}
		else if (distr == distribution::uniform)
			return unif_distr(generator);
		else if (distr == distribution::lognormal)
			return logn_distr(generator);
		else
		{
			cerr << "random number generator " << distr << " not supported" << endl;
			return 0;
		}
	}

	bool init_distr() //uses distr and param
	{
		string params_split(params);
		replace(params_split.begin(),params_split.end(),'_',' ');
		stringstream params_stream(params_split);

		if (distr==distribution::normal)
		{
			double mean;
			double sigma;
			params_stream >> mean >> sigma;
			normal_distribution<double> tmp(mean,sigma); // for some old compilers
			norm_distr.param(tmp.param());
		}
		
		else if (distr==distribution::exponential)
		{
			double lambda;
			params_stream >> lambda;
			exponential_distribution<double> tmp(lambda);
			exp_distr.param(tmp.param());
		}

		else if (distr == distribution::uniform)
		{
			double a;
			double b;
			params_stream >> a >> b;
			uniform_real_distribution<double> tmp(a, b);
			unif_distr.param(tmp.param());
		}

		else if (distr==distribution::weibull)
		{
			double a;
			double b;
			params_stream >> a >> b;
			weibull_distribution<double> tmp(a,b);
			weib_distr.param(tmp.param());
		}

		else if (distr == distribution::weibullAB)
		{
			double a;
			double b;
			params_stream >> a >> b >> multiA >> multiB;
			weibull_distribution<double> tmp(a, b);
			weib_distr.param(tmp.param());
			weight = (1 - multiB) / (multiA - multiB);
		}

		else if (distr == distribution::lognormal)
		{
			double m;
			double s;
			params_stream >> m >> s;
			lognormal_distribution<double> tmp(m, s);
			logn_distr.param(tmp.param());
		}
		else
		{
			cerr << "Cannot create random number generator with distribution " << distr
				 << " (with parameters " << params << ") "
				 << " and seed " << seed;
			return false;
		}
		return true;
	}


public:
	universal_rnd(distribution type, string params, int seed) : distr(type), params(params), multiA(0), multiB(0), weight(1), seed(seed), times_advanced(0), good(false)
	{
		generator.seed(seed);
		if (init_distr())
			good = true;
	}

	mt19937& getEngine(int increaseCounter = 1)
	{
		times_advanced += increaseCounter;
		return generator;
	}

	universal_rnd(string type_with_params, int seed) : multiA(0), multiB(0), weight(1), seed(seed), times_advanced(0), good(false) //constructor delegating is not allowed in vs2012
	{
		string type_name(type_with_params);
		replace(type_name.begin(),type_name.end(),'_',' ');
		stringstream type_name_stream(type_name);
		type_name_stream >> distr;

		stringstream param_stream(type_with_params);
		param_stream.ignore(numeric_limits<streamsize>::max(),'_');
		param_stream >> params;

		generator.seed(seed);
		if (init_distr())
			good = true;
	}

	double operator()()
	{
		++times_advanced;
		return rand_gen();
	}

	void discard(int plus_times)
	{
		times_advanced += plus_times;
		generator.discard(plus_times);
	}

	int getSeed() const
		{return seed;}

	distribution getType() const
		{return distr;}

	string getParams() const
		{return params;}

	int getTimesAdvanced() const
		{return times_advanced;};

	void increase_times_advanced(int number_of_occasion)
	    {times_advanced += number_of_occasion;}

	bool isGood() const
		{return good;};

	double getPar(int which) const
	{
		string params_split(params);
		replace(params_split.begin(), params_split.end(), '_', ' ');
		stringstream params_stream(params_split);
		double param = 0;
		for (int i = 0; i < which; ++i)
			params_stream >> param;
		
		return param;
	}
};
