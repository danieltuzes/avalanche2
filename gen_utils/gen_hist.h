// gen_hist.h : define basic tools for data analysis and histograms
//
#define GEN_HIST_VERSION 0.02

// the member counter in class accum is changed to unsigned long int

#include "stdafx.h"
#include "gen_math.h"
#include "gen_param.h"
#pragma once
using namespace std;


#pragma region scale

enum scale
{
	linear,      //equally distributed
	logarithmic, //equally distributed
	quasi_log   //distributed in such a way to comfortly store integer data, the delimiters are always half-integers, e.g.: 0.5-1.5, 1.5-2.5, 2.5-4.5, ... etc
};

ostream& operator<<(ostream& o, scale s)
{
	if (s==scale::linear)
		o << "linear";
	else if(s==scale::logarithmic)
		o << "logarithmic";
	else
		o << "quasi_log";
	return o;
}

istream& operator>>(istream& i, scale& s)
{
	string ss;
	i >> ss;
	if (ss.compare("linear") == 0)
		s = scale::linear;
	else if (ss.compare("logarithmic") == 0)
		s = scale::logarithmic;
	else if (ss.compare("quasi_log") == 0)
		s = scale::quasi_log;
	return i;
}

template <> string param<scale>::getValType() {return "scale";}

#pragma endregion


template <typename T> class accum // how many times have changed, what is the sum
{
protected:
	unsigned long int counter;
	T data;
public:
	accum(T data):
		counter(0), data(data) {};
	void add(const T& toAdd) // include another point in the histogram
	{
		data += toAdd;
		++counter;
	}
	void add(T&& toAdd)
	{
		this->add(toAdd);
	}
	T intensity() const
	{
		return data;
	}
	unsigned long int getCounter() const
	{
		return counter;
	}
	T averaged() const // returns the averaged data; 
	{
		return (data / counter);
	}
	friend ostream& operator<<(ostream& o, const accum<T>& m_accum)
	{
		o << m_accum.counter << "\t" << m_accum.data;
		return o;
	}
};


template <typename T> class bin : public accum<T>
{
protected:
	double lb;
	double ub;
public:
	bin(double lb, double ub, T data):
		accum<T>(data), lb(lb), ub(ub) {};

	double lowerb() const
	{
		return lb;
	}

	double upperb() const
	{
		return ub;
	}

	linvec<double> range() const // returns a 2 element long linvec with lb and ub
	{
		vector<double> basic_ret;
		basic_ret.push_back(lb);
		basic_ret.push_back(ub);
		return linvec<double>(basic_ret);
	}

	double width() const
	{
		return ub - lb;
	}

	bin& operator+=(const bin<T>& add)
	{
		this->counter += add.counter;
		this->data += add.data;
		return *this;
	}

	friend ostream& operator<<(ostream& o, const bin<T>& mybin)
	{
		o << mybin.lb << "\t" << mybin.ub << "\t" << accum<T>(mybin);
		return o;
	}
};


template <typename T> class hist // stores 1D histograms of values with type T
{
protected:
	double lb; //the lower boundary of the histogram domain
	double incr; // the width of the bins in case of linear scale
	double ub; //the upper boundary of the histogram domain
	int size; // the number of bins
	vector<bin<T>> m_bins; // the array of the bins; a bin contain its lb, ub, how many times it has been modified and the total value
	scale m_scale; // linear, logarithmic or quasi-logarithmic
	double logIncr; // the log of the quotient
	double logLb; // to calculate only once and store log(lb)
	bool printAv; // to print out the average in case the values are streamed out
	bool cumulative; //
	bool is_prop;

public:
	hist() : lb(0), incr(0), ub(0), size(0), m_scale(scale::linear), logIncr(0), logLb(0), printAv(false), cumulative(false), is_prop(false) {};
	hist(double lb, double ub, scale m_scale, int size, T data):
		lb(lb), ub(ub), size(size), m_scale(m_scale), logIncr(0), logLb(0), printAv(false), cumulative(false), is_prop(false)
	{
		if (m_scale == scale::linear)
		{
			incr = (ub-lb) / size;
			for (int i=0; i<size; ++i)
				m_bins.push_back(bin<T>(lb + i * incr, lb + (i + 1) * incr, data));
		}
		else //m_scale == scale::logarithmic || m_scale == scale::quasi_log must hold
		{
			if (lb <= 0)
				throw range_error("Bad lb value for scale::logarithmic");
			if (m_scale == scale::quasi_log)
				ub += 0.5; //otherwise value at the upper border will be out of range

			incr = pow(ub/lb,1./size); // the quotient
			logIncr = log(incr);
			logLb = log(lb);
			if (m_scale == scale::logarithmic)
			{
				for (int i=0; i<size; ++i)
					m_bins.push_back(bin<T>(lb * pow(incr,i), lb * pow(incr,i+1), data));
			}
			else if (m_scale == scale::quasi_log)
			{
				for (int i=0; i<size; ++i)
					m_bins.push_back(bin<T>(static_cast<int>(lb * pow(incr,i)) + 0.5, static_cast<int>(lb * pow(incr,i+1)) + 0.5, data));
			}

		}

	}
	
	hist(double lb, double ub, const vector<T>& data):
		lb(lb), incr((ub - lb) / data.size()), ub(ub), size(data.size()), m_scale(linear), printAv(false), cumulative(false), is_prop(false)
	{
		for (int i = 0; i < size; ++i)
			m_bins.push_back(bin<T>(lb + i*incr, lb + (i + 1)*incr, data[i]));
	}

	friend ostream& operator<< (ostream& o, const hist<T>& m_hist)
	{
		for (bin<T> X : m_hist.m_bins)
		{
			if (!m_hist.printAv)
				o << X << endl;
			else
				o << X << "\t" << X.averaged() << endl;
		}

		return o;
	}
	
	bool fWrite(string ofname,string header, int precision = 17)
	{
		ofstream f(ofname);
		if (!f)
		{
			cerr << "Can't create " << ofname << endl;
			return false;
		}
		f.precision(precision);
		f << header << *this;
		cout << "Writing " << ofname << endl;
		return true;
	}
	
	linvec<double> range()
	{
		vector<double> basic_ret;
		basic_ret.push_back(lb);
		basic_ret.push_back(ub);
		return linvec<double>(basic_ret);
	}
	
	hist& cumul()
	{
		cumulative = true;
		for (auto it = next(m_bins.begin()); it != m_bins.end(); ++it)
			*it += *prev(it);
		return *this;
	} // will output cumulated value on output
	
	hist& use_averaged()
	{
		printAv = true;
		return *this;
	}; // will include averaged data on output
	
	int getIndex(double position) const
	{
		int binid;
		if (m_scale == scale::linear)
			binid = static_cast<int>((position-lb)/incr);
		else if (m_scale == scale::logarithmic) //m_scale == scale::logarithmic must hold
			binid = static_cast<int>((log(position)-logLb)/logIncr);
		else // m_scale == scale::quasi_logarithmic must hold
			binid = static_cast<int>((log(static_cast<int>(position+0.5))-logLb)/logIncr);

		return binid;
	}
	
	bin<T>& operator[](double position)
	{
		return m_bins[getIndex(position)];
	};
	
	bin<T>& at(double position)
	{
		int binid = getIndex(position);
		if (binid >= size)
		{
			stringstream errormsg;
			errormsg << "histogram at position "
				<< position
				<< " is requested, but upper boundary is "
				<< ub
				<< "\nLower boundary: "
				<< lb
				<< "\nRequested bin: "
				<< binid
				<< "\nBins available: 0 to " << size
				<< "\nThe scaling is " << m_scale << endl;
			throw out_of_range(errormsg.str().c_str());
		}
		if (binid < 0)
		{
			stringstream errormsg;
			errormsg << "histogram at position "
				     << position
					 << " is requested, but lower boundary is "
					 << lb
					 << "\nUpper boundary: "
					 << ub << "\n"
					 << "\nRequested bin: "
					 << binid
					 << "\nBins available: 0 to " << size
					 << "\nThe scaling is " << m_scale << endl;
			throw out_of_range(errormsg.str().c_str());
		}
		return m_bins[binid];
	}
	
	bool ifIncrAt(double position, T value) //increment at position if the position is inside the range of the histogram
	{
		if (position < lb || position >= ub)
			return false;

		int binid = getIndex(position);

		m_bins[binid].add(value);
		return true;
	}
	
	vector<bin<T>>& bins()
	{
		return m_bins;
	};

	vector<bin<T>> bins() const
	{
		return m_bins;
	}

	long unsigned int binCounter(int binID) const
	{
		return m_bins[binID].getCounter();
	}
	
	double binIntensity(int binID) const
	{
		return m_bins[binID].intensity();
	}
	
	int getSize() const
	{
		return size;
	}
	
	double getIncr() const
	{
		return incr;
	}
	
	unsigned long int sumCounter() const // to normalisation
	{
		int ret = 0;
		for_each(m_bins.begin(), m_bins.end(), [&ret](bin<T> element){ret += element.getCounter(); }); //is it abusing lambda functions?
		return ret;
	}

	double totalCounterArea() const // to represent propability
	{
		double ret = 0;
		for_each(m_bins.begin(), m_bins.end(), [&ret](bin<T> element){ret += element.getCounter() * element.width(); }); //is it abusing lambda functions?
		return ret;
	}

	T sumIntensity() const // to normalisation
	{
		T ret = 0;
		for_each(m_bins.begin(), m_bins.end(), [&ret](bin<T> element){ret += element.total(); }); //is it abusing lambda functions?
		return ret;
	}

	double totalIntensityArea() const // to represent propability
	{
		double ret = 0;
		for_each(m_bins.begin(), m_bins.end(), [&ret](bin<T> element){ret += element.intensity() * element.width(); }); //is it abusing lambda functions?
		return ret;
	}

};
