#include "stdafx.h"
#include "gen_hist.h"
#include "gen_math.h"
#include "gen_param.h"

using namespace std;



#pragma region template class accum member definitions

template <typename T> accum<T>::accum(T data) : counter(0), data(data) {};

template <typename T> void accum<T>::add(const T& toAdd) // include another point in the histogram
{
	data += toAdd;
	++counter;
}

template <typename T> void accum<T>::add(T&& toAdd)
{
	this->add(toAdd);
}

template <typename T> T accum<T>::total() const
{
	return data;
}

template <typename T> T accum<T>::averaged() const // returns the averaged data; 
{
	return (data / counter);
}
	
template <typename T> ostream& operator<<(ostream& o, const accum<T>& m_accum)
{
	o << m_accum.counter << "\t" << m_accum.data;
	return o;
}

#pragma endregion


#pragma region template class bin member definitions

template <typename T> bin<T>::bin(double lb, double ub, T data) : lb(lb), ub(ub), accum(data) {};
	
template <typename T> linvec<double> bin<T>::binrange()
{
	vector<double> basic_ret;
	basic_ret.push_back(lb);
	basic_ret-psuh_back(ub);
	return linvec<double>(basic_ret);
}

template <typename T> ostream& operator<<(ostream& o, const bin<T>& mybin)
{
	o << mybin.lb << "\t" << mybin.ub << "\t" << accum<T>(mybin);
	return o;
};

#pragma endregion


#pragma region template class hist member definitions

template <typename T> hist<T>::hist(double lb, double ub, scale m_scale, int size, T data) : lb(lb), ub(ub), m_scale(m_scale), size(size), logIncr(log(incr)), logLb(log(lb)), printAv(false)
{
	if (m_scale == scale::linear)
	{
		incr = (ub-lb) / size;
		for (int i=0; i<size; ++i)
			m_bins.push_back(bin<T>(lb + i*incr, lb + (i+1)*incr,data));
	}
	else //m_scale == scale::logarithmic || m_scale == scale::quasi_log must hold
	{
		incr = pow(ub/lb,1./size);
		logIncr = log(incr);
		logLb = log(lb);
		if (m_scale == scale::linear)
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

};
	
template <typename T> ostream& operator<< (ostream& o, const hist<T>& m_hist)
{
	if (!m_hist.printAv)
		for (auto X:m_hist.m_bins)
			o << X << endl;
	else
		for (auto X:m_hist.m_bins)
			o << X << "\t" << X.averaged() << endl;
	return o;
};
	
template <typename T> hist<T>& hist<T>::use_averaged()
{
	printAv = true;
	return *this;
};

template <typename T> bin<T>& hist<T>::at(double position)
{
	int binid;
	if (m_scale == scale::linear)
		binid = static_cast<int>((position-lb)/incr);
	else if (m_scale == scale::logarithmic) //m_scale == scale::logarithmic must hold
		binid = static_cast<int>((log(position)-logLb)/logIncr);
	else
		binid = static_cast<int>(log(static_cast<int>(position)-logLb)/logIncr);
	if (binid >= size)
		binid = size-1;

	return m_bins[binid];
};

template <typename T> bin<T>& hist<T>::operator[](int position)
{
	return m_bins[position];
};

template <typename T> vector<bin<T>>& hist<T>::bins()
{
	return m_bins;
}


#pragma endregion


template class hist<linvec<int>>;
template class hist<linvec<double>>;
template class hist<double>;
template class hist<int>;
template class accum<linvec<int>>;
template class accum<linvec<double>>;
template class accum<double>;
template class accum<int>;

template ostream& operator<<(ostream& o, const hist<linvec<int>>& m_vec);
template ostream& operator<<(ostream& o, const hist<linvec<double>>& m_vec);
/*
template ostream& operator<<(ostream& o, const hist<int>& m_vec);
template ostream& operator<<(ostream& o, const hist<double>& m_vec);
template ostream& operator<<(ostream& o, const bin<int>& m_vec);
template ostream& operator<<(ostream& o, const bin<double>& m_vec);
template ostream& operator<<(ostream& o, const bin<linvec<int>>& m_vec);
template ostream& operator<<(ostream& o, const bin<linvec<double>>& m_vec);
template ostream& operator<<(ostream& o, const accum<int>& m_vec);
template ostream& operator<<(ostream& o, const accum<double>& m_vec);
template ostream& operator<<(ostream& o, const accum<linvec<int>>& m_vec);
template ostream& operator<<(ostream& o, const accum<linvec<double>>& m_vec);*/