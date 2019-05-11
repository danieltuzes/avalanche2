
#include "stdafx.h"
#include "for_gnuplot_utils.h"
#include "../gen_utils/gen_math.h"

using namespace std;


g_block_data::g_block_data() : block_data() {}  //constructor inheritance is not allowed in VS2012


g_block_data::g_block_data(const g_block_data& tocopy) : block_data(tocopy) {}  //constructor inheritance is not allowed in VS2012
g_block_data::g_block_data(int s) : block_data(s) {}  //constructor inheritance is not allowed in VS2012
g_block_data::g_block_data(int s, double val) : block_data(s, val) {}  //constructor inheritance is not allowed in VS2012
template <typename U> g_block_data::g_block_data(const block_data<U>& tocopy) : block_data(tocopy.getSize())
{
	for (int i=0; i<size*size; ++i)
		data[i] = static_cast<double>(tocopy[i]);
}
g_block_data::g_block_data(string fname, dataformat f) : block_data(fname, f, 0) {}  //constructor inheritance is not allowed in VS2012
g_block_data::g_block_data(char * fname, dataformat f) : block_data(fname, f, 0) {}  //constructor inheritance is not allowed in VS2012
g_block_data::g_block_data(vector<vector<double>> m_vector) : block_data(m_vector) {}


template <typename U> g_block_data& g_block_data::operator=(const block_data<U>& tocopy)
{
	for (int i=0; i<size*size; ++i)
		data[i] = static_cast<double>(tocopy[i]);
	return *this;
}

g_block_data& g_block_data::modulateXSin(double lambda)
{
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			this->operator()(i, j) *= sin(i * 2 * static_cast<double>(PI) / lambda);
		}
	}
	return *this;
}

g_block_data g_block_data::tilize() const
{
	g_block_data Ti(size);
	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
			Ti(x,y) = data[(((x + size / 2) & (size - 1)) << power) + ((y + size / 2) & (size-1))];
	return Ti;
}

g_block_data& g_block_data::origin_imported(string fname, int skip)
{
	ifstream ifile(fname);
	vector<vector<double>> filearray;
	for (int i = 0; i < skip; ++i, ifile.ignore(numeric_limits<streamsize>::max(), '\n'));

	for (string line; getline(ifile, line);)
	{
		istringstream oneline(line);
		for (int j = 0; j < skip; ++j, oneline.ignore(numeric_limits<streamsize>::max(), '\t'));
		filearray.push_back(vector<double>());
		for (double tmp; oneline >> tmp; filearray.back().push_back(tmp));
	}
	size = filearray.size();
	power = ::getPower(size);

	data = new double[size*size];
	*this = block_data(filearray);
	return *this;
}

/*
string g_block_data::stats(string commentChar) const // moved to block_data<T>
{
	ostringstream text;
	auto max = getMaxWhere();
	auto min = getMinWhere();
	text << commentChar << "Maximum value: " << get<0>(max) << " at (x;y): (" << get<1>(max) << ";" << get<2>(max) << ")" << endl;
	text << commentChar << "Minimum value: " << get<0>(min) << " at (x;y): (" << get<1>(min) << ";" << get<2>(min) << ")" << endl;
	text << commentChar << "Sum: " << getSum() << endl;
	text << commentChar << "Abs sum: " << Abs().getSum() << endl;
	text << commentChar << "Assymetry: " << getAssymetry() << endl;
	text << commentChar << "Standard deviation: " << getStdDev() << endl;
	
	return text.str();
}
*/

template g_block_data::g_block_data(const block_data<long double>&);
template g_block_data& g_block_data::operator=(const block_data<long double>&);
template g_block_data::g_block_data(const block_data<double>&);
template g_block_data& g_block_data::operator=(const block_data<double>&);