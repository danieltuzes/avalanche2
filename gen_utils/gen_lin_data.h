
#pragma once
#include "stdafx.h"
#include "gen_siomanip.h"
#include "gen_math.h"

using namespace std;


class lin_data
{
public:
	int size;
	int power;
	int mask;
	vector<double> elements;

	double operator[](int x) const
	{
		return elements[x];
	}
	double& operator[](int x)
	{
		return elements[x];
	}

	lin_data& operator += (double add)
	{
		for (int i = 0; i < size; ++i)
			elements[i] += add;
		return *this;
	}
	lin_data operator +(double add) const
	{
		lin_data ret(*this);
		ret += add;
		return ret;
	}
	lin_data& operator += (const lin_data& add)
	{
		for (int i = 0; i < size; ++i)
			elements[i] += add[i];
		return *this;
	}
	lin_data operator +(const lin_data& add) const
	{
		lin_data ret(*this);
		ret += add;
		return ret;
	}
	lin_data& operator -= (double subs)
	{
		for (int i = 0; i < size; ++i)
			elements[i] -= subs;
		return *this;
	}
	lin_data operator -(double subs) const
	{
		lin_data ret(*this);
		ret -= subs;
		return ret;
	}
	lin_data& operator -= (const lin_data& subs)
	{
		for (int i = 0; i < size; ++i)
			elements[i] -= subs[i];
		return *this;
	}
	lin_data operator -(const lin_data& subs) const
	{
		lin_data ret(*this);
		ret -= subs;
		return ret;
	}
	lin_data& operator *= (double factor)
	{
		for (int i = 0; i < size; ++i)
			elements[i] *= factor;
		return *this;
	}
	lin_data operator *(double factor) const
	{
		lin_data ret(*this);
		ret *= factor;
		return ret;
	}
	lin_data& operator *=(const lin_data& multi)
	{
		for (int i = 0; i < size; ++i)
			elements[i] *= multi[i];
		return *this;
	}
	lin_data operator *(const lin_data& multi)
	{
		vector<double> ret;
		for (int i = 0; i < size; ++i)
			ret.push_back(multi[i] * elements[i]);
		return ret;
	}
	lin_data& operator /= (double divisor)
	{
		return *this *= 1 / divisor;
	}
	lin_data operator / (double divisor) const
	{
		lin_data ret(*this);
		ret /= divisor;
		return ret;
	}

	double periodic(int x) const
	{
		return elements[x & mask];
	}
	double& periodic(int x)
	{
		return elements[x & mask];
	}

	int minplace()
	{
		int minplace = 0;
		for (int i = 1; i < size; ++i)
		if (elements[i] < elements[minplace])
			minplace = i;
		return minplace;

	}
	double laplace(int x) const
	{
		return (periodic(x - 1) + periodic(x + 1) - 2 * elements[x]);
	}
	vector<double> laplaced() const
	{
		vector<double> ret(size, 0);
		for (int i = 0; i < size; ++i)
			ret[i] = laplace(i);

		return ret;
	}
	vector<double> mhd() //mean high difference
	{
		vector<double> ret(size, 0);
		for (int i = 0; i < size; ++i)
		{
			for (int j = 1; j < size; ++j)
			{
				ret[i] += abs(elements[j] - periodic(i + j));
			}
			ret[i] /= size;
		}
		return ret;
	}
	vector<double> mhdsq() //mean high difference-square
	{
		vector<double> ret(size, 0);
		for (int i = 0; i < size; ++i)
		{
			for (int j = 1; j < size; ++j)
			{
				ret[i] += (elements[j] - periodic(i + j)) * (elements[j] - periodic(i + j));
			}
			ret[i] /= size;
		}
		return ret;
	}
	vector<double> autocorr_bruteforce()
	{
		vector<double> ret(size, 0);
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				ret[i] += elements[j] * periodic(j + i);
			}
			ret[i] /= size;
		}
		return ret;
	}

	lin_data reMap() const
	{
		vector<double> ret;
		for (int i = 0; i < size / 2; ++i)
			ret.push_back((elements[2 * i] + elements[2 * i + 1]) / 2);
		return ret;
	}
	lin_data reSample(int blockSize = 2, int shift = 0) const
	{
		vector<double> ret;
		for (int i = 0; i < size / blockSize; ++i)
			ret.push_back(elements[blockSize * i + shift]);
		return ret;
	}
	double length(double distance = 1) const
	{
		double sum = 0;
		for (int i = 0; i < size; ++i)
			sum += sqrt((elements[i] - periodic(i + 1)) * (elements[i] - periodic(i + 1)) + distance*distance);

		return sum;
	}
	int boxPointCount(int boxWidth, double verticalRescale = 1) //boxWidth is in units of distance of neighboring points 
	{
		vector<int> height;
		int boxPower = getPower(boxWidth);
		int sum = 0; //number of boxes
		for (int I = 0; I < size / boxWidth; ++I)
		{
			height.clear();
			for (int i = 0; i < boxWidth; ++i)
				height.push_back(fint(elements[(I << boxPower) + i] / boxWidth * size * verticalRescale));

			sort(height.begin(), height.end());
			sum += distance(height.begin(), unique(height.begin(), height.end()));
		}
		return sum;
	}
	int boxCurveCount(int boxWidth, double verticalRescale = 1) //boxWidth is in units of distance of neighboring points 
	{
		vector<int> height;
		int sum = 0; //number of boxes
		for (int i = 0; i < size;)
		{
			height.clear();
			for (int j = 0; j < boxWidth; ++j)
			{
				int thisbox = fint(elements[i] / boxWidth * size * verticalRescale);
				for (int heighint = fint((elements[i] + periodic(i - 1)) / 2 / boxWidth * size * verticalRescale); heighint != thisbox;)
				{
					height.push_back(heighint);
					if ((elements[i] + periodic(i - 1)) / 2 < elements[i])
						++heighint;
					else
						--heighint;
				}
				height.push_back(thisbox);
				for (int heighint = fint((elements[i] + periodic(i + 1)) / boxWidth * size / 2 * verticalRescale); heighint != thisbox;)
				{
					height.push_back(heighint);
					if ((elements[i] + periodic(i + 1)) / 2 < elements[i])
						++heighint;
					else
						--heighint;
				}
				++i;
			}
			sort(height.begin(), height.end());
			sum += distance(height.begin(), unique(height.begin(), height.end()));
		}
		return sum;
	}

	bool fWrite(string fname, string header = "") const
	{
		ofstream ofile(fname);
		if (!ofile)
		{
			cerr << "Cant't create " << fname << endl;
			return false;
		}
		cout << "Writing " << fname << endl;
		ofile.precision(17);
		ofile << header;
		for (int i = 0; i < size; ++i)
			ofile << i << "\t" << elements[i] << "\n";
		return true;
	}
	bool fRead(string fname)
	{
		ifstream ifile(fname);
		if (!ifile)
		{
			cerr << "Can't open " << fname << endl;
			return false;
		}
		skip_comment(ifile);
		int i;
		for (double b; ifile >> i >> b;)
			elements[i] = b;
		if (i + 1 != size)
		{
			cerr << "Size mismach: " << fname << " contains only " << i + 1 << " number of elements but " << size << " is expected" << endl;
			return false;
		}
		return true;
	}

	lin_data(int size, double value) : size(size), power(getPower(size)), mask(size - 1), elements(size, value) {}
	lin_data(vector<double> elements) : size(elements.size()), power(getPower(elements.size())), mask(elements.size() - 1), elements(elements) {}
	lin_data(string fname)
	{
		ifstream ifile(fname);
		if (!ifile)
			cerr << "Can't open fname." << endl;
		skip_comment(ifile);
		int tmp_size;
		for (double val; ifile >> tmp_size >> val;)
			elements.push_back(val);

		if (1 << getPower(tmp_size) != tmp_size + 1)
			cerr << "Size is not power of two.";

		size = tmp_size + 1;
		power = getPower(size);
		mask = size - 1;
	}
};
