#include "stdafx.h"
#include "gen_fmanip.h"

using namespace std;

class fractal2D : public block_data<>
{
public:

	fractal2D() : block_data<>() {};
	fractal2D(int size) : block_data<>(size) {};
	fractal2D(const block_data<>& tocopy) : block_data<>(tocopy) {};
	fractal2D(string fname, dataformat format) : block_data<>(fname, format) {};

	double getMinInRange(int leftbottomx, int leftbottomy, int righttopx, int righttopy) const
	{
		leftbottomx &= (size - 1);
		leftbottomy &= (size - 1);
		righttopy &= (size - 1);
		righttopx &= (size - 1);
		double minVal = data[(leftbottomx << power) + leftbottomy];

		for (int x = leftbottomx;;)
		{
			for (int y = leftbottomy;;)
			{
				if (data[(x << power) + y] < minVal)
					minVal = data[(x << power) + y];
				periodicIncr(y);
				if (y == righttopy)
					break;
			}
			periodicIncr(x);
			if (x == righttopx)
				break;
		}
		return minVal;
	}
	double getMaxInRange(int leftbottomx, int leftbottomy, int righttopx, int righttopy) const
	{
		leftbottomx &= (size - 1);
		leftbottomy &= (size - 1); 
		righttopy &= (size - 1);
		righttopx &= (size - 1);
		double maxVal = data[(leftbottomx << power) + leftbottomy];

		for (int x = leftbottomx;;)
		{
			for (int y = leftbottomy;;)
			{
				if (data[(x << power) + y] > maxVal)
					maxVal = data[(x << power) + y];
				periodicIncr(y);
				if (y == righttopy)
					break;
			}
			periodicIncr(x);
			if (x == righttopx)
				break;
		}
		return maxVal;
	}

	double area(double distance) const
	{
		double area = 0;
		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size; ++j)
			{
				// A: (0,0); B: (1,0); C: (1,1); D: (0,1)
				// a: AB; b: BC; c: CD; d: DA; e: AC
				// T1: ABC
				// T2: ACD
				double A = static_cast<double>(periodic(i, j));
				double B = static_cast<double>(periodic(i + 1, j));
				double C = static_cast<double>(periodic(i + 1, j + 1));
				double D = static_cast<double>(periodic(i, j + 1));

				double a = static_cast<double>(hypot(B - A, distance));
				double b = static_cast<double>(hypot(C - B, distance));
				double c = static_cast<double>(hypot(D - C, distance));
				double d = static_cast<double>(hypot(A - D, distance));
				double e = static_cast<double>(hypot(C - A, distance*SQRT2));

				double s = (a + b + e) / 2;
				area += sqrt(s*(s - a)*(s - b)*(s - e));
				s = (c + d + e) / 2;
				area += sqrt(s*(s - c)*(s - d)*(s - e));
			}
		return area;
	};
	int boxPointCount(int boxWidth) const
	{
		int power = ::getPower(boxWidth);
		int sum = 0;
		vector<bool> heights;
		
		//Tracker my_tracker(maxNofElements);
		for (int X = 0; X < size / boxWidth; ++X)
			for (int Y = 0; Y < size / boxWidth; ++Y)
			{
				//my_tracker.reset();
				double minval = getMinInRange(X * boxWidth, Y * boxWidth, (X + 1) * boxWidth, (Y + 1) * boxWidth);
				double maxval = getMaxInRange(X * boxWidth, Y * boxWidth, (X + 1) * boxWidth, (Y + 1) * boxWidth);
				int maxNofElements = static_cast<int>((maxval - minval) / boxWidth * size) + 1;
				heights.assign(maxNofElements, false);
				for (int x = 0; x < boxWidth; ++x)
					for (int y = 0; y < boxWidth; ++y)
					{
						//heights.push_back(fint(static_cast<double>((this->operator()(X*boxWidth + x, Y*boxWidth + y)) / boxWidth * size)));
						//int actval = fint(static_cast<double>((this->operator()(X*boxWidth + x, Y*boxWidth + y)) / boxWidth * size));
						int actvalpos = fint((static_cast<double>(this->operator()((X << power) + x, (Y << power) + y)) - minval) / boxWidth * size);
						if (!heights[actvalpos])
						{
							heights[actvalpos] = true;
							++sum;
						}
						//my_tracker.set(actvalpos);
					}

				//sum += my_tracker.count();
				//sum += heights.size();
			}
		return sum;
	};

	double mhdx(int distance) const
	{
		double diff = 0;
		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size; ++j)
			{
				diff += abs(periodic(i, j) - periodic(i + distance, j));
			}
		return diff;
	};
	double mhdy(int distance) const
	{
		double diff = 0;
		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size; ++j)
			{
				diff += abs(periodic(i, j) - periodic(i, j + distance));
			}
		return diff;
	};
	double mhdsqx(int distance) const
	{
		double diff = 0;
		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size; ++j)
			{
				diff += (periodic(i, j) - periodic(i + distance, j)) *
					(periodic(i, j) - periodic(i + distance, j));
			}
		return diff;
	};
	double mhdsqy(int distance) const
	{
		double diff = 0;
		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size; ++j)
			{
				diff += (periodic(i, j) - periodic(i, j + distance)) *
					(periodic(i, j) - periodic(i, j + distance));
			}
		return diff;
	};
	vector<double> mhdx() const
	{
		vector<double> ret(size, 0);
		for (int distance = 0; distance < size; ++distance)
			ret[distance] = mhdx(distance);
		return ret;
	}; 
	vector<double> mhdy() const
	{
		vector<double> ret(size, 0);
		for (int distance = 0; distance < size; ++distance)
			ret[distance] = mhdy(distance);
		return ret;
	};
	vector<double> mhdsqx() const
	{
		vector<double> ret(size, 0);
		for (int distance = 0; distance < size; ++distance)
			ret[distance] = mhdsqx(distance);
		return ret;
	};
	vector<double> mhdsqy() const
	{
		vector<double> ret(size, 0);
		for (int distance = 0; distance < size; ++distance)
			ret[distance] = mhdsqy(distance);
		return ret;
	};

	vector<double> pIntegral_y(int x)
	{
		vector<double> ret;
		for (int y = 0; y < size; ++y)
			ret.push_back(operator()(x, y));

		double sumy = accumulate(ret.begin(), ret.end(), 0.);

		double deltaFy = sumy / size;

		for (int y = 1; y < size; ++y)
			ret[y] += ret[y - 1] - deltaFy;

		return ret;
	}
	vector<double> pIntegral_x(int y)
	{
		vector<double> ret(size, 0);
		for (int x = 0; x < size; ++x)
			ret[x] = operator()(x, y);

		double sumx = 0;
		for (int x = 0; x<size; ++x)
			sumx += operator()(x, y);

		double deltaFx = sumx / size;

		for (int x = 1; x< size; ++x)
			ret[x] += ret[x - 1] + deltaFx;

		return ret;
	}
};