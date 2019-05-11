#include "stdafx.h"
#include "../avalanche2/make_simulation.h"
#include "tau_g_stat.h"

using namespace std;

class avalanche : public tau_g_stat
{
private:
	bool diam_calced;
	bool CoG_calced;
public:
	vector<pNd> m_pNd;

	avalanche(tau_g_line tau_g, vector<pNd>::const_iterator first, vector<pNd>::const_iterator last) :
		tau_g_stat(tau_g_line(tau_g), distance(first, last), first->posx, first->posy, prev(last)->posx, prev(last)->posy),
		diam_calced(false),
		CoG_calced(false),
		m_pNd(first, last) {};

	void calcDiam(int sysSize) // "calculates the diameter" does not, faulty
	{
		diam_calced = true;
		if (m_pNd.size() == 1)
		{
			diamx = 0;
			diamy = 0;
			return;
		}

		vector<bool> xaffectedline(sysSize, false);
		vector<bool> yaffectedline(sysSize, false);

		for (const auto& av : m_pNd)
		{
			xaffectedline[av.posx] = true;
			yaffectedline[av.posy] = true;
		}

		int xmaxgap = 0;
		int xfirst_place = 0;
		for (; !xaffectedline[xfirst_place]; ++xfirst_place);
		int xlast_place = xfirst_place;
		for (int line_place = xfirst_place + 1; line_place < sysSize; ++line_place)
		{
			if (xaffectedline[line_place])
			{
				if (line_place - xlast_place > xmaxgap)
					xmaxgap = line_place;
				xlast_place = line_place;
			}
		}
		if (sysSize - xlast_place + xfirst_place > xmaxgap)
			xmaxgap = sysSize - xlast_place + xfirst_place;

		int ymaxgap = 0;
		int yfirst_place = 0;
		for (; !yaffectedline[yfirst_place]; ++yfirst_place);
		int last_place = yfirst_place;
		for (int line_place = yfirst_place + 1; line_place < sysSize; ++line_place)
		{
			if (yaffectedline[line_place])
			{
				if (line_place - last_place > ymaxgap)
					ymaxgap = line_place;
				last_place = line_place;
			}
		}
		if (sysSize - last_place + yfirst_place > ymaxgap)
			ymaxgap = sysSize - last_place + yfirst_place;


		// an n log n complexity calculation for the maxgap, where n is the number of affected cells during the avalanche

		//vector<pNd> m_pNd_sorted_x(m_pNd);
		//vector<pNd> m_pNd_sorted_y(m_pNd);
		//sort(m_pNd_sorted_x.begin(), m_pNd_sorted_x.end(), [](pNd a, pNd b){return a.posx<b.posx; }); // complexity: n log n
		//sort(m_pNd_sorted_y.begin(), m_pNd_sorted_y.end(), [](pNd a, pNd b){return a.posy<b.posy; });

		//int maxgapx = 0;
		//for (auto X = m_pNd_sorted_x.begin(); X != prev(m_pNd_sorted_x.end()); ++X) // complexity: n
		//	maxgapx = max(maxgapx, next(X)->posx - X->posx);

		//maxgapx = max(maxgapx, m_pNd_sorted_x.front().posx + sysSize - m_pNd_sorted_x.back().posx);

		//int maxgapy = 0;
		//for (auto X = m_pNd_sorted_y.begin(); X != prev(m_pNd_sorted_y.end()); ++X)
		//	maxgapy = max(maxgapy, next(X)->posy - X->posy);

		//maxgapy = max(maxgapy, m_pNd_sorted_y.front().posy + sysSize - m_pNd_sorted_y.back().posy);

		diamx = sysSize - xmaxgap;
		diamy = sysSize - ymaxgap;
	};

	void calcDiam_v2(int sysSize) // O(n*n), still false output
	{
		diam_calced = true;
		diamx = 0;
		diamy = 0;

		for (vector<pNd>::const_iterator eventA = m_pNd.begin(); eventA != m_pNd.end(); ++eventA)
		{
			for (vector<pNd>::const_iterator eventB = next(eventA); eventB != m_pNd.end(); ++eventB) // this can be expensive in case of large avalanches; an order nlogn can be implemented by ordering the elements
			{
				diamx = max(diamx, dist(eventA->posx, eventB->posx, sysSize));
				diamy = max(diamy, dist(eventA->posy, eventB->posy, sysSize));
			}
		}
	}

	void calcDiam_v3(int sysSize) // O(n*log(n)), the good result
	{
		diam_calced = true;
		diamx = 0;
		diamy = 0;
		if (m_pNd.size() == 1) // the rest of the code works only if the size >= 2
			return;
		
		vector<pNd> sorted_pNd(m_pNd);

		sort(sorted_pNd.begin(), sorted_pNd.end(), [](const pNd& A, const pNd& B){return A.posx < B.posx; });
		int largestXgap = 0;
		for (vector<pNd>::const_iterator m_event = sorted_pNd.begin(); next(m_event) != sorted_pNd.end(); ++m_event)
			largestXgap = max(largestXgap, next(m_event)->posx - m_event->posx);
		largestXgap = max(largestXgap, sorted_pNd.front().posx + sysSize - sorted_pNd.back().posx);
		diamx = sysSize - largestXgap;

		sort(sorted_pNd.begin(), sorted_pNd.end(), [](const pNd& A, const pNd& B){return A.posy < B.posy; });
		int largestYgap = 0;
		for (vector<pNd>::const_iterator m_event = sorted_pNd.begin(); next(m_event) != sorted_pNd.end(); ++m_event)
			largestYgap = max(largestYgap, next(m_event)->posy - m_event->posy);
		largestYgap = max(largestYgap, sorted_pNd.front().posy + sysSize - sorted_pNd.back().posy);
		diamy = sysSize - largestYgap;
	}

	void calcCoG(int sysSize) //calculates the diameter
	{
		CoG_calced = true;
		double xuav = 0;
		double xpav = 0;
		double yuav = 0;
		double ypav = 0;
		for (const auto& X : m_pNd)
		{
			xuav += cos(X.posx * 2 * PI / sysSize) * abs(X.deform);
			xpav += sin(X.posx * 2 * PI / sysSize) * abs(X.deform);
			yuav += cos(X.posy * 2 * PI / sysSize) * abs(X.deform);
			ypav += sin(X.posy * 2 * PI / sysSize) * abs(X.deform);
		}

		CoGx = sysSize * (atan2(-xpav, -xuav) + PI) / (2 * PI); //atan2 returns value in range [-PI,PI)
		CoGy = sysSize * (atan2(-ypav, -yuav) + PI) / (2 * PI);
	};

	template <typename T> void addShapeto(block_data<T>& shape) const
	{
		const int x0 = m_pNd.front().posx;
		const int y0 = m_pNd.front().posy;
		for (const auto& x : m_pNd)
			shape.periodic(x.posx - x0, x.posy - y0) += x.deform;
	}

	friend ostream& operator<<(ostream& o, const avalanche& av)
	{
		o << tau_g_line(av) << "\t"
			<< av.avSize << "\t";
		if (av.diam_calced)
			o << av.diamx << "\t"
			<< av.diamy << "\t";
		else
			o << "n/a" << "\t"
			<< "n/a" << "\t";
		o << av.firstposx << "\t"
			<< av.firstposy << "\t"
			<< av.lastposx << "\t"
			<< av.lastposy << "\t";
		if (av.CoG_calced)
			o << linvec<>(av.CoG());
		else
			o << "n/a" << "\t"
			<< "n/a";
		return o;
	}
	//friend istream& operator>>(istream& i, avalanche& av);
};


bool fill_avalanches(string pNd_fname, bool readTillEnd, string tau_g_fname, vector<avalanche>& data, int sysSize, double factor)
{
	ifstream tau_g_file(tau_g_fname);
	if (!tau_g_file)
	{
		cerr << "Can't open " << tau_g_fname << endl;
		return false;
	}
	skip_comment(tau_g_file);


	vector<pNd> pNd_data;
	if (!get_pNd_content(pNd_fname, pNd_data, sysSize, factor))
	{
		cerr << "Can't get pNd content" << endl;
		return false;
	}


	vector<tau_g_line> tau_g_content;
	for (tau_g_line ttmp; tau_g_file >> ttmp; tau_g_content.push_back(ttmp*factor));

	if (tau_g_content.empty())
	{
		cerr << "Cannot get line from " << tau_g_fname << ". Maybe EOL error?" << endl;
		return false;
	}
	else if (!readTillEnd)
		tau_g_content.pop_back(); //the last line describe a non-stopped avalacnhe OR the same one as the previous one


	int initRnd = tau_g_content.begin()->random_times;
	if (initRnd != sysSize*sysSize * 2 && initRnd != sysSize * sysSize)
		cout << "Warning: the first rand_gen.getTimesAdvanced() value is not " << sysSize << "*" << sysSize << ", nor *=2" << endl; //tau_g_line.random_times first value is nonzero due to initialization

	for (auto it = next(tau_g_content.begin()); it != tau_g_content.end(); ++it)
	{
		auto m_begin = next(pNd_data.begin(), prev(it)->random_times - initRnd);
		auto m_end = next(pNd_data.begin(), it->random_times - initRnd);
		data.emplace_back(*it, m_begin, m_end);
	}

	return true;
}

bool fill_avalanches(string pNd_fname, string tau_g_fname, vector<avalanche>& data, int sysSize, double factor)
{
	return fill_avalanches(pNd_fname, false, tau_g_fname, data, sysSize, factor);
}