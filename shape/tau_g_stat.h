// tau_g_stat.h : stores tau_g information with additional information on avalanche shape

#include "stdafx.h"
#include "../gen_utils/tau_g.h"
#include "../gen_utils/gen_math.h"
#include "../gen_utils/gen_hist.h"
using namespace std;

#pragma once



class tau_g_stat : public tau_g_line
{
protected:
	int avSize;
	int diamx;
	int diamy;
	int firstposx;
	int firstposy;
	int lastposx;
	int lastposy;
	double CoGx;
	double CoGy;
public:
	tau_g_stat() {}

	tau_g_stat(tau_g_line tgl) : tau_g_line(tgl) {}
	
	tau_g_stat(tau_g_line tgl, int avSize) :
		tau_g_line(tgl),
		avSize(avSize) {}

	tau_g_stat(tau_g_line tgl, int avSize, int firstposx, int firstposy, int lastposx, int lastposy) :
		tau_g_line(tgl),
		avSize(avSize), 
		firstposx(firstposx), firstposy(firstposy),
		lastposx(lastposx), lastposy(lastposy) {}

	tau_g_stat(tau_g_line tgl, int avSize, int diamx, int diamy, int firstposx, int firstposy, int lastposx, int lastposy, double CoGx, double CoGy):
		tau_g_line(tgl),
		avSize(avSize),
		diamx(diamx), diamy(diamy), 
		firstposx(firstposx), firstposy(firstposy),
		lastposx(lastposx), lastposy(lastposy),
		CoGx(CoGx), CoGy(CoGy) {}

	vector<double> diam()
	{
		vector<double> ret;
		ret.push_back(diamx);
		ret.push_back(diamy);
		return ret;
	}

	vector<double> CoG() const
	{
		vector<double> ret;
		ret.push_back(CoGx);
		ret.push_back(CoGy);
		return ret;
	}

	friend int distX(const tau_g_stat&, const tau_g_stat&);
	friend int distY(const tau_g_stat&, const tau_g_stat&);

	friend ostream& operator<<(ostream&o, const tau_g_stat& m_tgs);
	friend istream& operator>>(istream&i, tau_g_stat& m_tgs);
};

int distX(const tau_g_stat& a, const tau_g_stat& b)
{
	return (a.firstposx - b.firstposx);
}

int distY(const tau_g_stat& a, const tau_g_stat& b)
{
	return (a.firstposy - b.firstposy);
}


ostream& operator<<(ostream& o, const tau_g_stat& m_tgs)
{
	o << tau_g_line(m_tgs) << "\t"
		<< m_tgs.avSize << "\t"
		<< m_tgs.diamx << "\t"
		<< m_tgs.diamy << "\t"
		<< m_tgs.firstposx << "\t"
		<< m_tgs.firstposy << "\t"
		<< m_tgs.lastposx << "\t"
		<< m_tgs.lastposy << "\t"
		<< m_tgs.CoGx << "\t"
		<< m_tgs.CoGy;

	return o;
}


istream& operator>>(istream& i, tau_g_stat& m_tgs)
{
	m_tgs = tau_g_line(i);
	i >> m_tgs.avSize
		>> m_tgs.diamx
		>> m_tgs.diamy
		>> m_tgs.firstposx
		>> m_tgs.firstposy
		>> m_tgs.lastposx
		>> m_tgs.lastposy
		>> m_tgs.CoGx
		>> m_tgs.CoGy;

	return i;
}