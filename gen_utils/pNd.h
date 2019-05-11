#include "stdafx.h"

typedef double deftype;

#pragma once
using namespace std;

class place
{
public:
	int x;
	int y;
	deftype tau_l;

	place();
	place(int x, int y, deftype tau_l);

	void setVal(int x, int y, deftype tau_l);
};

class pNd
{
public:
	int posx;
	int posy;
	deftype deform;

	pNd(int posx, int posy, deftype deform);
	pNd(int linPos, deftype deform, int sysSize, int power);
	pNd(const place& minPoint, deftype deform);

	int getLinPos(int power) const;
};

ostream& operator<< (ostream& o, const pNd& m_pNd);

bool get_pNd_content(string fname, vector<pNd>& data, int sysSize, double factor);