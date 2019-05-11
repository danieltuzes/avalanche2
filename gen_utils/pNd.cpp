#include "stdafx.h"
#include "pNd.h"
#include "../gen_utils/gen_math.h"


place::place() {}
place::place(int x, int y, deftype tau_l) : x(x), y(y), tau_l(tau_l) {}
void place::setVal(int x, int y, deftype tau_l)
{
	*this = place(x, y, tau_l);
}

pNd::pNd(int posx, int posy, deftype deform) : posx(posx), posy(posy), deform(deform) {}
pNd::pNd(int linPos, deftype deform, int sysSize, int power) : posx(getX(linPos,power)), posy(getY(linPos,sysSize - 1)), deform(deform) {}
pNd::pNd(const place& minPoint, deftype deform) : posx(minPoint.x), posy(minPoint.y), deform(deform) {}

int pNd::getLinPos(int power) const
{
	return ::getLinPos(posx, posy, power);
}

ostream& operator<< (ostream& o, const pNd& m_pNd)
{
	o << m_pNd.posx << "\t" << m_pNd.posy << "\t" << m_pNd.deform;
	return o;
}


bool get_pNd_content(string fname, vector<pNd>& data, int sysSize, double factor)
{
	int power = getPower(sysSize);

	FILE * ifile;
	if ((ifile = fopen(fname.c_str(), "rb")) == NULL)
	{
		cerr << "Unable to open " << fname << endl;
		return false;
	}

	int linPos_tmp;
	double def_tmp;
	while (fread(&linPos_tmp, sizeof(int), 1, ifile) && fread(&def_tmp, sizeof(deftype), 1, ifile))
		data.emplace_back(linPos_tmp, def_tmp*factor, sysSize, power);

	fclose(ifile);
	return true;
}