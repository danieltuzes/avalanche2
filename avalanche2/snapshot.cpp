// snapshot.cpp : defines snapshot savings and loadings
//


#include "stdafx.h"
#include "snapshot.h"
using namespace std;


snapshot::snapshot() : marked(false) {}
snapshot::snapshot(deftype sG, deftype sGn, deftype av_size, deftype av_size_n, deftype tau_ext, int rnd_state, int nominal_def)
	: marked(false), sG(sG), sGn(sGn), av_size(av_size), av_size_n(av_size_n), tau_ext(tau_ext), rnd_state(rnd_state), nominal_def(nominal_def) {}
snapshot::snapshot(const simVars& sim, int nominal_def)
	: marked(false), sG(sim.sG), sGn(sim.sGn), av_size(sim.av_size), av_size_n(sim.av_size_n), tau_ext(sim.tau_ext), rnd_state(sim.rand_gen.getTimesAdvanced()), nominal_def(nominal_def)
{
	m_sim = &sim;
}

	
istream& operator >> (istream& ist, snapshot& sh)
{
	if (ist.get() == '#')
		sh.marked = true;
	else
	{
		sh.marked = false; //ez meg miért kell, amikor a def ctor egy false-sal hozza létre?
		ist.unget();
	}
	//

	ist >> sh.sG
		>> sh.sGn
		>> sh.av_size
		>> sh.av_size_n
		>> sh.tau_ext
		>> sh.rnd_state
		>> sh.nominal_def;

	return ist;
}

ostream& operator << (ostream& ost, snapshot sh)
{
	if (sh.marked)
		ost << '#';

	streamsize orig_prec = ost.precision();
	ost.precision(18);
	ost << sh.sG << "\t"
		<< sh.sGn << "\t"
		<< sh.av_size << "\t"
		<< sh.av_size_n << "\t"
		<< sh.tau_ext << "\t"
		<< sh.rnd_state << "\t"
		<< sh.nominal_def << "\t";
	ost.precision(orig_prec);
	return ost;
}

bool operator==(const snapshot &lhs, const snapshot &rhs)
{
	if (lhs.av_size == rhs.av_size &&
		lhs.av_size_n == rhs.av_size_n &&
		lhs.nominal_def == rhs.nominal_def &&
		lhs.rnd_state == rhs.rnd_state &&
		lhs.sG == rhs.sG &&
		lhs.sGn == rhs.sGn &&
		lhs.tau_ext == rhs.tau_ext)
		return true;
	else return false;
}

bool operator!=(const snapshot &lhs, const snapshot &rhs)
{
	if (lhs == rhs)
		return false;
	else return true;
}

bool operator<(const snapshot &lhs, const snapshot &rhs)
{
	if (lhs.sG + lhs.av_size - lhs.sGn - lhs.av_size_n < rhs.sG + rhs.av_size - rhs.sGn - rhs.av_size_n)
		return true;
	else
		return false;
}

bool operator>(const snapshot &lhs, const snapshot &rhs)
{
	if (lhs.sG + lhs.av_size - lhs.sGn - lhs.av_size_n > rhs.sG + rhs.av_size - rhs.sGn - rhs.av_size_n)
		return true;
	else
		return false;
}

bool operator<=(const snapshot &lhs, const snapshot &rhs) //lhs == rhs || lhs < rhs
{
	if (lhs == rhs || lhs < rhs)
		return true;
	else
		return false;
}

bool operator>=(const snapshot &lhs, const snapshot &rhs) //lhs == rhs || lhs > rhs
{
	if (lhs == rhs || lhs > rhs)
		return true;
	else
		return false;
}

bool snapshot::mark(const simPars& pars) const
{
	vector<snapshot> snapshots = getSnapshots(pars);

	bool contiF = false;
	ofstream conti(fnameWp("conti",pars));
	for (auto x:snapshots)
	{
		if (x == *this)
		{
			if (contiF == true) //if another was already found
				cerr << "Error: found a contiLine equal to another." << endl;
			else
			{
				snapshot marked(*this);
				marked.setMarked();
				conti << marked;
				contiF = true;
			}
		}
		else if (x.getNominal_def() == getNominal_def())
			cerr << "Two different snapshots with the same nominal_def were found" << endl;
		else
			conti << x;
	}

	if (!contiF)
	{
		cerr << "contiLine " << *this << " not found in " << fnameWp("conti",pars) << " to mark it" << endl;
		return false;
	}

	string names[]={"Gamma","tau_l","tau_n","tau_p"};
	for (auto name:names)
	{
		string fname = fnameWp(name,pars,nominal_def);
		if (rename(fname.c_str(),(fname + "__mark").c_str()) != 0)
		{
			perror(("Error marking file " + fname).c_str());
			cerr << "File " << fname << " cannot marked" << endl;
			return false;
		}
	}

	cout << "Snapshot at nominal def " << this->nominal_def << " was marked" << endl;
	
	return true;
}


//bool contiLine::delSnapshot(const simPars& pars) const // delete a snapshot that correspont to this contiLine
//{
//	vector<contiLine> conti_content = getContiContent(pars);
//
//	bool contiF = false;
//	ofstream conti(fnameWp("conti",pars));
//	for (auto x:conti_content)
//	{
//		if (x == *this)
//		{
//			if (contiF = true) //if another was already found
//				cerr << "Error: found a contiLine equal to another." << endl;
//			else
//				contiF = true;
//		}
//		else if (x.getNominal_def() == getNominal_def())
//			cerr << "Two different snapshots with the same nominal_def were found" << endl;
//		else
//			conti << x;
//	}
//
//	if (!contiF)
//	{
//		cerr << "contiLine " << *this << " not found in " << fnameWp("conti",pars) << " to delete it" << endl;
//		return false;
//	}
//
//	string names[]={"Gamma","tau_l","tau_n","tau_p"};
//	for (auto name:names)
//		if (remove(fnameWp(name,pars).c_str()) != 0)
//		{
//			perror(string("Error deleting file " + fnameWp(name,pars) + "\n").c_str());
//			cerr << "File " << fnameWp(name,pars) << " cannot deleted" << endl;
//			return false;
//		}
//
//	cout << "Snapshot at nominal def " << this->nominal_def << " was deleted" << endl;
//	
//	return true;
//}

bool delMarkedSnapshots(const simPars& pars)
{
	vector <snapshot> snapshots = getSnapshots(pars);
	vector <snapshot> snapsRemained;

	for (auto x:snapshots)
		if (x.isMarked())
		{
			string names[]={"Gamma","tau_l","tau_n","tau_p"};
			for (auto name:names)
			{
				string fname = fnameWp(name,pars,x.getNominal_def());
				if (remove((fname + "__mark").c_str()) != 0)
				{
					perror(string("Error deleting file " + fname).c_str());
					cerr << "File " << fname << " cannot deleted" << endl;
					return false;
				}
			}
		}
		else
			snapsRemained.push_back(x);

	return writeContiFile(pars,snapsRemained);
}


bool snapshot::make(const simPars& pars) const
{
	m_sim->Gamma.writeToFile(fnameWp("Gamma",pars,nominal_def));
	m_sim->tau_l.writeToFile(fnameWp("tau_l",pars,nominal_def));
	m_sim->tau_n.writeToFile(fnameWp("tau_n",pars,nominal_def));
	m_sim->tau_p.writeToFile(fnameWp("tau_p",pars,nominal_def));

	vector<snapshot> snapshots;
	snapshots = getSnapshots(pars);

	snapshots.emplace_back(*m_sim,nominal_def);
	sort(snapshots.begin(),snapshots.end(),[](snapshot a, snapshot b) {return (a.getNominal_def() < b.getNominal_def());});
	writeContiFile(pars,snapshots);

	return true;
}

bool snapshot::check(const simPars& pars) const
{
	cerr << "Warning: checkSnapshot not implemented" << endl;
	return true;
}

deftype snapshot::getDeform() const
{return (sG - sGn + av_size - av_size_n);}


int snapshot::getNominal_def() const
{return nominal_def;}


bool snapshot::isMarked() const
{return marked;}


void snapshot::setMarked()
{marked = true;}



bool snapshot::load(simVars& sim, simPars pars, int nominalDef) const
{
	cerr << "loadSnapshot is not implemented" << endl;
	return false;

	// tau_g zárját el kell távolítani
}

//get the content of the conti file
vector<snapshot> getSnapshots(const simPars& pars)
{
	ifstream conti(fnameWp("conti",pars));
	if (!conti)
		cerr << fnameWp("conti",pars) << " was not found" << endl;

	vector<snapshot> snapshots;
	for (snapshot toadd;
		conti >> toadd;)
		snapshots.push_back(toadd);

	return snapshots;
}

//write the content to the conti file
bool writeContiFile(const simPars& pars, const vector<snapshot>& conti_content)
{
	ofstream conti(fnameWp("conti",pars));
	if (!conti)
	{
		cerr << "Cannot (over)write " << fnameWp("conti",pars) << endl;
		return false;
	}

	for (auto x : conti_content)
		conti << x << endl;

	return true;
}
