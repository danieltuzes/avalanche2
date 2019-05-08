// snapshot.h : defines snapshot savings and loadings
//
#define SNAPSHOT_VERSION 0.01

#pragma once
#include "stdafx.h"
#include "make_simulation.h"
using namespace std;

class snapshot
{
protected:
	bool marked;
	deftype sG;
	deftype sGn;
	deftype av_size;
	deftype av_size_n;
	deftype tau_ext;
	int rnd_state;
	int nominal_def;
	const simVars * m_sim;

public:
	snapshot();
	snapshot(deftype sG, deftype sGn, deftype av_size, deftype av_size_n, deftype tau_ext, int rnd_state, int nominal_def);
	snapshot(const simVars& sim, int nominal_def);

	deftype getDeform() const;
	int getNominal_def() const;
	bool isMarked() const;
	void setMarked();
	
	friend bool operator==(const snapshot &lhs, const snapshot &rhs);
	friend bool operator!=(const snapshot &lhs, const snapshot &rhs);
	friend bool operator<(const snapshot &lhs, const snapshot &rhs);
	friend bool operator>(const snapshot &lhs, const snapshot &rhs);
	friend bool operator<=(const snapshot &lhs, const snapshot &rhs);
	friend bool operator>=(const snapshot &lhs, const snapshot &rhs);
	friend istream& operator >> (istream& ist, snapshot& sh);
	friend ostream& operator << (ostream& ost, snapshot sh);

	bool make(const simPars& pars) const;
	bool mark(const simPars& pars) const;
	bool check(const simPars& pars) const;
	bool load(simVars& sim, simPars pars, int nominalDef) const;
};

bool operator==(const snapshot &lhs, const snapshot &rhs);
bool operator!=(const snapshot &lhs, const snapshot &rhs);
bool operator<(const snapshot &lhs, const snapshot &rhs);
bool operator>(const snapshot &lhs, const snapshot &rhs);
bool operator<=(const snapshot &lhs, const snapshot &rhs);
bool operator>=(const snapshot &lhs, const snapshot &rhs);

bool writeContiFile(const simPars& pars, const vector<snapshot>& snapshots); //write the content to the conti file
bool delMarkedSnapshots(const simPars& pars);
vector<snapshot> getSnapshots(const simPars& pars); //get the content of the conti file

