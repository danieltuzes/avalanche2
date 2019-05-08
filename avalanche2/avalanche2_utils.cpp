// avalanche2_utils.cpp : stores user defined utilities for avalanche2.cpp
//


#include "stdafx.h"
#include "avalanche2_utils.h"
#include "snapshot.h"
#include "../gen_utils/tau_g.h"


#pragma region defined enum type activity

activity::activity(int nominalDef, activity::type dol) : MyType(dol), nominalDef(nominalDef) {}
activity::activity(int nominalDef, string fname, activity::type dol) : MyType(dol), kernelName(fname) {}

ostream& operator<< (ostream& o, activity a)
{
	if (a.MyType == activity::deform)
		o <<  "deform until nominal def ";
	else if (a.MyType == activity::import)
		o << "import snapshot from nominal def ";
	else if (a.MyType == activity::deform_newKernel)
		o << "load new kernel " << a.kernelName << " after reached " << a.nominalDef;
	else
		cerr << "Unsopported activity type" << endl;

	o << a.nominalDef;
	return o;
}

#pragma endregion
#pragma region defined enum type infoMode

infoMode::infoMode() : MyType(undefined) {}
infoMode::infoMode(infoMode::type MyType) : MyType(MyType) {}


ostream& operator<< (ostream& ostr, infoMode info)
{
	if (info.MyType == infoMode::minimal)
		ostr << "minimal";
	else if (info.MyType == infoMode::progressBar)
		ostr << "progressBar";
	else if (info.MyType == infoMode::everyStep)
		ostr << "everyStep";
	
	return ostr;
}

istream& operator>> (istream& istr, infoMode& info)
{
	string name;
	istr >> name;
	if (name == "minimal")
		info.MyType = infoMode::minimal;
	else if (name == "progressBar")
		info.MyType = infoMode::progressBar;
	else if (name == "everyStep")
		info.MyType = infoMode::everyStep;
	
	return istr;
}

istream& operator>> (istream& istr, infoMode& info);

bool operator==(infoMode base, infoMode::type compare)
{
	return base.MyType == compare;
}

bool operator!=(infoMode base, infoMode::type compare)
{
	return !(base==compare);
}

template <> string param<infoMode>::getValType() {return "infoMode";}

#pragma endregion
#pragma region defined enum type DGamma

DGamma::DGamma() {}
DGamma::DGamma(deftype val) : MyType(definit), val(val) {}
DGamma::DGamma(type MyType) : MyType(MyType), val(0) {}
DGamma::DGamma(type MyType, deftype val) : MyType(MyType), val(val) {}

ostream& operator<< (ostream& o, const DGamma& DG)
{
	if (DG.MyType == DGamma::max)
		o << "max";
	else if (DG.MyType == DGamma::limited)
		o << "limited " << DG.val;
	else if (DG.MyType == DGamma::definit)
		o << "definit " << DG.val;
	return o;
}

istream& operator>> (istream& i, DGamma& DG)
{
	string tmp;
	i >> tmp;
	if (tmp.compare("max") == 0)
		DG = DGamma(DGamma::max);
	else if (tmp.substr(0,string("limited_").size()).compare("limited_") == 0)
	{
		stringstream tmp_stream(tmp.substr(string("limited_").size(), string::npos));
		deftype defval;
		tmp_stream >> defval;
		DG = DGamma(DGamma::limited,defval);
	}
	else // tmp treated as a number
	{
		stringstream tmp_stream(tmp);
		deftype defval;
		tmp_stream >> defval;
		DG = DGamma(defval);
	}

	return i;
}

template <> string param<DGamma>::getValType() {return "DGamma";}

#pragma endregion
#pragma region defined enum type leftFlowstress

leftFlowStress::leftFlowStress() : MyType(undefined) {}
leftFlowStress::leftFlowStress(type MyType) : MyType(MyType) {}

ostream& operator<< (ostream& o, const leftFlowStress& lfs)
{
	if (lfs.MyType == leftFlowStress::independent)
		o << "independent";
	else if (lfs.MyType == leftFlowStress::infinity)
		o << "infinity";
	else if (lfs.MyType == leftFlowStress::same)
		o << "same";
	else
		cerr << "unsupported leftFlowStress type";
	return o;
}

istream& operator>> (istream& i, leftFlowStress& lfs)
{
	string tmp;
	i >> tmp;
	if (tmp.compare("independent") == 0)
		lfs = leftFlowStress::independent;
	else if (tmp.compare("infinity") == 0)
		lfs = leftFlowStress::infinity;
	else if (tmp.compare("same") == 0)
		lfs = leftFlowStress::same;
	else
		cerr << "unsupported leftFlowStress type";
	return i;
}

template <> string param<leftFlowStress>::getValType() { return "leftFlowStress"; }

#pragma endregion
#pragma region defined enum type flowStressReGeneration

flowStressReGeneration::flowStressReGeneration() : MyType(undefined), measure(0) {}
flowStressReGeneration::flowStressReGeneration(type MyType) : MyType(MyType), measure(0) {}

ostream& operator<< (ostream& o, const flowStressReGeneration& lfs)
{
	if (lfs.MyType == flowStressReGeneration::independent)
		o << "independent";
	else if (lfs.MyType == flowStressReGeneration::infinity)
		o << "infinity";
	else if (lfs.MyType == flowStressReGeneration::no)
		o << "no";
	else if (lfs.MyType == flowStressReGeneration::exponential)
		o << "exponential_" << lfs.measure;
	else if (lfs.MyType == flowStressReGeneration::linear)
		o << "linear_" << lfs.measure;
	else
		cerr << "unsupported flowStressReGeneration type";
	return o;
}


istream& operator>> (istream& i, flowStressReGeneration& fsrg)
{
	string tmp;
	i >> tmp;
	if (tmp.compare("independent") == 0)
		fsrg = flowStressReGeneration::independent;
	else if (tmp.compare("infinity") == 0)
		fsrg = flowStressReGeneration::infinity;
	else if (tmp.compare("no") == 0)
		fsrg = flowStressReGeneration::no;
	else if (tmp.compare(0, strlen("linear"), "linear") == 0)
	{
		fsrg.MyType = flowStressReGeneration::linear;
		fsrg.measure = stod(tmp.substr(strlen("linear_")));
	}
	else if (tmp.compare(0, strlen("exponential"), "exponential") == 0)
	{
		fsrg.MyType = flowStressReGeneration::exponential;
		fsrg.measure = stod(tmp.substr(strlen("exponential_")));
	}
	else
		cerr << "unsupported leftFlowStress type";
	return i;
}

template <> string param<flowStressReGeneration>::getValType() { return "flowStressReGeneration"; }

#pragma endregion

#pragma region parameter settings

simPars::simPars(int s,
	int seed,
	string fsd,
	leftFlowStress lfs,
	DGamma DG,
	bool createTau_g,
	int fg,
	string fp,
	bool makeSN,
	int tbs,
	int mrt,
	int mna,
	bool pNd,
	flowStressReGeneration fsrg,
	bool dr,
	double sc,
	double maxLocG,
	bool mdp,
	string sfk,
	int randKernel,
	infoMode info,
	bool searchMinPoint,
	double tau_extMin,
	bool thermalActivation,
	double beta)
: s(s), seed(seed), fsd(fsd), lfs(lfs), DG(DG), createTau_g(createTau_g), fg(fg), fp(fp), makeSN(makeSN), tbs(tbs), mrt(mrt), mna(mna), pNd(pNd),
fsrg(fsrg), dr(dr), sc(sc), maxLocG(maxLocG), mdp(mdp), sfk(sfk), randKernel(randKernel), info(info), searchMinPoint(searchMinPoint),
tau_extMin(tau_extMin), thermalActivation(thermalActivation), beta(beta) {}

#pragma endregion

bool fillUpToDoList(param<int> DsGm, param<string> sFn, param<bool> ls, vector<activity>& loa)
{
	//fill loa with deformations
	ifstream sF(sFn());
	if (!sF)
	{
		cerr << "Cannot open " << sFn() << endl;
		if (sFn.setByUsr())
		{
			cerr << sFn.get_ext_name() << " set manually so file must exists" << endl;
			return false;
		}
	}

	vector<string> sw;
	for (string aWord; comment_skipped(sF) >> aWord; sw.push_back(aWord));
	for (auto it = sw.begin(); it != sw.end(); ++it)
	{
		if (next(it) == sw.end() || (next(it) != sw.end() && it->compare("load") != 0)) // if it is an integer
		{
			istringstream num(*it);
			int theNum;
			num >> theNum;
			loa.emplace_back(theNum, activity::type::deform);
		}
		else
		{
			istringstream num(*it);
			int theNum;
			num >> theNum;
			++it;
			if (!ifstream(*it))
			{
				cerr << "Cannot open " << *it << endl;
				return false;
			}
			loa.emplace_back(theNum, *it, activity::type::deform_newKernel);
		}
	}

	//fill up loa with loading snapshots
	vector<snapshot> snapshots;
	
	//if (ls())
	//{
	//	snapshots = getSnapshots(pars);
	//	for (auto x : snapshots)
	//	{
	//		if (x.check(pars))
	//			loa.emplace_back(x.getNominal_def(), activity::type::import);
	//		else
	//			x.mark(pars);
	//	}
	//}
	
		
	if (!loa.empty())
	{
		sort(loa.begin(), loa.end(), [](activity a, activity b) {return (a.nominalDef < b.nominalDef); });

		//truncate tau_g if necessary

		//if (!snapshots.empty() && snapshots.back().getDeform() != tau_g_line(lastNonEmptyLine(fnameWp("tau_g", pars))).DsG())
		//	cerr << "tau_g needs to be truncated but not implemented" << endl;


		//don't load snapshots if there's no snapshot to make
		//don't load or make snapshot if its nominalDef > DsGm()

		for (auto x = prev(loa.end()); x->MyType == activity::import || x->nominalDef > DsGm();)
		{
			--x;
			loa.pop_back();
		}
	}

	

	//if last snapshot's nominal def must be the DsGm
	if (loa.empty() || (loa.back().MyType != activity::deform && loa.back().nominalDef != DsGm()))
	{
		loa.emplace_back(DsGm(),activity::deform);
		cout << "Snapshot with deformation value " << DsGm << " added to the list of activity" << endl;
	}

	return true;
}


void intro()
{
	stringstream version;
	version << "VERSION:\n"
		<< "    " << "avalanche2:            " << XSTR(AVALANCHE2_VERSION) << "\n"
		<< "    " << "avalanche2_utils:      " << XSTR(AVALANCHE2_UTILS_VERSION) << "\n"
		<< "    " << "make_simulation:       " << XSTR(MAKE_SIMULATION) << "\n"
		<< "    " << "make_deform_serial:    " << XSTR(MAKE_DEFORM_SERIAL) << "\n"
		<< "    " << "make_deform_parallel:  " << XSTR(MAKE_DEFORM_PARALLEL) << "\n"
		<< "    " << "rnd_gen:               " << XSTR(RND_GEN_VERSION) << "\n"
		<< "    " << "snapshot:              " << XSTR(SNAPSHOT_VERSION) << "\n"
		<< "    " << "gen_siomanip:          " << XSTR(GEN_SIOMANIP_VERSION) << "\n"
		<< "    " << "gen_fmanip:            " << XSTR(GEN_FMANIP_VERSION) << "\n"
		<< "    " << "gen_version:           " << XSTR(GEN_VERSION_VERSION) << "\n"
		<< "    " << "compiler:              " << XSTR(COMPILER_VERSION) << "\n"
		<< "    " << "computer:              " << XSTR(MACHINE_INFO) << "\n";

	intro("This is an avalanche simulation program. ", version.str());
}

string fnameWp(string fname, simPars pars) { return fnameWp(fname, pars.s, pars.seed, pars.fg, pars.fp, ""); }


string fnameWp(string fname, simPars pars, int nominalDef) { return fnameWp(fname, pars.s, pars.seed, pars.fg, pars.fp, "_" + to_string(nominalDef)); }