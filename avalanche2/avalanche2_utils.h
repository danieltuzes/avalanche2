// avalanche2_utils.h : stores user defined utilities for avalanche2.cpp
//

#define AVALANCHE2_UTILS_VERSION 0.07
// 0.07: tau_g can contain time
// 0.06: tau_g contains the start positions of the avalanches
// 0.05: makeSN

#define AVALANCHE2_VERSION 1.0
// 1.0: thermal activation implemented
// 0.05: removed some templates
// 0.04: makeSN
// with flowStressReGeneration
// with springConstant
// stop simulation if the external stress reach inf

#define MAKE_DEFORM_PARALLEL 0.01
#define MAKE_DEFORM_SERIAL 1.0
// 1.0: thermal activation implemented
// 0.09: print out the start position of the avalanche
// 0.08: some templates are removed, it doesn't fasten up the simulations, but easier to read
// 0.07b: prepare the program for full_reconstruct
// 0.06: makeSN
// 0.05b: restricted randKernel implemented
// 0.05: randKernel implemented
// 0.04: DGamma::limited implemented
// 0.03: flowStressReGeneration::no case handled
// 0.02: create pNd after every deformation: no need to store the data, creates pNd info in the last avalanche as well; later, vector<pNd> placeNdeform can be removed completly


#pragma once
#include "../gen_utils/precision.h"
#include "../gen_utils/gen_version.h"
#include "../gen_utils/gen_param.h"
#include "../gen_utils/gen_math.h"
#include "../gen_utils/gen_fmanip.h"
#include "../gen_utils/teestream.h"
#include "../gen_utils/gen_interrupt.h"
#include "rnd_gen.h"
#include "../gen_utils/gen_fnameWp.h"


#pragma region defined enum types

class activity
{
public:
	enum type
	{
		import, //load a snapshot
		deform, //make a deformation and then create a snapshot
		deform_newKernel //load a new kernel
	};
	activity(int nominalDef, type dol);
	activity(int nominalDef, string fname, type dol);
	type MyType;
	int nominalDef;
	string kernelName;
	double tau_ext;
};

ostream& operator<< (ostream& o, activity a);


class infoMode
{
public:
	enum type
	{
	minimal, //show minimal informations
	progressBar, //show a progress bar, when it reach 100%, it creates a snapshot and returns to 0 again
	everyStep, //print infos from every deformations
	undefined //default constructor
	};
	infoMode();
	infoMode(type mode);
	type MyType;
};


bool operator==(infoMode base, infoMode::type compare);
bool operator!=(infoMode base, infoMode::type compare);

ostream& operator<< (ostream& ostr, infoMode info);
istream& operator>> (istream& istr, infoMode& info);


class DGamma
{
public:
	enum type
	{
		definit, // value is given by a number
		limited, // deformation is choosen to eliminate local stress but not bigger than the given value
		max // deformation is choosen to eliminate local stress
	};
	type MyType;
	deftype val;
	DGamma();
	DGamma(type MyType);
	DGamma(deftype val);
	DGamma(type MyType, deftype val);
};

ostream& operator<< (ostream& o, const DGamma& DG);
istream& operator>> (istream& i, DGamma& DG);

class leftFlowStress
{
public:
	enum type
	{
		independent, // the left and right flowstresses are independent, and basicly are from the same distribution
		infinity, // pratically this means that the sample can only deform to the right direction
		same, // the left and right flowstress values are the same as suggested by Michael
		undefined //for default constructor
	};
	type MyType;
	leftFlowStress();
	leftFlowStress(type MyType);
};

ostream& operator<< (ostream& o, const leftFlowStress& lfs);
istream& operator>> (istream& i, leftFlowStress& lfs);


class flowStressReGeneration
{
public:
	enum type
	{
		independent, // the new flow stress is generated from the very same distribution
		infinity,    // after a deformation occurs, the flow stress became infinitely large
		no,          // no flow stress regeneration, the flow stress remain the same
		linear,      // flowstress is increased by a factor of DsG * measure
		exponential, // flowstress is multiplied by a factor of DsG * measure
		undefined    // for default constructor
		
	};
	type MyType;
	double measure;
	flowStressReGeneration();
	flowStressReGeneration(type MyType);
};

ostream& operator<< (ostream& o, const flowStressReGeneration& lfs);
istream& operator>> (istream& i, flowStressReGeneration& lfs);


#pragma endregion

#pragma region parameter settings


class simPars // parameters needed during the simulations
{
public:
	simPars(int s,
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
		double beta);
	int s;                        // system size
	int seed;                     // seed number for random number generation
	string fsd;                   // flow stress distribution
	leftFlowStress lfs;           // allow avalanches only in one direction
	DGamma DG;                    // value to calculate deformation
	bool createTau_g;             // if the program should create tau_g file
	int fg;                       // filename grouping
	string fp;                    // filename prefix
	bool makeSN;                  // the program should make snapshots after a successful deformation
	int tbs;                      // time between snapshots
	int mrt;                      // max runtime
	int mna;                      // max number of avalanches
	bool pNd;                     // place and deform
	flowStressReGeneration fsrg;  // flow stress regeneration
	bool dr;                      // to increase the external stress after a deformation imitating disruption
	double sc;             // to choose the decreasement of the external stress after a deformation
	double maxLocG;        // to define the maximum amount of local strain, program stop simulation if the value is reached
	bool mdp;                     // make deformation parallel
	string sfk;                   // stress field kernel
	int randKernel;               // randomize the kernel except the (0,0) point
	infoMode info;                // setting verbosity
	bool searchMinPoint;                 // if program should search or load minpoint. For further development to make it possible to recreate whole simulations
	double tau_extMin; // when no more cell is activated, tau_ext = max(tau_p[minPoint] - local_s + EPSILON, tau_extMin). If sc is set, tau_ext may decrease below tau_extMin
	bool thermalActivation; // sets the simulation mode to thermal activation, i.e. at a given stress value, the propability that a cell will be deformed is proportional to e^-beta*t_th 
	double beta; // sets the value beta 

};

#pragma endregion

void intro();

bool fillUpToDoList(param<int> DsGm, param<string> sFn, param<bool> ls, vector<activity>& loa);

string fnameWp(string fname, simPars pars);
string fnameWp(string fname, simPars pars, int nominalDef);

bool makeSimulation(simPars pars, const vector<activity>& loa, const vector<double>& extraSnaps);

int nofdigit(int number);