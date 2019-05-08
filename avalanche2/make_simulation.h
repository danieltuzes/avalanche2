// make_simulation.h : creates avalanche with the given parameters
//

//uncomment the following line to create code that generates code with cuda
//#define IS_GPU 

#ifdef IS_GPU
#define MAKE_SIMUATION_GPU 0.01
#else
#define MAKE_SIMULATION 0.12

// 0.12: MonteCarloChoseNwait return value depends on tau_ext + tau_l
// 0.11: simVars contains Z, placeLock_tau_g uses printTau_g
// 0.1: simVars contains time and threshold, printTau_g get the start position from member value instead of call argument
// 0.09: simVars contains avStartPos, tau_g file can contain the start position of avalanches
// 0.08: prepare the program for full_reconstruct
// 0.07: makeSN
// 0.06b: restricted randKernel implemented
// 0.06 randKernel implemented
// 0.05 maxLocG implemented
// 0.04 DGamma::limited implemented
#endif
#define TAUG_LOCK_TEXT " # exited normally"
#pragma once

#include "avalanche2_utils.h"
#include "../gen_utils/pNd.h"
using namespace std;

enum direction
{
	none,
	left,
	right,
};

ostream& operator<< (ostream& o, const direction& dir);

deftype calcDeform(int minPoint, deftype local_s, deftype tau_ext, deftype stress_value, DGamma DG, direction D);


class simVars
{
private:
	int size; //the size of the simulation
	int mask;
	int power;
public:
	simVars(simPars pars);

	block_data<deftype> Gamma;
	block_data<deftype> tau_l, tau_n, tau_p, sf;
	deftype sG, sGn, av_size, av_size_n;
	int avStartPos = -1;
	deftype tau_ext;
	double simTime = 0;
	int nofAvs = 0;
	double threshold;
	double Z = 0; // canonical partition function
	double Zstartval;
	universal_rnd rand_gen;
	void initialize(leftFlowStress lfs); //fills up the maps with values
#ifdef IS_GPU
	bool makeDeformSerial_GPU(const simPars& pars, int nominalDef, bool writeToTau_g);
	bool makeDeformParallel_GPU(const simPars& pars, int nominalDef);
#else
	bool makeDeformSerial(const simPars& pars, int nominalDef, vector<double>& extraSnaps, bool printToTau_g);
	bool makeDeformParallel(const simPars& pars, int nominalDef);
#endif
	// bool reconstr2Stress(const simPars& pars, const vector<double>& extstress);
	int findMin(bool thermalActivation, double beta);
	int refresh(int minPoint, deftype deform, int randKernel, bool searchMinPoint, bool thermalActivation, double beta);
	direction MonteCarloChoseNwait(double beta, int& minPoint, double& threshold); // chose a cell for thermal activation and increase the time

	void printTau_g(ofstream& tau_g, bool thermalActivation, string note = "") const;

	int getSize() const;
	int getPower() const;
	int getLinPos(int x, int y) const;
	int getX(int pos) const;
	int getY(int pos) const;

	deftype getDeform() const; //sG - sGn + av_size - av_size_n
};


direction willYield(const simVars& sim, int minPoint, deftype local_s);

void printInfHeader(const simPars& pars);

void printInf(int minPoint, const simVars& sim, const simPars& pars);

deftype getLargestDeformInTau_g(const simPars pars);

bool makeSimulation(simPars pars, vector<activity>& loa, const vector<double>& extraSnaps, bool searchMinPoint);

bool placeLock_tau_g(const simVars& sim, const simPars& pars, string note = TAUG_LOCK_TEXT);