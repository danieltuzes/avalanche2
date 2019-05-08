// make_deform_parallel.h : makes deformation parallel until a given nominal def
//
// version number is in avalanche2_utils.h

#include "stdafx.h"
#include "make_simulation.h"
#include "snapshot.h"

using namespace std;

#ifndef IS_GPU
bool simVars::makeDeformParallel(const simPars& pars, int nominalDef)
{
	cerr << "Make deform parallel is empty" << endl;
	return false;

}
#endif