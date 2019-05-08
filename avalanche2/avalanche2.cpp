// avalanche2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "avalanche2_utils.h"

int main(int argc, char *argv[])
{
	time_t time_start = time(NULL);

#pragma region Defines parameters and set the values of them
	paramcontainer p,p_temp;
	param<int> s("sysSize",p,rq::must); //system size
	param<int> st("seedStart",p,rq::must); //seed start
	param<int> se("seedEnd",p,rq::must); //seed end
	param<int> DsGm("DsGm",p,rq::must); //Delta sum Gamma max
	param<bool> ls("loadSnapshot",p,rq::must); // loadsnapshot
	param<string> fsd("flowStressDistribution",p,rq::must); // the positive and negative flow stress distribution
	param<leftFlowStress> lfs("leftFlowStress", p, rq::optional, leftFlowStress::independent); // the formal oneway paramter
	param<DGamma> DG("DGamma",p,rq::optional,1); // value to calculate deformation
	param<string> extraSnapsFn("extraSnapsFn", p, rq::optional, "extrasnaps.ini"); // at which external stress values should the program make additional snapshots
	param<bool> searchMinPoint("searchMinPoint", p, rq::optional, true); // if the program should search for the minpoint or work from existing pNd file
	param<bool> createTau_g("createTau_g", p, rq::optional, true); // if the program should create tau_g files
	param<int> fg("fileGrouping",p,rq::optional,1000); //file grouping
	param<string> fp("fnamePrefix",p,rq::optional,""); //fname prefix
	param<string> sFn("snapshotFname",p,rq::optional,"snapshot.ini"); //snapshotfname
	param<bool> makeSN("makeSN", p, rq::optional, true); // if the program should make snapshots after a successful deformation
	param<int> tbs("timeBetweenSnapshots", p, rq::optional, 3600); //time between snapshots
	param<int> mrt("maxRunTime", p, rq::optional, numeric_limits<int>::max()); // maximal running time
	param<int> mna("maxNumberofAvalanches", p, rq::optional, numeric_limits<int>::max()); // the maximal number of avalanches
	param<bool> pNd("placeNdeform",p,rq::optional,false); //place and deform, create or not data file for deformation order
	param<flowStressReGeneration> fsrg("flowStressReGeneration", p, rq::optional, flowStressReGeneration::independent); //regenerate the flow stress after a deformation occures
	param<bool> dr("disruption", p, rq::optional, false); // to increase the external stress after a deformation imitating disruption
	param<double> sc("springConstant", p, rq::optional, 0); // to choose the decreasement of the external stress after a deformation
	param<double> maxLocG("maxLocG", p, rq::optional, 0); // to define the maximum amount of local strain, program stop simulation if the value is reached
	param<bool> mdp("makeDeformParallel", p, rq::optional, false); //make the deformation parallel or serial
	param<string> sfk("stressFieldKernelFname", p, rq::optional); // the filename of the stress field that generates the local stress field
	param<int> randKernel("randKernel", p, rq::optional, 0); // 0: false; 1: yes; 2: restricted, i.e. G(r) = G(-r)
	param<double> tau_extMin("tau_extMin", p, rq::optional, 0); // when no more cell is activated, tau_ext = max(tau_p[minPoint] - local_s + EPSILON, tau_extMin). If sc is set, tau_ext may decrease below tau_extMin
	param<bool> thermalActivation("thermalActivation", p, rq::optional, false); // sets the simulation mode to thermal activation, i.e. at a given stress value, the propability that a cell will be deformed is proportional to e^-beta*t_th 
	param<double> beta("beta", p, rq::optional); // sets the value beta 
	param<infoMode> info("infoMode",p,rq::optional,infoMode::minimal); //setting verbosity

	param<string> paramFname("paramFname",p_temp,rq::optional,"av_param.ini",argc,argv); //the name of the ini file

	intro();
	p.setFromIf(paramFname());
	if ((!paramFname.setByUsr() && p.setFromIf(argc,argv) * 2 != argc - 1) || (paramFname.setByUsr() && p.setFromIf(argc,argv) * 2 != argc - 3))
	{
		cerr << "Unknow parameter or missing value was found in program call, program terminates." << endl;
		p.description();
		return 0;
	}

	if (!p.summerize())
	{
		p.description();
		return 0;
	}
#pragma endregion


#pragma region copy cerr and cout to history; history is a file that contains previous series of simulations
	string logFname = fnameWp("log",s(),fg(),fp(),st(),se(),""); //the file name of history
	ofstream his(logFname,ofstream::app);
	if (!his) {
		cerr << "Error: unable to append history file " << logFname << ". Program terminates." << endl;
		return 0;
	}
	teebuf tee_cout(cout.rdbuf(), his.rdbuf());
	cout.rdbuf(&tee_cout);
	teebuf tee_err(cerr.rdbuf(), his.rdbuf());
	cerr.rdbuf(&tee_err);
#pragma endregion

#pragma region Checking parameters
	int i; //check the simulation size
	for (i = 16384; i != 1; i /= 2)
	{
		if (i == s)
			break;
	}
	if (i == 1)
	{
		cerr << s.get_ext_name() << " must be a power of two, bigger than 1 but not bigger than 16384. Program terminates." << endl;
		return 1;
	}

	if (se <= st) // check the seed start and seed end values
	{
		cerr << se.get_ext_name() << " must be greater than " << st.get_ext_name() << "= " << st << "  Program terminates." << endl;
		return 1;
	}

	universal_rnd rnd_test(fsd, st); // check if the random number generation is working correctly
	if (!rnd_test.isGood())
	{
		cerr << "Random number generator is not in a good state. Program terminates." << endl;
		return 1;
	}

	if (pNd && !searchMinPoint())
	{
		cerr << "This version of the program can only create pNd data if the minPoints are searched." << endl;
		return 1;
	}

	vector<activity> loa; //list of activity, defines the snapshot loading and creating
	if (!fillUpToDoList(DsGm, sFn, ls, loa))
	{
		cerr << "Unable to fill up toDoList." << endl;
		return false;
	}

	vector<double> extraSnaps;
	ifstream extraSnapsF(extraSnapsFn());
	if (!extraSnapsF && extraSnapsFn.setByUsr())
	{
		cerr << "Cannot open " << extraSnapsFn() << ". Program terminates." << endl;
		return false;
	}
	for (double tmp; comment_skipped(extraSnapsF) >> tmp; extraSnaps.push_back(tmp));

	if (thermalActivation && !beta.setByUsr())
	{
		cerr << "If thermalActivation is true, but beta is not set. Program terminates." << endl;
		return 1;
	}

	if (beta.setByUsr() && !thermalActivation.setByUsr())
	{
		cout << "beta is set, but thermalActivation was not set accordingly; program set thermalActivation to true" << endl;
		thermalActivation.modify(true);
	}

	if (beta.setByUsr() && thermalActivation.setByUsr() && !thermalActivation)
	{
		cerr << "beta is set and thermalActivation is to false. If beta is set, thermalActivation must be true. Program terminates." << endl;
		return 1;
	}


	if (!ifstream(sfk))
	{
		cerr << sfk << " not found" << endl;
		return false;
	}

#pragma endregion

	simPars pars(s, st, fsd, lfs, DG, createTau_g, fg, fp, makeSN, tbs, mrt, mna, pNd, fsrg, dr, sc, maxLocG, mdp, sfk, randKernel, info, searchMinPoint, tau_extMin, thermalActivation, beta);



	for (; pars.seed != se; ++pars.seed)
	{
		cout << "\nSeed: " << pars.seed << endl;
		if (interrupt("AvBegin",pars.seed)) {cerr << "Program terminates." << endl; return false;}
		
		time_t time_sim_start = time(NULL);
			
		////////////////////////////////////////////////////////////
		//// Simulation generating
		////////////////////////////////////////////////////////////

#pragma region tee cout and cerr to log file
		string logFname = fnameWp("log",pars);

		ofstream his(logFname,ofstream::app);
		if (!his)
		{
			cerr << "Error: unable to append log file " << logFname << endl;
			return false;
		}
		streambuf * coutbackup = cout.rdbuf();
		teebuf tee_cout(cout.rdbuf(), his.rdbuf());
		cout.rdbuf(&tee_cout);

		streambuf * cerrbackup = cerr.rdbuf();
		teebuf tee_err(cerr.rdbuf(), his.rdbuf());
		cerr.rdbuf(&tee_err);

		cout << currentDateTime() << " Start a simulation with\n";
		p.summerize(cout);

#pragma endregion

		bool result = makeSimulation(pars, loa, extraSnaps);

		if (result)
			cout << "Simulation with seed " << pars.seed << " was successful. "
				 << "Finished: " << currentDateTime() << endl
				 << "\tTime: " << time(NULL) - time_sim_start << endl;
		else
		{
			cerr << "Simulation with seed " << pars.seed << " was failed. "
				 << "Finished: " << currentDateTime() << endl
				 << "\tTime: " << time(NULL) - time_sim_start << endl;
		}

		cout.rdbuf(coutbackup);
		cerr.rdbuf(cerrbackup);
	}
	
	
	time_t time_end = time(NULL);
	cout << "\nBatch finished: " << currentDateTime()
		 << "\tTime: " << time_end - time_start << "s" << endl;


	return 0;
}