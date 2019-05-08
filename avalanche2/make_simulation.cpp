
#include "stdafx.h"
#include "make_simulation.h"
#include "snapshot.h"
#include "make_deform_parallel.h"
#include "../gen_utils/tau_g.h"

ostream& operator<< (ostream& o,const direction& dir)
{
	if (dir == direction::left)
		o << "<";
	else if (dir == direction::right)
		o << ">";
	else
		o << "-";
	return o;
}

direction willYield(const simVars& sim, int minPoint, deftype local_s)
{
	if (sim.tau_p[minPoint] < sim.tau_ext + local_s)
	{
		if (sim.tau_p[minPoint] - sim.tau_ext - local_s < sim.tau_n[minPoint] + sim.tau_ext + local_s)
			return direction::right;
		else return direction::left;
	}
	else if (sim.tau_n[minPoint] < -sim.tau_ext - local_s)
		return direction::left;
	
	return direction::none;
}

direction simVars::MonteCarloChoseNwait(double beta, int& minPoint, double& threshold) // chose a cell for thermal activation and increase the time
{
	uniform_real_distribution<double> ud(0, 1);
	double rndNum = ud(rand_gen.getEngine());
	simTime += -log(rndNum) / Z;
	rndNum *= Z;
	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
		{
			int linPos = getLinPos(x, y);
			double tau_th = tau_p[linPos] - tau_l[linPos] - tau_ext;
			rndNum -= nextafter(exp(-beta*tau_th), numeric_limits<double>::infinity());
			if (rndNum < 0)
			{
				minPoint = linPos;
				threshold = tau_th;
				if (tau_l[linPos] + tau_ext > 0)
					return direction::right;
				else
					return direction::left;
			}
			tau_th = tau_n[linPos] + tau_l[linPos] + tau_ext;
			rndNum -= nextafter(exp(-beta*tau_th), numeric_limits<double>::infinity());
			if (rndNum < 0)
			{
				minPoint = linPos;
				threshold = tau_th;
				if (tau_l[linPos] + tau_ext > 0)
					return direction::right;
				else
					return direction::left;
			}
		}
	cerr << "Didn't get to Z. Program temrinates." << endl;
	exit(1);
}

void printInfHeader(const simPars& pars)
{
	if (pars.info != infoMode::everyStep)
		return;

	cout << "tau_ext\t"
		<< "sG\t" << "sGn\t"
		<< "min.x\t" << "min.y\t"
		<< "min.dir\t"
		<< "tau_l[min]\t"
		<< "tau_n[min]\t"
		<< "tau_p[min]" << endl;
}

void printInf(int minPoint, const simVars& sim, const simPars& pars)
{
	if (pars.info != infoMode::everyStep)
		return;

	if (willYield(sim, minPoint, sim.tau_l[minPoint]) == direction::none) cout << endl;

	int x = sim.getX(minPoint);
	int y = sim.getY(minPoint);

	cout << sim.tau_ext << "\t" 
		 << sim.sG << "\t" << sim.sGn << "\t" 
		 << x << "\t" 
		 << y << "\t"
		 << willYield(sim, minPoint, sim.tau_l[minPoint]) << "\t"
		 << sim.tau_l[minPoint] << "\t"
		 << sim.tau_n[minPoint] << "\t"
		 << sim.tau_p[minPoint] << "\t" << endl;
}

simVars::simVars(simPars pars) : size(pars.s), mask(size-1), power(::getPower(size)), Gamma(size,0.), tau_l(size,0.), tau_n(size,0.), tau_p(size,0.), sf(size,0.), rand_gen(pars.fsd,pars.seed) {}

deftype simVars::getDeform() const {return (sG - sGn);}


//fills up the maps with useful values
void simVars::initialize(leftFlowStress lfs)
{
	for (int i=0; i<size; ++i)
		for (int j=0; j<size; ++j)
		{
			tau_p(i, j) = rand_gen();
			if (lfs.MyType == leftFlowStress::infinity)
			{
				rand_gen(); // to make the random generator consitent
				tau_n(i,j) = numeric_limits<double>::infinity();
			}
			else if (lfs.MyType == leftFlowStress::independent)
				tau_n(i, j) = rand_gen();
			else if (lfs.MyType == leftFlowStress::same)
				tau_n(i, j) = tau_p(i, j);
			tau_l(i, j) = 0;
			Gamma(i, j) = 0;
		};
	
	tau_ext = 0;
	sG = 0;
	sGn = 0;
	av_size = 0;
	av_size_n = 0;
}

int simVars::getSize() const{ return size; }

int simVars::getPower() const{ return power; }

int simVars::getLinPos(int x, int y) const{ return (x << power) + y; }

int simVars::getX(int pos) const{ return ::getX(pos, power); }
int simVars::getY(int pos) const{ return ::getY(pos, mask); }


bool init_tau_g(const simVars& sim, const simPars& pars)
{
	string fname(fnameWp("tau_g", pars));

	ofstream tau_g(fname);
	if (!tau_g)
	{
		cerr << "Can't initialize tau_g file" << fname << endl;
		return false;
	}
	tau_g << "# tau_ext(1)\t"
		<< "sG(2)\t"
		<< "sGn(3)\t"
		<< "av_size(4)\t"
		<< "av_size_n(5)\t"
		<< "rand_gen.getTimesAdvanced()(6)\t"
		<< "start_posx(7)\t"
		<< "start_posy(8)";
	if (pars.thermalActivation)
		tau_g << "\t"
		<< "time(9)\t"
		<< "threshold(10)\t"
		<< "Z(11)";
	tau_g << endl;
	tau_g.precision(17);
	tau_g_line firstLine = tau_g_line(sim.tau_ext,sim.sG,sim.sGn,sim.av_size,sim.av_size_n,sim.rand_gen.getTimesAdvanced());

	tau_g << firstLine << endl;
	return true;
}


void simVars::printTau_g(ofstream& tau_g, bool thermalActivation, string note) const
{
	tau_g << tau_ext << "\t"
		<< sG << "\t"
		<< sGn << "\t"
		<< av_size << "\t"
		<< av_size_n << "\t"
		<< rand_gen.getTimesAdvanced() << "\t"
		<< getX(avStartPos) << "\t"
		<< getY(avStartPos);
	if (thermalActivation)
		tau_g << "\t"
		<< simTime << "\t"
		<< threshold << "\t"
		<< Zstartval;
	tau_g << note << "\n";
}

deftype getLargestDeformInTau_g(const simPars pars/*, const vector<activity>& loa*/)
{
	ifstream file(fnameWp("tau_g",pars));
	if (!file)
		return 0;

	stringstream tau_g_lastLine(lastLine(file));
	deftype tau_ext = 0, sG = 0, sGn = 0;
	tau_g_lastLine >> tau_ext >> sG >> sGn;

	return (sG - sGn);
}

bool placeLock_tau_g(const simVars& sim, const simPars& pars, string note)
{
	ofstream tau_g(fnameWp("tau_g",pars),fstream::app);
	if (!tau_g)
	{
		cerr << "Cannot open " << fnameWp("tau_g", pars) << " to place the lock" << endl;
		return false;
	}

	tau_g.precision(17);
	sim.printTau_g(tau_g, pars.thermalActivation, note);
	return true;
}

bool removeLock_tau_g(const simVars& sim, const simPars& pars) //remove lock from tau_g if exists
{
	if (is_file_exist(fnameWp("tau_g",pars)))
	{
		ifstream tau_g_in(fnameWp("tau_g",pars));
		vector<string> tau_g_content;
		for (string tmp; getline(tau_g_in,tmp);)
			tau_g_content.push_back(tmp);

		tau_g_in.close();

		if (tau_g_content.end()->find(TAUG_LOCK_TEXT) == tau_g_content.end()->npos) //lock not found
			return false;

		ofstream tau_g_out(fnameWp("tau_g",pars)); //replace content without the last line
		for (auto x = tau_g_content.begin(); x != prev(tau_g_content.end()); ++x)
			tau_g_out << *x << endl;
		
		return true;
	}
	else
		return false;
}



//calculate deformation based on tau_l
deftype calcDeform(int minPoint, deftype local_s, deftype tau_ext, deftype stress_field_at_0, DGamma DG, direction D)
{
	if (DG.MyType == DGamma::definit)
		return D == direction::left ? -DG.val : DG.val;
	else if (DG.MyType == DGamma::limited)
	{
		return D == direction::left ?
			max(-DG.val, (-(tau_ext + local_s) / stress_field_at_0)) :
			min(DG.val, (-(tau_ext + local_s) / stress_field_at_0));
	}
	else // (DG.MyType == DGamma::max)
		return (-(tau_ext + local_s) / stress_field_at_0);
}






bool makeSimulation(simPars pars, const vector<activity>& loa, const vector<double>& extrasnaps)
{

	vector<double> extrasnaps_this_time = extrasnaps;

	cout << "List of activity:\n";
	for (auto m_activity:loa)
		cout << m_activity << endl;
	for (auto extrasnap : extrasnaps_this_time)
		cout << "\t(extra) snapshot at external stress " << extrasnap << endl;


	//// 
	//// create maps and initialize
	//// 

	simVars sim(pars);
	sim.sf.readFromFile(pars.sfk);
	sim.initialize(pars.lfs);
	sim.tau_ext = pars.tau_extMin;


	////
	//// import snapshots, make_deformation
	//// 

	int largestImport = 0;
	for (auto x:loa)
		if (x.MyType == activity::import)
			largestImport = x.nominalDef; //available snapshots are in ascending order

	if (largestImport != 0) //if a snapshot is available, than an already checked tau_g file is also available
	{
		if (!removeLock_tau_g(sim,pars))
		{
			cerr << "Can't remove lock from tau_g" << endl;
			return false;
		}
	}
	else
	{
		if (pars.createTau_g && !init_tau_g(sim,pars))
		{
			cerr << "Can't initialize tau_g" << endl;
			return false;
		}
	}
	if (pars.pNd)
	{
		ofstream ofs;
		string pndfname = fnameWp("pNd", pars);
		ofs.open(pndfname, std::ofstream::out | std::ofstream::trunc);
		if (!ofs)
		{
			cerr << "Cannot open " << pndfname << " for truncating. Program terminates;";
			exit(1);
		}
		ofs.close();
	}

		


	cout << "Start" << endl;
	for (auto x:loa)
	{
		if (x.MyType == activity::import) // import
		{
			cout << currentDateTime() << " Loading snapshot with nominal def " << x.nominalDef << " ... ";
			if (snapshot(pars,x.nominalDef).load(sim,pars,x.nominalDef))
				cout << " successful" << endl;
			else
			{
				cout << " failed\n" << endl;
				cerr << currentDateTime() << " Loading snapshot with nominal def " << x.nominalDef << " failed." << endl;
				return false;
			}
		}
		else if (x.MyType == activity::deform) //make_deformation 
		{
			if (interrupt("DefBegin",pars.seed)) {cerr << "Program skips this deformation" << endl; return false;}

			cout << currentDateTime() << " Making deformation until nominal def " << x.nominalDef << " ... " << endl;

			bool printToTau_g = false;
			if (x.nominalDef > largestImport && pars.createTau_g)
				printToTau_g = true;
				

			bool partial_result;

#ifdef IS_GPU
			if (!pars.mdp())
				partial_result = sim.makeDeformSerial_GPU(pars,x.nominalDef);
			else
			{
				cerr << "Parallel deformation mode is not implemented yet on GPU." << endl;
				return false;
			}
#else
			if (pars.mdp)
				partial_result = sim.makeDeformParallel(pars, x.nominalDef);
			else
				partial_result = sim.makeDeformSerial(pars, x.nominalDef, extrasnaps_this_time, printToTau_g);
#endif
			if (partial_result == true)
				cout << currentDateTime() << " Making deformation until nominal def " << x.nominalDef << " was successful." << endl;

			else
			{
				cerr << currentDateTime() << " Making deformation until nominal def " << x.nominalDef << " was failed." << endl;
				return false;
			}
		}
		else
		{
			cerr << "Unsupported activity type" << endl;
			return false;
		}
	}


	////
	//// Finishing
	////

	cout << "Simulation finished with the last activity." << endl;
	if (pars.createTau_g)
		placeLock_tau_g(sim, pars);

	return true;
}