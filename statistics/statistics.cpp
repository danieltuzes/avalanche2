// statistics.cpp : Defines the entry point for the console application.
// This program requires only information from tau_g folder

#include "stdafx.h"
#include "../for_gnuplot/for_gnuplot_utils.h"
#include "../gen_utils/gen_math.h"
#include "../gen_utils/gen_hist.h"
#include "../gen_utils/gen_fnameWp.h"
#include "../gen_utils/tau_g.h"
#include "../gen_utils/teestream.h"


template<typename T> bool size_stat_tofile(const evalPars& pars, const vector< hist<T> >& data, double tauIncr, string header, string extra_fs) // writes out the avalanche statistic results to files, handles only statistics in tau intervals
{
	ofstream av_file(fnameWp("tau_g", pars, extra_fs));
	if (!av_file)
	{
		cerr << "Can't create " << fnameWp("tau_g", pars, extra_fs);
		return false;
	}
	cout << "Writing " << fnameWp("tau_g", pars, extra_fs) << endl;

	av_file << "# This file contains information on avalanche size. Avalanches in each size intervals are counted at different external stress (tau) or DsG values." << endl
		<< "# The first non-commented line is the total count in that stress or deformation range." << endl
		<< header
		<< "# "
		<< "size: bin_min(1)\t" // header
		<< "bin_max(2)\t"
		<< "count in range:";
	for (int j = 0; j < static_cast<int>(data.size()); ++j)
		av_file << j*tauIncr << "-" << (j + 1)*tauIncr << "(" << j + 3 << ")\t";
	av_file << endl;

	av_file << 0 << "\t" << 0 << "\t";
	for (auto x : data)
		av_file << x.sumCounter() << "\t";
	av_file << endl;

	if (data.size() == 0)
	{
		cerr << "No statistics can be found, file " << fnameWp("tau_g", pars, extra_fs) << " will be empty." << endl;
		return false;
	}

	av_file.precision(av_file.precision() + 2);
	for (int j = 0; j < data[0].getSize(); ++j)
	{
		av_file << data[0].bins()[j].range() << "\t";
		for (size_t i = 0; i != data.size(); ++i)
			av_file << data[i].bins()[j].intensity() << "\t";
		av_file << endl;
	}
	return true;
}


int main(int argc, char* argv[])
{
	stringstream introtext;
	introtext << "This program calculates the\n"
		<< "* averaged stress - strain curve\n"
		<< "* the distribution of the negative and positive part of the avalanches, and the difference of them\n"
		<< "* the average stress at each of the first n avalanche\n"
		<< "* the tau_ext lists at given deformation values"
		<< "* the number of required avalanches to reach those given deformation values" << endl;

	intro(introtext.str(), "---");

#pragma region setting parameters
	paramcontainer p, p_temp;
	param<int> s("sysSize", p, rq::must); //system size, not needed in general
	param<int> st("seedStart", p, rq::must); //seed start
	param<int> se("seedEnd", p, rq::must); //seed end
	param<string> fp("fnamePrefix", p, rq::optional, ""); //fname prefix
	param<string> fs("fnameSuffix", p, rq::optional, ""); //fname suffix
	param<int> fg("fileGrouping", p, rq::optional, 1000); //file grouping
	param<double> sc("springConstant", p, rq::optional, 0); // to chose the decreasement of the external stress after a deformation
	param<bool> totStrain("totStrain", p, rq::optional, false); // if the program should calculate and plot the total strain (elastic+plastic): get_tau_ext() / springConstant is added to every sG-sGn values when reading in tau_g file files
	param<int> tauIntervalsAvStat("tauIntervalsAvStat", p, rq::optional, 10); // number of tau intervals to make avalanche statistics
	param<int> DsGIntervalsAvStat("DsGIntervalsAvStat", p, rq::optional, 10); // number of DsG intervals to make avalanche statistics
	param<double> tauIncrAvStat("tauIncrAvStat", p, rq::optional, 1); // how big a tau_ext range should be while sorting the avalanches for avalanche size distribution
	param<int> tauIntervalsAvAvSize("tauIntervalsAvAvSize", p, rq::optional, 10); // the number of x coord points of the average (av_size - av_size_n) - tau_ext plot
	param<double> tauIncrAvAvSize("tauIncrAvAvSize", p, rq::optional, 1); // the x coord increment of the average (av_size - av_size_n) - tau_ext plot
	param<int> sizeIntervals("sizeIntervals", p, rq::optional, 100); // number of avalanche-size intervals to make statistics
	param<bool> allowStressDecreasing("allowStressDecreasing", p, rq::optional, false); // skip checking increasing external stress
	param<string> DsGfname("DG_list_fname", p, rq::optional, "DG_list.txt"); // for each DsG values stored in the named file, the program creates a list of tau_ext values measured right after the reached DsG values, and the index of the avalanche with which it has been reached
	param<scale> sizeScale("sizeScale", p, rq::optional, scale::quasi_log); //scale of the histograms of the avalances sizes (av_size, av_size_n, av_size - av_size_n for positive values, av_size - av_size_n for negative values)
	param<double> factor("avFactor", p, rq::optional, 1); // multiply the avalanche size with this value. ***This modifies the statistical results!***
	param<double> minAvSize("minAvSize", p, rq::optional, 0.5); // for the size distribution what should be the lower limit of sizes; for logarithmic scaling it could be important
	param<double> maxAvSize("maxAvSize", p, rq::optional, 1e6); // end of xrange of the histogram of the size distribution in the positive case
	param<double> maxAvSize_n("maxAvSize_n", p, rq::optional, 1e6); // end of xrange of the histogram of the size distribution in the negative case
	param<int> sscRes("sscRes", p, rq::optional, 1000); // number of points in the stress-strain curve
	param<int> sscTauRes("sscTauRes", p, rq::optional, 1000); // number of points in the stress-strain curve
	param<double> sscIncr("sscIncr", p, rq::optional, 2); // the incrementation of the bins of the abscissa values of the stress-strain curve
	param<scale> sscScale("sscScale", p, rq::optional, scale::linear); // to be able to modify the abscissa of the stress-strain curve
	param<double> sscLB("sscLB", p, rq::optional, 0.5); // the lower boundary for logscale
	param<int> workUnit("workUnit", p, rq::optional, 1000); // number of files that represents the data range for the evaluation; too big value may cause unsufficient memory error. The program pre-reads this amount of files to calculate the correct values for the histograms - if requested via parameters
	param<int> iav("ithAvalanche", p, rq::optional, 0); //if set, the program makes 2 x iav number of columns that stores the external stress list at the ith avalanche and the avalanche size (ies and iDG_list, the first is in simulation order, the second is in ascending order
	param<int> onlyFirstN("onlyFirstN", p, rq::optional, 0); //if bigger than 0 is set, only the first N will be considered. Note: the size of the first avalanche is always 0.
	param<int> skipFirstN("skipFirstN", p, rq::optional, 0); // the first N lines will be ignored. Useful if you'd like to skip the 0 sized avalanche or in thermal activation
	param<bool> createPozStat("createPozStat", p, rq::optional, false); //if the program should create or not tau_n, tau_neg_part, tau_poz_part
	param<bool> createExtraStat("createExtraStat", p, rq::optional, false); //if the program should create or not tau_n, tau_neg_part, tau_poz_part
	param<bool> createStressStrainCurve("createStressStrainCurve", p, rq::optional, true); // if the program should create stress strain curves at gamma and tau_ext vals
	param<bool> createAvAvsize("createAvAvsize", p, rq::optional, true); // if the program should create average avalanche size - external stress vals
	param<string> paramFname("paramFname", p_temp, rq::optional, "statistics_param.ini", argc, argv); //the name of the ini file

	p.setFromIf(paramFname());
	if ((!paramFname.setByUsr() && p.setFromIf(argc, argv) * 2 != argc - 1) || (paramFname.setByUsr() && p.setFromIf(argc, argv) * 2 != argc - 3))
	{
		cerr << "Unknow parameter or missing value was found in program call, program terminates." << endl;
		p.description();
		return 0;
	}

	if (!p.checkRequired())
	{
		p.printOut(cerr);
		return 0;
	}

	evalPars pars(s, st(), fg, fp, st, se, fs);
#pragma endregion


#pragma region logfile
	string logFname = fnameWp("tau_g", pars, "_log");
	ofstream his(logFname, ofstream::app);
	if (!his) {
		cerr << "Error: unable to append history file " << logFname << ". Program temrinates." << endl;
		return 0;
	}
	teebuf tee_cout(cout.rdbuf(), his.rdbuf());
	cout.rdbuf(&tee_cout);
	teebuf tee_err(cerr.rdbuf(), his.rdbuf());
	cerr.rdbuf(&tee_err);
	p.printOut(cout);
#pragma endregion


#pragma region define variables
	
	vector<vector<tau_g_line>> tau_g_data; //every element stores the whole content of a tau_g file
	tau_g_line mixed_largest; //to store the largest values of the tau_g files (both tau_ext, sG, sGn, av_size, av_size_n, rand_gen.getTimesAdvanced()
	accum<double> av_size_av(0); //to store the average avalanche size; to create histograms in the good range
	accum<double> av_size_n_av(0); 

	hist<tau_g_line> tau_g_sampl_DsG; //sampling the tau_g files at DsG values
	hist<tau_g_line> tau_g_sampl_tau; //sampling the tau_g files at external stress values

	vector<double> DsGlist; 
	if (DsGfname.setByUsr()) // reading in the DsG values at which the external stress values and the the index of the avalanche will be investigated
	{
		ifstream DsGtelf(DsGfname());
		if (!DsGtelf)
		{
			cerr << "Cannot open " << DsGfname << ". Program terminates." << endl;
			return 0;
		}
		for (double DsGtmp; comment_skipped(DsGtelf) >> DsGtmp; DsGlist.push_back(DsGtmp));
	}
	vector<vector<double>> tau_exts(DsGlist.size(), vector<double>()); // to store the external stress values at the given DsGs
	vector<vector<int>> numof_avs(DsGlist.size(), vector<int>()); // to store the number of needed avalanches to reach the given DsGs
	vector<vector<double>> ith_tau_ext(se() - st(), vector<double>(iav() + 1, 0)); // to store the external stress values at the ith and the last avalanche; 1st index(from left to right) of ith_tau_ext : simulation(file) index; 2st index: ith avalanche
	vector<vector<bool>> ith_av_set(se() - st(), vector<bool>(iav() + 1, false)); // to store if the ith avalnche is exist or not for every simulation
	vector<vector<double>> ith_DG(se() - st(), vector<double>(iav() + 1, 0)); // to store the size of avalanche at the ith and at the last avalanche; 1st index(from left to right) of ith_DG : simulation(file) index; 2st index: ith avalanche

	
	hist<double> ssc_DsG, //sampling the external stress curve at DsG; takes into account sc and totStrain
		// ssc_tau, //sampling the stress-strain curve at external stress values, does not take into account sc and totStrain
		ssc_DsGSq; //to calculate the deviation of the external stress
	vector<hist<double>> av_p, //to store the statistics of the avalanche sizes if it is positive, for each tau_ext intervals
		av_p_DsG, //to store the statistics of the avalanche sizes if it is positive, for each DsG intervals
		av_n, //if it is negative, for each tau_ext intervals
		av_poz_part, //to store only the positive part, for each tau_ext intervals
		av_neg_part; //to store only the negative part, for each tau_ext intervals
	vector<accum<double>> average_av_size; //in different stress-intervals

	double DsGIncrAvStat; // how wide a DsG interval should be; will be calculated as the largest DsG / number of intervals

	int i_counter = 0; // the number of files read in
	int p_counter = 0; // the number of processed files; good to know while debugging and it is also printed out at the end
	bool preprocessed = false; // when true, the range of the measured quantities (sG, sGn, tau_ext ...) is reflected on the containers' size, while false, the range of the measured quantities are not used
	int fileindex = 0;

#pragma endregion


#pragma region read in some data, initialize vars, process, read and process the rest
	for (int seed = st(); seed < se();)
	{
		//Read in the first workUnit amount of realisation
		for (int seed_i = 0; seed_i < workUnit() && seed < se(); ++seed, ++seed_i) 
		{
			string tau_g_fname = fnameWp("tau_g", seed, fg(), fp(), "");
			ifstream if_tau_g(tau_g_fname, ios::binary);
			if (!if_tau_g)
			{
				cerr << "Can't open " << tau_g_fname << " to load data. Program skips this seed value." << endl;
				continue;
			}
			cout << "Reading " << tau_g_fname << endl;
			tau_g_data.push_back(vector<tau_g_line>());
			skip_comment(if_tau_g);
			bool incorrect_line = false;
			int linecounter = 0;
			for (int i = 0; i < skipFirstN; ++i)
				if_tau_g.ignore(numeric_limits<streamsize>::max(), '\n');
			for (tau_g_line ttmp; if_tau_g >> ttmp && (onlyFirstN() == 0 || linecounter < onlyFirstN()); linecounter++)
			{
				if_tau_g.ignore(numeric_limits<streamsize>::max(), '\n'); // a line in the tau_g files may conatin further element over tau_g_line
				ttmp *= factor(); // multiplies the sG, sGn, av_size, av_size_n values by factor()
				if (!tau_g_data.back().empty()) // check the consistency of the read-in line
				{
					if (tau_g_data.back().back().sG > ttmp.sG || // sG must not decrease
						tau_g_data.back().back().sGn > ttmp.sGn || // sGn must not decrease
						(tau_g_data.back().back().tau_ext > ttmp.tau_ext && !allowStressDecreasing())) // tau_ext must not decrease except it is allowed
					{
						incorrect_line = true;
						cerr << "File " << tau_g_fname << " contains at least 1 bad line at\n"
							<< ttmp << "\n"
							<< "Program skips this file." << endl;
						break;
					}
					
				}
				tau_g_data.back().push_back(ttmp);
				if (!preprocessed)// update the range of the data
				{
					mixed_largest.extend(ttmp);
					av_size_av.add(ttmp.av_size);
					av_size_n_av.add(ttmp.av_size_n);
				}
			}
			if (incorrect_line)
				continue;
			++i_counter;
		}

		// set the vector container, plastic strain is taken into account
		if (!preprocessed)
		{
			average_av_size = vector<accum<double>>(tauIntervalsAvAvSize(), 0);
			
			if (i_counter == 0)
			{
				cerr << "Cannot read at least one file for initial parameter setting." << endl;
				break;
			}
			if (tauIncrAvStat.setByUsr())
			{
				tauIntervalsAvStat.modify(static_cast<int>(mixed_largest.tau_ext / tauIncrAvStat()) + 1);
				if (tauIntervalsAvStat.setByUsr())
					cerr << "Both tauIncrAvStat and tauIntervalsAvStat are set by user. Value of tauIntervalsAvStat is now set to " << tauIntervalsAvStat.getVal() << endl;
			}
			else
				tauIncrAvStat.setDefVal(nextafter(mixed_largest.tau_ext / tauIntervalsAvStat(), DINFTY));
			
			DsGIncrAvStat = nextafter(mixed_largest.DsG() / DsGIntervalsAvStat(), DINFTY);

			if (tauIncrAvAvSize.setByUsr())
			{
				tauIntervalsAvAvSize.modify(static_cast<int>(mixed_largest.tau_ext / tauIncrAvAvSize()) + 1);
				cout << "Def val of tauIntervalsAvAvSize has been set to " << tauIntervalsAvAvSize() << endl;
				if (tauIntervalsAvAvSize.setByUsr())
					cerr << "Both tauIncrAvAvSize and tauIntervalsAvAvSize are set by user. Value of tauIntervals is now set to " << tauIntervalsAvAvSize.getVal() << endl;
			}
			else
			{
				tauIncrAvAvSize.setDefVal(nextafter(mixed_largest.tau_ext / tauIntervalsAvAvSize(), DINFTY));
				cout << "Def val of tauIncrAvAvSize has been set to " << tauIncrAvAvSize() << endl;
			}

			if (!maxAvSize.setByUsr())
			{
				maxAvSize.setDefVal(av_size_av.averaged() * 1000);
				cout << "Def val of maxAvSize has been set to " << maxAvSize() << endl;
			}

			if (!maxAvSize_n.setByUsr())
			{
				maxAvSize_n.setDefVal(av_size_n_av.averaged() * 1000);
				cout << "Def val of maxAvSize_n has been set to " << maxAvSize_n() << endl;
			}

			if (sscScale() == linear)
			{
				if (sscIncr.setByUsr() && !sscRes.setByUsr())
				{
					sscRes.setDefVal(static_cast<int>(nextafter(!totStrain() ? mixed_largest.DsG() : mixed_largest.tDsG(sc()), DINFTY) / sscIncr() + 1));
					cout << "Def val of sscRes has been set to " << sscRes() << endl;
				}
				if (!sscIncr.setByUsr())
				{
					sscIncr.setDefVal(nextafter(!totStrain() ? mixed_largest.DsG() : mixed_largest.tDsG(sc()), DINFTY) / sscRes());
					cout << "Def val of sscIncr has been set to " << sscIncr() << endl;
				}
				
				tau_g_sampl_DsG = hist<tau_g_line>(0, sscRes()*sscIncr(), scale::linear, sscRes(), tau_g_line());
				ssc_DsG = hist<double>(0, sscRes()*sscIncr(), scale::linear, sscRes(), 0);
				ssc_DsGSq = hist<double>(0, sscRes()*sscIncr(), scale::linear, sscRes(), 0);
			}
			else
			{
				if (sscIncr.setByUsr() && !sscRes.setByUsr())
				{
					sscRes.setDefVal(static_cast<int>(log(nextafter(!totStrain() ? mixed_largest.DsG() : mixed_largest.tDsG(sc()), DINFTY) / sscLB()) / log(sscIncr())) + 1);
					cout << "Def val of sscRes has been set to " << sscRes() << endl;
				}
				if (!sscIncr.setByUsr())
				{
					sscIncr.setDefVal(pow(nextafter(!totStrain() ? mixed_largest.DsG() : mixed_largest.tDsG(sc()), DINFTY) / sscLB(), 1. / sscRes()));
					cout << "Def val of sscIncr has been set to " << sscIncr() << endl;
				}
				tau_g_sampl_DsG = hist<tau_g_line>(sscLB(), sscLB()*pow(sscIncr(), sscRes()), sscScale(), sscRes(), tau_g_line());
				ssc_DsG = hist<double>(sscLB(), sscLB()*pow(sscIncr(), sscRes()), sscScale(), sscRes(), 0);
				ssc_DsGSq = hist<double>(sscLB(), sscLB()*pow(sscIncr(), sscRes()), sscScale(), sscRes(), 0);
			}

			tau_g_sampl_tau = hist<tau_g_line>(0, nextafter(mixed_largest.tau_ext, DINFTY), scale::linear, sscTauRes(), tau_g_line());
			// ssc_tau = hist<double>(0, nextafter(mixed_largest.tau_ext, DINFTY), scale::linear, sscTauRes(), 0);
			
			av_p = vector<hist<double>>(tauIntervalsAvStat, hist<double>(minAvSize, maxAvSize, sizeScale, sizeIntervals, 0));
			av_p_DsG = vector<hist<double>>(DsGIntervalsAvStat, hist<double>(minAvSize, maxAvSize, sizeScale, sizeIntervals, 0));
			av_n = vector<hist<double>>(tauIntervalsAvStat, hist<double>(minAvSize, maxAvSize_n, sizeScale, sizeIntervals, 0));
			av_poz_part = vector<hist<double>>(tauIntervalsAvStat, hist<double>(minAvSize, maxAvSize, sizeScale, sizeIntervals, 0));
			av_neg_part = vector<hist<double>>(tauIntervalsAvStat, hist<double>(minAvSize, maxAvSize_n, sizeScale, sizeIntervals, 0));
			average_av_size = vector<accum<double>>(tauIntervalsAvAvSize, 0);

			preprocessed = true;

		}

		for (auto tau_g_list : tau_g_data)
		{
			int tauStatIntervalId = 0, // in which tau_ext interval the iterator is for the avalanche size statistics
				tauSizeIntervalId = 0, // in which tau_ext interval the iterator is for the average avalanche size
				DsGStatIntervalId = 0; // in which DsG interval the iterator is for the avalanche size statistics
			int DsG_index = 0;
			int tau_index = 0;
			int which_DsG4tau_extList = 0;
			for (auto tau_g_it = tau_g_list.begin(); tau_g_it != tau_g_list.end(); ++tau_g_it)
			{
				if (DsGfname.setByUsr())
				{
					while (which_DsG4tau_extList < static_cast<int>(DsGlist.size()) && tau_g_it->DsG() > DsGlist[which_DsG4tau_extList]) // if an avalanche is large, a tau_ext value can belong to multiple required DsG values
					{
						tau_exts[which_DsG4tau_extList].push_back(tau_g_it->tau_ext); // will added to which_DsG4tau_extList, which_DsG4tau_extList + 1, as long as tau_g_it->DsG() > DsGlist[which_DsG4tau_extList]
						numof_avs[which_DsG4tau_extList].push_back(distance(tau_g_list.begin(), tau_g_it));
						++which_DsG4tau_extList;
					}
				}
				
				if (iav.setByUsr())
				{
					int av_number = distance(tau_g_list.begin(), tau_g_it);
					if (av_number < iav())
					{
						ith_tau_ext[fileindex][av_number] = tau_g_it->tau_ext;
						ith_DG[fileindex][av_number] = tau_g_it->DG();
						ith_av_set[fileindex][av_number] = true;
					}
					if (tau_g_it == prev(tau_g_list.end()))
					{
						ith_tau_ext[fileindex][iav()] = tau_g_it->tau_ext;
						ith_DG[fileindex][iav()] = tau_g_it->DG();
						ith_av_set[fileindex][iav()] = true;
					}
				}
				
				if (createStressStrainCurve)
				{
					// sampling the stress-strain curve at gamma values
					if (!totStrain()) // only plastic strain
					{
						while (DsG_index < tau_g_sampl_DsG.getSize() && tau_g_sampl_DsG.bins()[DsG_index].lowerb() <= tau_g_it->DsG())
						{
							tau_g_sampl_DsG.bins()[DsG_index].add(*tau_g_it);
							if (sc.setByUsr())
								ssc_DsG.bins()[DsG_index].add(tau_g_it->tau_ext + (tau_g_it->DsG() - tau_g_sampl_DsG.bins()[DsG_index].lowerb()) * sc()); // that's true, double checked
							else
								ssc_DsG.bins()[DsG_index].add(tau_g_it->tau_ext);
							ssc_DsGSq.bins()[DsG_index].add(pow(tau_g_it->tau_ext, 2));
							++DsG_index;
						}
					}
					else
					{
						while (DsG_index < tau_g_sampl_DsG.getSize() && tau_g_sampl_DsG.bins()[DsG_index].lowerb() <= tau_g_it->tDsG(sc()))
						{
							tau_g_sampl_DsG.bins()[DsG_index].add(*tau_g_it);
							if (tau_g_it == tau_g_list.begin())
								ssc_DsG.bins()[DsG_index].add(tau_g_it->tDsG(sc())); // in the add(), it is 0 usually
							else
								ssc_DsG.bins()[DsG_index].add(prev(tau_g_it)->tau_ext + (tau_g_sampl_DsG.bins()[DsG_index].lowerb() - prev(tau_g_it)->tDsG(sc()))*sc()); // that's true, checked a bit aain
							ssc_DsGSq.bins()[DsG_index].add(pow(tau_g_it->tau_ext, 2));
							++DsG_index;
						}
					}

					// sampling the stress-strain curve at external values
					while (tau_index < tau_g_sampl_tau.getSize() && tau_g_sampl_tau.bins()[tau_index].lowerb() < tau_g_it->tau_ext)
					{
						tau_g_sampl_tau.bins()[tau_index].add(*prev(tau_g_it));
						// ssc_tau.bins()[tau_index].add(prev(tau_g_it)->DsG());
						++tau_index;
					}
				}

				if (createAvAvsize)
				{
					for (; (tauSizeIntervalId + 1) * tauIncrAvAvSize() < tau_g_it->tau_ext; ++tauSizeIntervalId);
					if (tauSizeIntervalId < tauIntervalsAvAvSize())
						average_av_size[tauSizeIntervalId].add(tau_g_it->DG());
				}

				if (createPozStat || createExtraStat)
				{
					for (; (tauStatIntervalId + 1) * tauIncrAvStat() < tau_g_it->tau_ext; ++tauStatIntervalId);
					for (; (DsGStatIntervalId + 1) * DsGIncrAvStat   < tau_g_it->DsG();   ++DsGStatIntervalId);
					if (tauStatIntervalId < tauIntervalsAvStat())
					{
						if (tau_g_it->DG() >= 0)
							av_p[tauStatIntervalId].ifIncrAt(tau_g_it->DG(), 1);
						else
							av_n[tauStatIntervalId].ifIncrAt(tau_g_it->DG(), 1);

						av_poz_part[tauStatIntervalId].ifIncrAt(tau_g_it->av_size, 1);
						av_neg_part[tauStatIntervalId].ifIncrAt(tau_g_it->av_size_n, 1);
					}
					if (DsGStatIntervalId < DsGIntervalsAvStat())
					{
						if (tau_g_it->DG() >= 0)
							av_p_DsG[DsGStatIntervalId].ifIncrAt(tau_g_it->DG(), 1);
					}
				}

			}
			++p_counter;
			++fileindex;
		}

		tau_g_data.clear();
	}
	
#pragma endregion



#pragma region create output files
	if (p_counter == 0)
	{
		cerr << "No files were processed, no output will be created. Check log file " << logFname;
		return 0;
	}

	if (createStressStrainCurve)
	{
		stringstream tau_g_sampl_DsG_header;
		tau_g_sampl_DsG_header << "# Deformation resampled at " << (totStrain() ? string("total strain values") : string("plastic strain values")) << " .In each bin tau_g_line elements are collected.\n"
			<< "# In every simulation, each bin is enriched with a tau_g_line element that has equal or bigger DsG than the lower boundary of the bin and that tau_g_line element has the smallest DsG of all.\n"
			<< "# Stress decreasement due to the drop-off after plastic deformation is not taken into account.\n"
			<< "# bin DsG(1)\t"
			<< "next bin DsG(2)\t"
			<< "nof sims at DsG(3)\t"
			<< "tot tau_ext(4)\t"
			<< "tot sG(5)\t"
			<< "tot sGn(6)\t"
			<< "tot av_size(7)\t"
			<< "tot av_size_n(8)\t"
			<< "tot rand_gen(9)\t"
			<< "averaged tau_ext(10)\t"
			<< "averaged sG(11)\t"
			<< "averaged sGn(12)\t"
			<< "averaged av_size(13)\t"
			<< "averaged av_size_n(14)\t"
			<< "averaged rand_gen(15)\t"
			<< endl;
		tau_g_sampl_DsG.use_averaged().fWrite(fnameWp("tau_g", pars, "_tau_g_sampl_DsG"), tau_g_sampl_DsG_header.str());

		stringstream tau_g_sampl_tau_header;
		tau_g_sampl_tau_header << "# Deformation resampled at tau_ext. In each bin tau_g_line elements are collected.\n"
			<< "# In every simulation, each bin is enriched with a tau_g_line element that has equal or bigger tau_ext than the lower boundary of the bin and that tau_g_line element has the smallest tau_ext of all.\n"
			<< "# Stress decreasement due to the drop-off after plastic deformation is not taken into account.\n"
			<< "# bin tau_ext(1)\t"
			<< "next bin tau_ext(2)\t"
			<< "nof sims at DsG(3)\t"
			<< "tot tau_ext(4)\t"
			<< "tot sG(5)\t"
			<< "tot sGn(6)\t"
			<< "tot av_size(7)\t"
			<< "tot av_size_n(8)\t"
			<< "tot rand_gen(9)\t"
			<< "averaged tau_ext(10)\t"
			<< "averaged sG(11)\t"
			<< "averaged sGn(12)\t"
			<< "averaged av_size(13)\t"
			<< "averaged av_size_n(14)\t"
			<< "averaged rand_gen(15)\t"
			<< endl;
		tau_g_sampl_tau.use_averaged().fWrite(fnameWp("tau_g", pars, "_tau_g_sampl_tau"), tau_g_sampl_tau_header.str());

		stringstream ssc_DsG_header;
		ssc_DsG_header << "# Deformation resampled at DsG. In each bin tau_ext values are collected. "
			<< "In every simulation, each bin is enriched with tau_ext values with its corresponding DsG* equal or bigger than the lower boundary of the bin and that DsG* is the smallest DsG of all.\n"
			<< "# bin DsG(1)\tnext bin DsG(2)\tnof sims at DsG(3)\ttot tau_ext(4)\taverage tau_ext(5)" << endl;
		ssc_DsG.use_averaged().fWrite(fnameWp("tau_g", pars, "_ssc_DsG"), ssc_DsG_header.str());

		stringstream ssc_DsGdev_header;
		ssc_DsGdev_header << "# standard devition of tau at DsG."
			<< "In every simulation, each bin is enriched with tau_ext^2 values with its corresponding DsG* equal or bigger than the lower boundary of the bin and that DsG* is the smallest DsG of all. Before writing the results to a file, tau_ext average square is substracted and then square root calculated based on file " << ssc_DsGdev_header.str() << ".\n"
			<< "# bin DsG(1)\tnext bin DsG(2)\ttau_ext deviation(3)" << endl;
		ofstream ssc_DsGdev_file(fnameWp("tau_g", pars, "_ssc_DsGdev"));
		ssc_DsGdev_file << ssc_DsGdev_header.str();
		for (int i = 0; i < ssc_DsGSq.getSize(); ++i)
			ssc_DsGdev_file << ssc_DsGSq.bins()[i].lowerb() << "\t" << ssc_DsGSq.bins()[i].upperb() << "\t" << pow(ssc_DsGSq.bins()[i].averaged() - pow(ssc_DsG.bins()[i].averaged(), 2), 0.5) << "\n";

		// sampling ssc at tau_ext is useless because tau_g_sampl_tau contains more, and none of them considers sc and totStrain
		/*stringstream ssc_tau_header;
		ssc_tau_header << "# Deformation resampled at tau. In each bin DsG values are collected. "
			<< "In every simulation, each bin is enriched with DsG values with its corresponding tau_ext* equal or bigger than the lower boundary of the bin and that tau_ext* is the smallest tau_ext of all.\n"
			<< "# bin tau_ext(1)\tnext bin tau_ext(2)\tnof sims at tau_ext(3)\ttot DsG(4)\taverage DsG(5)" << endl;
		ssc_tau.use_averaged().fWrite(fnameWp("tau_g", pars, "_ssc_tau"), ssc_tau_header.str());*/ 
	}

	if (createPozStat)
	{
		size_stat_tofile(pars, av_p, tauIncrAvStat(), "# This file is a DG value histogram with positive values at tau_ext intervals.\n", "_av_p");
		size_stat_tofile(pars, av_p_DsG, DsGIncrAvStat, "# This file is a DG value histogram with positive values at DsG intervals.\n", "_av_p_DsG");
	}

	if (createExtraStat())
	{
		size_stat_tofile(pars, av_n, tauIncrAvStat(), "# This file is a DG value histogram with negative values at tau_ext intervals.\n", "_av_n");
		size_stat_tofile(pars, av_poz_part, tauIncrAvStat(), "# This file is a av_size value histogram with positive values at tau_ext intervals.\n", "_av_poz_part");
		size_stat_tofile(pars, av_neg_part, tauIncrAvStat(), "# This file is a av_size_n value histogram with positive values at tau_ext intervals.\n", "_av_neg_part");
	}

	if (createAvAvsize)
	{
		ofstream av_s_av(fnameWp("tau_g", pars, "_av_s_av"));
		cout << "Writing " << fnameWp("tau_g", pars, "_av_s_av") << endl;
		av_s_av << "# The averaged avalanche size in each tau intervals." << endl
			<< "# lb(1)\t" << "ub(2)\t" << "total count(3)\t" << "total value(4)\t" << "averaged value(5)" << endl;
		for (int i = 0; i < tauIntervalsAvAvSize(); ++i)
			av_s_av << i*tauIncrAvAvSize() << "\t" << (i + 1)*tauIncrAvAvSize() << "\t" << average_av_size[i] << "\t" << average_av_size[i].averaged() << endl;

	}
	
	if (!DsGlist.empty())
	{
		// Creating tau_ext_list
		{
			cout << "Writing " << fnameWp("tau_g", pars, "_tau_ext_list") << endl;
			ofstream tau_ext_list(fnameWp("tau_g", pars, "_tau_ext_list"));
			tau_ext_list << "# List of tau_ext at given DsG values.\n# ";
			for (int i = 0; i < static_cast<int>(DsGlist.size()); ++i)
				tau_ext_list << DsGlist[i] << " (" << i * 2 + 1 << ")\t" << DsGlist[i] << " (" << i * 2 + 2 << ") sorted" << "\t";
			tau_ext_list << endl;

			vector<vector<double>> tau_ext_inorder(tau_exts);
			for (int i = 0; i < static_cast<int>(DsGlist.size()); ++i)
				sort(tau_ext_inorder[i].begin(), tau_ext_inorder[i].end());

			for (int j = 0; j < static_cast<int>(tau_exts[0].size()); ++j)
			{
				for (int i = 0; i < static_cast<int>(tau_exts.size()) && j < static_cast<int>(tau_exts[i].size()); ++i)
					tau_ext_list << tau_exts[i][j] << "\t" << tau_ext_inorder[i][j] << "\t";

				tau_ext_list << "\n";
			}
		}
		
		
		//creating numof_avs
		{
			cout << "Writing " << fnameWp("tau_g", pars, "_numof_avs") << endl;
			ofstream numof_avs_list(fnameWp("tau_g", pars, "_numof_avs"));
			numof_avs_list << "# List of the ordinary number of avalanches that reached the given DsG value.\n# ";
			for (int i = 0; i < static_cast<int>(DsGlist.size()); ++i)
				numof_avs_list << DsGlist[i] << " (" << i * 2 + 1 << ")\t" << DsGlist[i] << " (" << i * 2 + 2 << ") sorted" << "\t";
			numof_avs_list << endl;

			vector<vector<int>> numof_avs_inorder(numof_avs);
			for (int i = 0; i < static_cast<int>(DsGlist.size()); ++i)
				sort(numof_avs_inorder[i].begin(), numof_avs_inorder[i].end());

			for (int j = 0; j < static_cast<int>(numof_avs[0].size()); ++j)
			{
				for (int i = 0; i < static_cast<int>(numof_avs.size()) && j < static_cast<int>(numof_avs[i].size()); ++i)
					numof_avs_list << numof_avs[i][j] << "\t" << numof_avs_inorder[i][j] << "\t";

				numof_avs_list << "\n";
			}
		}
		
	}

	if (iav.setByUsr())
	{
		ofstream ies_file(fnameWp("tau_g", pars, "_ies"));
		ofstream iDG_file(fnameWp("tau_g", pars, "_iDG"));
		ofstream itotDG_file(fnameWp("tau_g", pars, "_itotDG"));
		ofstream ies_av_file(fnameWp("tau_g", pars, "_ies_av"));
		ofstream iDG_av_file(fnameWp("tau_g", pars, "_iDG_av"));
		ofstream itotDG_av_file(fnameWp("tau_g", pars, "_itotDG_av"));
		cout << "Writing " << fnameWp("tau_g", pars, "_ies") << endl;
		cout << "Writing " << fnameWp("tau_g", pars, "_ies_av") << endl;
		cout << "Writing " << fnameWp("tau_g", pars, "_iDG") << endl;
		cout << "Writing " << fnameWp("tau_g", pars, "_iDG_av") << endl;
		cout << "Writing " << fnameWp("tau_g", pars, "_itotDG") << endl;
		cout << "Writing " << fnameWp("tau_g", pars, "_itotDG_av") << endl;
		ies_file << "# The external stress values at the ith avalanche. "
			<< "External stresses at the ith avalanche in simulation order can be found at the i*2-1th column and in ascending order in the i*2th column.\n# ";
		iDG_file << "# The strain size of the ith avalanche. "
			<< "Strain sizes at the ith avalanche in simulation order can be found at the i*2-1th column and in ascending order in the i*2th column.\n# ";
		itotDG_file << "# The total strain size at the ith avalanche. "
			<< "Total strain sizes at the ith avalanche in simulation order can be found at the i*2-1th column and in ascending order in the i*2th column.\n# ";
		ies_av_file << "# The average of external stress at the ith avalanche. "
			<< "The first column stores the number of simulations that contains at least i number of sims and second column contains the sum of these external stress values. The third line contains the standard deviation\n"
			<< "# The ith line corresponds to the ith avalanche.\n# Nof sims\tsum of ext_stress\tstandard dev." << endl;
		iDG_av_file << "# The average of the size of the strain of the ith avalanche, and its deviation. "
			<< "The first column stores the number of simulations that contains at least i number of sims and "
			<< "the second column contains the sum of the ith strain over these multiple sims. "
			<< "The third column contains the standard deviation of the ith strain."
			<< "# The ith line corresponds to the ith avalanche.\n# Nof sims\tsum of strains\tstandard dev." << endl;
		itotDG_av_file << "# The average of the total size of the strain at the ith avalanche, and its deviation. "
			<< "The first column stores the number of simulations that contains at least i number of sims and "
			<< "the second column contains the sum of the total strain at the ith avalanche. "
			<< "The third column contains the standard deviation of the tot strain.\n";

		for (int i = 0; i < iav(); ++i)
		{
			ies_file << (i + 1) * 2 - 1 << " (" << i + 1 << "th)\t" << (i + 1) * 2 << " (" << i + 1 << "th in order)" << "\t";
			iDG_file << (i + 1) * 2 - 1 << " (" << i + 1 << "th)\t" << (i + 1) * 2 << " (" << i + 1 << "th in order)" << "\t";
			itotDG_file << (i + 1) * 2 - 1 << " (" << i + 1 << "th)\t" << (i + 1) * 2 << " (" << i + 1 << "th in order)" << "\t";
		}
		ies_file << "last" << "\tlast (in order)\t" << endl;
		iDG_file << "last" << "\tlast (in order)\t" << endl;
		itotDG_file << "last" << "\tlast (in order)\t" << endl;

		vector<vector<double>> ith_tau_ext_ordered(iav() + 1, vector<double>()); //+1 for the last external stress, not for the i+1th
		vector<vector<double>> ith_DG_ordered(iav() + 1, vector<double>()); //+1 for the last external stress, not for the i+1th
		vector<vector<double>> ith_totDG(iav() + 1, vector<double>()); //+1 for the last external stress, not for the i+1th
		vector<vector<double>> ith_totDG_ordered(iav() + 1, vector<double>()); //+1 for the last external stress, not for the i+1th
		int nof_avs = 0;
		for (int ith = 0; ith < iav() + 1; ++ith)
		{
			for (int fileindex = 0; fileindex < se() - st(); ++fileindex)
			{
				if (ith_av_set[fileindex][ith])
				{
					ith_tau_ext_ordered[ith].push_back(ith_tau_ext[fileindex][ith]);
					ith_DG_ordered[ith].push_back(ith_DG[fileindex][ith]);
					ith_totDG[ith].push_back(accumulate(ith_DG[fileindex].begin(), ith_DG[fileindex].begin() + ith, 0.0));
					ith_totDG_ordered[ith].push_back(accumulate(ith_DG[fileindex].begin(), ith_DG[fileindex].begin() + ith, 0.0));
				}
			}
			
			nof_avs = ith_tau_ext_ordered[ith].size();
			sort(ith_tau_ext_ordered[ith].begin(), ith_tau_ext_ordered[ith].end());
			sort(ith_DG_ordered[ith].begin(), ith_DG_ordered[ith].end());
			sort(ith_totDG_ordered[ith].begin(), ith_totDG_ordered[ith].end());

			double sum_tau_ext = accumulate(ith_tau_ext_ordered[ith].begin(), ith_tau_ext_ordered[ith].end(), 0.0);
			double sumsq_tau_ext = accumulate(ith_tau_ext_ordered[ith].begin(), ith_tau_ext_ordered[ith].end(), 0.0, [](double s, double x){return s + x*x; });
			double sum_DG = accumulate(ith_DG_ordered[ith].begin(), ith_DG_ordered[ith].end(), 0.0);
			double sumsq_DG = accumulate(ith_DG_ordered[ith].begin(), ith_DG_ordered[ith].end(), 0.0, [](double s, double x){return s + x*x; });
			double sum_totDG = accumulate(ith_totDG_ordered[ith].begin(), ith_totDG_ordered[ith].end(), 0.0);
			double sumsq_totDG = accumulate(ith_totDG_ordered[ith].begin(), ith_totDG_ordered[ith].end(), 0.0, [](double s, double x){return s + x*x; });

			ies_av_file << nof_avs << "\t"
				<< sum_tau_ext << "\t"
				<< sqrt(sumsq_tau_ext / nof_avs - sum_tau_ext*sum_tau_ext / nof_avs / nof_avs) << "\n";

			iDG_av_file << nof_avs << "\t"
				<< sum_DG << "\t"
				<< sqrt(sumsq_DG / nof_avs - sum_DG*sum_DG / nof_avs / nof_avs) << "\n";

			itotDG_av_file << nof_avs << "\t"
				<< sum_totDG << "\t"
				<< sqrt(sumsq_totDG / nof_avs - sum_totDG*sum_totDG / nof_avs / nof_avs) << "\n";
		}

		for (int fileindex = 0; fileindex < static_cast<int>(ith_tau_ext.size()); ++fileindex)
		{
			for (int ith = 0; ith < iav() + 1; ++ith)
			{
				if (ith_av_set[fileindex][ith])
				{
					ies_file << ith_tau_ext[fileindex][ith] << "\t";
					iDG_file << ith_DG[fileindex][ith] << "\t";
					itotDG_file << ith_totDG[ith][fileindex] << "\t";
				}
				else
				{
					ies_file << "over" << "\t";
					iDG_file << "over" << "\t";
					itotDG_file << "over" << "\t";
				}

				if (fileindex < static_cast<int>(ith_tau_ext_ordered[ith].size()))
				{
					ies_file << ith_tau_ext_ordered[ith][fileindex] << "\t";
					iDG_file << ith_DG_ordered[ith][fileindex] << "\t";
					itotDG_file << ith_totDG_ordered[ith][fileindex] << "\t";
				}
				else
				{
					ies_file << "\t";
					iDG_file << "\t";
					itotDG_file << "\t";
				}
			}
			ies_file << "\n";
			iDG_file << "\n";
			itotDG_file << "\n";
		}
	}
	
#pragma endregion

	cout << "\nTotal of " << i_counter << " files have been read and a total of " << p_counter << " files have been processed." << endl;

	return 0;
}

