// correlation_integral.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "../gen_utils/gen_version.h"
#include "../gen_utils/gen_param.h"
#include "../gen_utils/gen_fnameWp.h"
#include "../gen_utils/gen_fmanip.h"
#include "../gen_utils/teestream.h"
#include "../gen_utils/tau_g.h"
#include "../gen_utils/gen_hist.h"

int main(int argc, char* argv[])
{
	intro("This program calculates correlation integral based on the file tau_g or tau_g_pNd. This file must contain the start positions of the avalanches, which are used to calculate the integral.\n", "***");
#pragma region setting parameters
	paramcontainer p, p_temp;
	param<int> s("sysSize", p, rq::must);                              // system size
	param<int> st("seedStart", p, rq::must); //seed start
	param<int> se("seedEnd", p, rq::must); //seed end
	param<string> fp("fnamePrefix", p, rq::optional, ""); //fname prefix
	param<string> fs("fnameSuffix", p, rq::optional, ""); //fname suffix
	param<int> fg("fileGrouping", p, rq::optional, 1000); //file grouping
	param<string> folderName("folderName", p, rq::optional, "tau_g"); // to read from which directory, maybe from tau_g_pNd
	param<int> posCol("posCol", p, rq::optional, 6); // index of the column that stores the x value of the position of the avalanche, next column must be the y position
	param<int> intervals("intervals", p, rq::optional, 10); // number of tau or gamma intervals to make avalanche statistics
	param<double> incr("incr", p, rq::optional, 1); // number of tau intervals to make avalanche statistics
	param<bool> gammaInterval("gammaInterval", p, rq::optional, false); // use gamma interval instead of tau
	param<int> cint_res("cint_res", p, rq::optional); //correlation integral resolution
	param<string> corrPatternListFname("corrPatternListFname", p, rq::optional); // program creates patterns at selected files for the intervals described by the prev param intervals

	param<string> paramFname("paramFname", p_temp, rq::optional, "correlation_integral.ini", argc, argv); //the name of the ini file

	p.setFromIf(paramFname());
	if ((!paramFname.setByUsr() && p.setFromIf(argc, argv) * 2 != argc - 1) || (paramFname.setByUsr() && p.setFromIf(argc, argv) * 2 != argc - 3))
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

#pragma region logfile
	string logFname = fnameWp(folderName(), fp(), st(), se(), fs() + "_corrint_log");
	ofstream his(logFname, ofstream::app);
	if (!his)
	{
		cerr << "Error: unable to append history file " << logFname << ". Program temrinates." << endl;
		return 0;
	}
	teebuf tee_cout(cout.rdbuf(), his.rdbuf());
	cout.rdbuf(&tee_cout);
	teebuf tee_err(cerr.rdbuf(), his.rdbuf());
	cerr.rdbuf(&tee_err);
#pragma endregion

#pragma region define and initialize variables and read in data

	if (!cint_res.is_set())
	{
		cint_res.setDefVal(s());
		cout << "Def val of cint_res has been set to " << cint_res() << endl;
	}

	vector<int> corrPatternListWish;
	vector<pair<int,int>> corrPatternList; // .first stores the desired seed number, .second stores the file ID at which the desired seed number has been reached
	if (corrPatternListFname.is_set())
	{
		ifstream corrPatternListStream(corrPatternListFname());
		if (!corrPatternListStream)
		{
			cerr << "Cannot open " << corrPatternListFname() << ". Program terminates." << endl;
			return 0;
		}
		for (int tmp; comment_skipped(corrPatternListStream) >> tmp; corrPatternListWish.push_back(tmp));
	}

	vector<vector<tau_g_line_pos>> data;
	double tau_ext_max = 0;
	double DsGmax = 0;

	int nof_processed_file = 0;
	for (int seed = st(); seed < se(); ++seed)
	{
		string tau_g_fname = fnameWp(folderName(), seed, fg(), fp(), "");
		ifstream if_tau_g(tau_g_fname);
		if (!if_tau_g)
		{
			cerr << "Can't open " << tau_g_fname << " to load data" << endl;
			continue;
		}
		else if (find(corrPatternListWish.begin(), corrPatternListWish.end(), seed) != corrPatternListWish.end())
			corrPatternList.push_back(pair<int, int>(seed, nof_processed_file));

		cout << "Reading " << tau_g_fname << endl;
		data.push_back(vector<tau_g_line_pos>());

		skip_comment(if_tau_g);
		if_tau_g.ignore(numeric_limits<streamsize>::max(), '\n'); // to skip the first noncommented line, which is without position
		for (tau_g_line_pos tmp(posCol()); if_tau_g >> tmp; data.back().push_back(tmp))
		{
			if_tau_g.ignore(numeric_limits<streamsize>::max(), '\n');
			tau_ext_max = max(tau_ext_max, tmp.tau_ext);
			DsGmax = max(DsGmax, tmp.DsG());
		}

		++nof_processed_file;
	}

	if (nof_processed_file == 0)
	{
		cerr << "No processable file was found. Program terminates." << endl;
		return 0;
	}

	if (incr.setByUsr() && !intervals.setByUsr())
	{
		intervals.unset();
		if (gammaInterval())
			intervals.setDefVal(static_cast<int>(DsGmax / incr()) + 1);
		else
			intervals.setDefVal(static_cast<int>(tau_ext_max / incr()) + 1);
		cout << "Value of intervals is now set to " << intervals.getVal() << endl;

		if (intervals() > 1e6)
		{
			cerr << "Number of intervals is too big: " << intervals() << ". Consider using larger increment value. Program terminates." << endl;
			return 0;
		}
	}
	else if (gammaInterval())
	{
		incr.setDefVal(nextafter(DsGmax / intervals(), numeric_limits<double>::infinity()));
		cout << "Value of incr is now set to " << incr.getVal();
	}
	else
	{
		incr.setDefVal(nextafter(tau_ext_max / intervals(), numeric_limits<double>::infinity()));
		cout << "Value of incr is now set to " << incr.getVal();
	}


	vector<hist<double>> corrint(intervals(), hist<double>(0, s(), scale::linear, cint_res(), 0));

	if (corrPatternListFname.is_set())
	{
		cout << "Program will create patterns at seed\n";
		for (const auto& x : corrPatternList)
			cout << x.first << " (as the " << ord_num(x.second) << " file read)\n";
	}
	

#pragma endregion

#pragma region calculate the integral

	int fileindex = 0;
	for (const auto& tau_g_line_pos_series : data)
	{
		cout << "Processing file " << fileindex << "\n";
		auto where_if_pattern = find_if(corrPatternList.begin(), corrPatternList.end(), [fileindex](const pair<int, int>& val) {return (val.second == fileindex); }); // points to the <seed,ID> pair if patterns should created at this simulation
		bool thisPattern = (where_if_pattern != corrPatternList.end()); // if patterns should printed at this realisation
		int seedPattern = -1; // the seed value of the realisation, defined if thisPattern is true
		if (thisPattern)
			seedPattern = where_if_pattern->first;
		vector<block_data<>> patterns(intervals(), block_data<>(s(), 0)); // stores the real space distribution of startpoint of the events
		vector<block_data<>> w_patterns(intervals(), block_data<>(s(), 0)); // stores the real space distribution of startpoint of the events weighted by the size of the avalanche
		vector<block_data<>> distances(intervals(), block_data<>(s(), 0)); // stores only the distances of the startpoints of the consequent events
		vector<block_data<>> w_distances(intervals(), block_data<>(s(), 0)); // stores only the distances of the startpoints of the consequent events weighted by the size of the avalanche
		int nof_dist_rec = 0; // number of points added to the 2D map
		int nof_patt_rec = 0; // number of points added to the 2D map

		int intervalIndex = 0;
		for (auto line_it = tau_g_line_pos_series.begin(); line_it != tau_g_line_pos_series.end(); ++line_it)
		{
			while (( gammaInterval() && (intervalIndex + 1) * incr() <= line_it->DsG()) ||
				   (!gammaInterval() && (intervalIndex + 1) * incr() <= line_it->tau_ext))
				   ++intervalIndex; // to identify which interval are we in


			for (auto line_it2 = next(line_it); line_it2 != tau_g_line_pos_series.end() &&
				(( gammaInterval() && (intervalIndex + 1) * incr() >= line_it2->DsG()) ||
				 (!gammaInterval() && (intervalIndex + 1) * incr() >= line_it2->tau_ext));
			++line_it2)
			{
				double d = dist(line_it->posx, line_it->posy, line_it2->posx, line_it2->posy, s);
				corrint[intervalIndex].at(d).add(line_it->DG() + line_it2->DG()); // the weight is stored as well, but counts will be printed out as well
				if (thisPattern)
				{
					int posx = round_int(dist(line_it->posx, line_it2->posx, s()));
					int posy = round_int(dist(line_it->posy, line_it2->posy, s()));
					distances[intervalIndex](posx,posy) += 1;
					w_distances[intervalIndex](posx, posy) += line_it->DG() + line_it2->DG();
					++nof_dist_rec;
				}
			}

			if (thisPattern)
			{
				patterns[intervalIndex](round_int(line_it->posx), round_int(line_it->posy + 0.5)) += 1;
				w_patterns[intervalIndex](round_int(line_it->posx), round_int(line_it->posy)) += line_it->DG();
				++nof_patt_rec;
			}
		}

		if (thisPattern)
		{
			for (int i = 0; i < intervals(); ++i)
			{
				distances[i].fWriteToFile(fnameWp("tau_g_pNd", seedPattern, fg(), fp(), fs() + "_" + to_string (i) + "_distances_map.txt"));
				w_distances[i].fWriteToFile(fnameWp("tau_g_pNd", seedPattern, fg(), fp(), fs() + "_" + to_string(i) + "_wighted_distances_map.txt"));
				patterns[i].fWriteToFile(fnameWp("tau_g_pNd", seedPattern, fg(), fp(), fs() + "_" + to_string(i) + "_patterns_map.txt"));
				w_patterns[i].fWriteToFile(fnameWp("tau_g_pNd", seedPattern, fg(), fp(), fs() + "_" + to_string(i) + "_weighted_patterns_map.txt"));
			}
		}
		++fileindex;
	}
#pragma endregion

#pragma region write out data
	stringstream header;
	header << "# Correlation integral values. Each count of a bin gives the \"possibility\" that the start point of two randomly selected avalanche is closer than that value. The sum of the deformation (i.e. sum of av_size - av_size_n) caused by those avalanches is also indicated.\n"
		<< "#bin left(1)\t"
		<< "bin right(2)\t";
	for (size_t i = 0; i < corrint.size(); ++i)
	{
		header << i*incr() << "-" << (i + 1)*incr() << "|nof av(" << i * 4 + 3 << ")\t"
			<< "propability(" << i * 4 + 4 << ")\t"
			<< "total DG(" << i * 4 + 5 << ")\t"
			<< "propability(" << i * 4 + 6 << ")\t";
		corrint[i].cumul();
	}
	string fname;
	if (gammaInterval())
		fname = fnameWp(folderName(), fp(), st(), se(), fs() + "_corrint_gammaIntervals");
	else
		fname = fnameWp(folderName(), fp(), st(), se(), fs() + "_corrint_tauIntervals");

	ofstream corrint_file(fname);
	cout << "Writing " << fname << endl;
	corrint_file << header.str() << endl;
	for (int i = 0; i < corrint[0].getSize(); ++i)
	{
		corrint_file << corrint[0].bins()[i].range() << "\t";
		for (size_t j = 0; j != corrint.size(); ++j)
			corrint_file << corrint[j].binCounter(i) << "\t"
			<< corrint[j].binCounter(i) / corrint[j].totalCounterArea() << "\t"
			<< corrint[j].binIntensity(i) << "\t"
			<< corrint[j].binIntensity(i) / corrint[j].totalIntensityArea() << "\t";
		corrint_file << endl;
	}

#pragma endregion

	return 0;
}

