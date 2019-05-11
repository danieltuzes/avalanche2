// shape.cpp : Defines the entry point for the console application.
// This program requires information from tau_g and pNd folders only

#include "stdafx.h"
#include "../for_gnuplot/for_gnuplot_utils.h"
#include "../gen_utils/gen_math.h"
#include "../gen_utils/gen_hist.h"
#include "../avalanche2/make_simulation.h"
#include "tau_g_stat.h"
#include "tau_g_pNd_parser.h"

using namespace std;

bool init_tau_g_pNd(ofstream& ofile)
{
	if (!ofile)
	{
		cerr << "Can't initilaize tau_g_pNd file" << endl;
		return false;
	}

	ofile << "# tau_ext(1)\t"
		  << "sG(2)\t"
		  << "sGn(3)\t"
		  << "av_size(4)\t"
		  << "av_size_n(5)\t"
		  << "times advanced rnd gen(6)\t"
		  << "tot nof events(7)\t"
		  << "diamx(8)\t"
		  << "diamy(9)\t"
		  << "first posx(10)\t"
		  << "first posy(11)\t"
		  << "last posx(12)\t"
		  << "last posy(13)\t"
		  << "CoGx(14)\t"
		  << "CoGy(15)" << endl;
	
	return true;
}

int getIndex(vector<double> list, double compare) // gives the index of the bin corresponding to compare
{
	return distance(list.begin(), upper_bound(list.begin(), list.end(), compare)) - 1;

}

int main(int argc, char * argv[])
{	
	stringstream introstream;
	introstream << "This program calculates various maps and others from pNd infos, such as\n"
		<< "* the shape of avalanches\n"
		<< "* the nth action distance from the first in the avalanche\n"
		<< "* enriched (by diameter, CoG, first point, last point ...) tau_g files at tau_g_pNd\n"
		<< "* Distance-map, where the cell value I(x,y) stores the number of cons. events with distance (x,y)\n"
		<< "Use shape_stat on tau_g_pNd to make histograms of shape-related infos.";
	intro(introstream.str(), "---");

#pragma region setting parameters
	paramcontainer p,p_temp;
	param<int> s("sysSize",p,rq::must); //system size
	param<int> st("seedStart",p,rq::must); //seed start
	param<int> se("seedEnd",p,rq::must); //seed end
	param<string> fp("fnamePrefix",p,rq::optional,""); //fname prefix
	param<string> ifs("ifnameSuffix", p, rq::optional, ""); //fname suffix
	param<string> ofs("ofnameSuffix", p, rq::optional, ""); //fname suffix
	param<int> fg("fileGrouping",p,rq::optional,1000); //file grouping
	param<double> factor("avFactor", p, rq::optional, 1); // multiply the avalanche size with this value. ***This modifies the statistics result!***
	param<double> tauIncr("tauIncr",p,rq::optional,1); // tau increment in shape; overflowing data will be ignored
	param<int> tauIntervals("tauIntervals", p, rq::optional, 10); // number of tau intervals to make avalanche statistics; overflowing data will be ignored
	param<string> sizeIntervalsListFname("sizeIntervalsListFname", p, rq::optional, "shape_size_list.ini"); // to store that at what size intervals should the program generate average avalanche shape
	param<int> actionTime("actionTime", p, rq::optional, 0); // how many steps should be counted in the nth event distance map (the distance from the first event)
	param<bool> calc_diam("calc_diam", p, rq::optional, true); // if the program should calculate diamter infos, such as avalanche diameter and center of gravity
	param<bool> calc_cog("calc_cog", p, rq::optional, true); // if the program should calculate center of gravity infos
	param<bool> calc_distances("calc_distances", p, rq::optional, true); // if the program should calculate the distances between consecutive deformation events, even during avalanches
	param<bool> calc_shapes("calc_shapes", p, rq::optional, true); // if the program should calculate shape infos, the shape maps
	param<double> allEventsTauExtMin("allEventsTauExtMin", p, rq::optional, -numeric_limits<double>::infinity()); // to print out tau_g_pNd records only at tau_ext larger than allEventsTauExtMin
	param<double> allEventsTauExtMax("allEventsTauExtMax", p, rq::optional, numeric_limits<double>::infinity()); // to print out tau_g_pNd records only at tau_ext smaller than allEventstauExtMax
	param<bool> allEvent("allEvent", p, rq::optional, false); //if the program should print out
	param<string> paramFname("paramFname", p_temp, rq::optional, "shape_param.ini", argc, argv); //the name of the ini file

	p.setFromIf(paramFname());
	if ((!paramFname.setByUsr() && p.setFromIf(argc,argv) * 2 != argc - 1) || (paramFname.setByUsr() && p.setFromIf(argc,argv) * 2 != argc - 3))
	{
		cerr << "Unknow parameter or missing value was found in program call, program terminates." << endl;
		p.description();
		return 1;
	}

	if (!p.summerize())
	{
		p.description();
		return 1;
	}
#pragma endregion

#pragma region define and initilaize variables
	vector<avalanche> avalanches;
	if (!fill_avalanches(fnameWp("pNd",st(),fg(),fp(),ifs()),fnameWp("tau_g",st(),fg(),fp(),ifs()),avalanches,s(),factor()))
	{
		cerr << "Can't fill_avalanches at seed " << st() << endl;
		return 1;
	}

	if (!tauIncr.setByUsr())
	{
		double set_tauIncr = nextafter(avalanches.back().tau_ext / tauIntervals(), numeric_limits<double>::infinity());
		tauIncr.setDefVal(set_tauIncr);
		cout << "tauIncr.setDefVal " << set_tauIncr << endl;
	}
	else if (!tauIntervals.setByUsr())
	{
		int set_tauIntervals = static_cast<int>(avalanches.back().tau_ext / tauIncr()) + 1;
		tauIntervals.setDefVal(set_tauIntervals);
		cout << "tauIncr.tauIntervals " << set_tauIntervals << endl;
	}

	int sizeIntervalsCount = 0;
	vector<double> sizeIntervalsList;
	if (calc_shapes())
	{
		ifstream sizeIntervalsListF(sizeIntervalsListFname());
		if (!sizeIntervalsListF)
		{
			cerr << "Can't open " << sizeIntervalsListFname() << endl;
			return 1;
		}
		for (double boundary; comment_skipped(sizeIntervalsListF) >> boundary; sizeIntervalsList.push_back(boundary));
		if (!is_sorted(sizeIntervalsList.begin(), sizeIntervalsList.end()))
		{
			cerr << sizeIntervalsListFname() << " must containss values in a monoton ascending order." << endl;
			return 1;
		}
		sizeIntervalsCount = sizeIntervalsList.size() - 1;
	}
	else if (sizeIntervalsListFname.setByUsr())
		cout << "Warning: sizeIntervalsListFname is set by user but calc_shapes is false." << endl;




	hist<linvec<double>> avDiams; // to store avalanche diameter histogram
	if (calc_diam())
	{
		avDiams = hist<linvec<double>>(0, tauIncr()*tauIntervals(), scale::linear, tauIntervals(), linvec<double>(vector<double>(2, 0)));
	}

	vector<block_data<double>> avShapesTot; //to store the stress-unbinned average avalanche shape
	vector<vector<block_data<double>>> avShapes; //to store stress-binned average avalanche shape
	vector<vector<int>> avShapes_counter; //to store how many avalanches were coung in each stress bin
	if (calc_shapes())
	{
		avShapesTot = vector<block_data<double>>(tauIntervals(), block_data<double>(s(), 0.)); 
		avShapes = vector<vector<block_data<double>>>(tauIntervals(), vector<block_data<double>>(sizeIntervalsCount, block_data<double>(s(), 0.)));
		avShapes_counter = vector<vector<int>>(tauIntervals(), vector<int>(sizeIntervalsCount, 0));
	}

	vector<vector<block_data<double>>> nthDistance; //for the stress-binned maps that store the distance of the nth activation from the first one
	vector<block_data<double>> dists; //to store the distance-map of consecutive deformation events
	vector<block_data<double>> alldists; //to store the distance-map of all of the deformations in a stress-range
	if (calc_distances)
	{
		nthDistance = vector<vector<block_data<double>>>(tauIntervals(), vector<block_data<double>>(actionTime(), block_data<double>(s(), 0.)));
		dists = vector<block_data<double>>(tauIntervals(), block_data<double>(s(), 0.));
		alldists = vector<block_data<double>>(tauIntervals(), block_data<double>(s(), 0.));
	}
	else if (actionTime.setByUsr())
		cout << "Wanting: actionTime is set by user, but calc_distances is false. " << endl;


#pragma endregion 

#pragma region read files in and do the analysis and print tau_g_pNd output
	for (int seed = st(); seed != se(); ++seed) // yes, the first file is already read in, but it makes the code shorter
	{
		string pNdFname = fnameWp("pNd", seed, fg(), fp(), ifs());
		string tau_gFname = fnameWp("tau_g", seed, fg(), fp(), ifs());
		cout << "Using files" << endl
			<< "\t" << pNdFname << endl
			<< "\t" << tau_gFname << endl;
		avalanches.clear();
		if (!fill_avalanches(pNdFname, tau_gFname, avalanches, s(), factor()))
		{
			cerr << "Can't fill_avalanches at seed " << seed << endl;
			continue;
		}

		string tau_g_pNd_fname = fnameWp("tau_g_pNd", seed, fg(), fp(), ofs());
		ofstream tau_g_pNd(tau_g_pNd_fname);
		if (!init_tau_g_pNd(tau_g_pNd))
		{
			cerr << "Can't initialize " << tau_g_pNd_fname << endl;
			continue;
		}

		int tau_c = 0;
		bool firstpNd = true;
		pNd lastdef(0, 0, 0);
		for (auto& m_avalanche : avalanches)
		{
			if (calc_diam())
				m_avalanche.calcDiam_v3(s());

			if (calc_cog())
				m_avalanche.calcCoG(s());

			tau_g_pNd << m_avalanche << endl;
			if (allEvent() && allEventsTauExtMin() < m_avalanche.tau_ext && m_avalanche.tau_ext < allEventsTauExtMax())
				for (auto& m_pNd : m_avalanche.m_pNd)
					tau_g_pNd << m_pNd << "\n";

			while ((tau_c + 1) * tauIncr() < m_avalanche.tau_ext)
				++tau_c;

			if (tau_c >= tauIntervals())
				continue;

			if (calc_diam())
				avDiams.bins()[tau_c].add(m_avalanche.diam());

			if (calc_shapes())
			{
				m_avalanche.addShapeto(avShapesTot[tau_c]);
				int sizebin = getIndex(sizeIntervalsList, m_avalanche.DG());
				if (0 <= sizebin && sizebin < sizeIntervalsCount)
				{
					m_avalanche.addShapeto(avShapes[tau_c][sizebin]);
					++avShapes_counter[tau_c][sizebin];
				}
			}
			
			if (calc_distances())
			{
				for (int i = 1; i <= actionTime() && i < static_cast<int>(m_avalanche.m_pNd.size()); ++i)
					nthDistance[tau_c][i - 1](
					dist(m_avalanche.m_pNd[i].posx,m_avalanche.m_pNd[0].posx, s()),
					dist(m_avalanche.m_pNd[i].posy, m_avalanche.m_pNd[0].posy, s())) += m_avalanche.m_pNd[i].deform;

				for (int i = 0; i < static_cast<int>(m_avalanche.m_pNd.size()); ++i) // to calculate distance maps
				{
					// to calculate only the consecutive deformation distances
					if (firstpNd)
					{
						lastdef = m_avalanche.m_pNd[0];
						firstpNd = false;
						continue;
					}
					int distx = dist(m_avalanche.m_pNd[i].posx, lastdef.posx, s());
					int disty = dist(m_avalanche.m_pNd[i].posy, lastdef.posy, s());
					dists[tau_c](distx, disty) += 1;
					lastdef = m_avalanche.m_pNd[i];
				}
			}
			
		}
		if (calc_distances())
		{
			// to calculate every avalanche in a stress range
			tau_c = 0;
			for (int i = 0; i != static_cast<int>(avalanches.size()); ++i)
			{
				while ((tau_c + 1) * tauIncr() < avalanches[i].tau_ext)
					++tau_c;
				if (tau_c >= tauIntervals())
					continue;

				for (int j = 0; j != static_cast<int>(avalanches[i].m_pNd.size()); ++j)
				{
					for (int i2 = 0; i2 <= i; ++i2)
					{
						for (int j2 = 0; j2 < static_cast<int>(avalanches[i2].m_pNd.size()); ++j2)
						{
							if (i == i2 && j == j2)
								break;
							int distx = dist(avalanches[i].m_pNd[j].posx, avalanches[i2].m_pNd[j2].posx, s());
							int disty = dist(avalanches[i].m_pNd[j].posy, avalanches[i2].m_pNd[j2].posy, s());
							alldists[tau_c](distx,disty) += 1;
						}
					}
				}
			}
		}

	}
	cout << endl;
#pragma endregion

#pragma region creating avalanche diameter histogram file
	if (calc_diam())
	{
		string diamFname = fnameWp("pNd", fg(), fp(), st(), se(), "_diam-stress" + ofs());
	
		ofstream diam(diamFname);
		diam << "# tau_bin_left(1)\t"
			 << "tau_bin_right(2)\t"
			  << "nof avalances(3)\t"
			 << "tot diamx(4)\t"
			 << "tot diamy(5)\t"
			 << "average diamx(6)\t"
			 << "average diamy(7)" << endl;
		cout << "Writing " << diamFname << endl;
		diam << avDiams.use_averaged();
	}
	
#pragma endregion

#pragma region creating avalanche shape histogram files
	if (calc_shapes())
	{
		for (int i = 0; i < tauIntervals(); ++i)
		{
			for (int j = 0; j < sizeIntervalsCount; ++j)
			{
				stringstream shapefnamesuffix;
				shapefnamesuffix << "_"
					<< i*tauIncr() << "-"
					<< (i + 1)*tauIncr() << "_"
					<< sizeIntervalsList[j] << "-"
					<< sizeIntervalsList[j+1] << "_"
					<< "shape"
					<< ofs();
				string shapefname = fnameWp("pNd/average_av_shape", fg(), fp(), st(), se(), shapefnamesuffix.str());
				avShapes[i][j] /= avShapes_counter[i][j];
				avShapes[i][j].fWriteToFile(shapefname);
			}

			stringstream avShapesSizeAveragedFnameSuffix;
			avShapesSizeAveragedFnameSuffix << "_"
				<< i*tauIncr() << "-"
				<< (i + 1)*tauIncr() << "_"
				<< "shape"
				<< ofs();
			string avShapesSizeAveragedFname = fnameWp("pNd/average_av_shape", fg(), fp(), st(), se(), avShapesSizeAveragedFnameSuffix.str());
			cout << "Writing " << avShapesSizeAveragedFname << endl;
			avShapesTot[i].fWriteToFile(avShapesSizeAveragedFname);
		}
	}
	

#pragma endregion

#pragma region creating nth activity distance mapfiles
	if (calc_distances())
	{
		for (int i = 0; i < tauIntervals(); ++i)
		{
			for (int j = 1; j <= actionTime(); ++j)
			{
				stringstream nthactdstfnamesuffix;
				nthactdstfnamesuffix << "_"
					<< i*tauIncr() << "-"
					<< (i + 1)*tauIncr() << "_"
					<< j << "th_"
					<< "shape"
					<< ofs();
				string nthactdstfname = fnameWp("pNd/activity_distance", fg(), fp(), st(), se(), nthactdstfnamesuffix.str());
				nthDistance[i][j-1].fWriteToFile(nthactdstfname);
			}
		}
	}
	

#pragma endregion

#pragma region creating distance-map files
	if (calc_distances())
	{
		for (int i = 0; i < tauIntervals(); ++i)
		{
			stringstream cdmfns;
			cdmfns << "_"
				<< i*tauIncr() << "-"
				<< (i + 1)*tauIncr() << "_"
				<< "consecDistmap"
				<< ofs();
			string cdmfn = fnameWp("pNd/distance-maps", fg(), fp(), st(), se(), cdmfns.str());
			dists[i].fWriteToFile(cdmfn, "# I(i,j) equals the number of consecutive deformations with distance i in x and j in y direction in a specific stress region.\n");

			stringstream edmfns;
			edmfns << "_"
				<< i*tauIncr() << "-"
				<< (i + 1)*tauIncr() << "_"
				<< "everyDistmap"
				<< ofs();
			string edmfn = fnameWp("pNd/distance-maps", fg(), fp(), st(), se(), edmfns.str());
			dists[i].fWriteToFile(edmfn, "# I(i,j) equals the number of total deformations with distance i in x and j in y direction in a specific stress region.\n");
		}
	}
	
#pragma endregion
	
	cout << "Done." << endl;
	return 0;
}

