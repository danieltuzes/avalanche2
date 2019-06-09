// corr.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "corr_utils.h" //to be compatible with precompiled headers, version for this cpp file can be found in that header file
#include "../for_gnuplot/for_gnuplot_utils.h"
#include "../gen_utils/gen_fnameWp.h"

using namespace std;


//defining the used parameters

int main(int argc, char * argv[])
{
#pragma region setting parameters
	intro();
	paramcontainer p, p_temp;
	param<corr_type> corr("correlation",p,rq::must);  // enum corr_type
	//X, //cross-correlation
	//A, //auto-correlation
	//F,  //frequencies, fourier space (abs square)
	//DF // frequency differences, files must be given in pairs by fs_o or fs_file_o
	param<bool> print_maps("printOutEveryMaps",p,rq::optional,false);
	param<int> st("seedStart",p,rq::must);
	param<int> se("seedEnd",p,rq::must);
	param<bool> sa("substractAverage", p, rq::optional, false);  // substract average
	param<int> fg("file_grouping",p,rq::optional,1000);
	
	param<string> fp("fnamePrefix",p,rq::optional,"");
	param<string> fp_file("fnamePrefixFname", p, rq::optional); 
	param<string> fs("fnameSuffix", p, rq::optional, "");
	param<string> fs_file("fnameSuffixFname", p, rq::optional);
	param<string> fs_o("fnameSuffixOther", p, rq::optional, "");
	param<string> fs_file_o("fnameSuffixFnameOther", p, rq::optional);

	param<string> mn("mapName",p,rq::must); //map name
	param<string> mn_o("mapNameOther",p,rq::optional,"");
	param<bool> mapOutput("mapOutput", p, rq::optional, false);
	param<bool> projectedOutput("projectedOutput", p, rq::optional, false);
	param<double> treatNaNas("treatNaNas", p, rq::optional, 1); // in case nan occurs, print it out as a treatNaNas number
	param<int> reMap("reMap",p,rq::optional,32); // remapping the field until reMap x reMap size

	param<string> paramfname("paramfname",p_temp,rq::optional,"corr_param.ini",argc,argv);

	p.setFromIf(paramfname());
	if (p.setFromIf(argc,argv) * 2 != argc - 1)
	{
		p.summerize();
		cerr << "Error!" << endl;
		p.description();
		cerr << "Unknow parameter or missing value was found in program call, program terminates." << endl;
		return 0;
	}

	if (!p.summerize())
	{
		p.description();
		return 0;
	}

	vector<string> fnamePrefixes;
	if (fp_file.setByUsr())
	{
		ifstream ifile(fp_file());
		for (string tmp; comment_skipped(ifile) >> tmp; fnamePrefixes.push_back(tmp));
		if (fnamePrefixes.empty())
		{
			cerr << "Cannot extract fnamePrefix values from file " << fp_file() << "\n. Program terminates.\n";
			return 0;
		}
	}
	else
		fnamePrefixes.push_back(fp());

	vector<string> fnameSuffixes;
	if (fs_file.setByUsr())
	{
		ifstream ifile(fs_file());
		for (string tmp; comment_skipped(ifile) >> tmp; fnameSuffixes.push_back(tmp));
		if (fnameSuffixes.empty())
		{
			cerr << "Cannot extract fnameSuffix values from file " << fs_file() << "\n. Program terminates.\n";
			return 0;
		}
	}
	else
		fnameSuffixes.push_back(fs());

	vector<string> fnameSuffixes_o;
	if (corr == corr_type::D)
	{
		if (fs_file_o.setByUsr())
		{
			ifstream ifile(fs_file_o());
			for (string tmp; comment_skipped(ifile) >> tmp; fnameSuffixes_o.push_back(tmp));
			if (fnameSuffixes_o.empty())
			{
				cerr << "Cannot extract fnameSuffixOther values from file " << fs_file_o() << "\n. Program terminates.\n";
				return 0;
			}
		}
		else
			fnameSuffixes_o.push_back(fs_o());

		if (fnameSuffixes.size() != fnameSuffixes.size())
		{
			cerr << "size of fnameSuffixes and fnameSuffixes_o must be equal, but " << fnameSuffixes.size() << " != " << fnameSuffixes_o.size() << " Program terminates." << endl;
			return 0;
		}
	}

#pragma endregion

	for (auto fnamePrefix : fnamePrefixes)
	{
		fp.modify(fnamePrefix);

		for (int i = 0; i < static_cast<int>(fnameSuffixes.size()); ++i)
		{
			string fnameSuffix = fnameSuffixes[i];
			fs.modify(fnameSuffix);

			block_data<> sum, sum_o, diff, ratio, ratio_c;
			FT_block_data FT_sum, FT_sum_o, FT_diff;
			int counter = 0; //number of elements added to block_data sum

			if (corr() == corr_type::A) //if it has to calculate autocorrelation
			{
				// read in the files and calculate autocorrelations
				for (int seed = st(); seed < se(); ++seed)
				{
					string ifname = fnameWp(mn(), seed, fg(), fp(), fs());

					cout << "Using " << ifname << " for autocorrelation ..." << endl;
					if (!is_file_exist(ifname))
					{
						cerr << ifname << " doesn't exist. Program skips this seed number." << endl;
						continue;
					}

					ft_block_data d(ifname, dataformat::binary);

					if (sa())
						d -= static_cast<double>(d.getAverage());

					d.replaceAutoCorrelation();

					if (print_maps())
						d.fWriteToFile(ifname + "_autocorr" + fs() + ".txt");

					if (!sum.hasData()) //sum is still not initialized until an autocorrelation is calculated
						sum.init(d.getSize(), 0);

					sum += d;
					++counter;
				}

				string ofname = fnameWp(mn(), fp(), st(), se(), fs() + "_autocorr" + fs());
				string ofname_tilized = fnameWp(mn(), fp(), st(), se(), fs() + fs() + "_autocorr_tilized" + fs());

				sum /= counter;
				//cout << "Write results to " << ofname << endl;
				sum.fWriteToFile(ofname);

				g_block_data(sum).tilize().fWriteToFile(fnameWp(mn(), fp(), st(), se(), fs() + "_autocorr_tilized" + fs()));

				if (reMap.setByUsr())
				{
					for (block_data<> reMapper(sum); reMapper.getSize() >= reMap(); reMapper = reMapper.reMap())
					{
						string reMapOfname = fnameWp(mn(), fp(), st(), se(), fs() + "_autocorr_reMap_" + to_string(reMapper.getSize()) + fs());
						cout << "Write reMap at size " << reMapper.getSize() << " to " << reMapOfname << endl;
						reMapper.fWriteToFile(reMapOfname);
					}
				}

			};

#if 0
			if (corr() == corr_type::F) //if it has to calculate the frequencies
			{
				// read in the files and calculate autocorrelations
				for (int seed = st(); seed < se(); ++seed)
				{
					string ifname = fnameWp(mn(), seed, fg(), fp(), fs());

					cout << "Using " << ifname << " for fourier space abs val square ..." << endl;
					if (!is_file_exist(ifname))
					{
						cerr << ifname << " doesn't exist. Program skips this seed number." << endl;
						continue;
					}

					ft_block_data d(ifname.c_str(), dataformat::binary);

					if (sa())
						d -= static_cast<double>(d.getAverage());

					d.replaceFourierAbsValSq();
					//d.sumNormalize();

					if (print_maps())
						d.fWriteToFile(ifname + "_freqs" + fs() + ".txt");

					if (!sum.hasData()) //sum is still not initialized until an autocorrelation is calculated
					{
						sum.init(d.getSize(), 0);
					}

					sum += d;
					++counter;
				}
				if (counter == 0)
				{
					cerr << "No files were found for " << fp() << ". Program skips this file." << endl;
					continue;
				}
				sum /= counter;


				string ofname = fnameWp(mn(), fp(), st(), se(), fs() + "_freq" + fs());

				if (!onlyXslice())
					sum.fWriteToFile(ofname);

				ofname += "_projectx.dat";
				ofstream(ofname) << linvec<>(sum.projectx()).printAsCol(true) << "\n";
			};
#endif

#if 1
			if (corr() == corr_type::F) //if it has to calculate the frequencies
			{
				// read in the files and calculate autocorrelations
				for (int seed = st(); seed < se(); ++seed)
				{
					string ifname = fnameWp(mn(), seed, fg(), fp(), fs());

					cout << "Using " << ifname << " for fourier space abs val square ..." << endl;
					if (!is_file_exist(ifname))
					{
						cerr << ifname << " doesn't exist. Program skips this seed number." << endl;
						continue;
					}

					block_data<double> signal(ifname, dataformat::binary);

					if (sa())
						signal -= signal.getAverage();

					FT_block_data FT_signal = signal;
					FT_signal.forward_transform();

					if (!sum.hasData()) //sum is still not initialized until an autocorrelation is calculated
					{
						sum.init(signal.getSize(), 0);
					}

					FT_signal.apply_complex2real_function(abs);
					signal = FT_signal.getReal_block_data();

					if (print_maps())
						signal.fWriteToFile(ifname + "_freqs.txt");

					sum += signal;
					++counter;
				}
				if (counter == 0)
				{
					cerr << "No files with seed[" << st() << "," << se() << "] found for fnamePrefix" << fp() << " and fnameSuffix" << fs() << ". No output was created." << endl;
					continue;
				}
				sum /= counter;
				fillUpperHalf(sum);

				string ofname = fnameWp(mn(), fp(), st(), se(), fs() + "_freq");

				if (mapOutput)
					sum.fWriteToFile(ofname);
				if (projectedOutput)
					ofstream(ofname + "_projected_to_x.dat") << linvec<>(sum.projectx()).printAsCol(true) << "\n";

			};
#endif
			if (corr() == corr_type::D)
			{
				fs_o.modify(fnameSuffixes_o[i]);
				// read in the files and calculate autocorrelations
				for (int seed = st(); seed < se(); ++seed)
				{
					string ifname = fnameWp(mn(), seed, fg(), fp(), fs());
					string ifname_o = fnameWp(mn(), seed, fg(), fp(), fs_o());

					cout << "Using " << ifname << " and " << ifname_o << " for fourier space abs val square ..." << endl;
					if (!is_file_exist(ifname))
					{
						cerr << ifname << " doesn't exist. Program skips this seed number." << endl;
						continue;
					}
					if (!is_file_exist(ifname_o))
					{
						cerr << ifname_o << " doesn't exist. Program skips this seed number." << endl;
						continue;
					}

					block_data<double> signal(ifname, dataformat::binary);
					block_data<double> signal_o(ifname_o, dataformat::binary);
					block_data<double> diff_sig(signal_o.getSize(), 0);
					block_data<double> ratio_sig(signal_o.getSize(), 0);
					block_data<double> ratio_sigC(signal_o.getSize(), 0);

					if (sa())
					{
						signal -= signal.getAverage();
						signal_o -= signal_o.getAverage();
					}

					FT_block_data FT_signal = signal;
					FT_block_data FT_signal_o = signal_o;
					FT_signal.forward_transform();
					FT_signal_o.forward_transform();

					if (!sum.hasData()) //sum is still not initialized until an autocorrelation is calculated
					{
						sum.init(signal.getSize(), 0);
						sum_o.init(signal.getSize(), 0);
						diff.init(signal.getSize(), 0);
						ratio.init(signal.getSize(), 0);
						ratio_c.init(signal.getSize(), 0);
					}

					diff_sig = (FT_signal_o - FT_signal).apply_complex2real_function(abs).getReal_block_data();
					ratio_sig = (FT_signal_o / FT_signal).apply_complex2real_function(abs).getReal_block_data();

					FT_signal.apply_complex2real_function(abs);
					FT_signal_o.apply_complex2real_function(abs);

					ratio_sigC = (FT_signal_o / FT_signal).getReal_block_data();

					signal = FT_signal.getReal_block_data();
					signal_o = FT_signal_o.getReal_block_data();

					if (print_maps())
						signal.fWriteToFile(ifname + "_freqs.txt");

					sum += signal;
					sum_o += signal_o;
					diff += diff_sig;
					ratio += ratio_sig;
					ratio_c += ratio_sigC;
					++counter;
				}
				if (counter == 0)
				{
					cerr << "No files with seed[" << st() << "," << se() << "] found for fnamePrefix" << fp() << " and fnameSuffix" << fs() << ". No output was created." << endl;
					continue;
				}
				sum /= counter;
				sum_o /= counter;
				diff /= counter;
				ratio /= counter;
				ratio_c /= counter;

				string ofname = fnameWp(mn(), fp(), st(), se(), fs() + "_freq");
				string ofname_o = fnameWp(mn(), fp(), st(), se(), fs_o() + "_freq");
				string ofname_diff = fnameWp(mn(), fp(), st(), se(), fs() + "_diff" + fs_o() + "_freq");
				string ofname_diff_b = fnameWp(mn(), fp(), st(), se(), fs() + "_diffB" + fs_o() + "_freq");
				string ofname_ratio = fnameWp(mn(), fp(), st(), se(), fs() + "_ratio" + fs_o() + "_freq");
				string ofname_ratio_b = fnameWp(mn(), fp(), st(), se(), fs() + "_ratioB" + fs_o() + "_freq");
				string ofname_ratio_c = fnameWp(mn(), fp(), st(), se(), fs() + "_ratioC" + fs_o() + "_freq");

				if (mapOutput)
				{
					sum.fWriteToFileTreatNN(ofname, true, treatNaNas);
					sum_o.fWriteToFileTreatNN(ofname_o, true, treatNaNas);
					diff.fWriteToFileTreatNN(ofname_diff, true, treatNaNas);
					(sum_o - sum).fWriteToFileTreatNN(ofname_diff_b, true, treatNaNas);
					ratio.fWriteToFileTreatNN(ofname_ratio, true, treatNaNas);
					ratio_c.fWriteToFileTreatNN(ofname_ratio_c, true, treatNaNas);
					(sum_o / sum).fWriteToFileTreatNN(ofname_ratio_b, true, treatNaNas);
				}
				if (projectedOutput)
				{
					ofstream(ofname + "_projected_to_x.dat") << linvec<>(sum.projectx()).printAsCol(true) << "\n";
					ofstream(ofname_o + "_projected_to_x.dat") << linvec<>(sum_o.projectx()).printAsCol(true) << "\n";
					ofstream(ofname_diff + "_projected_to_x.dat") << linvec<>(diff.projectx()).printAsCol(true) << "\n";
					ofstream(ofname_diff_b + "_projected_to_x.dat") << linvec<>((sum_o - sum).projectx()).printAsCol(true) << "\n";
					ofstream(ofname_ratio + "_projected_to_x.dat") << linvec<>(ratio.projectx()).printAsCol(true) << "\n";
					ofstream(ofname_ratio_c + "_projected_to_x.dat") << linvec<>(ratio_c.projectx()).printAsCol(true) << "\n";
					ofstream(ofname_ratio_b + "_projected_to_x.dat") << linvec<>((sum_o / sum).projectx()).printAsCol(true) << "\n";
				}

			}

			if (corr() == corr_type::X) //if it has to calculate cross-correlation
			{
				for (int seed = st(); seed < se(); ++seed)
				{
					string ifname_f = fnameWp(mn(), seed, fg(), fp(), fs());
					string ifname_o = fnameWp(mn_o(), seed, fg(), fp(), fs());

					cout << "Using " << ifname_f << " and " << ifname_o << " for cross-correlation ..." << endl;
					if (!is_file_exist(ifname_f) || !is_file_exist(ifname_o))
					{
						cerr << ifname_f << " or " << ifname_o << " is not exist. Program skips this seed number." << endl;
						continue;
					}

					ft_block_data d_f(ifname_f, dataformat::binary);
					ft_block_data d_o(ifname_o, dataformat::binary);

					d_f.replaceXCorrelation(d_o);

					if (print_maps())
						d_f.fWriteToFile(ifname_f + "_X_" + mn_o() + "_map" + fs() + ".txt");

					if (!sum.hasData()) //sum is still not initialized until an autocorrelation is calculated
						sum = ft_block_data(d_f.getSize(), 0);

					sum += d_f;
					++counter;
				}

				string ofname = fnameWp(mn() + "_X_" + mn_o(), fp(), st(), se(), fs() + "_corr" + fs() + ".dat");

				sum /= counter;

				sum.fWriteToFile(ofname);
			}
		}



	}
	return 0;
}

