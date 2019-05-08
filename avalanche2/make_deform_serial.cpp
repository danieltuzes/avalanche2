// make_deform_serial.h : makes deformation serial until a given nominal def
//
// version number is in avalanche2_utils.h
// https://www.youtube.com/watch?v=egoY0wzheaM#t=1m13s Apa kezdõdik!

#include "stdafx.h"
#include "make_simulation.h"
#include "snapshot.h"

#define linPosXY getLinPos(x, y) 

using namespace std;

// find the minimum value, where the driving force to deform a cell is the biggest
int simVars::findMin(bool thermalActivation, double beta)
{
	deftype tau_th_min = numeric_limits<deftype>::infinity();
	int minx=0, miny=0;
	Z = 0;

	for (int x=0; x < size; ++x)
	{
		for (int y=0; y < size; ++y)
		{
			deftype tau_th;
			if (((tau_th = tau_p[linPosXY] - tau_ext - tau_l[linPosXY]) < tau_th_min))
			{
				minx = x;
				miny = y;
				tau_th_min = tau_th;
			}
			if (thermalActivation && tau_th_min > 0)
				Z += exp(-beta*tau_th);
			if ((tau_th = tau_n[linPosXY] + tau_ext + tau_l[linPosXY]) < 0 && tau_th < tau_th_min)
			{
				minx = x;
				miny = y;
				tau_th_min = tau_th;
			}
			if (thermalActivation && tau_th_min > 0)
				Z += exp(-beta*tau_th);

		}
	}

	return ((minx << power) + miny);
}

//refresh the tau_l field with the convolution of the stress field and the previous deformation placed at min with and findMin
int simVars::refresh(int minPoint, deftype deform, int randKernel, bool searchMinPoint, bool thermalActivation, double beta)
{
	deftype tau_th_min = numeric_limits<deftype>::infinity();
	int minx=0;
	int miny=0;
	Z = 0;

	int minPointx = getX(minPoint);
	int minPointy = getY(minPoint);

	if (randKernel == 1)
		sf.shuffle(rand_gen.getEngine(size*size - 2), true);
	else if (randKernel == 2)
		rand_gen.increase_times_advanced(sf.shuffleRestricted(rand_gen.getEngine(0)));
	
	for (int x = 0; x < size; ++x)
	{
		for (int y = 0; y < size; ++y)
		{
			tau_l[linPosXY] += sf.periodic(x - minPointx, y - minPointy) * deform; // [(((x - minPointx + size) & (size - 1)) << power) + ((y - minPointy + size) & (size - 1))] * deform);

			if (searchMinPoint)
			{
				deftype tau_th;
				if ((tau_th = tau_p[linPosXY] - tau_ext - tau_l[linPosXY]) < tau_th_min)
				{
					minx = x;
					miny = y;
					tau_th_min = tau_th;
				}
				if (thermalActivation && tau_th_min > 0)
					Z += exp(-beta*tau_th);

				if ((tau_th = tau_n[linPosXY] + tau_ext + tau_l[linPosXY]) < 0 && tau_th < tau_th_min)
				{
					minx = x;
					miny = y;
					tau_th_min = tau_th;
				}
				if (thermalActivation && tau_th_min > 0)
					Z += exp(-beta*tau_th);
			}

		}
	}

	return getLinPos(minx, miny);
}

// ez a template-ezés is rohadt felesleges, Árpi mondta ezt is, hiszek neki
bool simVars::makeDeformSerial(const simPars& pars, int nominalDef, vector<double>& extraSnaps, bool printToTau_g)
{
	time_t lastWrite = time(NULL);

	FILE * pNd_file;
	vector<pNd> pNdlist; // the whole pNd list is loaded at the beginning
	vector<pNd>::const_iterator nextMin;
	if (pars.pNd)
	{
		if ((pNd_file = fopen(fnameWp("pNd", pars).c_str(), "ab")) == NULL)
		{
			cerr << "Can't append " << fnameWp("pNd", pars) << endl;
			return false;
		}
	}
	else
		pNd_file = 0; // now it is initialized, compiler will not complain
	
	if (!pars.searchMinPoint) // minden egyes rohadt nominalDefnél beolvassa az egész rohadt pNd filet és elmegy addig, amíg kell
	{
		if (!get_pNd_content(fnameWp("pNd", pars), pNdlist, pars.s, 1))
		{
			cerr << "searchMinPoint is set to true but no pNd file was found, program skips this realisation." << endl;
			return false;
		}
		nextMin = pNdlist.begin();
		for (deftype def = 0; (def += nextMin->deform) - EPSILON < getDeform(); ++nextMin);
	}
	
	int minPoint;
	if (pars.searchMinPoint)
	{
		minPoint = findMin(pars.thermalActivation, pars.beta);
		if (avStartPos == -1)
			avStartPos = minPoint;

	}
	else
	{
		if (nextMin != pNdlist.end())
		{
			minPoint = nextMin->getLinPos(power);
			++nextMin;
		}
		else
		{
			cerr << "End of pNd file has been reached." << endl;
			return false;
		}
	}

	deftype local_s = tau_l[minPoint]; // the local stress at the minPoint

	printInfHeader(pars);

	ofstream tau_g;
	if (printToTau_g)
	{
		tau_g.open(fnameWp("tau_g", pars), ios_base::app);
		if (!tau_g)
		{
			cerr << "Cannot append to " << fnameWp("tau_g", pars) << " Program skips the rest of this simulation." << endl;
			return false;
		}
		tau_g.precision(17); //http://en.wikipedia.org/wiki/Double-precision_floating-point_format#IEEE_754_double-precision_binary_floating-point_format:_binary64
	}
	
	while (getDeform() < nominalDef && tau_ext < numeric_limits<double>::infinity())
	{
		if (DEBUG_TEST) printInf(minPoint, *this, pars);

		if (willYield(*this, minPoint, local_s) == direction::none) //azért szükséges ellenõrizni, mert könnyen lehet, hogy épp lavina van folyamatban, mert snapshotonként megáll és újraindul ez a függvény
		{
			if (!pars.thermalActivation)
				tau_ext = max(tau_p[minPoint] - local_s + EPSILON, pars.tau_extMin);
			else
				tau_ext = pars.tau_extMin;
			av_size = 0;
			av_size_n = 0;
			avStartPos = minPoint;
			nofAvs++;
			if (!extraSnaps.empty())
				for (; !extraSnaps.empty() && extraSnaps.front() < tau_ext; extraSnaps.erase(extraSnaps.begin()))
				{
					tau_l.writeToFile(fnameWp("tau_l", pars.s, pars.seed, pars.fg, pars.fp, "tau_ext_" + to_string(extraSnaps.front())));
					tau_p.writeToFile(fnameWp("tau_p", pars.s, pars.seed, pars.fg, pars.fp, "tau_ext_" + to_string(extraSnaps.front())));
				}
		}

		do
		{
			if (static_cast<int>(time(0) - lastWrite) > pars.tbs) //Gamma and tau_l are in consistent state, not as like at the end of this loop
			{
				snapshot bak(*this,static_cast<int>(getDeform()));
				bak.make(pars);
				delMarkedSnapshots(pars);
				bak.mark(pars);
				lastWrite = time(0);
			}
			if (static_cast<int>(time(0) - lastWrite) > pars.mrt)
			{
				cerr << "Simulation reached the maximum simulation time, program skips this realisation." << endl;
				return false;
			}
			if (nofAvs > pars.mna)
			{
				cerr << "Simulation reached the maximum number of avalanches, program skips this realisation." << endl;
				return false;
			}

			direction nowFlow = willYield(*this, minPoint, local_s);
			if (pars.thermalActivation && nowFlow == direction::none)
			{
				
				nowFlow = MonteCarloChoseNwait(pars.beta, minPoint, threshold);
				avStartPos = minPoint;
				Zstartval = Z;
				local_s = tau_l[minPoint];
				if (DEBUG_TEST && pars.info == infoMode::everyStep)
				{
					cout << "x: " << getX(minPoint) << "\t"
						<< "y: " << getY(minPoint) << "\t"
						<< "th: " << tau_p[minPoint] - tau_l[minPoint] - tau_ext << "\t"
						<< "exp(-beta*th): " << exp(-pars.beta*(tau_p[minPoint] - tau_l[minPoint] - tau_ext)) << "\t"
						<< "Z: " << Z << endl;

					cout << "x: " << getX(minPoint) << "\t"
						<< "y: " << getY(minPoint) << "\t"
						<< "th: " << threshold << "\t"
						<< "exp(-beta*th): " << exp(-pars.beta*threshold) << "\t"
						<< "Z: " << Z << "\t"
						<< "dir: " << nowFlow << endl;
				}
			}
				
			deftype deform = calcDeform(minPoint, local_s, tau_ext, sf[0], pars.DG, nowFlow);
			Gamma[minPoint] += deform;
			tau_ext -= min(deform * pars.sc, tau_ext);
			if (pars.dr)
				tau_ext += tau_ext / (pars.s*pars.s - (sG + sGn) / pars.DG.val);


			if (pars.pNd)
			{
				fwrite(&minPoint, sizeof(int), 1, pNd_file);
				fwrite(&deform, sizeof(deftype), 1, pNd_file);
			}
				
			if (nowFlow == direction::right)
			{
				av_size += deform;
				sG += deform;
			}
			else
			{
				av_size_n -= deform;
				sGn -= deform;
			}


			if (pars.fsrg.MyType != flowStressReGeneration::no) //regenerate flow stress
			{
				if (nowFlow == direction::right)
				{
					switch (pars.fsrg.MyType)
					{
					case flowStressReGeneration::independent:
						tau_p[minPoint] = rand_gen();
						break;
					case flowStressReGeneration::infinity:
						tau_p[minPoint] = numeric_limits<deftype>::infinity();
						break;
					case flowStressReGeneration::linear:
						tau_p[minPoint] = rand_gen() * (1 + abs(Gamma[minPoint]) * pars.fsrg.measure);
						if (tau_p[minPoint] <= 0)
						{
							cerr << "Cell (" << (minPoint >> power) << "," << (minPoint & (size - 1)) << ") reached zero flowstress before end of simulation." << endl;
							if (printToTau_g)
								tau_g << tau_ext << "\t" << sG << "\t" << sGn << "\t" << av_size << "\t" << av_size_n << "\t" << rand_gen.getTimesAdvanced() << " # reached zero flowstress before end of simulation" << endl;
							fclose(pNd_file);
							return false;
						}
						break;
					case flowStressReGeneration::exponential:
						tau_p[minPoint] = rand_gen() * pow(pars.fsrg.measure, abs(Gamma[minPoint]));
						break;
					case flowStressReGeneration::no: //take no action
						break;
					default:
						cerr << "Not handled case.";
					}
					if (pars.lfs.MyType == leftFlowStress::same)
						tau_n[minPoint] = tau_p[minPoint];
				}
				else // nowFlow == direction::left
				{
					switch (pars.fsrg.MyType)
					{
					case flowStressReGeneration::independent:
						tau_n[minPoint] = rand_gen();
						break;
					case flowStressReGeneration::infinity:
						tau_n[minPoint] = numeric_limits<deftype>::infinity();
						break;
					case flowStressReGeneration::linear:
						tau_n[minPoint] = rand_gen() * (1 + abs(Gamma[minPoint]) * pars.fsrg.measure);
						if (tau_n[minPoint] <= 0)
						{
							cerr << "Cell (" << (minPoint >> power) << "," << (minPoint & (size - 1)) << ") reached zero flowstress before end of simulation." << endl;
							if (printToTau_g)
								tau_g << tau_ext << "\t" << sG << "\t" << sGn << "\t" << av_size << "\t" << av_size_n << "\t" << rand_gen.getTimesAdvanced() << " # reached zero flowstress before end of simulation" << endl;
							fclose(pNd_file);
							return false;
						}
						break;
					case flowStressReGeneration::exponential:
						tau_n[minPoint] = rand_gen() * pow(pars.fsrg.measure, abs(Gamma[minPoint]));
						break;
					case flowStressReGeneration::no: //take no action
						break;
					default:
						cerr << "Not handled case.";
					}
					if (pars.lfs.MyType == leftFlowStress::same)
						tau_p[minPoint] = tau_n[minPoint];
				}
			}

			if (pars.maxLocG > 0 && Gamma[minPoint] > pars.maxLocG)
			{
				cerr << "Cell (" << (minPoint >> power) << "," << (minPoint & (size - 1)) << ") reached the maximum amount of allowed strain: " << pars.maxLocG << " Program stops this simulation." << endl;
				if (printToTau_g)
					placeLock_tau_g(*this, pars, " # reached maxLocG local strain.");
				if (pars.pNd)
					fclose(pNd_file);
				return false;
			}

			if (pars.searchMinPoint)
				minPoint = refresh(minPoint, deform, pars.randKernel, true, pars.thermalActivation, pars.beta);
			else
			{
				if (nextMin != pNdlist.end())
				{
					refresh(minPoint, deform, pars.randKernel, false,pars.thermalActivation, pars.beta);
					minPoint = nextMin->getLinPos(power);
					++nextMin;
				}
				else if (getDeform() < nominalDef)
				{
					cerr << "End of pNd file has been reached." << endl;
					return false;
				}
				else
					minPoint = refresh(minPoint, deform, pars.randKernel, true, pars.thermalActivation, pars.beta);
			}
			local_s = tau_l[minPoint];

			if (DEBUG_TEST) printInf(minPoint,*this,pars);
		} while (willYield(*this, minPoint,local_s) != direction::none && getDeform() < nominalDef);

		if (willYield(*this, minPoint, local_s) == direction::none && printToTau_g)
			printTau_g(tau_g, pars.thermalActivation);
	}
	
	if (pars.makeSN && ( !snapshot(*this,nominalDef).make(pars) || !delMarkedSnapshots(pars)))
		return false;
	
	if (pNd_file != 0)
		fclose(pNd_file);

	return true;
}