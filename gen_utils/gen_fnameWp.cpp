
#include "stdafx.h"
#include "gen_fnameWp.h"

basic_simVars::basic_simVars(int s, int seed, int fg, string fp) : s(s), seed(seed), fg(fg), fp(fp) {}

evalPars::evalPars(param<int> s, int seed, param<int> fg, param<string> fp, param<int> st, param<int> se, param<string> fs) :
	basic_simVars(s(), seed, fg(), fp()), st(st()), se(se()), fs(fs()) {}

string fnameWp(string fname, int s, int seed, int fg, string fp, string fs)
{
	return fnameWp(fname, s, fg, fp, seed, seed, fs);
}


string fnameWp(string fname, int s, int fg, string fp, int seedstart, int seedend, string fs)
{
	stringstream fname_wp_stream;
	fname_wp_stream << "sim_" << s << "/"
		<< fnameWp(fname, fg, fp, seedstart, seedend, fs);

	return fname_wp_stream.str();
}


string fnameWp(string fname, int seed, int fg, string fp, string fs)
{
	return fnameWp(fname, fg, fp, seed, seed, fs);
}

string fnameWp(string fname, string fp, int seedstart, int seedend, string fs)
{
	return fnameWp(fname, 0, fp, seedstart, seedend, fs);
}

string fnameWp(string fname, int fg, string fp, int seedstart, int seedend, string fs)
{
	stringstream fname_wp_stream;
	fname_wp_stream << fname << "/";

	if (seedstart == seedend)
		fname_wp_stream << static_cast<int>(seedstart / fg) * fg << "/"
		<< fp
		<< seedstart;
	else
		fname_wp_stream << fp
		<< seedstart << "-"
		<< seedend;

	fname_wp_stream << fs;
	if (!fname.compare("tau_l") ||
		!fname.compare("tau_n") ||
		!fname.compare("tau_p") ||
		!fname.compare("rho") ||
		!fname.compare("kappa") ||
		!fname.compare("tau_sc") ||
		!fname.compare("tau_b") ||
		!fname.compare("tau_d") ||
		!fname.compare("Gamma") ||
		!fname.compare("gamma") ||
		!fname.compare("tau_y_m") ||
		!fname.compare("tau_y_p") ||
		!fname.compare("tau_y_m_r") ||
		!fname.compare("tau_y_p_r") ||
		!fname.compare("pNd") ||
		!fname.compare("pNe")) //result of compare is 0 is they are equal
		fname_wp_stream << ".dat";
	else if (!fname.compare("log"))
		fname_wp_stream << ".log";
	else
		fname_wp_stream << ".txt";

	return fname_wp_stream.str();
}

string fnameWp(string fname, evalPars pars, string extra_fs)
{
	return fnameWp(fname, pars.fg, pars.fp, pars.st, pars.se, pars.fs + extra_fs);
}