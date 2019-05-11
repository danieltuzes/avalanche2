#define GEN_fnameWp_VERSION 0.03

// GEN_fnameWp_VERSION 0.03: evalPars is moved here from evaluate.h, so statistics has no dependence on shape
// GEN_fnameWp_VERSION 0.02: pNe is added to dat file extension

#pragma once
#include "stdafx.h"
#include "gen_param.h"
using namespace std;

class basic_simVars
{
public:
	basic_simVars(int s, int seed, int fg, string fp);
	int s;      // system size
	int seed;          // seed value
	int fg;     // filename grouping
	string fp;  // filename prefix
};

class evalPars : public basic_simVars
{
public:
	evalPars(param<int> s, int seed, param<int> fg, param<string> fp, param<int> st, param<int> se, param<string> fs);
	int st;
	int se;
	string fs;
};

string fnameWp(string fname, evalPars pars, string extra_fs);
string fnameWp(string fname, int s, int seed, int fg, string fp, string fs);
string fnameWp(string fname, int s, int fg, string fp, int seedstart, int seedend, string fs);
string fnameWp(string fname, int fg, string fp, int seedstart, int seedend, string fs);
string fnameWp(string fname, string fp, int seedstart, int seedend, string fs);
string fnameWp(string fname, int seed, int fg, string fp, string fs);