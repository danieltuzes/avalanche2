#define GEN_fnameWp_VERSION 0.02

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

string fnameWp(string fname, int s, int seed, int fg, string fp, string fs);
string fnameWp(string fname, int s, int fg, string fp, int seedstart, int seedend, string fs);
string fnameWp(string fname, int fg, string fp, int seedstart, int seedend, string fs);
string fnameWp(string fname, string fp, int seedstart, int seedend, string fs);
string fnameWp(string fname, int seed, int fg, string fp, string fs);