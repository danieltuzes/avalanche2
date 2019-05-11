// gen_siomanip.h : stores user defined general function relatad to stream manipulation
//
#define GEN_SIOMANIP_VERSION 0.01

#pragma once
#include "stdafx.h"
using namespace std;




//skips empty and commented lines
istream& comment_skipped(istream& stream, const char comment_symbol = '#');


//skips empty and commented lines
void skip_comment(istream& stream, const char comment_symbol = '#');

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime();

string ord_num(int a);