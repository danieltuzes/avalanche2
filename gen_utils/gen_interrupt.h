// to stop applications in a smarter way than sending sigterm signal to them

#pragma once
#include "stdafx.h"
#include "gen_fmanip.h"
#include "gen_siomanip.h"

using namespace std;

bool interrupt(string fnamePart);

bool interrupt(string fnamePart, int seed);
