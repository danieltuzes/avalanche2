// version.h : stores user defined utilities for version handling
//
#pragma once
#include "stdafx.h"
#define GEN_VERSION_VERSION 0.06b

// GEN_VERSION_VERSION 0.06b: bug corrected not showing extra info
// GEN_VERSION_VERSION 0.06: intro can skip extraInfo
// GEN_VERSION_VERSION 0.05: suggested to use at least Visual Studio 2013

#define STR(x) #x
#define XSTR(x) STR(x)

#ifdef __GNUC__
#define COMPILER_VERSION GNUC __GNUC__.__GNUC_MINOR__.__GNUC_PATCHLEVEL__
#define MACHINE_INFO GCC_MACHINE_INFO
#endif

#ifdef _MSC_VER
#define COMPILER_VERSION MSVSC++  _MSC_VER
#define MACHINE_INFO A windows machine
#endif


#ifdef _DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

using namespace std;

void intro(string header, string footer, bool printExtraInfo = true, char indent_char = '#', int indent_size = 4, int wide = 75);