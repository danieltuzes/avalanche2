
#pragma once

#include "stdafx.h"
#include "../gen_utils/gen_version.h"
#include "../gen_utils/gen_fmanip.h"
#include "../gen_utils/gen_param.h"


class g_block_data : public block_data<double>
{
public:
	g_block_data();  //constructor inheritance is not allowed in VS2012
	g_block_data(const g_block_data& tocopy);  //constructor inheritance is not allowed in VS2012
	g_block_data(int s);  //constructor inheritance is not allowed in VS2012
	g_block_data(int s, double val);  //constructor inheritance is not allowed in VS2012
	template<typename U> g_block_data(const block_data<U>& tocopy);
	g_block_data(string fname, dataformat f);  //constructor inheritance is not allowed in VS2012
	g_block_data(char * fname, dataformat f);  //constructor inheritance is not allowed in VS2012
	g_block_data(vector<vector<double>> m_vector);

	template <typename U> g_block_data& operator=(const block_data<U>& tocopy);

	g_block_data& modulateXSin(double lambda);

	g_block_data& origin_imported(string fname, int skip = 0);

	g_block_data tilize() const;

	//string stats() const;
};