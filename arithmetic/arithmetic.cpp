// test.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "../avalanche2/rnd_gen.h"
#include "../gen_utils/gen_fmanip.h"
#include "../for_gnuplot/for_gnuplot_utils.h"

using namespace std;

int main(int argc, char* argv[]) {

	cout << "This program calculates new field based on the given ones, the size must be 2^x time 2^x, x integer.\nThe program can do the followings:\n"
		<< "\tadd(field, field): matrix addition of two fields\n"
		<< "\tmultiply(field, number): pointwise product of the field and the number*identity\n"
		<< "\tsubstract(field, field): substract the second field from the first\n\n"
		<< "1st input: type of operator (add, multiply, substract)\n"
		<< "2nd input: fname of the 1st argument"
		<< "3rd input: fname of the 2nd argument" << endl;

	if (argc != 4)
		return 0;

	string type(argv[1]);
	string fnameA(argv[2]);
	g_block_data marrayA(fnameA, bin_text);

	g_block_data C(marrayA.getSize());
	string ofname;

	if (type.compare("add") == 0)
	{
		string fnameB(argv[3]);
		g_block_data marrayB(fnameB, bin_text);
		C = marrayA + marrayB;
		ofname = fnameA + "+" + fnameB;
	}

	if (type.compare("multiply") == 0)
	{
		istringstream factorstream(argv[3]);
		double factor;
		factorstream >> factor;
		C = factor * marrayA;
		ofname += fnameA + "_multiplied" + string(argv[3]);
	}

	if (type.compare("substract") == 0)
	{
		string fnameB(argv[3]);
		g_block_data marrayB(fnameB, bin_text);
		C = marrayA - marrayB;
		ofname = fnameA + "-" + fnameB;
	}

	C.writeToFile(ofname + "_bin.dat");
	C.fWriteToFile(ofname + ".txt");
	C.tilize().fWriteToFile(ofname + "_tilized.txt");

	cout << C.stats();

	return 0;
}