#include "stdafx.h"

// From http://rosettacode.org/wiki/Sorting_algorithms/Radix_sort
// Note: the LSD radix sort uses the standard library std::stable_partition algorithm.
// This algorithm is guaranteed to preserve relative order and has a higher runtime cost.
// The MSD radix sort uses std::partition and can be significantly faster. 

// Radix sort comparator for 32-bit two's complement integers
class radix_test
{
	const int bit; // bit position [0..31] to examine
public:
	radix_test(int offset) : bit(offset) {} // constructor

	bool operator()(int value) const // function call operator
	{
		if (bit == 31) // sign bit
			return value < 0; // negative int to left partition
		else
			return !(value & (1 << bit)); // 0 bit to left partition
	}
};

// Least significant digit radix sort
void lsd_radix_sort(int *first, int *last)
{
	for (int lsb = 0; lsb < 32; ++lsb) // least-significant-bit
	{
		std::stable_partition(first, last, radix_test(lsb));
	}
}

// Most significant digit radix sort (recursive)
void msd_radix_sort(int *first, int *last, int msb = 31)
{
	if (first != last && msb >= 0)
	{
		int *mid = std::partition(first, last, radix_test(msb));
		msb--; // decrement most-significant-bit
		msd_radix_sort(first, mid, msb); // sort left partition
		msd_radix_sort(mid, last, msb); // sort right partition
	}
}