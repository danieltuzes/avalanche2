#pragma once
#include "stdafx.h"
#include "gen_fmanip.h"
#include "../for_gnuplot/for_gnuplot_utils.h"
#include "../fftw/fftw3.h"

#define FT_BLOCK_DATA_VERSION 1.1

// FT_BLOCK_DATA_VERSION 1.1: FT_block_data operator/(const FT_block_data& divisor) const; and FT_block_data& operator/=(const FT_block_data& divisor); are added
// FT_BLOCK_DATA_VERSION 1: multiplication does not imply normalisation, it has to be asked by multiply_normalised

enum class qunatity_type
{
	real,
	half_complex,
	complex
};

class FT_block_data
{
protected:
	qunatity_type data_type = qunatity_type::real;
	int msize;
	double * mdata;
	void in_place_forward_transform();
	void in_place_backward_transform();
public:
	FT_block_data();
	FT_block_data(int size);
	FT_block_data(int size, double val);
	FT_block_data(int size, const complex<double>& val);
	FT_block_data(const block_data<double>& in_data);
	~FT_block_data();

	/* Copy constructor */
	FT_block_data(const FT_block_data& other);

	/** Move constructor */
	FT_block_data(FT_block_data&& other);

	/** Copy assignment operator */
	FT_block_data& operator= (const FT_block_data& other);

	/** Move assignment operator */
	FT_block_data& operator= (FT_block_data&& other);

	FT_block_data operator+(const FT_block_data& add) const; // add the complex values
	FT_block_data& operator+=(const FT_block_data& add); // add the complex values

	FT_block_data operator-(const FT_block_data& subs) const; // substract the complex values
	FT_block_data& operator-=(const FT_block_data& subs); // substract the complex values

	FT_block_data operator*(const FT_block_data& multiplier) const; // multiply the complex values pointwise
	FT_block_data& operator*=(const FT_block_data& multiplier); // multiply the complex values pointwise

	FT_block_data operator/(const FT_block_data& divisor) const; // divide the complex values pointwise
	FT_block_data& operator/=(const FT_block_data& divisor); // divide the complex values pointwise


	FT_block_data& forward_transform(); // apply forward, unnormalised transformation
	FT_block_data& backward_transform(); // apply backward, unnormalised transformation

	// FT_block_data& sine_transform(const block_data<double>& odd_data); // odd around -0.5 and odd around size/2 - 0.5; FFTW_RODFT10 (DST-II)
	// not implemented, because it's rare that one need to apply sine_transform multiply times
	// but it could be useful because half the memory is needed, and it could fasten up complex multiplication because of the better memory-cache usage
	// please note that multiply_pure_imag is available, i.e. one could apply forward_transform() and then benefit a bit from the the fact that real part is 0

	FT_block_data& normalise(); // apply size*size normalisation

	FT_block_data& multiply_unnormalised(const FT_block_data& multiply); // multiply the complex values pointwise; operator*= redirects here
	FT_block_data& multiply_normalised(const FT_block_data& multiply); // multiply the complex values pointwise and apply size*size normalisation

	FT_block_data& multiply_pure_imag(const FT_block_data& multiply); // multiply the complex values pointwise if it is known that all conmplex values of multiply is imaginary

	block_data<double> get_block_data(); // returns a new block_data object
	block_data<double> getReal_block_data(); // returns a new block_data object
	block_data<double> getImag_block_data(); // returns a new block_data object

	FT_block_data& apply_complex2real_function(double(*function)(const complex<double>& x)); // apply the desired real valued complex function on each element; the data is still stored as complex numbers with imaginary part = 0
	
	FT_block_data& apply_complex2complex_function(complex<double>(*function)(const complex<double>& x)); // apply the desire complex valued complex function on each element
};
