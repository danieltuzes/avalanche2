#include "stdafx.h"
#include "gen_fmanip.h"
#include "FT_block_data.h"
#include "../for_gnuplot/for_gnuplot_utils.h"
#include "../fftw/fftw3.h"


void FT_block_data::in_place_forward_transform()
{
	if (data_type != qunatity_type::real)
	{
		cerr << "Cannot fourier transform forward a non-real quantity. Program terminates." << endl;
		exit(0);
	}
	fftw_plan forward_plan = fftw_plan_dft_r2c_2d(msize, msize, mdata, reinterpret_cast<fftw_complex*>(mdata), FFTW_ESTIMATE);
	fftw_execute(forward_plan);
	data_type = qunatity_type::complex;
}
void FT_block_data::in_place_backward_transform()
{
	if (data_type != qunatity_type::complex)
	{
		cerr << "Cannot fourier transform back a non-complex quantity. Program terminates." << endl;
		exit(0);
	}
	fftw_plan backward_plan = fftw_plan_dft_c2r_2d(msize, msize, reinterpret_cast<fftw_complex*>(mdata), mdata, FFTW_ESTIMATE);
	fftw_execute(backward_plan);
	data_type = qunatity_type::real;
}

FT_block_data::FT_block_data() : msize(0), mdata(NULL) {}
FT_block_data::FT_block_data(int size) : msize(size), mdata(fftw_alloc_real(msize*(msize + 2))) {}
FT_block_data::FT_block_data(int size, double val) : FT_block_data(size)
{
	fill(mdata, mdata + msize*(msize + 2), val);
}
FT_block_data::FT_block_data(int size, const complex<double>& val) : FT_block_data(size)
{
	data_type = qunatity_type::complex;
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
		reinterpret_cast<complex<double>*>(mdata)[i] = val;
}
FT_block_data::FT_block_data(const block_data<double>& in_data) : FT_block_data(in_data.getSize())
{
	for (int i = 0; i < msize; ++i)
	{
		memcpy(mdata + i*(msize + 2), in_data.getData() + i*msize, msize*sizeof(double));
	}
}
FT_block_data::~FT_block_data()
{
	fftw_free(mdata);
}

/* Copy constructor */
FT_block_data::FT_block_data(const FT_block_data& other) : FT_block_data(other.msize)
{
	data_type = other.data_type;
	memcpy(mdata, other.mdata, msize*(msize + 2)*sizeof(double));
}

/** Move constructor */
FT_block_data::FT_block_data(FT_block_data&& other) : data_type(other.data_type), msize(other.msize), mdata(other.mdata)
{
	other.mdata = NULL;
}

/** Copy assignment operator */
FT_block_data& FT_block_data::operator= (const FT_block_data& other)
{
	FT_block_data tmp(other);         // re-use copy-constructor
	*this = move(tmp); // re-use move-assignment
	return *this;
}

/** Move assignment operator */
FT_block_data& FT_block_data::operator= (FT_block_data&& other)
{
	msize = other.msize;
	delete[] mdata;
	mdata = other.mdata;
	data_type = other.data_type;
	other.mdata = NULL;
	return *this;
}


FT_block_data FT_block_data::operator+(const FT_block_data& add) const // pointwise complex mulitplication
{
	FT_block_data tmp(*this);
	return tmp = add;
}

FT_block_data& FT_block_data::operator+=(const FT_block_data& add) // pointwise mulitplication and division by size*size
{
	if (mdata == add.mdata)
	{
		cerr << "Operator += is not tested when two matrices are the same. Program terminates." << endl;
		exit(0);
	}
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
		reinterpret_cast<complex<double>*>(mdata)[i] += reinterpret_cast<complex<double>*>(add.mdata)[i];
	return *this;
}

FT_block_data FT_block_data::operator-(const FT_block_data& add) const // pointwise complex mulitplication
{
	FT_block_data tmp(*this);
	return tmp -= add;
}

FT_block_data& FT_block_data::operator-=(const FT_block_data& add) // pointwise mulitplication and division by size*size
{
	if (mdata == add.mdata)
	{
		cerr << "Operator += is not tested when two matrices are the same. Program terminates." << endl;
		exit(0);
	}
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
		reinterpret_cast<complex<double>*>(mdata)[i] -= reinterpret_cast<complex<double>*>(add.mdata)[i];
	return *this;
}

FT_block_data FT_block_data::operator*(const FT_block_data& multiplier) const // pointwise complex mulitplication
{
	FT_block_data tmp(*this);
	return tmp *= multiplier;
}

FT_block_data& FT_block_data::operator*=(const FT_block_data& multiplier) // pointwise mulitplication
{
	this->multiply_unnormalised(multiplier);
	return *this;
}

FT_block_data FT_block_data::operator/(const FT_block_data& divisor) const // pointwise complex division
{
	FT_block_data tmp(*this);
	return tmp /= divisor;
}

FT_block_data& FT_block_data::operator/=(const FT_block_data& divisor) // pointwise complex division
{
	if (mdata == divisor.mdata)
	{
		cerr << "Operator /= is not tested when two matrices are the same. Program terminates." << endl;
		exit(0);
	}
	if (data_type != qunatity_type::complex && divisor.data_type != qunatity_type::complex)
	{
		cerr << "Complex operator / is used on real data. Program terminates." << endl;
		exit(0);
	}
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
		reinterpret_cast<complex<double>*>(mdata)[i] /= reinterpret_cast<complex<double>*>(divisor.mdata)[i];
	return *this;
}




FT_block_data& FT_block_data::forward_transform()
{
	in_place_forward_transform();
	return *this;
}

FT_block_data& FT_block_data::backward_transform()
{
	in_place_backward_transform();
	return *this;
}

/*
FT_block_data& FT_block_data::sine_transform(const block_data<double>& odd_data) // odd around -0.5 and odd around size/2 - 0.5; FFTW_RODFT10 (DST-II)
{
	// not implemented, because it's rare that one need to apply sine_transform multiply times
	// but it could be useful because half the memory is needed, and it could fasten up complex multiplication because of the better memory-cache usage
	// please note that multiply_pure_imag is available, i.e. one could apply forward_transform() and then benefit a bit from the the fact that real part is 0

	return *this;
}
*/
FT_block_data& FT_block_data::normalise()
{
	for (int i = 0; i < msize*(msize + 2); ++i)
		mdata[i] /= (msize*msize);
	return *this;
}

FT_block_data& FT_block_data::multiply_unnormalised(const FT_block_data& multiply)
{
	if (data_type != qunatity_type::complex && multiply.data_type != qunatity_type::complex)
	{
		cerr << "Complex operator * is used on real data. Program terminates." << endl;
		exit(0);
	}
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
		reinterpret_cast<complex<double>*>(mdata)[i] *= reinterpret_cast<complex<double>*>(multiply.mdata)[i];
	return *this;
}

FT_block_data& FT_block_data::multiply_normalised(const FT_block_data& multiply)
{
	if (mdata == multiply.mdata)
	{
		cerr << "Operator *= is not tested when two matrices are the same. Program terminates." << endl;
		exit(0);
	}
	if (data_type != qunatity_type::complex && multiply.data_type != qunatity_type::complex)
	{
		cerr << "Complex operator * is used on real data. Program terminates." << endl;
		exit(0);
	}
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
		reinterpret_cast<complex<double>*>(mdata)[i] *= reinterpret_cast<complex<double>*>(multiply.mdata)[i] / (msize*msize + 0.);
	return *this;
}

FT_block_data& FT_block_data::multiply_pure_imag(const FT_block_data& multiply)
{
	if (data_type != qunatity_type::complex && multiply.data_type != qunatity_type::complex)
	{
		cerr << "Complex operator * is used on real data. Program terminates." << endl;
		exit(0);
	}
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
	{
		double tmp = mdata[2 * i];
		mdata[2 * i] = -mdata[2 * i + 1] * multiply.mdata[2 * i + 1] / msize / msize;
		mdata[2 * i + 1] = tmp * multiply.mdata[2 * i + 1] / msize / msize;
	}
	return *this;
}

block_data<double> FT_block_data::get_block_data()
{
	if (data_type != qunatity_type::real)
	{
		cerr << "Cannot get block_data of a complex FT_block_data" << endl;
		exit(0);
	}
	block_data<double> ret(msize);
	for (int i = 0; i < msize; ++i)
	{
		memcpy(ret.getData() + i*msize, mdata + i*(msize + 2), msize*sizeof(double));
	}
	return ret;
}

block_data<double> FT_block_data::getReal_block_data()
{
	if (data_type != qunatity_type::complex)
	{
		cerr << "Cannot use getReal_block_data of a non-complex FT_block_data" << endl;
		exit(0);
	}
	block_data<double> ret(msize, 0);
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
	{
		ret[i + i / (msize / 2 + 1)*(msize / 2 - 1)] = reinterpret_cast<complex<double>*>(mdata)[i].real();
	}
	return ret;
}

block_data<double> FT_block_data::getImag_block_data()
{
	if (data_type != qunatity_type::complex)
	{
		cerr << "Cannot use getReal_block_data of a non-complex FT_block_data" << endl;
		exit(0);
	}
	block_data<double> ret(msize, 0);
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
	{
		ret[i + i / (msize / 2 + 1)*(msize / 2 - 1)] = reinterpret_cast<complex<double>*>(mdata)[i].imag();
	}
	return ret;
}

FT_block_data& FT_block_data::apply_complex2real_function(double(*function)(const complex<double>& x))
{
	if (data_type != qunatity_type::complex)
	{
		cerr << "Complex operator is used on real data. Program terminates." << endl;
		exit(0);
	}
	// the matrix is in complex format
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
	{
		reinterpret_cast<complex<double>*>(mdata)[i] = function(reinterpret_cast<complex<double>*>(mdata)[i]);
	}
	return *this;
}

FT_block_data& FT_block_data::apply_complex2complex_function(complex<double>(*function)(const complex<double>& x))
{
	if (data_type != qunatity_type::complex)
	{
		cerr << "Complex operator is used on real data. Program terminates." << endl;
		exit(0);
	}
	for (int i = 0; i < msize*(msize / 2 + 1); ++i)
	{
		reinterpret_cast<complex<double>*>(mdata)[i] = function(reinterpret_cast<complex<double>*>(mdata)[i]);
	}
	return *this;

}