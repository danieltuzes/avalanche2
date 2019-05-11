#pragma once

#include "stdafx.h"
#include "gen_lin_data.h"
#include "gen_fmanip.h"
#include "../for_gnuplot/for_gnuplot_utils.h"
#include "../fftw/fftw3.h"

using namespace std;

void abs_val2(fftw_complex * c, int size) //absolute value square for complex numbers
{
	for (int j = 0; j<size; ++j)
	{
		c[j][0] = c[j][0] * c[j][0] + c[j][1] * c[j][1];
		c[j][1] = 0;
	}
}

vector<double> autoCorrelation1D(const vector<double>& data)
{
	static int size = data.size();
	vector<double> ret(size, 0);
	fftw_complex * tmp;
	tmp = new fftw_complex[size];

	fftw_plan a = fftw_plan_dft_r2c_1d(size, const_cast<double*>(&data.front()), tmp, FFTW_ESTIMATE);
	fftw_execute(a);
	abs_val2(tmp, size);

	fftw_plan b = fftw_plan_dft_c2r_1d(size, tmp, &ret[0], FFTW_ESTIMATE);
	fftw_execute(b);

	fftw_destroy_plan(b);
	fftw_destroy_plan(a);
	for (int i = 0; i < size; ++i)
		ret[i] /= size*size;

	delete[] tmp;
	tmp = nullptr;
	return ret;
}

block_data<double> autoCorrelation2D(const block_data<double>& data)
{
	const int size = data.getSize();
	const int power = data.getPower();

	block_data<double> ret(size);
	fftw_complex * fftw_data; //data for fftw in the good byte order, in the good type
	fftw_data = new fftw_complex[size*size];
	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
		{
			fftw_data[(x << power) + y][0] = data[(x << power) + y];
			fftw_data[(x << power) + y][1] = 0;
		}

	fftw_plan FT = fftw_plan_dft_2d(size, size, fftw_data, fftw_data, FFTW_FORWARD, FFTW_ESTIMATE); //fourier transform fftw_data
	fftw_execute(FT); //execute operation

	abs_val2(fftw_data, size*size); //calculate absolute value square

	fftw_plan BFT = fftw_plan_dft_2d(size, size, fftw_data, fftw_data, FFTW_BACKWARD, FFTW_ESTIMATE); //fourier transform back
	fftw_execute(BFT); //execute back-transformation

	for (int x = 0; x < size; ++x) //copy the result to double * data
		for (int y = 0; y < size; ++y)
			ret[(x << power) + y] = fftw_data[(x << power) + y][0] / (size*size); //to ensure normality

	delete[] fftw_data;
	fftw_data = nullptr;
	return ret;
}

bool fourier_decomp_incr(const block_data<>& in, block_data<double>& abs, block_data<double>& arg) // decomp the fourier transformation of in and add to abs and arg
{
	fftw_complex * fftw_data; //data for fftw in the good byte order, in the good type
	if (abs.getSize() != arg.getSize())
	{
		cerr << "the size of abs and arg must much" << endl;
		return false;
	}
	static int size = abs.getSize();
	static int power = abs.getPower();

	fftw_data = new fftw_complex[size*size];
	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
		{
			fftw_data[(x << power) + y][0] = in(x, y);
			fftw_data[(x << power) + y][1] = 0;
		}

	fftw_plan FT = fftw_plan_dft_2d(size, size, fftw_data, fftw_data, FFTW_FORWARD, FFTW_ESTIMATE); //fourier transform fftw_data
	fftw_execute(FT); //execute operation

	for (int i = 0; i < size*size; ++i)
	{
		abs[i] += sqrt(fftw_data[i][0] * fftw_data[i][0] + fftw_data[i][1] * fftw_data[i][1]);
		arg[i] += atan2(fftw_data[i][1], fftw_data[i][0]);
	}

	delete[] fftw_data;
	fftw_data = nullptr;
	return true;
}

bool fourier_decomp_incr(const vector<double>& signal, vector<double>& abs, vector<double>& arg) // decomp the fourier transformation of in and add to abs and arg
{
	static int size = signal.size();

	fftw_complex * fftw_data;
	fftw_data = new fftw_complex[size];
	for (int i = 0; i < size; ++i)
	{
		fftw_data[i][0] = signal[i];
		fftw_data[i][1] = 0;
	}

	fftw_plan a = fftw_plan_dft_1d(size, fftw_data, fftw_data, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(a);

	for (int i = 0; i < size; ++i)
	{
		abs[i] += sqrt(fftw_data[i][0] * fftw_data[i][0] + fftw_data[i][1] * fftw_data[i][1]);
		arg[i] += atan2(fftw_data[i][1], fftw_data[i][0]);
	}

	fftw_destroy_plan(a);
	delete[] fftw_data;
	fftw_data = nullptr;
	return true;
}

block_data<double> fourierAbsValSq(block_data<double>& in)
{
	const int size = in.getSize();
	const int power = in.getPower();

	fftw_complex * fftw_data; //data for fftw in the good byte order, in the good type
	fftw_data = new fftw_complex[size*size];
	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
		{
			fftw_data[(x << power) + y][0] = in[(x << power) + y];
			fftw_data[(x << power) + y][1] = 0;
		}

	fftw_plan FT = fftw_plan_dft_2d(size, size, fftw_data, fftw_data, FFTW_FORWARD, FFTW_ESTIMATE); //fourier transform fftw_data
	fftw_execute(FT); //execute operation

	abs_val2(fftw_data, size*size); //calculate absolute value square

	block_data<double> ret(size);
	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
			ret[(x << power) + y] = fftw_data[(x << power) + y][0];

	return ret;
}


class ft_block_data : public block_data<double> {
protected:

public:
	ft_block_data() { };
	ft_block_data(string fname, dataformat f) : block_data(fname, f, 0) { } //constructor inheritance is not allowed in VS2012
	ft_block_data(char * fname, dataformat f) : block_data(fname, f, 0) { }  //constructor inheritance is not allowed in VS2012
	ft_block_data(int s) : block_data(s) { }  //constructor inheritance is not allowed in VS2012
	ft_block_data(int s, double val) : block_data(s, val) { }  //constructor inheritance is not allowed in VS2012
	ft_block_data(block_data<double> cpy) : block_data(cpy) { } //copy block_data
	~ft_block_data() { };

	bool replaceAutoCorrelation()
	{
		fftw_complex * fftw_data; //data for fftw in the good byte order, in the good type
		fftw_data = new fftw_complex[size*size];
		for (int x = 0; x < size; ++x)
			for (int y = 0; y < size; ++y)
			{
				fftw_data[(x << power) + y][0] = data[(x << power) + y];
				fftw_data[(x << power) + y][1] = 0;
			}

		fftw_plan FT = fftw_plan_dft_2d(size, size, fftw_data, fftw_data, FFTW_FORWARD, FFTW_ESTIMATE); //fourier transform fftw_data
		fftw_execute(FT); //execute operation

		abs_val2(fftw_data, size*size); //calculate absolute value square

		fftw_plan BFT = fftw_plan_dft_2d(size, size, fftw_data, fftw_data, FFTW_BACKWARD, FFTW_ESTIMATE); //fourier transform back
		fftw_execute(BFT); //execute back-transformation

		for (int x = 0; x < size; ++x) //copy the result to double * data
			for (int y = 0; y < size; ++y)
				data[(x << power) + y] = fftw_data[(x << power) + y][0] / (size*size); //to ensure normality

		delete[] fftw_data;
		fftw_data = nullptr;
		return true;
	}

	bool replaceFourierAbsValSq()
	{
		fftw_complex * fftw_data; //data for fftw in the good byte order, in the good type
		fftw_data = new fftw_complex[size*size];
		for (int x = 0; x < size; ++x)
			for (int y = 0; y < size; ++y)
			{
				fftw_data[(x << power) + y][0] = data[(x << power) + y];
				fftw_data[(x << power) + y][1] = 0;
			}

		fftw_plan FT = fftw_plan_dft_2d(size, size, fftw_data, fftw_data, FFTW_FORWARD, FFTW_ESTIMATE); //fourier transform fftw_data
		fftw_execute(FT); //execute operation

		abs_val2(fftw_data, size*size); //calculate absolute value square

		for (int x = 0; x < size; ++x) //copy the result to double * data
			for (int y = 0; y < size; ++y)
				data[(x << power) + y] = fftw_data[(x << power) + y][0] / (size*size); //to ensure normality

		delete[] fftw_data;
		fftw_data = nullptr;
		return true;
	}

	ft_block_data bruteAutoCorr() const
	{
		ft_block_data ret(size);
		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size; ++j)
			{
				ret(i, j) = 0;
				for (int k = 0; k < size; ++k)
					for (int l = 0; l < size; ++l)
						ret(i, j) += periodic(k + i, l + j) * operator()(k, l);

			}

		return ret;
	} // this is a debugging tool only, use replaceAutoCorrelation() instead for faster calculation

	bool replaceXCorrelation(const ft_block_data& other) //replace data with the cross-correlation of the argument and itself; replace the argument with its un-normalized fourier transformed val
	{
		if (size != other.getSize())
		{
			cerr << "Program cannot calculate cross-correlation of different sized datas. Program temrinates." << endl;
			exit(1);
		}

		fftw_complex * fftw_data_first, *fftw_data_other; //data for fftw in the good byte order, in the good type
		fftw_data_first = new fftw_complex[size*size];
		fftw_data_other = new fftw_complex[size*size];
		for (int i = 0; i<size; ++i)
			for (int j = 0; j<size; ++j)
			{
				fftw_data_first[(i << power) + j][0] = data[(i << power) + j];
				fftw_data_first[(i << power) + j][1] = 0;
				fftw_data_other[(i << power) + j][0] = other[(i << power) + j];
				fftw_data_other[(i << power) + j][1] = 0;
			}

		fftw_plan FT_f = fftw_plan_dft_2d(size, size, fftw_data_first, fftw_data_first, FFTW_FORWARD, FFTW_ESTIMATE); //fourier transform fftw_data_first
		fftw_plan FT_o = fftw_plan_dft_2d(size, size, fftw_data_other, fftw_data_other, FFTW_FORWARD, FFTW_ESTIMATE); //fourier transform fftw_data_other
		fftw_execute(FT_f); //execute operation
		fftw_execute(FT_o); //execute operation


		for (int i = 0; i<size; ++i)
			for (int j = 0; j<size; ++j)
			{
				double tmp = fftw_data_first[(i << power) + j][0];
				fftw_data_first[(i << power) + j][0] = fftw_data_first[(i << power) + j][0] * fftw_data_other[(i << power) + j][0] + fftw_data_first[(i << power) + j][1] * fftw_data_other[(i << power) + j][1];
				fftw_data_first[(i << power) + j][1] = fftw_data_first[(i << power) + j][1] * fftw_data_other[(i << power) + j][0] - tmp * fftw_data_other[(i << power) + j][1];
			}

#pragma warning( push )
#pragma warning( disable : 4127 )
		if (DEBUG_TEST)
#pragma warning( pop )
		{
			block_data<double> rp(size, 0.);
			block_data<double> ip(size, 0.);
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < size; ++j)
				{
					if (j < size)
					{
						rp(i, j) = fftw_data_first[(i << power) + j][0];
						ip(i, j) = fftw_data_first[(i << power) + j][1];
					}
					else
					{
						rp(i, j) = 0;
						ip(i, j) = 0;
					}
				}
			}
			rp.fWriteToFile("Xreal_part.dat", "# real part, orig\n");
			ip.fWriteToFile("Ximag_part.dat", "# imag part\n");
		}

		fftw_plan BFT = fftw_plan_dft_2d(size, size, fftw_data_first, fftw_data_first, FFTW_BACKWARD, FFTW_ESTIMATE); //fourier transform back
		fftw_execute(BFT); //execute back-transformation

		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size; ++j)
				data[(i << power) + j] = fftw_data_first[(i << power) + j][0];


		*this /= size*size;

		return true;
	}

	bool replace_rrXCorrelation(const ft_block_data& other) //replace data with the cross-correlation of the argument and itself; replace the argument with its un-normalized fourier transformed val
	{
		if (size != other.getSize())
		{
			cerr << "Program cannot calculate cross-correlation of different sized datas. Program temrinates." << endl;
			exit(1);
		}

		fftw_complex *fftw_data_first_out, *fftw_data_other_out; //data for fftw in the good byte order, in the good type
		fftw_data_first_out = new fftw_complex[size*size];
		fftw_data_other_out = new fftw_complex[size*size];
		memset(*fftw_data_first_out, 0, (size*size)*sizeof(fftw_complex));
		memset(*fftw_data_other_out, 0, (size*size)*sizeof(fftw_complex));



		fftw_plan FT_f = fftw_plan_dft_r2c_2d(size, size, this->data, fftw_data_first_out, FFTW_ESTIMATE); //fourier transform fftw_data_first
		fftw_plan FT_o = fftw_plan_dft_r2c_2d(size, size, other.data, fftw_data_other_out, FFTW_ESTIMATE); //fourier transform fftw_data_other
		fftw_execute(FT_f); //execute operation
		fftw_execute(FT_o); //execute operation



		for (int i = 0; i < size; ++i)
			for (int j = 0; j < size / 2 + 1; ++j)
			{
				double tmp = fftw_data_first_out[(i * size / 2 + 1) + j][0];
				fftw_data_first_out[(i * size / 2 + 1) + j][0] =
					fftw_data_first_out[(i * size / 2 + 1) + j][0] * fftw_data_other_out[(i * size / 2 + 1) + j][0] + fftw_data_first_out[(i * size / 2 + 1) + j][1] * fftw_data_first_out[(i * size / 2 + 1) + j][1];
				fftw_data_first_out[(i * size / 2 + 1) + j][1] =
					fftw_data_first_out[(i * size / 2 + 1) + j][1] * fftw_data_other_out[(i * size / 2 + 1) + j][0] - tmp * fftw_data_other_out[(i * size / 2 + 1) + j][1];
			}


#pragma warning( push )
#pragma warning( disable : 4127 )
		if (DEBUG_TEST)
#pragma warning( pop )
		{
			block_data<double> rp(size, 0.);
			block_data<double> ip(size, 0.);
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < size; ++j)
				{
					if (j > size / 2)
					{
						rp(i, j) = 0;
					}
					else
					{
						rp(i, j) = fftw_data_first_out[i * (size / 2 + 1) + j][0];
						ip(i, j) = fftw_data_first_out[i * (size / 2 + 1) + j][1];
					}
				}
			}
			rp.fWriteToFile("rrX_real_part.dat", "# real part\n");
			ip.fWriteToFile("rrX_imag_part.dat", "# imag part\n");
		}

		fftw_plan BFT = fftw_plan_dft_c2r_2d(size, size / 2, fftw_data_first_out, this->data, FFTW_ESTIMATE); //fourier transform back
		fftw_execute(BFT); //execute back-transformation

		return true;
	}
};

#pragma region "1D convolution functions"

vector<double> convolution_1d_rr_bruteforce(const vector<double>& signal, const vector<double>& kernel)
{
	unsigned int s = signal.size();
	vector<double> ret(s);
	for (unsigned int i = 0; i < s; ++i)
	{
		for (unsigned int ii = 0; ii < s; ++ii)
		{
			ret[i] += signal[ii] * kernel[(i - ii + s) & (s - 1)];
			
		}
	}
	return ret;
}

vector<pair<double, double>> convolution_1d_rr_fftw_via_c(const vector<double>& signal, const vector<double>& kernel)
{
	unsigned int s = signal.size();
	vector<pair<double, double>> ret(s, pair<double, double>(0, 0));

	fftw_complex *in_signal, *in_kernel, *out_signal, *out_kernel, *ret_in, *ret_out;
	fftw_plan p_signal, p_kernel, p_ret;
	in_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s);
	out_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s);
	in_kernel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s);
	out_kernel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s);
	ret_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s);
	ret_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s);
	p_signal = fftw_plan_dft_1d(s, in_signal, out_signal, FFTW_FORWARD, FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_1d(s, in_kernel, out_kernel, FFTW_FORWARD, FFTW_ESTIMATE);

	for (unsigned int i = 0; i < s; ++i)
	{
		in_signal[i][0] = signal[i];
		in_signal[i][1] = 0;

		in_kernel[i][0] = kernel[i];
		in_kernel[i][1] = 0;
	}


	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (unsigned int i = 0; i < s; ++i)
	{
		ret_in[i][0] = out_signal[i][0] * out_kernel[i][0] - out_signal[i][1] * out_kernel[i][1];
		ret_in[i][1] = out_signal[i][0] * out_kernel[i][1] + out_signal[i][1] * out_kernel[i][0];
	}

	p_ret = fftw_plan_dft_1d(s, ret_in, ret_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p_ret);

	for (unsigned int i = 0; i < s; ++i)
	{
		ret[i].first = ret_out[i][0] / s;
		ret[i].second = ret_out[i][1] / s;
	}

	return ret;
}

vector<pair<double, double>> convolution_1d_rr_fftw_via_c_in_place(const vector<double>& signal, const vector<double>& kernel)
{
	unsigned int s = signal.size();
	vector<pair<double, double>> ret(s, pair<double, double>(0, 0));

	fftw_complex *f_signal, *f_kernel, *f_ret;
	fftw_plan p_signal, p_kernel, p_ret;
	f_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s);
	f_kernel= (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s);
	f_ret = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s);
	p_signal = fftw_plan_dft_1d(s, f_signal, f_signal, FFTW_FORWARD, FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_1d(s, f_kernel, f_kernel, FFTW_FORWARD, FFTW_ESTIMATE);

	for (unsigned int i = 0; i < s; ++i)
	{
		f_signal[i][0] = signal[i];
		f_signal[i][1] = 0;

		f_kernel[i][0] = kernel[i];
		f_kernel[i][1] = 0;
	}


	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (unsigned int i = 0; i < s; ++i)
	{
		//f_ret[i] = f_signal[i] * f_kernel[i];
		f_ret[i][1] = f_signal[i][0] * f_kernel[i][1] + f_signal[i][1] * f_kernel[i][0];
	}

	p_ret = fftw_plan_dft_1d(s, f_ret, f_ret, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p_ret);

	for (unsigned int i = 0; i < s; ++i)
	{
		ret[i].first = f_ret[i][0] / s;
		ret[i].second = f_ret[i][1] / s;
	}

	return ret;
}

vector<complex<double>> convolution_1d_rr_fftw_via_c_in_place_native(const vector<double>& signal, const vector<double>& kernel)
{
	unsigned int s = signal.size();
	vector<complex<double>> ret(s, complex<double>(0,0));

	complex<double> *f_signal, *f_kernel, *f_ret;
	fftw_plan p_signal, p_kernel, p_ret;
	f_signal = (complex<double>*)fftw_malloc(sizeof(complex<double>) * s);
	f_kernel = (complex<double>*)fftw_malloc(sizeof(complex<double>) * s);
	f_ret = (complex<double>*)fftw_malloc(sizeof(complex<double>) * s);
	p_signal = fftw_plan_dft_1d(s, reinterpret_cast<fftw_complex*>(f_signal), reinterpret_cast<fftw_complex*>(f_signal), FFTW_FORWARD, FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_1d(s, reinterpret_cast<fftw_complex*>(f_kernel), reinterpret_cast<fftw_complex*>(f_kernel), FFTW_FORWARD, FFTW_ESTIMATE);

	for (unsigned int i = 0; i < s; ++i)
	{
		f_signal[i] = signal[i];
		f_kernel[i] = kernel[i];
	}


	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (unsigned int i = 0; i < s; ++i)
	{
		f_ret[i] = f_signal[i] * f_kernel[i];
	}

	p_ret = fftw_plan_dft_1d(s, reinterpret_cast<fftw_complex*>(f_ret), reinterpret_cast<fftw_complex*>(f_ret), FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p_ret);

	for (unsigned int i = 0; i < s; ++i)
	{
		ret[i] = f_ret[i] / static_cast<double>(s);
	}

	return ret;
}

vector<complex<double>> convolution_1d_rr_fftw_via_c_in_place_native_pointedResult(const vector<double>& signal, const vector<double>& kernel)
{
	unsigned int s = signal.size();
	vector<complex<double>> ret(s, complex<double>(0, 0));

	complex<double> *f_signal, *f_kernel;
	fftw_plan p_signal, p_kernel, p_ret;
	f_signal = (complex<double>*)fftw_malloc(sizeof(complex<double>) * s);
	f_kernel = (complex<double>*)fftw_malloc(sizeof(complex<double>) * s);
	p_signal = fftw_plan_dft_1d(s, reinterpret_cast<fftw_complex*>(f_signal), reinterpret_cast<fftw_complex*>(f_signal), FFTW_FORWARD, FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_1d(s, reinterpret_cast<fftw_complex*>(f_kernel), reinterpret_cast<fftw_complex*>(f_kernel), FFTW_FORWARD, FFTW_ESTIMATE);

	for (unsigned int i = 0; i < s; ++i)
	{
		f_signal[i] = signal[i];
		f_kernel[i] = kernel[i];
	}


	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (unsigned int i = 0; i < s; ++i)
	{
		ret[i] = f_signal[i] * f_kernel[i];
	}

	p_ret = fftw_plan_dft_1d(s, reinterpret_cast<fftw_complex*>(ret.data()), reinterpret_cast<fftw_complex*>(ret.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p_ret);

	for (unsigned int i = 0; i < s; ++i)
	{
		ret[i] /= s;
	}

	return ret;
}

vector<double> convolution_1d_rr_fftw_via_r(const vector<double>& signal, const vector<double>& kernel)
{
	unsigned int s = signal.size();
	vector<double> ret(s);

	double *d_signal, *d_kernel, *d_ret;
	fftw_complex *f_signal, *f_kernel, *f_ret;
	fftw_plan p_signal, p_kernel, p_ret;

	d_signal = (double*)fftw_malloc(sizeof(double) * (s + 2));
	d_kernel = (double*)fftw_malloc(sizeof(double) * (s + 2));
	d_ret = (double*)fftw_malloc(sizeof(double) * (s + 2));
	f_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (s / 2 + 1));
	f_kernel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (s / 2 + 1));
	f_ret = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (s / 2 + 1));

	p_signal = fftw_plan_dft_r2c_1d(s, d_signal, f_signal, FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_r2c_1d(s, d_kernel, f_kernel, FFTW_ESTIMATE);

	for (unsigned int i = 0; i < s; ++i)
	{
		d_signal[i] = signal[i];
		d_kernel[i] = kernel[i];
	}


	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (unsigned int i = 0; i < s / 2 + 1; ++i)
	{
		f_ret[i][0] = f_signal[i][0] * f_kernel[i][0] - f_signal[i][1] * f_kernel[i][1];
		f_ret[i][1] = f_signal[i][0] * f_kernel[i][1] + f_signal[i][1] * f_kernel[i][0];
	}

	p_ret = fftw_plan_dft_c2r_1d(s, f_ret, d_ret, FFTW_ESTIMATE);
	fftw_execute(p_ret);

	for (unsigned int i = 0; i < s; ++i)
	{
		ret[i] = d_ret[i] / s;
	}

	return ret;
}

vector<double> convolution_1d_rr_fftw_via_r_in_place(const vector<double>& signal, const vector<double>& kernel)
{
	unsigned int s = signal.size();
	vector<double> ret(s);

	double *d_signal, *d_kernel, *d_ret;
	fftw_plan p_signal, p_kernel, p_ret;

	d_signal = (double*)fftw_malloc(sizeof(double) * (s + 2));
	d_kernel = (double*)fftw_malloc(sizeof(double) * (s + 2));
	d_ret = (double*)fftw_malloc(sizeof(double) * (s + 2));

	p_signal = fftw_plan_dft_r2c_1d(s, d_signal, reinterpret_cast<fftw_complex*>(d_signal), FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_r2c_1d(s, d_kernel, reinterpret_cast<fftw_complex*>(d_kernel), FFTW_ESTIMATE);

	for (unsigned int i = 0; i < s; ++i)
	{
		d_signal[i] = signal[i];
		d_kernel[i] = kernel[i];
	}


	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (unsigned int i = 0; i < s / 2 + 1; ++i)
	{
		reinterpret_cast<fftw_complex*>(d_ret)[i][0] = reinterpret_cast<fftw_complex*>(d_signal)[i][0] * reinterpret_cast<fftw_complex*>(d_kernel)[i][0] - reinterpret_cast<fftw_complex*>(d_signal)[i][1] * reinterpret_cast<fftw_complex*>(d_kernel)[i][1];
		reinterpret_cast<fftw_complex*>(d_ret)[i][1] = reinterpret_cast<fftw_complex*>(d_signal)[i][0] * reinterpret_cast<fftw_complex*>(d_kernel)[i][1] + reinterpret_cast<fftw_complex*>(d_signal)[i][1] * reinterpret_cast<fftw_complex*>(d_kernel)[i][0];
	}

	p_ret = fftw_plan_dft_c2r_1d(s, reinterpret_cast<fftw_complex*>(d_ret), d_ret, FFTW_ESTIMATE);
	fftw_execute(p_ret);

	for (unsigned int i = 0; i < s; ++i)
	{
		ret[i] = d_ret[i] / s;
	}

	return ret;
}

vector<double> convolution_1d_rr_fftw_via_r_in_place_native(const vector<double>& signal, const vector<double>& kernel)
{
	unsigned int s = signal.size();
	vector<double> ret(s);

	double *d_signal, *d_kernel, *d_ret;
	fftw_plan p_signal, p_kernel, p_ret;

	d_signal = (double*)fftw_malloc(sizeof(double) * (s + 2));
	d_kernel = (double*)fftw_malloc(sizeof(double) * (s + 2));
	d_ret = (double*)fftw_malloc(sizeof(double) * (s + 2));

	p_signal = fftw_plan_dft_r2c_1d(s, d_signal, reinterpret_cast<fftw_complex*>(d_signal), FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_r2c_1d(s, d_kernel, reinterpret_cast<fftw_complex*>(d_kernel), FFTW_ESTIMATE);

	for (unsigned int i = 0; i < s; ++i)
	{
		d_signal[i] = signal[i];
		d_kernel[i] = kernel[i];
	}


	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (unsigned int i = 0; i < s / 2 + 1; ++i)
	{
		reinterpret_cast<complex<double>*>(d_ret)[i] = reinterpret_cast<complex<double>*>(d_signal)[i] * reinterpret_cast<complex<double>*>(d_kernel)[i];
	}

	p_ret = fftw_plan_dft_c2r_1d(s, reinterpret_cast<fftw_complex*>(d_ret), d_ret, FFTW_ESTIMATE);
	fftw_execute(p_ret);

	for (unsigned int i = 0; i < s; ++i)
	{
		ret[i] = d_ret[i] / s;
	}

	return ret;
}

vector<double> convolution_1d_rr_fftw_via_r_in_place_native_pointed(vector<double>& signal, vector<double>& kernel)
{
	unsigned int s = signal.size();
	vector<double> ret(s + 2);
	signal.resize(signal.size() + 2);
	kernel.resize(signal.size() + 2);

	fftw_plan p_signal, p_kernel, p_ret;

	p_signal = fftw_plan_dft_r2c_1d(s, signal.data(), reinterpret_cast<fftw_complex*>(signal.data()), FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_r2c_1d(s, kernel.data(), reinterpret_cast<fftw_complex*>(kernel.data()), FFTW_ESTIMATE);
	
	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (unsigned int i = 0; i < s / 2 + 1; ++i)
	{
		reinterpret_cast<complex<double>*>(ret.data())[i] = reinterpret_cast<complex<double>*>(signal.data())[i] * reinterpret_cast<complex<double>*>(kernel.data())[i];
	}

	p_ret = fftw_plan_dft_c2r_1d(s, reinterpret_cast<fftw_complex*>(ret.data()), ret.data(), FFTW_ESTIMATE);
	fftw_execute(p_ret);

	for (unsigned int i = 0; i < s; ++i)
	{
		ret[i] /= s;
	}

	ret.resize(s);
	return ret;
}

#pragma endregion

#pragma region "2D convolution functions"

vector<vector<double>> convolution_2d_rr_bruteforce(const vector<vector<double>>& signal, const vector<vector<double>>& kernel)
{
	int s = static_cast<int>(signal.size());
	vector<vector<double>> ret(s, vector<double>(s));

	for (int i = 0; i < s; ++i)
	{
		for (int j = 0; j < s; ++j)
		{
			for (int k = 0; k < s; ++k)
			{
				for (int l = 0; l < s; ++l)
				{
					ret[i][j] += signal[k][l] * kernel[(i - k + s) & (s - 1)][(j - l + s) & (s - 1)];
				}
			}
		}
	}
	return ret;
}

pair<vector<vector<double>>, vector<vector<double>>> convolution_2d_rr_via_c(const vector<vector<double>>& signal, const vector<vector<double>>& kernel)
{
	int s = static_cast<int>(signal.size());
	pair<vector<vector<double>>, vector<vector<double>>> ret(vector<vector<double>>(s, vector<double>(s)), vector<vector<double>>(s, vector<double>(s)));

	fftw_complex *in_signal, *in_kernel, *out_signal, *out_kernel, *ret_in, *ret_out;
	fftw_plan p_signal, p_kernel, p_ret;
	in_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	out_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	in_kernel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	out_kernel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	ret_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	ret_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	p_signal = fftw_plan_dft_2d(s, s, in_signal, out_signal, FFTW_FORWARD, FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_2d(s, s, in_kernel, out_kernel, FFTW_FORWARD, FFTW_ESTIMATE);

	for (int i = 0; i < s; ++i)
	{
		for (int j = 0; j < s; ++j)
		{
			int index = i*s + j;
			in_signal[index][0] = signal[i][j];
			in_signal[index][1] = 0;
			in_kernel[index][0] = kernel[i][j];
			in_kernel[index][1] = 0;
		}
	}

	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (int i = 0; i < s; ++i)
	{
		for (int j = 0; j < s; ++j)
		{
			int index = i*s + j;
			ret_in[index][0] = out_signal[index][0] * out_kernel[index][0] - out_signal[index][1] * out_kernel[index][1];
			ret_in[index][1] = out_signal[index][0] * out_kernel[index][1] + out_signal[index][1] * out_kernel[index][0];
		}
	}


	p_ret = fftw_plan_dft_2d(s, s, ret_in, ret_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p_ret);


	for (int i = 0; i < s; ++i)
	{
		for (int j = 0; j < s; ++j)
		{
			int index = i*s + j;
			ret.first[i][j] = ret_out[index][0] / s / s;
			ret.second[i][j] = ret_out[index][1] / s / s;
		}
	}

	return ret;
}

vector<vector<double>> convolution_2d_rr_via_c_native(const vector<vector<double>>& signal, const vector<vector<double>>& kernel)
{
	int s = static_cast<int>(signal.size());
	vector<vector<double>> ret(vector<vector<double>>(s, vector<double>(s)));

	fftw_complex *in_signal, *in_kernel, *out_signal, *out_kernel, *ret_in, *ret_out;
	fftw_plan p_signal, p_kernel, p_ret;
	in_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	out_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	in_kernel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	out_kernel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	ret_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	ret_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * s);
	p_signal = fftw_plan_dft_2d(s, s, in_signal, out_signal, FFTW_FORWARD, FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_2d(s, s, in_kernel, out_kernel, FFTW_FORWARD, FFTW_ESTIMATE);

	for (int i = 0; i < s; ++i)
	{
		for (int j = 0; j < s; ++j)
		{
			int index = i*s + j;
			reinterpret_cast<complex<double>*>(in_signal)[index] = signal[i][j];
			reinterpret_cast<complex<double>*>(in_kernel)[index] = kernel[i][j];
		}
	}

	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (int i = 0; i < s; ++i)
	{
		for (int j = 0; j < s; ++j)
		{
			int index = i*s + j;
			reinterpret_cast<complex<double>*>(ret_in)[index] = reinterpret_cast<complex<double>*>(out_signal)[index] * reinterpret_cast<complex<double>*>(out_kernel)[index];
		}
	}


	p_ret = fftw_plan_dft_2d(s, s, ret_in, ret_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p_ret);


	for (int i = 0; i < s; ++i)
	{
		for (int j = 0; j < s; ++j)
		{
			int index = i*s + j;
			ret[i][j] = ret_out[index][0] / s / s;
		}
	}

	return ret;
}

vector<vector<double>> convolution_2d_rr_via_r_native(const vector<vector<double>>& signal, const vector<vector<double>>& kernel)
{
	int s = static_cast<int>(signal.size());
	vector<vector<double>> ret(vector<vector<double>>(s, vector<double>(s)));

	double *in_signal, *in_kernel, *ret_out;
	fftw_complex *out_signal, *out_kernel, *ret_in;
	fftw_plan p_signal, p_kernel, p_ret;
	in_signal = (double*)fftw_malloc(sizeof(double) * s * s);
	out_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * (s / 2 + 1));
	in_kernel = (double*)fftw_malloc(sizeof(double) * s * s);
	out_kernel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * (s / 2 + 1));
	ret_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * (s / 2 + 1));
	ret_out = (double*)fftw_malloc(sizeof(double) * s * s);
	p_signal = fftw_plan_dft_r2c_2d(s, s, in_signal, out_signal, FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_r2c_2d(s, s, in_kernel, out_kernel, FFTW_ESTIMATE);

	for (int i = 0; i < s; ++i)
	{
		for (int j = 0; j < s; ++j)
		{
			int index = i*s + j;
			in_signal[index] = signal[i][j];
			in_kernel[index] = kernel[i][j];
		}
	}

	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (int i = 0; i < s*(s / 2 + 1); ++i)
	{
		reinterpret_cast<complex<double>*>(ret_in)[i] = reinterpret_cast<complex<double>*>(out_signal)[i] * reinterpret_cast<complex<double>*>(out_kernel)[i];
	}


	p_ret = fftw_plan_dft_c2r_2d(s, s, ret_in, ret_out, FFTW_ESTIMATE);
	fftw_execute(p_ret);


	for (int i = 0; i < s; ++i)
	{
		for (int j = 0; j < s; ++j)
		{
			int index = i*s + j;
			ret[i][j] = ret_out[index] / s / s;
		}
	}

	return ret;
}

block_data<double> convolution_2d_rr_via_r_native(const block_data<double>& signal, const block_data<double>& kernel)
{
	int s = static_cast<int>(signal.getSize());
	block_data<double> ret(s);

	fftw_complex *out_signal, *out_kernel, *ret_in;
	fftw_plan p_signal, p_kernel, p_ret;
	out_signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * (s / 2 + 1));
	out_kernel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * (s / 2 + 1));
	ret_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * s * (s / 2 + 1));
	p_signal = fftw_plan_dft_r2c_2d(s, s, signal.getData(), out_signal, FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_r2c_2d(s, s, kernel.getData(), out_kernel, FFTW_ESTIMATE);

	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	block_data<> out_signal_real(s, 0), out_signal_imag(s, 0), out_kernel_real(s, 0), out_kernel_imag(s, 0);

	for (int i = 0; i < s*(s / 2 + 1); ++i)
	{
		out_signal_real[i] = out_signal[i][0];
		out_signal_imag[i] = out_signal[i][1];
		out_kernel_real[i] = out_kernel[i][0];
		out_kernel_imag[i] = out_kernel[i][1];
	}
	out_signal_real.fWriteToFile("out_signal_real.txt");
	out_signal_imag.fWriteToFile("out_signal_imag.txt");
	out_kernel_real.fWriteToFile("out_kernel_real.txt");
	out_kernel_imag.fWriteToFile("out_kernel_imag.txt");

	for (int i = 0; i < s*(s / 2 + 1); ++i)
	{
		reinterpret_cast<complex<double>*>(ret_in)[i] = reinterpret_cast<complex<double>*>(out_signal)[i] * reinterpret_cast<complex<double>*>(out_kernel)[i];
	}


	p_ret = fftw_plan_dft_c2r_2d(s, s, ret_in, ret.getData(), FFTW_ESTIMATE);
	fftw_execute(p_ret);

	return ret / s / s;
}

block_data<double> convolution_2d_rr_via_r_native_in_place(const block_data<double>& in_signal, const block_data<double>& in_kernel)
{
	int s = static_cast<int>(in_signal.getSize());

	double *signal, *kernel;
	signal = fftw_alloc_real(s*(s + 2));
	kernel = fftw_alloc_real(s*(s + 2));

	for (int i = 0; i < s*s; ++i)
	{
		signal[i + 2 * (i / s)] = in_signal[i];
		kernel[i + 2 * (i / s)] = in_kernel[i];
	}

	fftw_plan p_signal, p_kernel, p_ret;

	p_signal = fftw_plan_dft_r2c_2d(s, s, signal, reinterpret_cast<fftw_complex*>(signal), FFTW_ESTIMATE);
	p_kernel = fftw_plan_dft_r2c_2d(s, s, kernel, reinterpret_cast<fftw_complex*>(kernel), FFTW_ESTIMATE);

	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	//block_data<> out_signal_real(s, 0), out_signal_imag(s, 0), out_kernel_real(s, 0), out_kernel_imag(s, 0);

	//for (int i = 0; i < s*(s / 2 + 1); ++i)
	//{
	//	out_signal_real[i] = signal[2 * i];
	//	out_signal_imag[i] = signal[2 * i + 1];
	//	out_kernel_real[i] = kernel[2 * i];
	//	out_kernel_imag[i] = kernel[2 * i + 1];
	//}
	//out_signal_real.fWriteToFile("out_signal_real_ip.txt");
	//out_signal_imag.fWriteToFile("out_signal_imag_ip.txt");
	//out_kernel_real.fWriteToFile("out_kernel_real_ip.txt");
	//out_kernel_imag.fWriteToFile("out_kernel_imag_ip.txt");

	
	for (int i = 0; i < s*(s / 2 + 1); ++i)
		reinterpret_cast<complex<double>*>(signal)[i] *= reinterpret_cast<complex<double>*>(kernel)[i] /(s*s+0.);

	block_data<> ret(s);

	p_ret = fftw_plan_dft_c2r_2d(s, s, reinterpret_cast<fftw_complex*>(signal), ret.getData(), FFTW_ESTIMATE);
	fftw_execute(p_ret);

	return ret;
}

vector<double> convolution_2d_rr_via_r_native_in_place_vector(block_data<double>& signal_bd, block_data<double>& kernel_bd) // doesn't work because one has to use fftw_malloc
{
	int s = static_cast<int>(signal_bd.getSize());
	vector<double> ret(s * (s + 2)), signal(signal_bd.getData(), signal_bd.getData() + s*(s + 2)), kernel(kernel_bd.getData(), kernel_bd.getData() + s*(s + 2));

	fftw_plan p_signal, p_kernel, p_ret;
	int n[] = { s, s };
	p_signal = fftw_plan_many_dft_r2c(2, n, 1, signal.data(), NULL, 1, 0, reinterpret_cast<fftw_complex*>(signal.data()), NULL, 1, 0, FFTW_ESTIMATE);
	p_kernel = fftw_plan_many_dft_r2c(2, n, 1, kernel.data(), NULL, 1, 0, reinterpret_cast<fftw_complex*>(kernel.data()), NULL, 1, 0, FFTW_ESTIMATE);

	fftw_execute(p_signal);
	fftw_execute(p_kernel);

	for (int i = 0; i < s*(s / 2 + 1); ++i)
	{
		reinterpret_cast<complex<double>*>(ret.data())[i] = reinterpret_cast<complex<double>*>(signal.data())[i] * reinterpret_cast<complex<double>*>(kernel.data())[i];
	}

	p_ret = fftw_plan_dft_c2r_2d(s, s, reinterpret_cast<fftw_complex*>(ret.data()), ret.data(), FFTW_ESTIMATE);
	fftw_execute(p_ret);

	ret.resize(s*s);
	for (int i = 0; i < s*s; ++i)
	{
		ret[i] /= s*s;
	}
	return ret;
}

#pragma endregion
