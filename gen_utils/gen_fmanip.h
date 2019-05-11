// gen_fmanip.h : stores user defined general function relatad to file and data_block manipulation
//
#define GEN_FMANIP_VERSION 0.09

// What's new?

// 0.09: fWriteFile's specialNan applied to a copy of the original data
// 0.08: stats implemented, fWriteFile have autoheader and treatNanAs
// 0.07: projectx and projecty are implemented; block_data(string fname, dataformat f, int skip = 0) is corrected, constructors are nested
// 0.06: block_data.data is more safely deleted, block_data& operator=(const block_data& tocopy) is more safely working, block_data& operator=(block_data&& tocopy) is added
// GEN_FMANIP_VERSION 0.05: getMinWhere and getMaxWhere didn't give the correct values, now they do

#pragma once
#include "stdafx.h"

using namespace std;

enum dataformat
{
	binary,            // the data in binary format, without index values
	text,              // the data in ascii text format, with index values
	bin_text,          // the data format is unknown, but binary is should test first
	text_bin,          // the data format is unknown, but text is should test first
	origin,            // the data format is known, from origin
	bin_text_origin,   // the data format is unknown
	zip                // the data is compressed with zip algorithm
};

ifstream::pos_type filesize(string fname);
ifstream::pos_type filesize(const char * fname);
ifstream::pos_type filesize(string fname); //tell the size of a file, exit if doesn't exist
ifstream::pos_type filesize(const char * fname);


bool is_file_exist(string fname);
bool is_file_exist(char * fname);

void fstream_must_good(bool stream_is_good, string fname);

string lastLine(ifstream& file, int lineSize = 80);  //finds the last line of the stream; lineSize is just an estimation on line size
string lastNonEmptyLine(string fname, int lineSize = 80);


template <typename T = double> class block_data {
protected:
	int size;
	int power;
	T * data;
	void periodicIncr(int& x) const;
public:
	block_data();
	block_data(const block_data& tocopy);
	block_data(int s);
	block_data(int s, T val);
	block_data(vector<vector<T>> m_vector);
	block_data(vector<T> m_vector);
	
	block_data(string fname, dataformat f, int skip = 0);
	block_data(const char * fname, dataformat f, int skip = 0);
	~block_data();

	bool hasData();
	block_data& operator=(const block_data& tocopy);
	block_data& operator=(block_data&& tocopy);
	block_data& operator=(T initVal);

	bool init(int size, T val);
	bool init(const block_data<T>& newData);

	bool readFromFile(char * fname);
	bool readFromFile(string fname);
	bool originReadFromFile(string fname, int skip);

	bool fReadFromFile(char * fname);
	bool fReadFromFile(string fname);

	block_data<T> reMap() const;
	block_data<T> reSample(int blockSize = 2, int shiftx = 0, int shifty = 0) const;
	bool symmetrize(block_data<T>& dest) const;

	int getSize() const;
	int getPower() const;

	T operator[](int i) const;
	T& operator[](int i);
	T* getData();
	T* getData() const;
	T operator()(int x, int y) const;
	T& operator()(int x, int y);
	T periodic(int x, int y) const;
	T& periodic(int x, int y);
	T periodic_x(int x, int y) const;
	T& periodic_x(int x, int y);
	T periodic_y(int x, int y) const;
	T& periodic_y(int x, int y);

	T rightAv(int x, int y) const; // gives the average of the cell's and it's right neighbour's values
	T localAv(int x, int y, int radiX = 2, int radY = 1) const; // gives the local average around the border (x,y) and (x+1,y)
	T rightGrad(int x, int y, T dist = 1) const; // gives the average of the cell's and it's right neighbour's values
	
	vector<vector<T>> std2Dvector() const;

	vector<T> slicey(int x0) const;
	vector<T> slicex(int y0) const;
	vector<T> projecty() const;
	vector<T> projectx() const;
	vector<double> rdependence(double x0, double y0) const;
	vector<double> rabsdependence(double x0, double y0) const;
	vector<double> rcos4fidependence(double x0, double y0) const;

	bool fWriteToFile(string fname, string header = "", int precision = 9, bool swapOrder = false) const; //(over)write formatted data to fname; create file if not exist
	bool fWriteToFileTreatNN(string fname, bool speacialNan = true, double treatNanAs = 1, bool autoheader = true, string header = "", int precision = 9, bool swapOrder = false) const; //(over)write formatted data to fname; create file if not exist
	// bool fWriteToFile(char * fname, char * header = "", int precision = 9) const;  //(over)write formatted data to fname; create file if not exist

	bool writeToFile(string fname) const; //(over)write binary data to fname; create file if not exist
	// bool writeToFile(char * fname) const;  //(over)write data to fname; create file if not exist


	void sumNormalize(T n = 1); //devide every data value to make sum of them to the given n
	void maxNormalize(T n = 1); //devide every data value to make max value to the given n
	block_data Abs() const;

	T getMax() const;
	T getMin() const;
	T getSum() const;
	double getAverage() const;
	T getSumSquare() const;
	long double getStdDev() const;
	long double getStdDev(long double average) const;
	long double getCoeffofVar() const;
	long double getCoeffofVar(long double average) const;
	long double getAssymetry() const;

	tuple<T, int, int> getMaxWhere() const;
	tuple<T, int, int> getMinWhere() const;

	block_data operator+(const block_data& toadd) const;
	block_data& operator+=(const block_data& toadd);
	block_data operator-(const block_data& tusub) const;
	block_data& operator-=(const block_data& tusub);
	block_data operator+(T toadd) const;
	block_data& operator+=(T toadd);
	block_data operator-(T tusub) const;
	block_data& operator-=(T tusub);

	block_data operator*(const block_data& mutliply) const; // pointwise mulitplication
	block_data& operator*=(const block_data& mutliply); // pointwise mulitplication
	block_data operator/(const block_data& divisor) const; // pointwise division
	block_data& operator/=(const block_data& divisor); // pointwise division

	bool operator==(const block_data& tocmp) const;

	block_data operator/(T denom) const;
	block_data& operator/=(T denom);

	block_data& operator*=(T factor);

	void shuffle(mt19937& generator, bool leaveFirstUntouched);
	int shuffleRestricted(mt19937& generator);

	string stats(string commentChar = "#") const;
};

template <typename T> block_data<T> operator*(T const& scalar, block_data<T> rhs);

template <typename T> block_data<T> operator*(block_data<T> lhs, T const& scalar);/*

template <typename T> block_data<T> getSumSquare(T average);

template <typename T> block_data<T> getStdDev(T average);

template <typename T> block_data<T> getCoeffofVar(T average);*/
