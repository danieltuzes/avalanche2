


#include "stdafx.h"
#include "gen_fmanip.h"
#include "gen_version.h"
#include "gen_math.h"
#include "gen_siomanip.h"

using namespace std;

ifstream& operator>>(ifstream& i, dataformat& df)
{
	string name;
	i >> name;
	if (name.compare("binary"))
		df = binary;
	else if (name.compare("text"))
		df = text;
	else if (name.compare("zip"))
		df = zip;
	else if (name.compare("bin_text"))
		df = bin_text;
	else if (name.compare("text_bin"))
		df = text_bin;

	return i;
}

ofstream& operator<<(ofstream& o, const dataformat& df)
{
	if (df == binary)
		o << "binary";
	else if (df == text)
		o << "text";
	else if (df == zip)
		o << "zip";
	else if (df == bin_text)
		o << "bin_text";
	else if (df == text_bin)
		o << "text_bin";

	return o;
}

ifstream::pos_type filesize(string fname) //tell the size of a file, exit if doesn't exist
{
    ifstream in(fname, ifstream::in | ifstream::binary);
	fstream_must_good((bool)in, fname);
	in.seekg(0, ifstream::end);
    return in.tellg(); 
}
ifstream::pos_type filesize(const char * fname) {return filesize(string(fname));}





bool is_file_exist(string fname)
{
	ifstream in(fname);
	if (in)
		return true;
	else
		return false;
}
bool is_file_exist(char * fname) {return is_file_exist(string(fname));}

void fstream_must_good(bool stream_is_good, string fname)
{
	if (!stream_is_good)
	{
		cerr << "Error: Cannot open " << fname << endl
			<< "Program termiantes." << endl;
		exit(0);
	}
}

string lastNonEmptyLine(string fname, int lineSize)
{
	ifstream bin_ifile(fname, ios_base::binary);
	if (!bin_ifile)
	{
		cerr << "Can't open " << fname;
		return string();
	}

	int bufSize = lineSize; // because why not? tweak if need to
	char * buf;
	buf = new char[bufSize];
	string line;

	int seek, nloff;
	bool foundNonNewLine = false;
	// iterate over multiples of bufSize while file ok
	for (size_t n = 1; bin_ifile; ++n)
	{
		// next seek position will be a multiple of bufSize
		seek = -static_cast<int>(n * bufSize);
		bin_ifile.seekg(seek, bin_ifile.end);
		// read "bufSize" bytes into buffer
		bin_ifile.read(buf, bufSize);

		// in case no newline found, seek past eof
		nloff = -seek;
		// find offset of last newline in buffer
		for (int i = bufSize - 1; i >= 0; --i)
		{
			if (buf[i] == '\n' || buf[i] == '\r')
			{
				if (foundNonNewLine)
				{
					nloff = i;
					break;
				}
			}
			else
				foundNonNewLine = true;
		}
		seek += nloff + 1; // new seek position is one character after found newline
		if (seek >= 0) continue; // just kidding about the "past eof" part ;)

		// seek to after found newline and get line
		bin_ifile.seekg(seek, bin_ifile.end);
		getline(bin_ifile, line);
		if (!line.empty()) break; // have result, break and return
	}

	delete[] buf;

	if (bin_ifile) return line;
	else
		return string();
}

string lastLine(ifstream& bin_ifile, int lineSize) //finds the last line of the stream; lineSize is just an estimation on line size
{
    if (!bin_ifile)
	{
		cerr << ("Bad stream on input");
		return string();
	}

    size_t bufSize = lineSize; // because why not? tweak if need to
    char * buf;
	buf = new char [bufSize];
    string line;

    int seek, nloff;
    // iterate over multiples of bufSize while file ok
    for (size_t n = 1; bin_ifile; ++n)
    {
        // next seek position will be a multiple of bufSize
        seek = -static_cast<int>(n * bufSize);
        bin_ifile.seekg(seek, bin_ifile.end);
        // read "bufSize" bytes into buffer
        bin_ifile.read(buf, bufSize);

        // in case no newline found, seek past eof
        nloff = -seek;
        // find offset of last newline in buffer
        for (size_t i = 0; i < bufSize; ++i)
        {
            if (buf[i] == '\n') nloff = i;
        }
        seek += nloff + 1; // new seek position is one character after found newline
        if (seek >= 0) continue; // just kidding about the "past eof" part ;)

        // seek to after found newline and get line
        bin_ifile.seekg(seek, bin_ifile.end);
        getline(bin_ifile, line);
        if (!line.empty()) break; // have result, break and return
    }

	delete[] buf;

    if (bin_ifile) return line;
    else
		return string();
}

template <typename T> block_data<T>::block_data() : size(0), data(NULL) {}

template <typename T> block_data<T>::block_data(int size) : size(size), power(::getPower(size)), data(new T[size*size]) { }

template <typename T> block_data<T>::block_data(int size, T val) : block_data(size)
{
	fill(data, data + size*size, val);
}

template <typename T> block_data<T>::block_data(const block_data<T>& tocopy) : block_data(tocopy.getSize())
{
	memcpy(data, tocopy.data, size*size*sizeof(T));
}

template <typename T> block_data<T>::block_data(vector<vector<T>> m_vector) : block_data(m_vector.size())
{
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			this->operator()(i, j) = m_vector[i][j];
}

template <typename T> block_data<T>::block_data(vector<T> m_vector) : block_data(pow2(lb(m_vector.size()) / 2))
{
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			this->operator()(i, j) = m_vector[i*size + j];
}

template <typename T> block_data<T>::block_data(string fname, dataformat f, int skip) : block_data()
{
	//cout << "block_data string, dataformat, int constructor called" << endl;
	if (!is_file_exist(fname))
	{
		cerr << fname << " doesn't exist, cannot construct block_data." << endl;
		exit(0);
	}

	if (f==binary)
		readFromFile(fname);
	else if (f==text)
		fReadFromFile(fname);
	else if (f == bin_text)
	{
		if (!readFromFile(fname))
		{
			cerr << "Can't construct a block_data element from binary file " << fname << endl
				<< "Try to read it as a text file instead" << endl;
			fReadFromFile(fname);
		}
		cout << "block_data size: " << size << "x" << size << endl;
	}
	else if (f == text_bin)
	{
		if (!fReadFromFile(fname))
		{
			cerr << "Constructor can't construct a block_data element from text file " << fname << endl
				<< "Program tries to read it as a binary file instead" << endl;
			readFromFile(fname);
		}
		cout << "block_data size: " << size << "x" << size << endl;
	}
	else if (f == origin)
	{
		originReadFromFile(fname, skip);
	}
	else if (f == bin_text_origin)
	{
		if (!readFromFile(fname))
		{
			cerr << "Can't construct a block_data element from binary file " << fname << endl
				<< "Try to read it as a text file instead" << endl;
			if (!fReadFromFile(fname))
			{
				cerr << "Can't construct a block_data element from text file " << fname << endl
					<< "Try to read it as an origin file instead" << endl;
				originReadFromFile(fname, skip);
			}
		}
		cout << "block_data size: " << size << "x" << size << endl;
	}
	else
		cerr << "dataformat with ID " << f << " is not supported" << endl;
	//cout << "block_data string, dataformat, int constructor ended" << endl;

}

template <typename T> block_data<T>::block_data(const char * fname, dataformat f, int skip) : block_data(string(fname), f, skip) {}

template <typename T> block_data<T>::~block_data()
{
	delete[] data;
}

template <typename T> bool block_data<T>::hasData() {return data!=NULL;}

template <typename T> block_data<T>& block_data<T>::operator=(const block_data& tocopy)
{
	block_data<T> tmp(tocopy);
	*this = move(tmp);
	return *this;
}

template <typename T> block_data<T>& block_data<T>::operator=(block_data&& tocopy)
{
	delete[] data;
	data = tocopy.data;
	tocopy.data = nullptr;
	return *this;
}

template <typename T> block_data<T>& block_data<T>::operator=(T initVal)
{
	fill(data, data + size*size, initVal);
	return *this;
}

template <typename T> bool block_data<T>::init(int size, T val) 
{return init(block_data<T>(size,val));}


template <typename T> bool block_data<T>::init(const block_data<T>& newData) 
{
	if (hasData())
		delete[] data;

	block_data::size = newData.getSize();
	block_data::power = newData.getPower();
	
	data = new T[size*size];

	for (int i=0; i<size*size; ++i)
		data[i] = newData[i];

	return true;
}

template <typename T> bool block_data<T>::readFromFile(string fname)
{
	//cout << "block_data readFromFile called" << endl;

	double flinesize_d = sqrt(filesize(string(fname))/sizeof(T));
	int flinsize = static_cast<int>(flinesize_d);
	if (flinesize_d - flinsize > 1e-10) // flinesize_d is not an integer
	{
		cerr << "Warning: " << fname << " contains " << flinesize_d << " number of rows and cols" << endl;
		return false;
	}

	int powr_size = 1; 
	for (;powr_size<=flinsize; powr_size*=2);
	powr_size /=2;
	if (!(flinsize -1 == powr_size || flinsize == powr_size))
	{
		cout << "Warning: file contains not power of 2 sized data block, nor with a 1 extra border value." << endl;
		return false;
	}

	FILE * ifile;
	if ((ifile = fopen(fname.c_str(),"rb")) == NULL)
	{
		cerr << "Can't open " << fname << endl;
		return false;
	}
	//cout << "flinsize " << flinsize << ", size: " << size << ", powr_size:" << powr_size << endl;

	if (size != powr_size)
	{
		size = powr_size;
		power = ::getPower(size);
		delete[] data;
		data = new T[size*size];
	}
	//cout << "block_data readFromFile after mem alloc called" << endl;
	for (int i = 0; i < flinsize; ++i)
	{
		//cout << "i: " << i << endl;
		for (int j = 0; j < flinsize; ++j)
		{
			//cout << "j: " << j << endl;

			if (i >= size || j >= size)
			{
				T rubbish;
				if (fread(&rubbish, sizeof(T), 1, ifile) != 1)
				{
					cerr << "Can't read rubbish element " << i << " " << j << endl;
					return false;
				}

			}

			else if (fread(&data[i*size + j], sizeof(T), 1, ifile) != 1)
			{
				cerr << "Can't read element " << i << " " << j << endl;
				return false;
			}
		}
	}

		
	fclose(ifile); //a FILE * should be closed manually
	//cout << "block_data readFromFile fclose called" << endl;

	return true;
}
template <typename T> bool block_data<T>::readFromFile(char * fname) {return readFromFile(string(fname));}


template <typename T> bool block_data<T>::fReadFromFile(string fname)
{
	stringstream lastline(lastNonEmptyLine(fname, 80));
	int x, y;
	lastline >> x >> y;
	if (!lastline)
	{
		cerr << "Cannot read the last two indices." << endl;
		return false;
	}

	else if (x != y)
	{
		cerr << "Last line indices are not equal (" << x << ", " << y << ")." << endl;
		return false;
	}
	
	int tmp_power = ::getPower(x);
	if ((1 << tmp_power) != x && (1 << tmp_power) != x + 1)
	{
		cout << "Error: file contains not power of 2 sized data block, nor with a 1 extra border value." << endl;
		return false;
	}
	power = tmp_power;
	size = (1 << power);

	ifstream in(fname);
	skip_comment(in);

	data = new T[size*size];



	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
		{
			int checkx = 0, checky = 0;
			T probe;
			if(!(in >> checkx >> checky >> probe))
			{
				cerr << "Can't read element at " << x << " " << y << endl;
				return false;
			}
			if (checkx == size || checky == size)
			{
				--y;
				continue;
			}

			if (checkx > size || checky > size)
			{
				cerr << "Wrong index value at \"" << checkx << " " << checky << "\" in file " << fname << endl;
				return false;
			}
			data[(x << power) + y] = probe;
		}
		
	return true;
}
template <typename T> bool block_data<T>::fReadFromFile(char * fname) {return fReadFromFile(string(fname));}

template <typename T> bool block_data<T>::originReadFromFile(string fname, int skip)
{
	ifstream ifile(fname);
	if (!ifile)
	{
		cerr << "Cannot open " << fname;
		return false;
	}
	vector<vector<T>> filearray;
	for (int i = 0; i < skip; ++i, ifile.ignore(numeric_limits<streamsize>::max(), '\n'));

	for (string line; getline(ifile, line);)
	{
		istringstream oneline(line);
		for (int j = 0; j < skip; ++j, oneline.ignore(numeric_limits<streamsize>::max(), '\t'));
		filearray.push_back(vector<T>());
		for (T tmp; oneline >> tmp; filearray.back().push_back(tmp));
	}
	size = filearray.size();
	power = ::getPower(size);

	data = new T[size*size];

	*this = block_data(filearray);

	return true;
}



template <typename T> block_data<T> block_data<T>::reMap() const
{
	block_data<T> ret(size / 2);

	for (int i = 0; i < size / 2; ++i)
	{
		for (int j = 0; j < size / 2; ++j)
		{
			ret(i, j) = periodic((i * 2) - 1, (j * 2) - 1) / 16;
			ret(i, j) += periodic((i * 2) - 1, (j * 2)) / 8;
			ret(i, j) += periodic((i * 2) - 1, (j * 2) + 1) / 16;
			
			ret(i, j) += periodic((i * 2), (j * 2) - 1) / 8;
			ret(i, j) += periodic((i * 2), (j * 2)) / 4;
			ret(i, j) += periodic((i * 2), (j * 2) + 1) / 8;

			ret(i, j) += periodic((i * 2) + 1, (j * 2) - 1) / 16;
			ret(i, j) += periodic((i * 2) + 1, (j * 2)) / 8;
			ret(i, j) += periodic((i * 2) + 1, (j * 2) + 1) / 16;
		}
	}

	return ret;
}

template <typename T> block_data<T> block_data<T>::reSample(int blockSize, int shiftx, int shifty) const
{
	int blockPower = ::getPower(blockSize);
	if (1 << ::getPower(blockSize) != blockSize)
		cerr << "blockSize is not a power of 2" << endl;

	if (shiftx >= blockSize)
		cerr << "shiftx must be smaller than blocksize" << endl;
	if (shifty >= blockSize)
		cerr << "shifty must be smaller than blocksize" << endl;

	block_data<T> ret(size / blockSize);

	for (int i = 0; i < size / blockSize; ++i)
	for (int j = 0; j < size / blockSize; ++j)
		ret(i, j) = periodic((i << blockPower) + shiftx, (j << blockPower) + shifty);

	return ret;
}

template <typename T> bool block_data<T>::symmetrize(block_data<T>& dest) const
{
	dest.init(size,0);
	for (int i=0; i<size; ++i)
	{
		for (int j=0; j<size; ++j)
		{
			dest(i,j) = (operator()(i,j) + operator()(j,i)) / 2;
		}
	}
	return true;
}



template <typename T> int block_data<T>::getSize() const {return size;}
template <typename T> int block_data<T>::getPower() const {return power;}

template <typename T> T block_data<T>::operator[](int i) const {return data[i];}
template <typename T> T& block_data<T>::operator[](int i) {return data[i];}
template <typename T> T* block_data<T>::getData() { return data; }
template <typename T> T* block_data<T>::getData() const { return data; }
template <typename T> T block_data<T>::operator()(int x, int y) const {return data[(x << power) + y];}
template <typename T> T& block_data<T>::operator()(int x, int y) {return data[(x << power) + y];}
template <typename T> T block_data<T>::periodic(int x, int y) const { return operator()((x + size) & (size - 1), (y + size) & (size - 1)); }
template <typename T> T& block_data<T>::periodic(int x, int y) { return operator()((x + size) & (size - 1), (y + size) & (size - 1)); }

template <typename T> T block_data<T>::periodic_x(int x, int y) const { return operator()((x + size) & (size - 1), y); }
template <typename T> T& block_data<T>::periodic_x(int x, int y) { return operator()((x + size) & (size - 1), y); }
template <typename T> T block_data<T>::periodic_y(int x, int y) const { return operator()(x, (y + size) & (size - 1)); }
template <typename T> T& block_data<T>::periodic_y(int x, int y) { return operator()(x, (y + size) & (size - 1)); }

template <typename T> T block_data<T>::rightAv(int x, int y) const { return (operator()(x, y) + periodic_x(x + 1, y)) / 2; } // gives the average of the cell's and it's right neighbour's values
template <typename T> T block_data<T>::localAv(int x, int y, int radX, int radY) const // gives the local average around the border (x,y) and (x+1,y)
{
	T sum = 0;
	for (int xp = x - (radX - 1); xp != x + radX; ++xp)
	for (int yp = y - radY; yp != y + radY; ++yp)
		sum += periodic(xp, yp);

	return sum / ((radX + 1)*(radY + 1));
}
template <typename T> T block_data<T>::rightGrad(int x, int y, T dist) const { return ((periodic_x(x + 1, y) - operator()(x, y)))/dist; } // gives the average of the cell's and it's right neighbour's values

template <typename T> vector<vector<T>> block_data<T>::std2Dvector() const
{
	vector<vector<T>> ret(size, vector<T>(size));
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			ret[i][j] = this->operator()(i, j);
		}
	}
	return ret;
}

template <typename T> vector<T> block_data<T>::slicex(int y0) const
{
	vector<T> ret;
	for (int x = 0; x < size; ++x)
		ret.push_back(operator()(x, y0));

	return ret;
}

template <typename T> vector<T> block_data<T>::slicey(int x0) const
{
	vector<T> ret;
	for (int y = 0; y < size; ++y)
		ret.push_back(operator()(x0, y));

	return ret;
}

template <typename T> vector<T> block_data<T>::projectx() const
{
	vector<T> ret(size, 0);
	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
			ret[x] += operator()(x, y);

	return ret;
}

template <typename T> vector<T> block_data<T>::projecty() const
{
	vector<T> ret(size, 0);
	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
			ret[y] += operator()(x,y);

	return ret;
}

template <typename T> vector<double> block_data<T>::rdependence(double x0, double y0) const
{
	int length = static_cast<int>(sqrt(2)*size/2) + 2;
	vector<double> value_ret(length, 0);
	vector<double> counter(length, 0);

	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
		{
			double pos = dist(x, y, x0, y0, size);
			int a = static_cast<int>(pos);
			int b = a + 1;
			double weighta = (b - pos);
			double weightb = (pos - a);

			counter[a] += weighta;
			counter[b] += weightb;
			value_ret[a] += static_cast<double>(weighta*operator()(x, y));
			value_ret[b] += static_cast<double>(weightb*operator()(x, y));
		}

	for (int i = 0; i < length; ++i)
	{
		if (counter[i] != 0)
			value_ret[i] /= counter[i];
		else
			value_ret[i] = 0;

	}

	return value_ret;
}

template <typename T> vector<double> block_data<T>::rabsdependence(double x0, double y0) const
{
	int length = static_cast<int>(sqrt(2)*size / 2) + 2;
	vector<double> value_ret(length, 0);
	vector<double> counter(length, 0);

	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
		{
			double pos = dist(x, y, x0, y0, size);
			int a = static_cast<int>(pos);
			int b = a + 1;
			double weighta = (b - pos);
			double weightb = (pos - a);

			counter[a] += weighta;
			counter[b] += weightb;
			value_ret[a] += static_cast<double>(weighta*abs(operator()(x, y)));
			value_ret[b] += static_cast<double>(weightb*abs(operator()(x, y)));
		}

	for (int i = 0; i < length; ++i)
	{
		if (counter[i] != 0)
			value_ret[i] /= counter[i];
		else
			value_ret[i] = 0;

	}

	return value_ret;
}

template <typename T> vector<double> block_data<T>::rcos4fidependence(double x0, double y0) const
{
	int length = static_cast<int>(sqrt(2)*size / 2) + 2;
	vector<double> value_ret(length, 0);
	vector<double> counter(length, 0);

	for (int x = 0; x < size; ++x)
		for (int y = 0; y < size; ++y)
		{
			double pos = dist(x, y, x0, y0, size);
			int a = static_cast<int>(pos);
			int b = a + 1;
			double weighta = (b - pos);
			double weightb = (pos - a);

			counter[a] += weighta;
			counter[b] += weightb;
			value_ret[a] += static_cast<double>(weighta*operator()(x, y)*cos(4 * atan(dist(static_cast<double>(y), y0, size))));
			value_ret[b] += static_cast<double>(weightb*operator()(x, y)*cos(4 * atan(dist(static_cast<double>(y), y0, size))));
			// itt valamiért az alábbiak álltak
			//value_ret[a] += static_cast<double>(weighta*operator()(x, y)*cos(4 * atan(dist(static_cast<double>(y), y0, size) / dist(static_cast<double>(x), x0, size))));
			//value_ret[b] += static_cast<double>(weightb*operator()(x, y)*cos(4 * atan(dist(static_cast<double>(y), y0, size) / dist(static_cast<double>(x), x0, size))));
		}

	for (int i = 0; i < length; ++i)
	{
		if (counter[i] != 0)
			value_ret[i] /= counter[i];
		else
			value_ret[i] = 0;

	}

	return value_ret;
}

template <typename T> bool block_data<T>::fWriteToFile(string fname, string header, int precision, bool swapOrder) const //(over)write formatted data to fname; create file if not exist
{
	return fWriteToFileTreatNN(fname, false, 1, false, header, precision, swapOrder);
}

template <typename T> bool block_data<T>::fWriteToFileTreatNN(string fname, bool specialNan, double treatNanAs, bool autoheader, string header, int precision, bool swapOrder) const //(over)write formatted data to fname; create file if not exist
{
	ofstream o(fname);
	if (!o)
	{
		cerr << "Can't create " << fname << endl;
		return false;
	}
	cout << "Writing " << fname << endl;

	block_data<T> datacpy = *this;
	if (specialNan)
		for (int i = 0; i < size*size; ++i)
		{
			if (std::isnan(static_cast<double>(datacpy[i]))) // libstdc++ bug, https://gcc.gnu.org/bugzilla/show_bug.cgi?id=60407
				datacpy[i] = static_cast<T>(treatNanAs);
			if (datacpy[i] == numeric_limits<double>::infinity())
				datacpy[i] = 1;
			if (datacpy[i] == -numeric_limits<double>::infinity())
				datacpy[i] = -1;
		}

	o << header;
	if (autoheader)
		o << datacpy.stats();
	o << setprecision(precision);
	if (swapOrder)
	{
		for (int y = 0; y < size; ++y)
		{
			for (int x = 0; x < size; ++x)
				o << x << "\t" << y << "\t" << datacpy[(x << power) + y] << "\n";
			if (y != size - 1)
				o << "\n";
		}
	}
	else
	{
		for (int x = 0; x < size; ++x)
		{
			for (int y = 0; y < size; ++y)
				o << x << "\t" << y << "\t" << datacpy[(x << power) + y] << "\n";
			o << x << "\t" << size << "\t" << 0 << "\n\n";
		}

		for (int y = 0; y <= size; ++y)
			o << size << "\t" << y << "\t" << 0 << "\n";
	}


	return true;
}


template <typename T> bool block_data<T>::writeToFile(string fname) const//(over)write data to file fname; create file if not exist
{
	FILE * o;
	if ((o=fopen(fname.c_str(),"wb")) == NULL)
	{
		cerr << "Can't open " << fname << endl;
		return false;
	}

	cout << "Writing " << fname << endl;
	int s = 0;

	for (int i=0; i<size*size; ++i)
		s += fwrite(&data[i],sizeof(T),1,o);
		
	if (s != size * size)
	{
		cerr << "Writing " << size*size << " number of values (with typeof " << typeid(T).name() << ") to " << fname << " was not successful." << endl;
		return false;
	}
	fclose(o); //a FILE * should be closed manually
	return true;
}
//template <typename T> bool block_data<T>::writeToFile(char * fname) const
//	{return writeToFile(string(fname));} //(over)write data to fname; create file if not exist


template <typename T> void block_data<T>::sumNormalize(T n) //devide every data value to make sum of them to the given n
{
	T sum = getSum();

	for (int i = 0; i < size*size; ++i)
		data[i] *= n / sum;
}



template <typename T> void block_data<T>::maxNormalize(T n) //devide every data value to make max value to the given n
{
	T max = getMax();

	for (int i=0; i<size*size; ++i)
			data[i] *= n / max;
}

template <typename T> block_data<T> block_data<T>::Abs() const // return a new block_data with the absolute values of the original
{
	block_data tmp(*this);
	for (int i = 0; i < size*size; ++i)
		tmp[i] = abs(tmp[i]);

	return tmp;
}

template <typename T> T block_data<T>::getSum() const //get the sum of every value of the data
{
	T sum = 0;
	for (int i=0; i<size*size; ++i)
		sum += data[i];

	return sum;
}

template <typename T> double block_data<T>::getAverage() const //get the average of the data
{
	return static_cast<double>(this->getSum()) / size / size;
}


template <typename T> T block_data<T>::getSumSquare() const //get the sum of every value of the data
{
	T sumSquare = 0;
	for (int i = 0; i<size*size; ++i)
		sumSquare += data[i] * data[i];

	return sumSquare;
}

template <typename T> long double block_data<T>::getStdDev(long double average) const //get standard deviation if average is given
{
	long double sum = 0;
	for (int i = 0; i < size*size; ++i)
		sum += (data[i] - average) * (data[i] - average);

	return sqrt(sum / (size * size));
}

template <typename T> long double block_data<T>::getStdDev() const //get standard deviation if average is unknown
{
	return this->getStdDev(this->getAverage());
}


template <typename T> long double block_data<T>::getCoeffofVar(long double average) const //get the coefficient of variation
{
	return this->getStdDev(average) / average;
}

template <typename T> long double block_data<T>::getCoeffofVar() const //get the coefficient of variation
{
	return this->getStdDev(this->getAverage()) / this->getAverage();
}


template <typename T> T block_data<T>::getMax() const //get the sum of every value of the data
{
	T max = data[0];
	for (int i=1; i<size*size; ++i)
		if (data[i] > max)
			max = data[i];

	return max;
}


template <typename T> T block_data<T>::getMin() const //get the sum of every value of the data
{
	T min = data[0];
	for (int i=1; i<size*size; ++i)
		if (data[i] < min)
			min = data[i];

	return min;
}

template <typename T> long double block_data<T>::getAssymetry() const //get Sum data(r) - data(-r)
{
	T diff = 0;
	for (int i = -size/2; i < size / 2; ++i)
	for (int j = -size/2; j < 1; ++j)
	{
		if (i >= 0 && (j == 0 || j == -size / 2))
			continue;
		diff += abs(periodic(i, j) - periodic(-i, -j));
	}
	
	return diff;
}




template <typename T> tuple<T, int, int> block_data<T>::getMaxWhere() const
{
	T max = data[0];
	int coord = 0;
	for (int i = 0; i < size*size; ++i)
	{
		if (data[i] > max)
		{
			max = data[i];
			coord = i;
		}
	}
	return make_tuple(max, getX(coord,power), getY(coord,size-1));
}


template <typename T> tuple<T, int, int> block_data<T>::getMinWhere() const
{
	T min = data[0];
	int coord = 0;
	for (int i=0; i<size*size; ++i)
		if (data[i] < min)
		{
			min = data[i];
			coord = i;
		}
	
	return make_tuple(min, getX(coord, power), getY(coord, size - 1));
}

	
template <typename T> block_data<T> block_data<T>::operator+(const block_data<T>& toAdd) const
{
	block_data tmp(*this);
	return tmp+=toAdd;
}

template <typename T> block_data<T>& block_data<T>::operator+=(const block_data<T>& toAdd)
{
	for (int i=0; i<size*size; ++i)
		data[i] += toAdd[i];

	return *this;
}

template <typename T> block_data<T> block_data<T>::operator-(const block_data<T>& toSub) const
{
	block_data tmp(*this);
	return tmp-=toSub;
}

template <typename T> block_data<T>& block_data<T>::operator-=(const block_data<T>& toSub)
{
	for (int i=0; i<size*size; ++i)
		data[i] -= toSub[i];

	return *this;
}


template <typename T> block_data<T> block_data<T>::operator*(const block_data<T>& mutliply) const
{
	block_data tmp(*this);
	return tmp *= mutliply;
}

template <typename T> block_data<T>& block_data<T>::operator*=(const block_data<T>& mutliply)
{
	for (int i = 0; i<size*size; ++i)
		data[i] *= mutliply[i];

	return *this;
}

template <typename T> block_data<T> block_data<T>::operator/(const block_data<T>& divisor) const
{
	block_data tmp(*this);
	return tmp /= divisor;
}

template <typename T> block_data<T>& block_data<T>::operator/=(const block_data<T>& divisor)
{
	for (int i = 0; i<size*size; ++i)
		data[i] /= divisor[i];

	return *this;
}


template <typename T> block_data<T> block_data<T>::operator+(T toAdd) const
{
	block_data tmp(*this);
	return tmp += toAdd;
}

template <typename T> block_data<T>& block_data<T>::operator+=(T toAdd)
{
	for (int i = 0; i<size*size; ++i)
		data[i] += toAdd;

	return *this;
}

template <typename T> block_data<T> block_data<T>::operator-(T toSub) const
{
	block_data tmp(*this);
	return tmp -= toSub;
}

template <typename T> block_data<T>& block_data<T>::operator-=(T toSub)
{
	for (int i = 0; i < size*size; ++i)
		data[i] -= toSub;

	return *this;
}


template <typename T> bool block_data<T>::operator==(const block_data<T>& toCmp) const
{
	if (size != toCmp.getSize()) return false;
	for (int i=0; i<size*size; ++i)
		if (data[i] != toCmp[i]) return false;
	return true;
}



template <typename T> block_data<T> block_data<T>::operator/(T denom) const
{
	block_data tmp(*this);
	return tmp /= denom;
}



template <typename T> block_data<T>& block_data<T>::operator/=(T denom)
{
	for (int i=0; i<size*size; ++i)
		data[i] /= denom;

	return *this;
}



template <typename T> block_data<T>& block_data<T>::operator*=(T factor)
{
	for (int i=0; i<size*size; ++i)
		data[i] *= factor;

	return *this;
}

template <typename T> void block_data<T>::shuffle(mt19937& generator, bool leaveFirstUntouched)
{
	for (int i = size*size - 1; i > static_cast<int>(leaveFirstUntouched); --i)
	{
		uniform_int_distribution<int> distr(static_cast<int>(leaveFirstUntouched), i);
		swap(data[i], data[distr(generator)]);
	}
}

template <typename T> int block_data<T>::shuffleRestricted(mt19937& generator)
{
	int c = 0; //how many times had been the generator advanced
	for (int i = size * size / 2 + size / 2 - 1 - 2; i > 1; --i)
	{
		int sw = uniform_int_distribution<int>(1, i)(generator);
		++c;
		if (sw >= size / 2)
			++sw;
		if (sw >= size*size / 2)
			++sw;

		int linpos = i;
		if (linpos >= size / 2)
			++linpos;
		if (linpos >= size*size / 2)
			++linpos;
		swap(data[linpos], data[sw]);

		int x = linpos >> power;
		int y = linpos & (size - 1);

		int sw_x = sw >> power;
		int sw_y = sw & (size - 1);
		swap(periodic(-x, -y), periodic(-sw_x, -sw_y));
	}
	return c;
}

template <typename T> void block_data<T>::periodicIncr(int& x) const
{
	++x;
	x &= (size - 1);
}

template <typename T> block_data<T> operator*(T const& scalar, block_data<T> rhs)
{
	return rhs *= scalar;
}


template <typename T> block_data<T> operator*(block_data<T> lhs, T const& scalar)
{
	return lhs *= scalar;
}



template <typename T> string block_data<T>::stats(string commentChar) const
{
	ostringstream text;
	auto max = getMaxWhere();
	auto min = getMinWhere();
	text << commentChar << "Maximum value: " << get<0>(max) << " at (x;y): (" << get<1>(max) << ";" << get<2>(max) << ")" << endl;
	text << commentChar << "Minimum value: " << get<0>(min) << " at (x;y): (" << get<1>(min) << ";" << get<2>(min) << ")" << endl;
	text << commentChar << "Sum: " << getSum() << endl;
	text << commentChar << "Abs sum: " << Abs().getSum() << endl;
	text << commentChar << "Assymetry: " << getAssymetry() << endl;
	text << commentChar << "Standard deviation: " << getStdDev() << endl;

	return text.str();
}


template block_data<long double> operator*(long double const&, block_data<long double>);
template block_data<long double> operator*(block_data<long double>, long double const&);
template block_data<double> operator*(double const&, block_data<double>);
template block_data<double> operator*(block_data<double>, double const&);
template block_data<int> operator*(int const&, block_data<int>);
template block_data<int> operator*(block_data<int>, int const&);

template class block_data<int>;
template class block_data<float>;
template class block_data<double>;      //template Instantiation
template class block_data<long double>; //template Instantiation