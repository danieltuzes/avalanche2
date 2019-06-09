// gen_param.h : stores user defined general function relatad to class param
//
#define GEN_PARAM_VERSION 0.08b

// version notes
// 0.08b: bool setFromLine(string s_line) has argument called s_line to not be covered by temporal local variable s
// 0.08: if value is modified by program call, printOut will tell it
// 0.07: parameter value can be accessed via the name of the variable, i.e. the following is defined: operator paramtype() const { return getVal(); }
// 0.06b: parameter value is padded
// 0.06: at summary, all parameters are printed out, so the info if the parameter is set from line or file or used by def

#pragma once
#include "stdafx.h"
#include "gen_siomanip.h"
using namespace std;



enum rq //to set requirement of something
{
	useless,  //it does not take effect to set the value; set a param's rq to this if you get to know uselessness in runtime
	must,     //you must provide it, othervise the program doesn't start
	optional, //you don't have to provide it, but setting this value may moidifes the program's behaviour
	elective  //you must use a value of a group, not used yet
};

rq opt_if(bool statement);

class basic_param
{
protected:
	bool set;
	bool set_by_call;
	bool set_from_file;
	bool moded_by_prg;
	bool disabled;
	string ext_name;
	rq require;
	string Val_str;
public:
	basic_param()
		: set(false), set_by_call(false), set_from_file(false), moded_by_prg(false), disabled(false), require(rq::optional) {}

	basic_param(string ext_name) : basic_param()  //set the ext_name of the parameter, it is almost always used
	{
		this->ext_name = ext_name;
	}

	basic_param(string ext_name, rq require) : basic_param(ext_name)  //set the ext_name of the parameter and to note if it is required or not
	{
		this->require = require;
	} 

	/*virtual basic_param * copy() = 0;*/

	void set_rq(rq q) {require = q;}
	
	bool is_set() const {return set;}
	bool is_disabled() const {return disabled;}
	string get_ext_name() const {return ext_name;}
	rq get_require() const {return require;}

	virtual string getValType()                                  = 0; //returns the typeid(...).name()
	virtual bool setFromLine (string fname)                      = 0; //basic_param is just an interface for template <> class param
	virtual bool setFromProgramCall(int argc, char * argv[])     = 0; //set parameter from formats like program calls arguments
	virtual bool setFromStream(istream&)                         = 0; //virtual functions are for reaching derived class type dependent members
	virtual bool setFromFile(string fname)                       = 0;

	virtual ostream& operator << (ostream& o) {return o;}
	virtual istream& operator >> (istream& i) {return i;}

	string getVal_str() const {return Val_str;}

	void unset() {set = false;}
	bool setByUsr() const {return (set_by_call || set_from_file || moded_by_prg);}; //if value has been set by the user; program modifies only if user input vlaues
	bool isModifiedByPrg() const {return moded_by_prg; };
	bool is_setFromCall() const {return set_by_call;};
	bool is_setFromFile() const { return set_from_file; };
};


class paramcontainer : public vector<basic_param*>
{
protected:
	int nof_sfpc;  //number of set from program call
public:
	paramcontainer();
	int setFromIf(int argc, char * argv[]);        // set every parameter from program call if possible
	int setFromIf(string fname);                   // set every parameter from a file if possible
	bool checkRequired() const;                    // check for every rq::must paramters to set	
	int getLongestNameSize() const;                // gives the longest paramter name in order to fancy print
	void printOut(ostream& o) const;               // prints out every value if it has been set
	bool summerize(ostream& o = cout) const;       // runs checkRequired() and printOut()
	void setFrom(int argc, char * argv[]) const;   // set every parameter from program call
	void description(ostream& o = cout) const;
};


template <typename paramtype> class param : public basic_param
{
protected:
	paramtype Val;
	paramtype defVal;

	bool setVal(const paramtype& value)
	{
		if (require == rq::useless)
			return false;

		if (setByUsr() && !moded_by_prg)
			cerr << "Warning: assign a value to " << (*this).ext_name << ": " << *this << " again." << endl;
		
		Val = value;
		stringstream o;
		o << value;
		Val_str = o.str();
		set = true;
		return true;
	}

public:

	explicit param(string ext_name) :
		basic_param(ext_name) {}


	explicit param(string ext_name, paramcontainer& p) :
		basic_param(ext_name) {p.push_back(this);} //declarator inheritance is not implemented in msvs2012

	explicit param(string ext_name, paramcontainer& p, rq require) :
		basic_param(ext_name, require), Val(paramtype()), defVal(paramtype())
	{
		p.push_back(this);
	} //declarator inheritance is not implemented in msvs2012
	
	explicit param(string ext_name, paramcontainer& p, rq require, paramtype def_val) :
		basic_param(ext_name, require), Val(paramtype()), defVal(paramtype())
	{
		p.push_back(this);
		setDefVal(def_val);
	}

	explicit param(string ext_name, paramcontainer& p, rq require, paramtype def_val, string fname) :
		basic_param(ext_name, require), Val(paramtype()) //set the paramter from a file
	{
		p.push_back(this);
		if(!setFromFile(fname))
			setDefVal(def_val);
	}

	explicit param(string ext_name, paramcontainer& p, rq require, paramtype def_val, int argc, char * argv[]) :
		basic_param(ext_name, require) //set the parameter from program call
	{
		p.push_back(this);
		if(!setFromProgramCall(argc,argv))
			setDefVal(def_val);
	}

	bool setDefVal(const paramtype& value)
	{
		if (require == rq::useless)
			return false;

		if (setByUsr())
			cerr << "Warning: assign a default value to " << (*this).ext_name << ": " << *this << " that has a user-assigned value." << endl;
		
		defVal = value;
		stringstream o;
		o << value;
		Val_str = o.str();
		set = true;
		return true;
	}

	void modify(paramtype val)
	{
		moded_by_prg = true;
		setVal(val);
	}

	string getValType(); //special cases are defined later


	bool setFromLine(string s_line) //set parameter from a single line; if setting is failed, the parameter will be disabled
	{
		stringstream i(s_line);
		bool retval = false;
		for (string s; i >> s; )
		{
			if (s==ext_name) //if variable name found implicitely
			{
				if (require == rq::useless)
				{
					cerr << "Setting parameter " << s << " is useless in this environment." << endl;
					return false;
				}
				paramtype MyVal;
				if (i >> MyVal)
				{

					retval = (retval || setVal(MyVal)); //open the possibility to set a param value multiple times
				}
				else
				{
					retval = false; //but cannot read the value of it
					disabled = true;
				}
			}
		}
		return retval;
	}
	
	bool setFromProgramCall(int argc, char * argv[]) //set parameter from formats like program calls arguments
	{
		string s;
		for (int i=1; i<argc; ++i) //argv[0] is the program name itself
			s += string(argv[i]) + string(" ");

		if (setFromLine(s))
		{
			set_by_call = true;
			return true;
		}
		else
			return false;

	}

	bool setFromStream(istream& i) //set parameter from stream, skips comments
	{
		bool retval = false;
		for (string s; getline(comment_skipped(i),s); )
			retval |= setFromLine(s); //open the possibility to set a param value multiple times

		return retval;
	}

	bool setFromFile(string fname) //set parameter from file, skips comments
	{
		ifstream i(fname);
		if (!i)
			cerr << "Cannot open " << fname << " to set " << ext_name << endl;


		if (setFromStream(i))
		{
			set_from_file = true;
			return true;
		}
		else
			return false;

	}
	
	paramtype operator() () const
		{return getVal();}

	paramtype getVal() const
	{
		if (this->setByUsr())
			return Val;
		else
			return defVal;
	}

	void printOut(ostream& o = cout) const
		{o << ext_name << " " << getVal();}

	string getNameNVal() const
	{
		stringstream o;
		printOut(o);
		return o.str();
	}

	operator paramtype() const { return getVal(); }
};

template <typename T> ostream& operator << (ostream& o, param<T>& MyParam)
{
	if (!MyParam.is_set()) o << "NA";
	else o << MyParam();
	return o;
}



template <typename T> istream& operator >> (istream& i, param<T>& MyParam)
{
	T value;
	i >> value;
	MyParam.setVal(value);
	return i;
}