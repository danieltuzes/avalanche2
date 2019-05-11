// gen_param.h : stores user defined general function relatad to class param
//


#include "stdafx.h"
#include "gen_param.h"
using namespace std;


rq opt_if(bool statement)
{
	if (statement)
		return rq::optional;
	else
		return rq::useless;
}


paramcontainer::paramcontainer() : nof_sfpc(0) {}


bool paramcontainer::checkRequired() const //check for every rq::must paramters to set
{
	bool retval = true;
	for (unsigned int i = 0; i < size(); ++i)
		if (at(i)->get_require() == rq::must && !at(i)->is_set())
		{
			cerr << "Error: required paramter " << at(i)->get_ext_name() << " not set by the user" << endl;
			retval = false;
		}

	return retval;
}
	

int paramcontainer::getLongestNameSize() const //gives the longest paramter name in order to fancy print
{
	size_t retval = 0; //the size of the biggest name of a set param
	for (unsigned int i = 0; i < size(); ++i)
		if (at(i)->get_ext_name().size() > retval)
		//if (at(i)->setByUsr() && at(i)->get_ext_name().size() > retval)
			retval = at(i)->get_ext_name().size();
			
	return static_cast<int>(retval);
}


void paramcontainer::printOut(ostream& o) const //prints out every value if it has been set
{
	streamsize orig_widt = o.width();
	ios_base::fmtflags orig_flags = o.flags();

	o << left;

	for (unsigned int i=0; i<size(); ++i)
	{
		o.width(getLongestNameSize() + 2);
		o << at(i)->get_ext_name();
		o.width(0);
		if (at(i)->is_disabled())
			o << "disabled";
		else
			o << at(i)->getVal_str();
		if (at(i)->is_setFromCall())
			o << " (set by user at program call)";
		else if (at(i)->is_setFromFile())
			o << " (set by user from a file)";
		else if (at(i)->setByUsr())
			o << " (set by user)";
		else
			o << " (set by default value)";
		if (at(i)->isModifiedByPrg())
			o << " (modified by prg in runtime)";
		o << endl;
	}

	o.width(orig_widt);
	o.flags(orig_flags);
}

	
bool paramcontainer::summerize(ostream& o) const
{
	printOut(o);
	if (checkRequired())
	{
		o << "Every needed parameters have been set.\n" << endl;
		return true;
	}
	else
		return false;
}

void paramcontainer::setFrom(int argc, char * argv[]) const
{
	stringstream line;
	for (int j=1; j < argc; ++j)
		line << argv[j] << " ";

	for (unsigned int i = 0; i<this->size(); ++i)
		this->at(i)->setFromLine(line.str());
}

//set every parameter from program call if possible
int paramcontainer::setFromIf(int argc, char * argv[])
{
	for (auto x = this->begin(); x!=this->end(); ++x)
	{
		if ((*x)->setFromProgramCall(argc,argv))
			++nof_sfpc;
	}
	return nof_sfpc;
}

int paramcontainer::setFromIf(string fname)
{
	cout << "Using " << fname << " to set parameters" << endl;
	int nof_sff = 0; //number of parameter set from the file
	for (auto x = this->begin(); x!=this->end(); ++x)
	{
		if ((*x)->setFromFile(fname))
			++nof_sff;
	}

	return nof_sff;
}


void paramcontainer::description(ostream& o) const
{
	bool firsttime = true;
	for (auto it = this->begin(); it != this->end(); ++it)
	{
		if ((*it)->get_require() == rq::must)
		{
			if (firsttime)
			{
				o << "Set of parameters required:" << endl;
				firsttime = false;
			}
			o << "\t" << (*it)->getValType() << " " << (*it)->get_ext_name() << endl;
		}
	}

	firsttime = true;
	for (auto it = this->begin(); it != this->end(); ++it)
	{
		if ((*it)->get_require() == rq::optional)
		{
			if (firsttime)
			{
				o << "Set of parameters optional:" << endl;
				firsttime = false;
			}
			o << "\t" << (*it)->getValType() << ": " << (*it)->get_ext_name();
			if ((*it)->is_set())
				o << " (" << (*it)->getVal_str() << ")";
			o << endl;
		}
	}

	firsttime = true;
	for (auto it = this->begin(); it != this->end(); ++it)
	{
		if ((*it)->get_require() == rq::useless)
		{
			if (firsttime)
			{
				o << "Set of parameters not available with the parameter values provided:" << endl;
				firsttime = false;
			}
			o << "\t" << (*it)->getValType() << " " << (*it)->get_ext_name() << endl;
		}
	}
}



template <> string param<string>::getValType() {return "string";}
template <> string param<unsigned long int>::getValType() { return "unsigned long int"; }
template <> string param<int>::getValType() {return "int";}
template <> string param<bool>::getValType() {return "bool";}
template <> string param<double>::getValType() {return "double";}
template <> string param<long double>::getValType() { return "long double"; }
template <> string param<float>::getValType() {return "float";}

