// to stop applications in a smarter way than sending sigterm signal to them

#include "stdafx.h"
#include "gen_interrupt.h"

using namespace std;

bool interrupt(string fnamePart)
{
	if (is_file_exist("stopAllAt" + fnamePart + ".txt"))
	{
		cerr << "stopAllAt" + fnamePart + ".txt" << " found." << endl;
		return true;
	}
	return false;
}

bool interrupt(string fnamePart, int seed)
{
	ifstream ifile("stopSpecificAt" + fnamePart + ".txt");
	if (!ifile)
		return false;

	for (int s; comment_skipped(ifile) >> s;)
		if (seed == s)
		{
			cerr << "stopSpecificAt" + fnamePart + ".txt" << " found with value " << s << endl;
			return true;
		};
		
	return false;
		
}

