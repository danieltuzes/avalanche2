

#include "stdafx.h"
#include "gen_version.h"
using namespace std;

void intro(string header, string footer, bool printExtraInfo, char indent_char, int indent_size, int wide)
{
	stringstream moreindent;
	for (int i=0; i<indent_size; ++i)
		moreindent << indent_char;
	moreindent << " ";

	stringstream separator;
	for (int i=0; i<wide; ++i)
		separator << indent_char;

	stringstream content;
	content << header << endl;
	if (printExtraInfo)
	{
		content << "AUTHORS: Peter Dusan Ispanovity (ispanovity@metal.elte.hu), Daniel Tuzes (tuzes@metal.elte.hu)\n"
			<< "COMPATIBILITY: This program is written in c++,\ntested under MSVS2013 and gcc4.7,\nbut should work under every c++11 compatible compiler.\n"
			<< "COPYRIGHT: some part of the program contains external sources\nthat are copyrighted under various terms.\nThe rest of the program and the sources of it are copyrighted under\nAttribution-NonCommercial-ShareAlike 2.5 Hungary.\n"
			<< "MORE INFO: for more information,\nplease search for a doc folder and look up its contents.\n";
	}
	content << footer;

	cout << endl << separator.str() << endl;
	int linesize = 0;
	for (unsigned int i=0; i<content.str().size(); ++i)
	{
		if (linesize == wide - 1) {
			cout << "\n";
			linesize = 0;
		}

		if (linesize == 0)
		{
			 cout << moreindent.str();
			 linesize += moreindent.str().size();
		}

		cout << content.str()[i];
		++linesize;

		if (content.str()[i] == '\n') {
			cout << moreindent.str();
			linesize = 0;
			linesize += moreindent.str().size();
		}
	}
	cout << endl << separator.str() << endl << endl;
}