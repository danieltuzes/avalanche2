


#include "stdafx.h"
#include "gen_siomanip.h"
using namespace std;





istream& comment_skipped(istream& stream, const char comment_symbol) //move stream position to the first non-commented line TODO: check compatibility with non-binary ifstream
{
	streamoff pos_b_r = stream.tellg();
	for (string oneline; getline(stream,oneline); )
	{
		unsigned int i;
		for (i = 0; i < oneline.size() && (oneline[i] == '\t' || oneline[i] == ' '); ++i);
		if (i < oneline.size() && oneline[i] != comment_symbol) 
		{
			stream.seekg(pos_b_r);
			return stream;
		}
		pos_b_r = stream.tellg();
	}
	return stream;
}



void skip_comment(istream& stream, const char comment_symbol)
{
	streamoff pos_b_r = stream.tellg();
	for (string oneline; getline(stream,oneline); )
	{
		unsigned int i=0;
		for (i=0; i<oneline.size() && (oneline[i] == '\t' || oneline[i] == ' '); ++i); //\r ki kell hagyni, különben windows-unix sortördelés miatt gondok adódhatnak
		if (i < oneline.size() && oneline[i] != comment_symbol) //ha az elsõ nem szóköz és tabulátor a sorban van, és nem a komment jel, akkor
		{
			stream.seekg(pos_b_r);
			return;
		}
		pos_b_r = stream.tellg();
	}
	return;
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(NULL);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

string ord_num(int a)
{
	string ret;
	if (a == 1)
		ret = "1st";
	else if (a == 2)
		ret = "2nd";
	else if (a == 3)
		ret = "3rd";
	else
		ret = to_string(a) + "th";

	return ret;
}