#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <QDebug>

#include "event.h"

using namespace std;

// 0x09 - horizontal tab
// 0x0a - linefeed
// 0x0b - vertical tab
// 0x0c - form feed
// 0x0d - carriage return
// 0x20 - space
// \0   - end of string
char *as_trim(char *str) {
	char *p = &str[0];
	int run = 1;
	while(run) {
		run = 0;
		switch(p[0]) {
			case 0x09: case 0x0a: case 0x0b: case 0x0c: case 0x0d: case 0x20:
				run = 1;
				p++;
			break;
		}
	}
		return p;
}

// return true if line starts with '#' or '$'
int as_is_comment(char *str) {
	char *p = as_trim(str);
	return ((p[0] == '#')||(p[0] == '$'));
}

// checks if line only constists of spaces
int as_is_empty(char *str) {
	char *p = as_trim(str);
	return (strlen(p)==0);
}

void as_replace_detector(char *dest, const char *src, int det) {
	// replace "$d"
	string str(src);
	size_t pos = 0;
	char c_detstr[200];
	sprintf(c_detstr, "%d", det);
	while((pos=str.find("$", pos))!=string::npos) {
		unsigned int digits = 0;
		size_t startpos = pos;
		while(src[pos]=='$') {
			if(pos==strlen(src)) {
                qDebug() << "no end of $ expression in string " << src;
				exit(1);
			}
			digits++;
			pos++;
		}
		if(src[pos]!='d') {
            qDebug() << "d after $ expression not found in string " << src;
				exit(1);
		}
		string detstr(c_detstr);
		while(detstr.length()<digits) {
			detstr = "0" + detstr;
		}
		str.replace(startpos, digits+1, detstr, 0, detstr.length());
		//str.replace(pos, strlen(detstr), detstr);		
	}
	strcpy(dest, str.c_str());
}

void as_replace_wildcards(char *dest, const char *src, int run, int det) {
	// replace "$d"
	// TODO avoid conversion between cstring and c++string
	// replace "$d"
	char buf[2000];
	as_replace_detector(buf, src, det);
	string str(buf);


	// replace "#..#"
	char runstr[200];
	sprintf(runstr, "%d", run);
	size_t pos = 0;
	while((pos=str.find("#", pos))!=string::npos) {
		// determine number # in a row
		int length = 0;
		while((pos+length<str.length()) && str.at(pos+length)=='#') {
			length++;
		}

		// write leading zeros
		string buf(runstr);
		for(int i=0; i<length - (int)strlen(runstr); i++) {
			buf = "0" + buf;
		}

		str.replace(pos, length, buf);
	}
	
	strcpy(dest, str.c_str());
	
}
