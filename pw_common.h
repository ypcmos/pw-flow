/**
*	pw_common.h
*	About some based function and definition and containing common head files
*	Copyright(c) Peng Yao
*	Date: 2012.12
*	All rights reserved
*/
#pragma once

#ifndef _PW_COMMON_H
#define _PW_COMMON_H
#endif

#include <list>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stack>
using namespace std;

typedef stack<string> StrStack;
typedef vector<string> StrArray;
typedef vector<int> IntArray;
typedef vector<int>::iterator PInt;
typedef vector<int>::const_iterator CPInt;
typedef vector<string>::iterator PString;
typedef vector<string>::const_iterator CPString;
typedef string::iterator PChar;

#define _ZERO 1.0e-10 
#define PI 3.14159265357

#ifndef BYTE
#define BYTE unsigned char
#endif

#ifndef UINT
#define UINT unsigned
#endif

class _Exception
{
public:
	_Exception(string err);
	string Error() const;

private:
	string Message;
};

class ExString
{
public:
	static string DelSpace(string str);
	//Each string in strc can not be space!
	static StrArray Splite(string str, const StrArray& strc);
	//Each string in strc can not be space!
	static StrArray SpliteElement(string str, const StrArray& strc);
	static string DelExcept(string str, string exp);
	static string Trim(string str);
	static string ToString(const StrArray& strs);
	static void Merger(StrArray& des,const StrArray& src);
	static void UniqueMerger(StrArray& des,const StrArray& src);
	static StrArray ToStringArray(string src);
};

