/*
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <assert.h>
#include <cassert>
*/

using namespace std;

typedef double Doub;
typedef long double Ldoub;
typedef int Int;
typedef vector<double> VecDoub;

#ifdef _MSC_VER /* Microsoft C++ */
  typedef unsigned __int64 Ullong; // 64-bit unsigned integer
#else /* Macintosh, Linux */
  typedef unsigned long long Ullong;

  typedef bool Bool;
  typedef int Int;
#endif

//global functions
//VecDoub logfact(1024,-1);
//VecDoub logfact(1048576,-1);	//1048576 = 1024*1024

struct Ran{//Strictly Ranq2 in Numerical Recipes
	Ullong v,w;
        //constructor
	Ran(Ullong j) : v(4101842887655102017LL), w(1){
	v ^= j;
	w = int64();
	v = int64();
  	}

  	inline Ullong int64(){
    		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    		w = 429457665U*(w & 0xffffffff) + (w >> 32);
    		return v ^ w;
  	}
  	inline double doub() {return 5.42101086242752217E-20 * int64();}
  	inline unsigned int int32() {return (unsigned int) int64();}
};

