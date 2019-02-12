#ifndef GLOBAL_H
#define GLOBAL_H

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
#include <random>
#include <algorithm>
#include <functional>
#include <iterator>
#include <sys/stat.h>

#include <errno.h>    // errno, ENOENT, EEXIST
#if defined(_WIN32)
#include <direct.h>   // _mkdir
#endif

#include "rng.h"

/*SodaPop
Copyright (C) 2019 Louis Gauthier

    SodaPop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    SodaPop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with SodaPop.  If not, see <http://www.gnu.org/licenses/>.
 */

/******** CONTAINERS AND ITERATORS ********/

typedef std::vector<int> VectInt;
typedef std::vector<int>::iterator VectInt_iterator;
typedef std::vector<int>::const_iterator VectInt_citerator;

typedef std::vector<double> VectDouble;
typedef std::vector<double>::iterator VectDouble_iterator;
typedef std::vector<double>::const_iterator VectDouble_citerator;

typedef std::vector<std::string> VectStr;
typedef std::vector<std::string>::iterator VectStr_iterator;
typedef std::vector<std::string>::const_iterator VectStr_citerator;

extern VectStr PrimordialAASeq;

template <typename T>
T remove_at(std::vector<T>&v, typename std::vector<T>::size_type n)
{
    T ans = std::move_if_noexcept(v[n]);
    v[n] = std::move_if_noexcept(v.back());
    v.pop_back();
    return ans;
}

/******* CONSTANTS *********/

/*****
N.B. The physically allowed value for mutational DDG is DGG_min to DGG_max.
If the estimated energy is out of this range, the mutation is ignored.
*****/

const int POPSIZEMAX(1000000);
const int GENECOUNTMAX(10);
const int PBWIDTH(70);

// for pretty printing of progress
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"

extern double const ddG_low_bound;
extern double const ddG_high_bound;
extern double const CONC_MAX ;
extern double const kT; //defines the energy units
extern double const MISFOLDING_COST; // misfolding cost, see Geiler-Samerotte et al. 2011
extern double const fNS; //fraction of non-synonymous substitutions in a typical protein
extern double const PREFACTOR; // prefactor for growth rate fitness function

// exponent values are precalculated to be used readily
extern double const DDG_min;
extern double const DDG_max;

extern int const Bigbuffer_max;
extern double const PI;

// If the mutation is to a stop codon
// DG_mutant is set to 99 kcal/mol 
// -> all copies are effectively aggregated
extern double const DG_STOP;

// Create a 3D matrix for fitness landscape
const int gene_number(100);
const int res_number(1200);

extern double matrix[gene_number][res_number][20];
extern double matrix_supp[gene_number][res_number][20];

extern double fold_DG;
extern double bind_DG;

enum Matrix_Type {
    is_folding,
    is_binding
};

/******* FUNCTION DECLARATIONS *******/
int GetIndexFromAA(char);
int GetIndexFromCodon(std::string);
std::string GetProtFromNuc(std::string);
std::string n3_to_n3(std::string, std::string, int);
std::string getBarcode();
//double RandomNumber();
double Ran_Gaussian(double const, double const);
const char AdjacentBP(char, int);
void InitMatrix();
double ExtractDDGMatrix(std::string,Matrix_Type);
void ExtractDMSMatrix(std::string);
int LoadPrimordialGenes(const std::string&, const std::string&);
int StringDiff(const std::string&, const std::string&);
std::string trim(const std::string&);
bool isDirExist(const std::string&);
bool makePath(const std::string&);
void printProgress (double);

void qread_Cell(std::fstream&, std::fstream&);
void seqread_Cell(std::fstream&, std::fstream&);
void read_Parent(std::fstream&, std::fstream&);
void read_Cell(std::fstream&, std::fstream&, bool);
#endif