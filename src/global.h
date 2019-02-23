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

const int maxPopSize(1000000);
const int maxGeneCount(100);
const int PBWidth(70);

const double ddG_low_bound(-10);
const double ddG_high_bound(99);
const double maxConcentration(1e15);
const double kT(0.5922); //defines the energy units
const double misfoldingCost(1e-4); // misfolding cost, see Geiler-Samerotte et al. 2011
const double fNS(0.775956284); //fraction of non-synonymous substitutions in a typical protein
const double prefactor(16000);

// exponent values are precalculated to be used readily
double const DDG_min = exp(-1*(ddG_low_bound)/kT);
double const DDG_max = exp(-1*(ddG_high_bound)/kT);
extern char buffer[];

// If the mutation is to a stop codon
// DG_mutant is set to 99 kcal/mol 
// -> all copies are effectively aggregated
double const DG_stop = exp(-1*(99)/kT);

// for pretty printing of progress
#define PBstr "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"

// Create a 3D matrix for fitness landscape
const int gene_number(100);
const int res_number(1200);

extern double matrix[gene_number][res_number][20];
extern double matrix_supp[gene_number][res_number][20];

extern double fold_DG;
extern double bind_DG;

extern int Total_Cell_Count;
extern int dummy;
extern double frame_time;
extern std::string outPath;

enum Matrix_Type {
    is_folding,
    is_binding
};

enum Encoding_Type {
    full,
    by_default,
    no_sequence,
    other
};

enum Init_Pop {
    from_snapFile,
    from_cellFile
};

/******* FUNCTION DECLARATIONS *******/
Encoding_Type intToEncoding_Type(int);
Init_Pop intToPop_Type(int);

int GetIndexFromAA(char);
int GetIndexFromCodon(std::string);
std::string GetProtFromNuc(std::string);
std::string n3_to_n3(std::string, std::string, int);

std::string getBarcode();

//double RandomNumber();

double Ran_Gaussian(double const, double const);

const char AdjacentBP(char, int);

void openStartingPop(std::string, std::ifstream&);
void readSnapshotHeader(std::ifstream&);
void createOutputDir(std::string);

void openCommandLog(std::ofstream&, std::string, char *[], int);
void openMutationLog(std::ofstream&, std::string);

void InitMatrix();
double ExtractDDGMatrix(std::string,Matrix_Type);
void ExtractDMSMatrix(std::string);
int LoadPrimordialGenes(const std::string&, const std::string&);
int StringDiff(const std::string&, const std::string&);
std::string trim(const std::string&);
bool isDirExist(const std::string&);
bool makePath(const std::string&);
void printProgress (double);

void qread_Cell(std::ifstream&, std::ofstream&);
void seqread_Cell(std::ifstream&, std::ofstream&);
void read_Parent(std::ifstream&, std::ofstream&);
void read_Cell(std::ifstream&, std::ofstream&, bool);
#endif