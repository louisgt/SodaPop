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
#include "pcg_random.hpp"
#include <sys/stat.h>
#include <random/include_headers.h>
#include <random/gamma.h>
#include <random/deviates.h>
#include <errno.h>    // errno, ENOENT, EEXIST
#if defined(_WIN32)
#include <direct.h>   // _mkdir
#endif

/*AUTHOR: LOUIS GAUTHIER*/

#define POPSIZEMAX 	1000000
#define GENECOUNTMAX 	100

template<class T = std::mt19937, std::size_t N = T::state_size>
auto ProperlySeededRandomEngine () -> typename std::enable_if<!!N, T>::type {
    typename T::result_type random_data[N];
    std::random_device source;
    std::generate(std::begin(random_data), std::end(random_data), std::ref(source));
    std::seed_seq seeds(std::begin(random_data), std::end(random_data));
    T seededEngine (seeds);
    return seededEngine;
}

/******* RANDOM GENERATORS *********/

//Global declarations of deviate generators
Ran *uniformdevptr;
Normaldev *normaldevptr;
Poissondev *poissondevptr;
Binomialdev *binomialdevptr;

/******** PSEUDO RANDOM NUMBER GENERATOR

adapted from      PCG: A Family of Simple Fast Space-Efficient Statistically Good
                  Algorithms for Random Number Generation
                  by MELISSA E. O’NEILL

this RNG provides high quality uniform random floats
without a big memory or time penalty
*/

//Set seed source to rand device
pcg_extras::seed_seq_from<std::random_device> seed_source;
//Create rng
pcg32 rng(seed_source);
//use uniform real dist
std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

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

VectStr PrimordialAASeq;

/******* CONSTANTS *********/

/*****
N.B. The physically allowed value for mutational DDG is DGG_min to DGG_max.
If the estimated energy is out of this range, the mutation is ignored.
*****/
const double ddG_min = -7.5;
const double ddG_max = 99;
const double CONC_MAX = 1e15;
const double kT = 0.5922; //defines the energy units
const double COST = 1e-4; // see Geiler-Samerotte et al. 2011
const double A_FACTOR = 0.005754898; // see Serohijos & Shakhnovich 2013


// exponent values are precalculated to be used readily
const double DDG_min = exp(-1*(ddG_min)/kT);
const double DDG_max = exp(-1*(ddG_max)/kT);
const int Bigbuffer_max = 80;
constexpr double PI  = 3.141592653589793238463;


// If the mutation is to a stop codon
// DG_mutant is set to 99 kcal/mol 
// -> all copies are effectively aggregated
const double DG_STOP = exp(-1*(99)/kT);

// Create a 3D matrix for pDDG values
const int max_gene = 1200;
const int max_resi = 640;
double pDDG[max_gene][max_resi][20];

/******* FUNCTION DECLARATIONS *******/
int GetIndexFromAA(std::string);
int GetIndexFromAA(char);
int GetIndexFromCodon(std::string);
std::string GetProtFromNuc(std::string);
int CheckBP(std::string);
std::string n3_to_n3(std::string, int);
std::string getBarcode();
double RandomNumber();
void InitDDGMatrix();
void ExtractPDDGMatrix(std::string);
void LoadPrimordialGenes(const std::string&,const std::string&);
int StringDiff(const std::string&, const std::string&);
std::string trim(const std::string&);
bool isDirExist(const std::string&);
bool makePath(const std::string&);

/******* GENETIC CODE MAPPINGS *******/
// these const mappings are hard-coded and populated at compile-time

struct codon_to_num{
    static std::map<std::string,int> create_map()
        {
          std::map<std::string,int> m;
          m["GCA"] = 1;
          m["GCC"] = 1;
          m["GCG"] = 1;
          m["GCT"] = 1;
          m["TGC"] = 2;
          m["TGT"] = 2;
          m["GAC"] = 3;
          m["GAT"] = 3;
          m["GAA"] = 4;
          m["GAG"] = 4;
          m["TTC"] = 5;
          m["TTT"] = 5;
          m["GGA"] = 6;
          m["GGC"] = 6;
          m["GGG"] = 6;
          m["GGT"] = 6;
          m["CAC"] = 7;
          m["CAT"] = 7;
          m["ATA"] = 8;
          m["ATC"] = 8;
          m["ATT"] = 8;
          m["AAA"] = 9;
          m["AAG"] = 9;
          m["CTA"] = 10;
          m["CTC"] = 10;
          m["CTG"] = 10;
          m["CTT"] = 10;
          m["TTA"] = 10;
          m["TTG"] = 10;
          m["ATG"] = 11;
          m["AAC"] = 12;
          m["AAT"] = 12;
          m["CCA"] = 13;
          m["CCC"] = 13;
          m["CCG"] = 13;
          m["CCT"] = 13;
          m["CAA"] = 14;
          m["CAG"] = 14;
          m["AGA"] = 15;
          m["AGG"] = 15;
          m["CGA"] = 15;
          m["CGC"] = 15;
          m["CGG"] = 15;
          m["CGT"] = 15;
          m["AGC"] = 16;
          m["AGT"] = 16;
          m["TCA"] = 16;
          m["TCC"] = 16;
          m["TCG"] = 16;
          m["TCT"] = 16;
          m["ACA"] = 17;
          m["ACC"] = 17;
          m["ACG"] = 17;
          m["ACT"] = 17;
          m["GTA"] = 18;
          m["GTC"] = 18;
          m["GTG"] = 18;
          m["GTT"] = 18;
          m["TGG"] = 19;
          m["TAC"] = 20;
          m["TAT"] = 20;
          m["TAA"] = 21;
          m["TAG"] = 21;
          m["TGA"] = 21;
          return m;
        }
    static const std::map<std::string,int> cnum;
};

struct codon_to_prot{
    static std::map<std::string,std::string> create_map()
        {
          std::map<std::string,std::string> m;
          m["AAA"] = "K";
          m["AAC"] = "N";
          m["AAG"] = "K";
          m["AAT"] = "N";
          m["ACA"] = "T";
          m["ACC"] = "T";
          m["ACG"] = "T";
          m["ACT"] = "T";
          m["AGA"] = "R";
          m["AGC"] = "S";
          m["AGG"] = "R";
          m["AGT"] = "S";
          m["ATA"] = "I";
          m["ATC"] = "I";
          m["ATG"] = "M";
          m["ATT"] = "I";
          m["CAA"] = "Q";
          m["CAC"] = "H";
          m["CAG"] = "Q";
          m["CAT"] = "H";
          m["CCA"] = "P";
          m["CCC"] = "P";
          m["CCG"] = "P";
          m["CCT"] = "P";
          m["CGA"] = "R";
          m["CGC"] = "R";
          m["CGG"] = "R";
          m["CGT"] = "R";
          m["CTA"] = "L";
          m["CTC"] = "L";
          m["CTG"] = "L";
          m["CTT"] = "L";
          m["GAA"] = "E";
          m["GAC"] = "D";
          m["GAG"] = "E";
          m["GAT"] = "D";
          m["GCA"] = "A";
          m["GCC"] = "A";
          m["GCG"] = "A";
          m["GCT"] = "A";
          m["GGA"] = "G";
          m["GGC"] = "G";
          m["GGG"] = "G";
          m["GGT"] = "G";
          m["GTA"] = "V";
          m["GTC"] = "V";
          m["GTG"] = "V";
          m["GTT"] = "V";
          m["TAA"] = "X";
          m["TAC"] = "Y";
          m["TAG"] = "X";
          m["TAT"] = "Y";
          m["TCA"] = "S";
          m["TCC"] = "S";
          m["TCG"] = "S";
          m["TCT"] = "S";
          m["TGA"] = "X";
          m["TGC"] = "C";
          m["TGG"] = "W";
          m["TGT"] = "C";
          m["TTA"] = "L";
          m["TTC"] = "F";
          m["TTG"] = "L";
          m["TTT"] = "F";
          return m;
        }
    static const std::map<std::string,std::string> cprot;
};

struct prot_to_num{
    static std::map<std::string,int> create_map()
        {
          std::map<std::string,int> m;
          m["A"] = 1;
          m["C"] = 2;
          m["D"] = 3;
          m["E"] = 4;
          m["F"] = 5;
          m["G"] = 6;
          m["H"] = 7;
          m["I"] = 8;
          m["K"] = 9;
          m["L"] = 10;
          m["M"] = 11;
          m["N"] = 12;
          m["P"] = 13;
          m["Q"] = 14;
          m["R"] = 15;
          m["S"] = 16;
          m["T"] = 17;
          m["V"] = 18;
          m["W"] = 19;
          m["Y"] = 20;
          m["X"] = 21;//STOP
          return m;
        }
    static const std::map<std::string,int> pnum;
};

// populate genetic code mappings
const std::map<std::string,int> codon_to_num::cnum = codon_to_num::create_map();
const std::map<std::string,int> prot_to_num::pnum = prot_to_num::create_map();
const std::map<std::string,std::string> codon_to_prot::cprot = codon_to_prot::create_map();

/******* MAPPING FUNCTIONS *******/

// input:  3 nucleotide codon as string
// output: index number of codon
int GetIndexFromCodon(std::string in_codon)
{
    std::map <std::string, int>::const_iterator it;
    it = codon_to_num::cnum.find(in_codon);  
    if(it == codon_to_num::cnum.end()){
        std::cerr << "Invalid codon: "<< in_codon << std::endl;
        exit(2);
    }
    return it->second;
}

// input:  3 nucleotide codon sequence as string
// output: amino acid sequence as string
std::string GetProtFromNuc(std::string in_seq)
{
    int ln = in_seq.length();
    if((ln % 3) != 0)
    {
        std::cerr << "Invalid length for nucleotide sequence: " << ln << std::endl;
        exit(2);
    }
    int la=ln/3;   
    std::string AA="";
    for(int i=0; i<la;i++){
        std::string temp=in_seq.substr(i*3,3);
        
        //check for valid code
        std::map <std::string, std::string> :: const_iterator Iter;
        Iter = codon_to_prot::cprot.find(temp);
        if (Iter == codon_to_prot::cprot.end()){
          std::cerr << "Invalid codon: "<< temp << std::endl;
          exit(2);
        }   
        AA.append(codon_to_prot::cprot.at(temp));
    }
    return AA;
}

// input:  amino acid letter as string
// output: index number of amino acid
int GetIndexFromAA(std::string aa)
{  
    //check for valid amino acid
    std::map <std::string, int> :: const_iterator it;
    it = prot_to_num::pnum.find(aa);
    if( it == prot_to_num::pnum.end()){
        std::cerr << "Invalid amino acid: "<< aa << std::endl;
        exit(2);
    }
    return it->second;
}

//duplicate function, allows for a single char input
int GetIndexFromAA(char aa){
    std::map <char, int> m;
    m['A'] = 1;
    m['C'] = 2;
    m['D'] = 3;
    m['E'] = 4;
    m['F'] = 5;
    m['G'] = 6;
    m['H'] = 7;
    m['I'] = 8;
    m['K'] = 9;
    m['L'] = 10;
    m['M'] = 11;
    m['N'] = 12;
    m['P'] = 13;
    m['Q'] = 14;
    m['R'] = 15;
    m['S'] = 16;
    m['T'] = 17;
    m['V'] = 18;
    m['W'] = 19;
    m['Y'] = 20;
    m['X'] = 21;//default for missing amino acids

    //check for valid amino acid
    std::map <char, int> :: const_iterator Iter;
    Iter = m.find(aa);
    if ( Iter == m.end( ) ){
      std::cerr << "Invalid amino acid: "<< aa << std::endl;
      exit(2);
    }
    return m[aa];
}

// verifies validity of nucleotide
// input: a single nucleotide
// output: nucleotide number from index
int CheckBP(std::string a){
    std::map <std::string, int> A;
    A["A"] = 1;
    A["T"] = 2;
    A["G"] = 3;
    A["C"] = 4;

    std::map <std::string, int> :: const_iterator Iter;
    Iter = A.find(a);
    if (Iter == A.end()){
        std::cerr << "Invalid nucleotide: " << a << std::endl;
        exit(2);
    }
    return A[a];
}

// returns the next or second to next bp from a given nucleotide
std::string AdjacentBP(std::string a, int j){
 
    if ( j > 2 ){
        std::cerr << "Invalid bp distance. Error in AdjacentBP(). "<< j << std::endl;
        exit(2);
    }

    std::map <std::string, int> A;
    A["A"] = 0;
    A["T"] = 1;
    A["G"] = 2;
    A["C"] = 3;

    std::map <int, std::string> B;
    B[0] = "A";
    B[1] = "T";
    B[2] = "G";
    B[3] = "C";

    int x = (A[a] + j + 1) % 4; 	//get (j+1) nucleotide from a
    return B[x];
}

//THIS FUNCTION NEEDS AN EXTREME MAKEOVER

//Checks if codon b when mutated resulted in a STOP codon a.
//If b is a STOP codon, it is replaced by 
//Input: 	a, new codon
//		b, old codon
//		i, mutation site in codon (<3)
std::string n3_to_n3(std::string a, std::string b, int i){
  double r = RandomNumber();
  double l;

  if ( a == "TAA")
  {
      switch (i)
      {
          case 0:
              l =  (2*r);// std::cout << l << std::endl;

              if (b == "CAA" ){
                if  ( l<1 )	a = "GAA";
                else 		a = "AAA";
              }else if (b == "GAA" ){
                if  ( l<1 )	a = "CAA";
                else 		a = "AAA";
              }else if (b == "AAA" ){
                 if  ( l<1 )	a = "GAA";
                else 		a = "CAA";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

          case 1:
              if 	(b == "TTA") a = "TCA";
              else if (b == "TCA") a = "TTA";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

            break;

          case 2:
              if 	(b == "TAT") a = "TAC";
              else if (b == "TAC") a = "TAT";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

            break;

           default:
            std::cerr << "ERROR in translating STOP codon to NON-STOP codon." << std::endl;    
      }
  }
  else if (a == "TAG")
  {
      switch (i)
      {
          case 0:
              l = (2*r);

              if (b == "AAG" ){
                if  ( l<1 )	a = "GAG";
                else 		a = "CAG";
              }else if (b == "GAG" ){
                if  ( l<1 )	a = "AAG";
                else 		a = "CAG";
              }else if (b == "CAG" ){
                 if  ( l<1 )	a = "GAG";
                else 		a = "AAG";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

          case 1:
              l = (2*r);

              if (b == "TTG" ){
                if  ( l<1 )	a = "TGG";
                else 		a = "TGC";
              }else if (b == "TGG" ){
                if  ( l<1 )	a = "TTG";
                else 		a = "TGC";
              }else if (b == "TCG" ){
                 if  ( l<1 )	a = "TTG";
                else 		a = "TGG";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

          case 2:
              if 	(b == "TAT") a = "TAC";
              else if (b == "TAC") a = "TAT";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

              break;

             default:
             std::cerr << "ERROR in translating STOP codon to NON-STOP codon." << std::endl;    
       }
   }
   else if (a == "TGA")
   {
      switch (i)
      {
          case 0:
              l =  (2*r);

              if (b == "AGA" ){
                if  ( l<1 )	a = "GGA";
                else 		a = "CGA";
              }else if (b == "GGA" ){
                if  ( l<1 )	a = "AGA";
                else 		a = "CGA";
              }else if (b == "CGA" ){
                 if  ( l<1 )	a = "AGA";
                else 		a = "GGA";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

         case 1:
              if 	(b == "TTA") a = "TCA";
              else if (b == "TCA") a = "TTA";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

              break;

          case 2:
              l =  (2*r);

              if (b == "TGT" ){
                if  ( l<1 )	a = "TGG";
                else 		a = "TGC";
              }else if (b == "TGG" ){
                if  ( l<1 )	a = "TGT";
                else 		a = "TGC";
              }else if (b == "TGC" ){
                 if  ( l<1 )	a = "TGT";
                else 		a = "TGG";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

           default:
            std::cerr << "ERROR in translating STOP codon to NON-STOP codon." << std::endl;    
      }  
  }
  return a;
}

// generates a random, 15nt barcode
std::string getBarcode()
{
    char seq [16];
    for(int i = 0; i < 15; i++){
        if(RandomNumber()<0.5){
            if(RandomNumber()<0.5){
                seq[i] = 'G';
            }
            else seq[i] = 'C';
        }
        else{
            if(RandomNumber()<0.5){
                seq[i] = 'A';
            }
            else seq[i] = 'T';
        }
    }
    seq[15] = '\0';
    return std::string(seq);
}

// initializes the 3D matrix for DDG values
void InitDDGMatrix()
{
    for(int i = 0; i != max_gene; ++i)
      for(int j = 0; j != max_resi; ++j)
        for(int k = 0; k != 20; ++k)
          pDDG[i][j][k] = 1; //i.e., exp(-0/kT) = 1
}

// extracts values from the DDG file and stores them in the matrix
void ExtractPDDGMatrix(std::string filepath)
{
    std::fstream temp(filepath);
    if(!temp.is_open()){
        std::cerr << "File could not be open: "<< filepath << std::endl;
        exit(2);
    }
    std::string line;
    int gene_num = 0;
    while(!temp.eof()){
        std::string word;
        getline(temp,line);
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;
        if ( word == "DDG"){
            iss >> word;
            //residue index
            int i = atoi(word.c_str());
            for(int j = 0; iss>>word; j++){
                //extract DDG values
                double x = atof(word.c_str());
                pDDG[gene_num][i-1][j] = exp(-x/kT);
            }
        }
        else if ( word == "Gene_NUM"){
            iss >> word;
            gene_num = atoi(word.c_str());
        }
    }
    temp.close();
}

// returns uniformly distributed floating-point random number on interval [0.0,1.0]
double RandomNumber()
{
    return uniform_dist(rng);
}

double Ran_Gaussian(const double mean, const double sigma)
{
    double x, y, r2;
    do{
        // choose x,y in uniform square [-1,+1]
        x = -1 + 2 * RandomNumber();
        y = -1 + 2 * RandomNumber();
        // check if it is in the unit circle
        r2 = x * x + y * y;
    }while (r2 > 1.0 || r2 == 0); 
    // Box-Muller transform
    return mean + sigma * y * sqrt (-2.0 * log (r2) / r2);
}

// Loads primordial genes in a VectStr
void LoadPrimordialGenes(const std::string& genelistfile, const std::string& genesPath)
{  
    std::fstream genelistIN (genelistfile.c_str());
    if (!genelistIN.is_open()){
        std::cerr << "File could not be open: "<< genelistfile <<std::endl;
        exit(2);
    }
    int flag_AASeq = 0; 
    int gc = 0;
    while(!genelistIN.eof()){
        std::string word, line;
        getline(genelistIN,line);
        std::istringstream iss(line, std::istringstream::in);
        iss >> word; 
        if ( word=="G" ){
            iss >> word;
            word = genesPath + word;
            std::fstream genefileIN (word.c_str(), std::fstream::in | std::fstream::out);
            if (!genefileIN.is_open()){
                std::cerr << "File could not be open: "<< word <<std::endl;
                exit(2);
            }
            assert ( gc > 0);

            int gn = -1;
            while(!genefileIN.eof()){
                  std::string w, l;
                  getline(genefileIN,l);
                  std::istringstream iss(l, std::istringstream::in);
                  iss >> w;
                  if(w=="Gene_NUM"){
                      iss >> w;
                      gn = atoi(w.c_str());
                      assert( gn == flag_AASeq); // Gene list must be ordered from 0..(N-1);
                  }
                  else if ( w=="N_Seq" ){
                      assert( gn >= 0);
                      iss >> w;
                      std::string aaseq=GetProtFromNuc(w);
                      //check stop codons in midsequence
                      std::string::size_type loc = aaseq.find("X", 0);
                      assert( loc == std::string::npos ); // no match
                      VectStr_iterator iter = PrimordialAASeq.begin();
                      PrimordialAASeq.insert(iter+gn, aaseq); 
                      flag_AASeq += 1;
                  }
            }
        }
        else if ( word=="Gene_Count" ){
            iss >> word;
            gc = atoi(word.c_str());
            PrimordialAASeq.reserve(gc);
        }
    }
    assert( flag_AASeq==gc );
}

// Reads a unit cell stored in binary format using Cell::dumpShort()
void qread_Cell(std::fstream& IN, std::fstream& OUT)
{
    char buffer[140];
    int na, ns;
    double f;
    std::string barcode;

    int l;
    IN.read((char*)&l, sizeof(int));
    std::vector<char> buf(l);
    IN.read(&buf[0], l);
    barcode.assign(buf.begin(), buf.end());

    std::cout << barcode << std::endl;

    OUT << barcode;

    IN.read((char*)(&na),sizeof(int));
    IN.read((char*)(&ns),sizeof(int));
    IN.read((char*)(&f),sizeof(double));
    
    sprintf(buffer," C %6d %6d %12e", na, ns, f);
    OUT << buffer << std::endl;
}

// Reads a unit cell stored in binary format using Cell::dump()
void read_Cell(std::fstream& IN, std::fstream& OUT)
{
    char buffer[140];
    int cell_id, cell_index, gene_size;
    double m,m0;
    std::string barcode;

    IN.read((char*)(&cell_index),sizeof(int));
    IN.read((char*)(&cell_id),sizeof(int));

    int l;
    IN.read((char*)&l, sizeof(int));
    std::vector<char> buf(l);
    IN.read(&buf[0], l);  
    barcode.assign(buf.begin(), buf.end());

    OUT << barcode;

    IN.read((char*)(&m0),sizeof(double));
    IN.read((char*)(&m),sizeof(double));
    IN.read((char*)(&gene_size),sizeof(int));
    
    sprintf(buffer," C %6d %12e %12e", cell_index, m0, m);
    OUT << buffer << std::endl;

    for(int j=0; j<gene_size; j++){
        double e, c, dg;
        int gene_nid, Ns, Na;
        std::string DNAsequence;

        IN.read((char*)(&gene_nid),sizeof(int));   
        IN.read((char*)(&e),sizeof(double));
        IN.read((char*)(&c),sizeof(double));
        IN.read((char*)(&dg),sizeof(double));

        IN.read((char*)(&Ns),sizeof(int));
        IN.read((char*)(&Na),sizeof(int));
        
        //read DNA sequence
        int nl;
        IN.read((char*)&nl, sizeof(int));
        std::vector<char> buff(nl);
        IN.read(&buff[0], nl);  
        DNAsequence.assign(buff.begin(), buff.end());

        sprintf(buffer,"G %12e %12e %6d %6d ", c, dg, Na, Ns);
        OUT << buffer << std::endl;
    } 
}

// Counts the number of dissimilar characters in 2 strings
int StringDiff(const std::string& A, const std::string& B)
{
    unsigned int L = A.length();
    assert(L == B.length());
    int ctr = 0;  
    for(unsigned int i =0; i<L; i++){
        if( A.at(i) != B.at(i)) ctr+=1; 
    } 
    return ctr;
}

std::string trim(const std::string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first){
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

// checks if given directory exists
bool isDirExist(const std::string& path)
{
#if defined(_WIN32)
    struct _stat info;
    if (_stat(path.c_str(), &info) != 0){
        return false;
    }
    return (info.st_mode & _S_IFDIR) != 0;
#else 
    struct stat info;
    if (stat(path.c_str(), &info) != 0){
        return false;
    }
    return (info.st_mode & S_IFDIR) != 0;
#endif
}

// creates directory from specified path
bool makePath(const std::string& path)
{
#if defined(_WIN32)
    int ret = _mkdir(path.c_str());
#else
    mode_t mode = 0755;
    int ret = mkdir(path.c_str(), mode);
#endif
    if (ret == 0)
        return true;

    switch (errno){
    case ENOENT:
        // parent didn't exist, try to create it
        {
            std::size_t pos = path.find_last_of('/');
            if (pos == std::string::npos)
#if defined(_WIN32)
                pos = path.find_last_of('\\');
            if (pos == std::string::npos)
#endif
                return false;
            if (!makePath( path.substr(0, pos) ))
                return false;
        }
        // try to create again
#if defined(_WIN32)
        return 0 == _mkdir(path.c_str());
#else 
        return 0 == mkdir(path.c_str(), mode);
#endif
    case EEXIST:
        // done
        return isDirExist(path);

    default:
        return false;
    }
}




