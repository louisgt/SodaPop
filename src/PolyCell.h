#ifndef POLYCELL_H
#define POLYCELL_H
#include "Cell.h"
#include "rng.h"

/*SodaPop
Copyright (C) 2018 Louis Gauthier

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

class PolyCell: public Cell
{
    typedef double(PolyCell::*funcPtr)(void);

public:
    static int ff_;
    static bool useDist_;
    static bool fromS_;

    PolyCell();
    PolyCell(std::fstream&);                
    PolyCell(std::fstream&, const std::string&);

    void FillGene_L();

    // Fitness functions
    void selectFitness();
    double flux();
    double toxicity();
    double metabolicOutput();
    double multiplicative();
    double neutral();
    double noMut();
    double fold();
    void UpdateRates();

    void ranmut_Gene();
    void ranmut_Gene(std::ofstream&, int);
    void change_exprlevel();
    double normalizeFit(double);
    void dump(std::fstream&, int);
    void dumpShort(std::fstream&);
    void dumpSeq(std::fstream&, int);
    void dumpParent(std::fstream&);
    void PrintCell(int);
    
    int Na(){return Total_Na_;}
    int Ns(){return Total_Ns_;}
    void UpdateNsNa();

protected:
    int Total_Ns_;
    int Total_Na_;
    // function pointer to select fitness function
    funcPtr fit;
};

// By default the fitness function is set to neutral
int PolyCell::ff_ = 5;
bool PolyCell::useDist_ = false;
bool PolyCell::fromS_ = false;

PolyCell::PolyCell(){}
PolyCell::PolyCell(std::fstream& f) : Cell(f)
{
    selectFitness();
	// Update current rates
  	this->UpdateRates();  
  	// Fill gene length array
  	this->FillGene_L();
}    
PolyCell::PolyCell(std::fstream& f, const std::string& s) : Cell(f,s)
{
    selectFitness();
	// Update current rates
  	this->UpdateRates();  
  	// Fill gene length array
  	this->FillGene_L();
}

// Initialize the cumulative gene length array
void PolyCell::FillGene_L()
{
    int sum = 0;
    std::vector<Gene>::iterator i;
    for(i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        sum+= (*i).length();
        Gene_L_.push_back(sum);
    }
}

void PolyCell::selectFitness()
{
    switch(PolyCell::ff_){
        case 1: fit = &PolyCell::fold;
            break;
        case 2: fit = &PolyCell::flux;
            break;
        case 3: fit = &PolyCell::toxicity;
            break;
        case 4: fit = &PolyCell::metabolicOutput;
            break;
        case 5: fit = &PolyCell::multiplicative;
            break;
        case 6: fit = &PolyCell::neutral;
            break;
        case 7: fit = &PolyCell::noMut;
            break;
        default:;
    }
}

// FOLDING-STABILITY BASED FLUX FITNESS FUNCTION
double PolyCell::fold()
{
    double sum_func = 0;
    //sum (concentration*Pnat) over all genes
    for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        sum_func += gene_it->functional();
    }
    return sum_func;
}

// METABOLIC FLUX FITNESS FUNCTION
double PolyCell::flux()
{
    double sum_func = 0;
    double a = 0;
    for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        // efficiency * Pnat * abundance
        sum_func += 1/(gene_it->eff()*gene_it->functional());
        a += gene_it->A_factor();
    }
    return static_cast<double>(a)/sum_func;
}

// MISFOLDING TOXICITY FITNESS FUNCTION
double PolyCell::toxicity()
{
    double f = 0;
    //sum (concentration*(1-Pnat)) over all genes
    for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        f += gene_it->misfolded();
    }
    return exp(-(COST*f));
}

// COMBINED METABOLIC OUTPUT FITNESS FUNCTION
double PolyCell::metabolicOutput()
{
    double flux = 0;
    double toxicity = 0;
    double a = 0;
    for (auto& it : Gene_arr_) {
        flux += 1 / (it.eff()*it.functional());
        toxicity += it.misfolded();
        a += it.A_factor();
    }

    flux = a/flux;
    toxicity = COST * toxicity;

    double fitness = flux - toxicity;
    // if toxicity > flux, return 0
    return (fitness < 0) ? 0 : fitness;
}

// MULTIPLICATIVE FITNESS FUNCTION
double PolyCell::multiplicative()
{
    double fitness = 0;
    for (auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it)
    {
        fitness += gene_it->f()*gene_it->e();
    }
    return fitness/gene_count();
}

// NEUTRAL FITNESS FUNCTION
double PolyCell::neutral()
{
    return 1;
}

// NEUTRAL FITNESS FUNCTION
double PolyCell::noMut()
{
    return fitness();
}

void PolyCell::UpdateRates()
{
    fitness_ = (this->*fit)();
}

void PolyCell::ranmut_Gene(std::ofstream& log,int ctr)
{
    // get genome size
    int L = Gene_L_.back();

    // pick random site to mutate

    int site = (int) ( L * randomNumber());

    // find the corresponding gene
    std::vector<Gene>::iterator j = Gene_arr_.begin();
    VectInt_citerator k = Gene_L_.begin();

    if(site >= (*k)){
    // random number generated is greater than
    // the cumulative sum of genes
         for(k = Gene_L_.begin(); k != Gene_L_.end(); ++k){
             if( site< (*k) ) break;
             j++; 
         }        
         k--;
         site = site - (*k);        
    }

    std::string mutation = "";

    int bp = (int) (3 * randomNumber());

    double wi = fitness();
    if(fromS_)
    {
        if(useDist_)
        {
            (*j).Mutate_Select_Dist(site,bp);
        }
        else
        {
            mutation = (*j).Mutate_Select(site,bp);
        }
    }
    else
    {
        if(useDist_)
        {
            (*j).Mutate_Stabil_Gaussian(site,bp);
        }
        else
        {
            mutation = (*j).Mutate_Stabil(site,bp);
        }
    }

    UpdateRates();
    double wf = fitness();
    double s = wf - wi;

    // save beneficial mutations to log
    // we could save all mutations with abs(s) >= some value x
    log << barcode().c_str() << "\t";
    log << std::fixed;
    log << mutation << "\t";
    log << s << "\t";
    log << ctr << std::endl;
}

void PolyCell::ranmut_Gene()
{
    // get genome size
    int L = Gene_L_.back();
    // pick random site to mutate

    int site = (int) ( L * randomNumber());

    // find the corresponding gene
    std::vector<Gene>::iterator j = Gene_arr_.begin();
    VectInt_citerator k = Gene_L_.begin();

    if(site >= (*k)){
    // random number generated is greater than
    // the cumulative sum of genes
         for(k = Gene_L_.begin(); k != Gene_L_.end(); ++k){
             if(site< (*k) ) break;
             j++; 
         }        
         k--;
         site = site - (*k);        
    }

    int bp = (int) (3 * randomNumber());
    // what is the input type?
    if(fromS_)
    {
        if(useDist_)
        {
            (*j).Mutate_Select_Dist(site,bp);
        }
        else
        {
            (*j).Mutate_Select(site,bp);
        }
    }
    else
    {
        if(useDist_)
        {
            (*j).Mutate_Stabil_Gaussian(site,bp);
        }
        else
        {
            (*j).Mutate_Stabil(site,bp);
        }
    }
         
    UpdateRates();
}

double PolyCell::normalizeFit(double fittest){
    if(fittest <= 0) {
        std::cerr << "Population collapse, average fitness is null.\n";
        exit(1);
    }
    double newfit = (Gene_arr_.begin()->f())/fittest;
    Gene_arr_.begin()->ch_f(newfit);
    UpdateRates();
    return newfit;
}

// Dump cell information to binary file
void PolyCell::dump(std::fstream& OUT, int cell_index)
{
    int x;
    double y;
    
    OUT.write((char*)(&cell_index),sizeof(int));
    
    //cell ID
    OUT.write((char*)(&ID_),sizeof(int));

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    y = fitness();
    OUT.write((char*)(&y),sizeof(double));
    
    y = c_mrate_;		 	 
    OUT.write((char*)(&y),sizeof(double));

    x = (int)(Gene_arr_.size());		 	 
    OUT.write((char*)(&x),sizeof(int));

   for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        int gene_nid = gene_it->num();
        double s = gene_it->e();
        double c = gene_it->conc();
        double eff = gene_it->eff();
        double dg = -kT*log(gene_it->dg());
        double f = gene_it->f();
        //std::cout << f << std::endl;

        int Ns = gene_it->Ns();
        int Na = gene_it->Na();

        OUT.write((char*)(&gene_nid),sizeof(int));
        OUT.write((char*)(&s),sizeof(double));
        OUT.write((char*)(&c),sizeof(double));
        OUT.write((char*)(&eff),sizeof(double));
        OUT.write((char*)(&dg),sizeof(double));
        OUT.write((char*)(&f),sizeof(double));
        OUT.write((char*)(&Na),sizeof(int));
        OUT.write((char*)(&Ns),sizeof(int));

        //Save length of nucleo sequence
        std::string DNAsequence = gene_it->nseq();
        int nl = DNAsequence.length();
        OUT.write((char*)&nl, sizeof(int));
        OUT.write(DNAsequence.data(), nl);
    }
}

// Dump cell summary to binary file
void PolyCell::dumpShort(std::fstream& OUT)
{
    int x;
    double y;

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    x = Na();
    OUT.write((char*)(&x),sizeof(int));

    x = Ns();
    OUT.write((char*)(&x),sizeof(int));

    y = fitness();         
    OUT.write((char*)(&y),sizeof(double));
}

// Dump cell summary to binary file
void PolyCell::dumpSeq(std::fstream& OUT, int cell_index)
{
    int x;
    double y;

    OUT.write((char*)(&cell_index),sizeof(int));

    //cell ID
    OUT.write((char*)(&ID_),sizeof(int));

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    y = fitness();         
    OUT.write((char*)(&y),sizeof(double));

    y = c_mrate_;            
    OUT.write((char*)(&y),sizeof(double));

    x = (int)(Gene_arr_.size());             
    OUT.write((char*)(&x),sizeof(int));

    for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){

        int Ns = gene_it->Ns();
        int Na = gene_it->Na();

        OUT.write((char*)(&Na),sizeof(int));
        OUT.write((char*)(&Ns),sizeof(int));

        //Save length of nucleo sequence
        std::string DNAsequence = gene_it->nseq();
        int nl = DNAsequence.length();
        OUT.write((char*)&nl, sizeof(int));
        OUT.write(DNAsequence.data(), nl);
    }
}

// Dump cell parent to binary file
void PolyCell::dumpParent(std::fstream& OUT)
{
    uint32_t a;
    a = parent();
    OUT.write((char*)(&a),sizeof(a));
}

void PolyCell::UpdateNsNa()
{
    int new_Na = 0;
    int new_Ns = 0;
    for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        new_Na += gene_it->Na();
        new_Ns += gene_it->Ns();
    }
    Total_Na_ = new_Na;
    Total_Ns_ = new_Ns;
}

// Print cell information to stdout
void PolyCell::PrintCell(int cell_ndx)
{
      char buffer[140];
      sprintf(buffer,"C %6d %6d %12e %12e %d", cell_ndx, ID_, o_mrate_, mrate(), (int)Gene_arr_.size());  
      std::cout << buffer << std::endl;
      for(auto gene_it = Gene_arr_.begin(); gene_it != Gene_arr_.end(); ++gene_it){
        int gene_nid = gene_it->num();
        double e = gene_it->e();
        double c = gene_it->conc();
        double dg = -kT*log(gene_it->dg());
        int Ns = gene_it->Ns();
        int Na = gene_it->Na();
           
        sprintf(buffer,"G %d% 2.2f %10.8f %10.8f %d %d ", gene_nid, e, c, dg, Ns, Na);
        std::cout << buffer << std::endl;  
      }
      std::cout << std::endl;
}
#endif