#ifndef POLYCELL_H
#define POLYCELL_H

#include "Cell.h"
/*SodaPop
Copyright (C) 2017 Louis Gauthier

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
	protected:
        double fitness_;
        int Total_Ns_;
        int Total_Na_;

	public:
        static int ff_;
        static bool useGauss_;
		PolyCell();
	    PolyCell(std::fstream&);			    
	    PolyCell(std::fstream&, const std::string&);

	    void UpdateRates();
	    void ranmut_Gene(std::ofstream&, int);
	    void change_exprlevel();
	    void dump(std::fstream&, int);
        void dumpShort(std::fstream&);
	    void PrintCell(int);
	    void FillGene_L();
        void ch_Fitness(double f){fitness_ = f;}
        // Fitness functions
        double flux();
        double toxicity();
        double metabolicOutput();
        const double neutral();
        ///////////////////////
        const double fitness();
        int Na(){return Total_Na_;}
        int Ns(){return Total_Ns_;}
        void UpdateNsNa();
};

// By default the fitness function is set to neutral
int PolyCell::ff_ = 4;
bool PolyCell::useGauss_ = false;

PolyCell::PolyCell(){}
PolyCell::PolyCell(std::fstream& f) : Cell(f)
{
	// Update current rates
  	this->UpdateRates();  
  	// Fill gene length array
  	this->FillGene_L();
}    
PolyCell::PolyCell(std::fstream& f, const std::string& s) : Cell(f,s)
{
	// Update current rates
  	this->UpdateRates();  
  	// Fill gene length array
  	this->FillGene_L();
}

// Initialize the cummulative gene length array
void PolyCell::FillGene_L()
{
    int sum = 0;
    std::vector<Gene>::iterator i;
    for(i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        sum+= (*i).length();
        Gene_L_.push_back(sum);
    }
}

// FLUX FITNESS FUNCTION
double PolyCell::flux()
{
    double f = 0;
    for(std::vector<Gene>::iterator i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        f = f + 1/((*i).functional());
    }
    return A_FACTOR/f;
}

// TOXICITY FITNESS FUNCTION
double PolyCell::toxicity()
{
    double f = 0;
    for(std::vector<Gene>::iterator i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        f +=(*i).misfolded();
    }
    double w = 1-(COST*f);
    if(w<0) return 0;
    else return w;
}

// METABOLIC FLUX FITNESS FUNCTION
double PolyCell::metabolicOutput()
{
    double f = 0;
    double t = 0;
    for(std::vector<Gene>::iterator i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        f = f + 1/((*i).functional());
        t +=(*i).misfolded();
    }
    f = A_FACTOR/f;
    t = COST*f;

    double w = f - t;
    if(w<0) return 0;
    else return w;
}

// NEUTRAL FITNESS FUNCTION
const double PolyCell::neutral()
{
    return 1;
}

const double PolyCell::fitness()
{
    return fitness_;
}

void PolyCell::UpdateRates()
{
    switch(PolyCell::ff_){
        case 1: ch_Fitness(this->flux());
            break;
        case 2: ch_Fitness(this->toxicity());
            break;
        case 3: ch_Fitness(this->metabolicOutput());
            break;
        case 4: ch_Fitness(this->neutral());
            break;
        default:;
    }
}

void PolyCell::ranmut_Gene(std::ofstream& log,int ctr)
{
    // get genome size
    int L = Gene_L_.back();
    // pick random site to mutate
    int site = (int) ( L * (uniformdevptr->doub()));

    // find the corresponding gene
    std::vector<Gene>::iterator j = Gene_arr_.begin();
    VectInt_citerator k = Gene_L_.begin();

    if(site >= (*k)){
    // random number generated is greater than
    // the cummulative sum of genes
         for(k = Gene_L_.begin(); k != Gene_L_.end(); ++k){
             if( site<(*k) ) break;
             j++; 
         }        
         k--;
         site = site - (*k);        
    }

    std::string mutation = "";

    int bp = (int) (3 * (uniformdevptr->doub()));
    double wi = fitness();
    if(useGauss_){
        (*j).Mutate_BP_Gaussian(site,bp);
    }
    else{
        mutation = (*j).Mutate_BP(site,bp);
    }
            
    UpdateRates();
    double wf = fitness();
    double s = wf - wi;

    // save beneficial mutations to log
    // we could save all mutations with abs(s) >= some value x
    log << barcode().c_str() << "\t";
    log << fixed;
    log << mutation << "\t";
    log << s << "\t";
    log << ctr << endl;
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

    y = o_mrate_;
    OUT.write((char*)(&y),sizeof(double));
    
    y = c_mrate_;		 	 
    OUT.write((char*)(&y),sizeof(double));

    x = (int)(Gene_arr_.size());		 	 
    OUT.write((char*)(&x),sizeof(int));

   for(std::vector<Gene>::iterator i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        int gene_nid = (*i).num();
        double s = (*i).e;
        double c = (*i).conc;
        double dg = -kT*log((*i).dg());

        int Ns = i->Ns();
        int Na = i->Na();

        OUT.write((char*)(&gene_nid),sizeof(int));
        OUT.write((char*)(&s),sizeof(double));
        OUT.write((char*)(&c),sizeof(double));
        OUT.write((char*)(&dg),sizeof(double));
        OUT.write((char*)(&Na),sizeof(int));
        OUT.write((char*)(&Ns),sizeof(int));

        //Save length of nucleo sequence
        std::string AAsequence = GetProtFromNuc((*i).nseq());
        int nl = AAsequence.length();
        OUT.write((char*)&nl, sizeof(int));
        OUT.write(AAsequence.data(), nl);
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

void PolyCell::UpdateNsNa()
{
    int new_Na = 0;
    int new_Ns = 0;
    for(std::vector<Gene>::iterator i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        new_Na += (*i).Na();
        new_Ns += (*i).Ns();
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
      for(std::vector<Gene>::iterator i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        //cout << "X ";
        int gene_nid = (*i).num();
        double e = (*i).e;
        double c = (*i).conc;
        double dg = -kT*log((*i).dg());
        double ddg = -kT*log((*i).CheckDG());
        int Ns = (*i).Ns();
        int Na = (*i).Na();
           
        sprintf(buffer,"G %d% 2.2f %10.8f %10.8f %10.8f %d %d ", gene_nid, e, c, dg, ddg, Ns, Na);
        std::cout << buffer << std::endl;  
      }
      std::cout << std::endl;
}

#endif
