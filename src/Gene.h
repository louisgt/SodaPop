#ifndef GENE_H
#define GENE_H

#include "global.h"
#include "Cell.h"

class Cell;

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

class Gene 
{
public:
    Gene();
    Gene(std::fstream&,Cell *);
    Gene(const Gene&);
    ~Gene(); 
  
    bool operator==(Gene&);
    Gene& operator=(const Gene&);

    static void initGamma(double, double);
    static void initNormal(double, double);
    static double RandomGamma();
    static double RandomNormal();

    double Mutate_Stabil_Gaussian(int, int);
    std::string Mutate_Stabil(int, int);
    double Mutate_Select_Dist(int, int);
    std::string Mutate_Select(int, int);

    void Update_Sequences(std::string);

    int num() const {return gene_idx_;}
    int geneLength() const {return gene_len_;}
    int protLength() const {return prot_len_;}
    std::string geneSeq() const {return gene_seq_;}
    double dg() const {return dg_;}
    double eff() const {return eff_;}
    double f() const {return f_;}
    int Ns() const {return Ns_;}
    int Na() const {return Na_;}

    double conc() const {return conc_;}
    double e() const {return e_;}

    double functional() const;
    double misfolded() const;
    double Pnat() const;
    double A_factor() const;
    double DDG_mean() const;

    void ch_dg(const double a){dg_ = a;}
    void ch_f(const double a){f_ = a;}
    void ch_eff(const double a){eff_ = a;}
    void ch_conc(const double c){conc_ = c;}
    void ch_Na(const int a){Na_ = a;}
    void ch_Ns(const int a){Ns_ = a;}
    void ch_e(const double e){e_ = e;}

    Cell *GetCell() const;
    const void setCell(Cell*);

    private:
        int gene_idx_;     //numeric ID pointing to primordial gene
        int gene_len_;        //length nuc seq
        int prot_len_;        //length aa seq

        int Na_;        //number of non-synonymous substitutions
        int Ns_;        //number of sysnonymous substitutions
        
        std::string gene_seq_;    //nucleotide sequence
        
        double dg_;     //stability
        double f_;      //gene "fitness"
        double eff_;    //enzymatic efficiency

        double conc_;    //concentration
        double e_;       //essentiality: between 0 and 1, can be used as a coefficient

        Cell *myCell_;

        static std::gamma_distribution<> gamma_;
        static std::normal_distribution<> normal_;     
};

#endif