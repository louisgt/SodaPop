// Cell.h

#ifndef CELL_H
#define CELL_H

#include "Gene.h"
#include "global.h"

class Gene;

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

class Cell {
public:
    
    Cell();
    Cell(std::fstream & );
    Cell(std::fstream & ,const std::string & );

    virtual ~Cell() {};

    int total_mutations(const int & );
    void FillGene_L();
    void linkGenes();

    virtual void UpdateRates() = 0;
    virtual void dump(std::fstream & , int) = 0;
    virtual void PrintCell(int) = 0;

    const int ID() {return ID_;}
    uint32_t parent() {return parent_;}
    const double mrate() {return c_mrate_;}
    const int gene_count() {return Gene_arr_.size();}
    const int genome_size() {return Gene_L_.back();}
    const std::string barcode() {return barcode_;}

    void change_ID(int a) {ID_ = a;}
    void setParent(uint32_t a) {parent_ = a;}
    void ch_barcode(std::string s) {barcode_ = s;}

    void ch_Fitness(double f){fitness_ = f;}
    const double fitness();
    
protected:
    // organism barcode
    std::string barcode_;

    // organism ID
    int ID_;

    // parent ID
    int parent_;

    // initial mutation rate
    double o_mrate_;

    // current mutation rate
    double c_mrate_;

    // organismal fitness
    double fitness_;

    //Array of genes
    std::vector <Gene> Gene_arr_;

    //Cummulative sum of gene lengths (i.e. genome size)
    VectInt Gene_L_;
    
};
#endif