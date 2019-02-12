#ifndef POLYCELL_H
#define POLYCELL_H
#include "Cell.h"

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
    typedef double(PolyCell::*funcPtr)(void) const;

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
    double flux() const;
    double toxicity() const;
    double metabolicOutput() const;
    double multiplicative() const;
    double neutral() const;
    double noMut() const;
    double fold() const;
    double growthRate() const;
    void UpdateRates();

    void ranmut_Gene();
    void ranmut_Gene(std::ofstream&, int);
    void change_exprlevel();
    double normalizeFit(double);
    void dump(std::fstream&, int) const;
    void dumpShort(std::fstream&) const;
    void dumpSeq(std::fstream&, int) const;
    void dumpParent(std::fstream&) const;

    void PrintCell(int) const;
    
    int Na() const {return Total_Na_;}
    int Ns() const {return Total_Ns_;}
    void UpdateNsNa();

protected:
    int Total_Ns_;
    int Total_Na_;

    // function pointer to select fitness function
    funcPtr fit;
};

#endif