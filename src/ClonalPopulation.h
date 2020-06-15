#ifndef CLONALPOPULATION_H
#define CLONALPOPULATION_H

#include "Clone.h"

class ClonalPopulation {

public:
    ClonalPopulation();
    //constructor for c clones with sizes defined by arg1 and fitness arg2
    ClonalPopulation(std::vector<int>, std::vector<double>, double, int); 
    //constructor using a snap file
    ClonalPopulation(const std::string& );

    //dump population information
    void dump(std::ofstream&) const;

    void sortClonesByFitness();

protected:
    int Csize_;
    int Nsize_;
    double sumFitness_;
    int mutationRate_;
    int generation_;
    //std::vector<double> fitnessVector_; //fitness vector of all clones in the cell
    std::vector<Clone> clones_;
    
};

#endif