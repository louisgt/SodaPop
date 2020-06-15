#ifndef CLONE_H
#define CLONE_H

//#include "Cell.h"
#include "global.h"

class Clone{

public:
    Clone();
    Clone(double, int, int); //constructor for c clones with n cells each
    Clone(double , std::vector<unsigned long int>); //constructor for vector of cell IDs
    
    void dump(std::ofstream&, int) const;
    double sumFitness(){return sumFitness_ = cloneFitness_*(clonalCells_.size());}

    /*
    overloaded relational operators required for sorting 
    */
    bool operator< (const Clone &other) const {
        return sumFitness_ < other.sumFitness_;
    }
    bool operator> (const Clone &other) const {
        return sumFitness_ > other.sumFitness_;
    }
    /*
    bool operator== (const Clone &other) const {
        return sumFitness_ == other.sumFitness_;
    }
    bool operator<= (const Clone &other) const {
        return sumFitness_ <= other.sumFitness_;
    }
    bool operator>= (const Clone &other) const {
        return sumFitness_ >= other.sumFitness_;
*/

protected:
	//int cloneSize_;
	double cloneFitness_;
    double sumFitness_;
	std::vector<unsigned long int> clonalCells_; //numeric IDs for each cell in the clone   
};

#endif