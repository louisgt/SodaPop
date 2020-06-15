#include "Clone.h"

Clone::Clone():
    cloneFitness_(0),
    sumFitness_(0)
{
	clonalCells_.reserve(1000);
}

//fill clone with size cells labeled by same ID
Clone::Clone(double fitness, int size, int ID)
{
    cloneFitness_ = fitness;
    
    //generate cell IDs
    clonalCells_ = {};
    std::fill (clonalCells_.begin(),  clonalCells_.begin() + size, ID); 
    //std::shuffle(v.begin(), v.end(), g_rng);
    
    sumFitness_ = cloneFitness_*(clonalCells_.size());
}

Clone::Clone(double fitness, std::vector<unsigned long int> carray)
{
    cloneFitness_ = fitness;
    clonalCells_ = carray;
    sumFitness_ = cloneFitness_*(clonalCells_.size());
    std::cout << "create clone " << sumFitness_;
}
    
// Dump clone information to file
void Clone::dump(std::ofstream& OUT, int clone_index) const
{   
    OUT << "Clone\t" << clone_index << "\t" << clonalCells_.size() << "\t" << cloneFitness_ << "\t";
    //OUT << "Cells\t";
    for(const auto& x : clonalCells_){
       OUT << x <<" ";
    }
    OUT << std::endl;
}