#ifndef POPULATION_H
#define POPULATION_H

#include "Cell.h"

class Population {

public:
	Population();
	Population(int);
	Population(std::ifstream & ,const std::string &, int, Init_Pop);

	void initMonoclonal(std::ifstream & ,const std::string &, int);
	void initPolyclonal(std::ifstream & ,const std::string &, int);

	static void initLandscape(int, std::vector<std::string>,std::string,std::string);

	void saveSnapshot(std::ofstream&, std::string, int, Encoding_Type);
	void writeSnapshotHeader(std::ofstream&, Encoding_Type);
	void writePop(std::ofstream&, Encoding_Type);

	void divide(int, int, std::ofstream&);

	int getSize() const {return size_;}

	int getMutationCount() const {return mutationCounter_;}

	double getSumFitness() const {return sumFitness_;}

	int getGeneration() const {return generation_;}

	void incrementGeneration() {generation_++;}

	std::vector<Cell> getCells() const {return cells_;}

	void fill_n(int, const Cell&);

	void incrementMutationCount(int c) {mutationCounter_ += c;}

	void incrementSize(int c) {size_ += c;}
	void decrimentSize(int c) {size_ -= c;}


	void setSize(int s) {size_ = s;}

	void resetSumFitness() {sumFitness_ = 0;}

	double addSumFitness(double); //this adds only the fitness of one cell

	void reBarcode();

	static int numberOfGenes;

	static Input_Type simType;

	static bool noMut;
    
    //Adrian's added functions
    double calcTotalFitness(); //calculates the total fitness of all cells in the population
    void MoranBirth();
    void randomCellDeath();
    
	void sortPopulationByFitness();

    
    

protected:
	int size_;
	double sumFitness_;
	int mutationCounter_;
	int generation_;
	std::vector<Cell> cells_;

};

#endif