#ifndef POPULATION_H
#define POPULATION_H

#include "Cell.h"

class Population {

public:
	Population();
	Population(int);
	Population(Population&);
	Population(std::ifstream & ,const std::string &, int, Init_Pop);

	void initMonoclonal(std::ifstream & ,const std::string &, int);
	void initPolyclonal(std::ifstream & ,const std::string &, int);
	void initMicrobiota(Population&);

	static void initLandscape(int, std::vector<std::string>,std::string,std::string);

	static void initExponential(double);
	static int RandomExponential();

	void saveSnapshot(std::ofstream&, std::string, int, Encoding_Type);
	void writeSnapshotHeader(std::ofstream&, Encoding_Type);
	void writePop(std::ofstream&, Encoding_Type);

	void divide(int, int, std::ofstream&, bool);

	bool addPacket(Population&);

	int getSize() const {return cells_.size();}

	int getMutationCount() const {return mutationCounter_;}

	double getSumFitness() const {return sumFitness_;}

	int getGeneration() const {return generation_;}

	void incrementGeneration() {generation_++;}

	std::vector<Cell> getCells() const {return cells_;}

	void fill_n(int, const Cell&);

	void incrementMutationCount(int c) {mutationCounter_ += c;}

	void calculateFitness();

	void resetSumFitness() {sumFitness_ = 0;}

	double addSumFitness(double);

	void reBarcode();

	void shuffle(pcg32);

	static int numberOfGenes;

	static Input_Type simType;

	static bool noMut;

	static std::exponential_distribution<> packet_exponential_;

protected:
	double sumFitness_;
	int mutationCounter_;
	int generation_;
	std::vector<Cell> cells_;

};

#endif
