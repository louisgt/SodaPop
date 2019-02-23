// Population.h

#ifndef POPULATION_H
#define POPULATION_H

#include "Cell.h"

class Population {
public:
	Population();
	Population(std::ifstream & ,const std::string &, int, Init_Pop);

	void initMonoclonal(std::ifstream & ,const std::string &, int);
	void initPolyclonal(std::ifstream & ,const std::string &, int);

	void saveSnapshot(std::ofstream&, std::string, int, Encoding_Type);
	void writeSnapshotHeader(std::ofstream&, Encoding_Type);
	void writePop(std::ofstream&, Encoding_Type);

	int getSize() const {return size_;}

protected:

	int size_;

	std::vector<Cell> cells_;

};
#endif