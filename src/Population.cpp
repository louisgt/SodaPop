#include "Population.h"

Population::Population(){
	size_ = 0;
	cells_.reserve(1000);
}

Population::Population(std::ifstream& startFile,const std::string & genesPath, int targetSize, Init_Pop popType){
	if(popType==Init_Pop::from_snapFile){
		initPolyclonal(startFile, genesPath, targetSize);
	}
	else{
		initMonoclonal(startFile, genesPath, targetSize);
	}
}

void Population::initMonoclonal(std::ifstream& startFile,const std::string & genesPath, int targetSize){
	std::cout << "Creating a population of " << targetSize   << " cells ..." << std::endl;
    //CMDLOG << "Creating a population of " << targetSize   << " cells ..." << std::endl;
    Cell A(startFile, genesPath);
    cells_ = std::vector <Cell>(targetSize , A);
    size_ = targetSize;
    for (auto& cell : cells_) {
        cell.ch_barcode(getBarcode());
    }
    if (Cell::ff_ == 5){
        for (auto& cell : cells_) {
            cell.UpdateRates();
        }
    }
}

void Population::initPolyclonal(std::ifstream& startFile,const std::string & genesPath, int targetSize){
	// ELSE IT MUST BE POPULATED CELL BY CELL FROM SNAP FILE
    cells_.reserve(targetSize) ;
    int count = 0;
    while (count < Total_Cell_Count && !startFile.eof()){
        cells_.emplace_back(startFile, genesPath);
        ++count;  
    }
    if (Cell::ff_ == 5){
        for (auto& cell : cells_) {
            cell.UpdateRates();
        }
    }
    size_ = count;
}

void Population::saveSnapshot(std::ofstream& toSnapshot, std::string dirName, int currentGen, Encoding_Type encoding){
    sprintf(buffer,"%s/%s.gen%010d.snap",outPath.c_str(),dirName.c_str(), currentGen); 

    // Open snapshot file
    toSnapshot = std::ofstream(buffer, std::ios::out | std::ios::binary);
    if (toSnapshot.is_open()){
        //write
        writeSnapshotHeader(toSnapshot, encoding);
        writePop(toSnapshot, encoding);
    }
    else{
        // error opening file, throw exception
        throw std::runtime_error("Unable to open file for output.");
    }
}

void Population::writeSnapshotHeader(std::ofstream& toSnapshot, Encoding_Type encoding)
{
    toSnapshot.write((char*)(&frame_time),sizeof(double));
    toSnapshot.write((char*)(&Total_Cell_Count),sizeof(int));
    toSnapshot.write((char*)(&encoding),sizeof(int));
}

void Population::writePop(std::ofstream& toSnapshot, Encoding_Type encoding){
    int idx;
    switch (encoding){
        case Encoding_Type::by_default: //"normal" output format
        case Encoding_Type::full: 
                idx=1;
                for (const auto& cell : cells_) {
                    w_sum += cell.fitness();
                    cell.dump(toSnapshot,idx++);
                } 
            break;
        case Encoding_Type::no_sequence: //"short" output format
                for (const auto& cell : cells_) {
                    w_sum += cell.fitness();
                    cell.dumpShort(toSnapshot);
                }
            break;
        case Encoding_Type::other: //dump with parent data, to be implemented
            break;
    }   
}