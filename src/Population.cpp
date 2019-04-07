#include "Population.h"

int Population::numberOfGenes = 0;
Input_Type Population::simType = Input_Type::selection_coefficient;

Population::Population():
    size_(0),
    sumFitness_(0)
{
	cells_.reserve(1000);
}

Population::Population(int targetSize):
    size_(0),
    sumFitness_(0)
{
    size_ = 0;
    cells_.reserve(targetSize);
}

Population::Population(std::ifstream& startFile,const std::string & genesPath, int targetSize, Init_Pop popType):
    size_(0),
    sumFitness_(0)
{
	if(popType==Init_Pop::from_snapFile){
		initPolyclonal(startFile, genesPath, targetSize);
	}
	else{
		initMonoclonal(startFile, genesPath, targetSize);
	}
    for (auto& cell : cells_) {
        addSumFitness(cell.fitness());
    }
}

void Population::initLandscape(int fitArg, std::vector<std::string> matrixVec, std::string geneListFile, std::string genesPath){
    switch(Population::simType){
      case Input_Type::selection_coefficient:
          Cell::fromS_ = true;
          if (fitArg<5){
              // not allowed, should throw error
              std::cout << "Not a valid argument" << std::endl;
          }
          else Cell::ff_ = fitArg;

          InitMatrix();
          Population::numberOfGenes = LoadPrimordialGenes(geneListFile,genesPath);

          std::cout << "Gene count: " << numberOfGenes << std::endl;

          if(matrixVec.size()==1)
              ExtractDMSMatrix(matrixVec.front().c_str());
          else{
              Cell::useDist_ = true;
          }

        break;


      case Input_Type::stability:
          InitMatrix();
          Population::numberOfGenes = LoadPrimordialGenes(geneListFile,genesPath);
          Cell::ff_ = fitArg;
          // if DDG matrix is given
          if(matrixVec.size()){
              switch (matrixVec.size()){
                  case 2:
                      bind_DG = ExtractDDGMatrix(matrixVec.front().c_str(),Matrix_Type::is_binding);
                      std::cout << "-> Average ∆∆G_binding is " << bind_DG << " ..." << std::endl;
                  case 1:
                      fold_DG = ExtractDDGMatrix(matrixVec.front().c_str(),Matrix_Type::is_folding);
                      std::cout << "-> Average ∆∆G_folding is " << fold_DG << " ..." << std::endl;
                    break;
                  }
          }
          else{
              Cell::useDist_ = true;
          }

        break;
        
      default:
        break;
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
    if (Cell::ff_ == 5 || Cell::ff_ == 9){
        for (auto& cell : cells_) {
            cell.UpdateRates();
        }
    }
}

void Population::initPolyclonal(std::ifstream& startFile,const std::string & genesPath, int targetSize){
	// ELSE IT MUST BE POPULATED CELL BY CELL FROM SNAP FILE
    std::cout << "Creating population from file ..." << std::endl;
    cells_.reserve(targetSize) ;
    int count = 0;
    while (count < Total_Cell_Count && !startFile.eof()){
        cells_.emplace_back(startFile, genesPath);
        ++count;  
    }
    if (Cell::ff_ == 5 || Cell::ff_ == 9){
        for (auto& cell : cells_) {
            cell.UpdateRates();
        }
    }
    size_ = count;
}

void Population::divide(int targetBuffer, int targetSize, std::ofstream& LOG){
    // allocate space for temporary population
    Population newPopulation(targetBuffer);
    for (const auto& cell : cells_) {

        // fitness of cell j with respect to sum of population fitness
        double relative_fitness = cell.fitness()/getSumFitness();

        // probability parameter of binomial distribution
        std::binomial_distribution<> binCell(targetSize, relative_fitness);

        // number of progeny k is drawn from binomial distribution with N trials and mean w=relative_fitness
        int n_progeny = binCell(g_rng);
            
        // if nil, the cell will be wiped from the population
        if(n_progeny == 0) continue; 

        // iterator to current available position
        auto it = std::end(newPopulation.cells_);

        // iterator to end position of fill
        auto last = it + n_progeny;

        //cell_it->setParent(cell_it - cells_.begin());
        
        newPopulation.fill_n(n_progeny,cell);

        auto link = it;

        do{
            link->linkGenes();
            ++link;
        }while(link < last);

            //if (!noMut){
            // after filling with children, go through each one for mutation
                do{
                std::binomial_distribution<> binMut(it->genome_size(), it->mrate());
                int n_mutations = binMut(g_rng);
                    // attempt n mutations
                    for (int i=0;i<n_mutations;++i){
                        incrementMutationCount(1);
                        // change statement to switch
                        if (false){
                            // mutate and write mutation to file
                            it->ranmut_Gene(LOG,getGeneration());
                        }
                        else{
                            it->ranmut_Gene();
                        }       
                    }
                    ++it;
                }while(it < last);
            //}
        }

        // if the population is below N
        // randomly draw from progeny to pad
        while (newPopulation.getSize() < targetSize ){
            auto cell_it = newPopulation.cells_.begin();
            newPopulation.cells_.emplace_back(*(cell_it + randomNumber()*newPopulation.getSize()));
            newPopulation.incrementSize(1);
        }

        if (newPopulation.getSize() > targetSize ){
            std::shuffle(newPopulation.cells_.begin(), newPopulation.cells_.end(), g_rng);
            newPopulation.cells_.resize(targetSize);
            newPopulation.setSize(targetSize);
        }

        //alternative to shuffling
       /* while(v_size > N){
            int rand_idx = v_size*randomNumber();
            remove_at(newPopulation.cells_,rand_idx);
            v_size--;
        }*/

        Total_Cell_Count = newPopulation.getSize();
        assert (Total_Cell_Count == targetSize) ;
        
        // swap population with initial vector
        cells_.swap(newPopulation.cells_);

        // reset and update sumFitness_
        // update Ns and Na for each cell

        resetSumFitness();
        double fittest = 0;
        for (auto& cell : cells_) {
            double current = cell.fitness();
            addSumFitness(current);
            if (current > fittest) 
                fittest = current;
            cell.UpdateNsNa();
        }
        //normalize by fittest individual to prevent overflow
        if (Population::simType == Input_Type::selection_coefficient){
            sumFitness_ = 0;
            for (auto& cell : cells_) {
                sumFitness_ += cell.normalizeFit(fittest);
            }
        }
}

void Population::fill_n(int n_progeny, const Cell& c)
{
    std::fill_n(std::back_inserter(cells_),n_progeny,c);
    incrementSize(n_progeny);
}

void Population::saveSnapshot(std::ofstream& toSnapshot, std::string dirName, int currentGen, Encoding_Type encoding){
    sprintf(buffer,"%s/%s.gen%010d.snap",outPath.c_str(),dirName.c_str(), currentGen); 

    // Open snapshot file
    toSnapshot = std::ofstream(buffer, std::ios::out | std::ios::binary);
    if (toSnapshot.is_open()){
        //write
        writeSnapshotHeader(toSnapshot, encoding);
        writePop(toSnapshot, encoding);
        toSnapshot.close();
        // compress
        std::string command = "gzip -f ";
        command += buffer;
        const char *cmd = command.c_str();
        system(cmd);
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
    int idx = 1;
    switch (encoding){
        case Encoding_Type::by_default: //"normal" output format
        case Encoding_Type::full: 
                for (const auto& cell : cells_) {
                    addSumFitness(cell.fitness());
                    cell.dump(toSnapshot,idx++);
                } 
            break;
        case Encoding_Type::no_sequence: //"short" output format
                for (const auto& cell : cells_) {
                    addSumFitness(cell.fitness());
                    cell.dumpShort(toSnapshot);
                }
            break;
        case Encoding_Type::other: //dump with parent data, to be implemented
            break;
    }
}

double Population::addSumFitness(double w){
    sumFitness_ += w;
    return sumFitness_;
}