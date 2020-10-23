#include "Population.h"

int Population::numberOfGenes = 0;
Input_Type Population::simType = Input_Type::selection_coefficient;
bool Population::noMut = false;

Population::Population():
    sumFitness_(0),
    mutationCounter_(0)
{
	cells_.reserve(1000);
}

Population::Population(int targetSize):
    sumFitness_(0),
    mutationCounter_(0)
{
    cells_.reserve(targetSize);
}

Population::Population(int capacity, Population& inoculum):
    sumFitness_(0),
    mutationCounter_(0)
{
    initMicrobiota(capacity,inoculum);
}

Population::Population(std::ifstream& startFile,const std::string & genesPath, int targetSize, Init_Pop popType):
    sumFitness_(0),
    mutationCounter_(0)
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

          std::cout << "-> Gene count: " << numberOfGenes << std::endl;

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
    for (auto& cell : cells_) {
        cell.ch_barcode(getBarcode());
    }
    if (Cell::ff_ == 5){
        for (auto& cell : cells_) {
            cell.UpdateRates();
        }
    }
}

void Population::initMicrobiota(int carryingCapacity, Population& inoculum){
    // from the inoculum population, transfer N cells (packetSize) to microbiota population
    std::cout << "-> Creating microbiota from inoculum." << std::endl;
    
    addPacket(carryingCapacity, inoculum);

    if (Cell::ff_ == 5){
        for (auto& cell : cells_) {
        cell.propagateFitness();
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
	    cell.propagateFitness();
            cell.UpdateRates();  
        }
    }
    cells_.shrink_to_fit();
    std::cout << "-> Created population of size " << cells_.size() << " from file ..." << std::endl;
}

// core wright-fisher process
// targetSize: estimated size of the next generation
// capacity: maximum carrying capacity of the population
// these values may or may not be equal
// for instance, when microbiota is inoculated, the targetSize can be much lower than max capacity
// max capacity is only effective when a rescaling is applied after replication of the population
void Population::divide(int targetSize, int capacity, std::ofstream& LOG, bool rescale){
    // allocate space for temporary population
    Population newPopulation(targetSize);
    double relative_fitness(1);
    int n_progeny(0);

    std::cout << "*** Size of population before WF: " << cells_.size() << std::endl;

    for (const auto& cell : cells_) {

	    //std::cout << cell.fitness() << std::endl;

        // fitness of cell j with respect to sum of population fitness
        relative_fitness = cell.fitness()/getSumFitness();

	    //std::cout << relative_fitness << std::endl;

        // probability parameter of binomial distribution
        std::binomial_distribution<> binCell(cells_.size(), relative_fitness);

        // number of progeny k is drawn from binomial distribution with N trials and mean w=relative_fitness
        n_progeny = binCell(g_rng);
            
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

            if(!Population::noMut){
                // after filling with children, go through each one for mutation
                do{
		            // hardcode beneficial mutation rate here
                    std::binomial_distribution<> binMut(100, 10e-8);
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
            }
        }

        std::cout << "*** Size of new generation: " << newPopulation.getSize() << std::endl;

        // if the population is below N
        // randomly draw from progeny to pad
        // MICROBIOTA SIMULATION:
        // * do not rescale population if below N
        // * let cells grow to carrying capacity with the trickling in of inoculum

        if(rescale){
            while (newPopulation.getSize() < capacity ){
                auto cell_it = newPopulation.cells_.begin();
                newPopulation.cells_.emplace_back(*(cell_it + randomNumber()*newPopulation.getSize()));
            }

            if (newPopulation.getSize() > capacity ){
                std::shuffle(newPopulation.cells_.begin(), newPopulation.cells_.end(), g_rng);
                newPopulation.cells_.resize(capacity);
            }

            Total_Cell_Count = newPopulation.getSize();
            assert (Total_Cell_Count == capacity) ;
        }
        else
        {
            Total_Cell_Count = newPopulation.getSize();
        }

          
        // swap population with initial vector
        cells_.swap(newPopulation.cells_);

        calculateFitness();
}

bool Population::addPacket(int carryingCapacity, Population& inoculum){

    // calculate packetSize from inoculum
    int packetSize = Population::getPacketSize(carryingCapacity, inoculum);

    // reserve space for incoming packet
    cells_.reserve(packetSize);

    // add packet
    auto cell_it = *(inoculum.cells_.begin() + randomNumber()*inoculum.getSize());

    std::cout << "-> ... Adding packet of size " << packetSize << " from inoculum." << std::endl;
    fill_n(packetSize,cell_it);

    std::cout << "-> ... " << getSize() << " cells in microbiota." << std::endl;

    calculateFitness();

    // calculate size of remainder
    int remainder = inoculum.cells_.size() - packetSize;

    if(remainder <= 0){
        std::cout << "-> ... Inoculum has been exhausted" << std::endl;
        remainder = 0;

        // clear vector
        inoculum.cells_.clear();
    }
    else{
        // resize inoculum by trimming off the added cells
        inoculum.cells_.resize(remainder);
        std::cout << "-> ... " << inoculum.getSize() << " cells remaining in inoculum." << std::endl;
    }

    return remainder > 0;
}

int Population::getPacketSize(int carryingCapacity, Population& inoculum){
    return 10000;
}

void Population::fill_n(int n_progeny, const Cell& c)
{
    std::fill_n(std::back_inserter(cells_),n_progeny,c);
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
                    cell.dump(toSnapshot,idx++);
                } 
            break;
        case Encoding_Type::no_sequence: //"short" output format
                for (const auto& cell : cells_) {
                    cell.dumpShort(toSnapshot);
                }
            break;
        case Encoding_Type::other: //dump with parent data, to be implemented
            break;
    }
}

void Population::calculateFitness(){
    resetSumFitness();
    double fittest = 0;
    for (auto& cell : cells_) {
        double current = cell.fitness();
        addSumFitness(current);
        if (current > fittest) 
                fittest = current;
        cell.UpdateNsNa();
    }

    //std::cout << "Fittest individual at " << fittest << std::endl;
    //std::cout << "Sum of fitness " << getSumFitness() << std::endl;    

    //normalize by fittest individual to prevent overflow
    if (Population::simType == Input_Type::selection_coefficient){
        //sumFitness_ = 0;
        for (auto& cell : cells_) {
            cell.normalizeFit(fittest);
        }
    }
}

double Population::addSumFitness(double w){
    sumFitness_ += w;
    return sumFitness_;
}

void Population::shuffle(pcg32 engine){
    std::shuffle(this->cells_.begin(), this->cells_.end(), engine);
}

void Population::reBarcode(){
    // open stream to read population summary
    std::ifstream popf("files/start/initial_population_barcoded.dat");
    if (!popf.is_open()) {
        std::cerr << "File could not be opened";
        exit(1);
    }

    std::shuffle(cells_.begin(), cells_.end(), g_rng);

    int k = 0;
    int b = 0;

    std::string line;
    while (!popf.eof()) {
        getline(popf, line);
        std::string word;
        std::istringstream iss(line, std::istringstream:: in );
        iss >> word;
	if (word == "C") {
            iss >> word; //cell count
            int count = atoi(word.c_str());
            iss >> word; //cell files
	    iss >> word; // barcode
	    //std::string bc = getBarcode();
	    std::string bc = word.c_str();
	    for (int i = 0; i < count; i++) {
            	//std::cout << cells_.at(i+k).barcode() << std::endl;
		cells_.at(i+k).ch_barcode(bc);
            }
	    k += count;
            b++;
	}

    }
    std::cout << "Generated " << b << " barcodes" << std::endl;
    assert(k == getSize());
    std::cout << "Finished rebarcoding" << std::endl;
}
