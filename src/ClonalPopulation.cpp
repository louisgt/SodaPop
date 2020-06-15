#include "ClonalPopulation.h"


ClonalPopulation::ClonalPopulation():
    Csize_(0),
    Nsize_(0),
    sumFitness_(0),
    mutationRate_(0),
    generation_(0)
{
    //fitnessVector_.reserve(10);
    clones_.reserve(10);
}

//constructor for c clones with sizes defined by arg1 and fitness arg2
ClonalPopulation::ClonalPopulation(std::vector<int> C, std::vector<double> f, double m, int g)
{
    mutationRate_ = m;
    generation_ = g;
    
    if(f.size() == C.size()){
        Csize_ = C.size();
        std::cout << "Total clone count: " << Csize_ << std::endl;

    }
    else {
        throw std::runtime_error("Sizes of clone vector and fitness vector don't match.");
    }

    //total population size
    Nsize_ = 0;
    for(int i=0; i<Csize_; i++){
        Nsize_ += C[i];
    }
    std::cout << "Total cell count: " << Nsize_ << std::endl;

    //total population fitness
    sumFitness_ = 0;
    for(int i=0; i<Csize_; i++){
        sumFitness_ += f[i] * C[i];
    }
    std::cout << "Mean population fitness: " << (sumFitness_/Nsize_) << std::endl;

    
    //generate cell IDs
    std::vector<unsigned long int> v(Nsize_) ; 
    std::iota (std::begin(v), std::end(v), 0); 
    //std::shuffle(v.begin(), v.end(), g_rng);
    
    //populate Clone array
    unsigned int cs = 0;
    for(int i=0; i<Csize_; i++){
              
        std::vector<unsigned long int>::const_iterator first = v.begin() + cs;
        std::vector<unsigned long int>::const_iterator last = v.begin() + cs + C[i];
        std::vector<unsigned long int> subV(first, last);
        
        cs += C[i];

        Clone tempClone(f[i], subV);
        clones_.emplace_back(tempClone);
    }
    //sort clones based on fitness
}

//Construct population from a snapshot
ClonalPopulation::ClonalPopulation(const std::string& snapFile)
{  
    std::ifstream file (snapFile.c_str());
    if (!file.is_open()){
        std::cerr << "File could not be open: "<< snapFile <<std::endl;
        exit(2);
    }
    std::cout << "Loading the clonal population " << snapFile << std::endl;
    
    
    Csize_ = 0;
    Nsize_ = 0;
    sumFitness_ = 0;
    mutationRate_ = 0;
    generation_ = 0;
    
    clones_.reserve(10);    
    
    //parse the input snap file
    int clonectr = 0;
    int Nctr = 0;
    while(!file.eof()){
        std::string word, line;
        getline(file,line);
        
        std::istringstream iss(line, std::istringstream::in);
        iss >> word; 
        if ( word =="Population" ){
            iss >> word;
            Csize_ = atoi(word.c_str()); 
            assert ( Csize_ > 0); //number of clones must be at least 1

            iss >> word;
            Nsize_ = atoi(word.c_str()); 
            assert ( Nsize_ > 0); //Population size must be at least 1
            
            iss >> word;
            sumFitness_ = atof(word.c_str()); 
            assert ( sumFitness_ > 0); //Average fitness 
            
            iss >> word;
            generation_ = atoi(word.c_str()); 
            assert ( generation_ > 0);              
        }else if(word =="Clone"){
            
            iss >> word;
            int Cid = atoi(word.c_str()); 
            if (Cid != clonectr){
                std::cerr << "Clone IDs in snap file are not ordered." << std::endl;
                exit(2);
            }
            
            iss >> word;
            int cellNumber = atoi(word.c_str()); 
            assert ( cellNumber > 0); //Clone must at least have 1 cell
            
            iss >> word;
            double cloneFitness = atof(word.c_str()); 
            
            //parse cell IDs
            std::vector<unsigned long int> cells={};
            int i;
            for(i=0; i<cellNumber; i++){
                iss >> word;
                int ID = atoi(word.c_str()); assert(ID>=0);
                cells.emplace_back(ID); 
            }
            if (i != cellNumber){
                std::cerr << "Number of cell IDs do not match the indicated cell count in clone: " << Cid << std::endl;
                exit(2);
            }
            Nctr += cellNumber;
            
            //Create clone vector
            Clone tempClone(cloneFitness, cells);
            clones_.emplace_back(tempClone);
            
            clonectr++;
        }
    }
    
    if (Csize_ != clonectr){
        std::cerr << "Inconsistent clone size in snap file." << std::endl;
        exit(2);
    }
    
    if (Nsize_ != Nctr){
        std::cerr << "Inconsistent population size in snap file." << std::endl;
        exit(2);
    }
}


void ClonalPopulation::dump(std::ofstream& OUT) const
{ 
    OUT << "Population\t" << Csize_ << "\t" << Nsize_ << "\t" << sumFitness_/Nsize_ << "\t" << generation_ << std::endl;
    
    int i = 0;
    for (auto& clone : clones_){
        clone.dump(OUT,i);
        i++;
    }
}

void ClonalPopulation::sortClonesByFitness()
{
    std::sort(clones_.begin(), clones_.end(), std::greater<Clone>());
}

