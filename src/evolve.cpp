#include "PolyCell.h"
#include <tclap/CmdLine.h>
#include <unistd.h>
#include "rng.h"

/*SodaPop

Copyright (C) 2019 Louis Gauthier

    SodaPop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    SodaPop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with SodaPop.  If not, see <http://www.gnu.org/licenses/>.
 */

int main(int argc, char *argv[])
{
    // these variables will hold the parameters input (or not) by the user
    int CURRENT_GEN = 1;
    int MAX_GEN = CURRENT_GEN + 1;
    int MUTATION_CTR = 0;
    int GENE_COUNT = 0;
    int ENCODING = 0;
    unsigned int POP_SIZE=1;
    int STEP = 1;
    char buffer[200];
    bool enableAnalysis = false;
    bool trackMutations = false;
    bool createPop = false;
    bool noMut = false;

    std::string inputType;
    std::string geneListFile, genesPath;
    std::string outDir, startSnapFile;
    std::vector<std::string> matrixVec;

    // Wrap everything in a try block
    // errors in input are caught and explained to user
    try { 
        // Define the command line object
        TCLAP::CmdLine cmd("SodaPop: a multi-scale model of molecular evolution", ' ', "v1.0");

        // Define value arguments
        TCLAP::ValueArg<int> maxArg("m","maxgen","Number of generations",false,10,"int");
        TCLAP::ValueArg<int> popArg("n","size","Initial population size",false,1,"int");
        TCLAP::ValueArg<int> dtArg("t","dt","Time interval for snapshots",false,1,"int");

        //files
        TCLAP::ValueArg<std::string> prefixArg("o","prefix","Prefix to be used for snapshot files",false,"sim","filename");
        TCLAP::ValueArg<std::string> geneArg("g","gene-list","Gene list file",true,"null","filename");
        TCLAP::ValueArg<std::string> startArg("p","pop-desc","Population description file",true,"null","filename");
        TCLAP::ValueArg<std::string> libArg("l","gene-lib","Gene library directory",false,"files/genes/","filename");

        TCLAP::MultiArg<std::string> matrixArg("i","input","Input file(s) defining the fitness landscape",false,"filepath(s)");
        
        // fitness function
        TCLAP::ValueArg<int> fitArg("f","fitness","Fitness function",false,5,"integer ID");
        
        // boolean switch to use DDG as input type
        TCLAP::ValueArg<std::string> inputArg("","sim-type","Define simulation type\n<s> (from selection coefficient, DMS or otherwise)\n<stability> (from DDG matrix or distribution)", false,"s","string");

        //use gamma distribution to draw selection coefficients
        TCLAP::SwitchArg gammaArg("","gamma","Draw selection coefficients from gamma distribution", cmd, false);

        //use normal distribution to draw selection coefficients
        TCLAP::SwitchArg normalArg("","normal","Draw selection coefficients from normal distribution", cmd, false);

        //first parameter of distribution
        TCLAP::ValueArg<double> alphaArg("","alpha","Alpha parameter of distribution\nGamma -> shape\nNormal -> mean",false,1,"double");

        //second parameter of distribution
        TCLAP::ValueArg<double> betaArg("","beta","Beta parameter of distribution\nGamma -> scale\nNormal -> S.D.",false,1,"double");

        // boolean switch to create population from scratch
        TCLAP::SwitchArg initArg("c","create-single","Create initial population on the fly", cmd, false);
        
        // boolean switch to enable analysis
        TCLAP::SwitchArg analysisArg("a","analysis","Enable analysis scripts", cmd, false);
        
        // boolean switch to track mutations
        TCLAP::SwitchArg eventsArg("e","track-events","Track mutation events", cmd, false);

        // sequence output selection
        TCLAP::ValueArg<int> seqArg("s","seq-output","Sequence output format",false,0,"integer");

        // RNG seed
        TCLAP::ValueArg<unsigned long> seedArg("", "seed", "Seed value for RNG.", false, 0, "unsigned int (64-bit)");

        // Add the arguments to the CmdLine object.
        cmd.add(maxArg);
        cmd.add(popArg);
        cmd.add(dtArg);
        cmd.add(prefixArg);
        cmd.add(geneArg);
        cmd.add(startArg);
        cmd.add(libArg);
        cmd.add(fitArg);
        cmd.add(seqArg);
        cmd.add(matrixArg);
        cmd.add(alphaArg);
        cmd.add(betaArg);
        cmd.add(inputArg);

        // Parse the argv array.
        cmd.parse(argc, argv);

        // Get values from args. 
        MAX_GEN = maxArg.getValue();
        POP_SIZE = popArg.getValue();
        STEP = dtArg.getValue();

        geneListFile = geneArg.getValue();
        outDir = prefixArg.getValue();
        startSnapFile = startArg.getValue();
        genesPath = libArg.getValue();

        inputType = inputArg.getValue();

        if (seedArg.isSet())
            setRngSeed(seedArg.getValue());

        std::cout << "Begin ... " << std::endl;
        if (inputType == "s"){
            PolyCell::fromS_ = true;
            if (fitArg.getValue()<5){
                PolyCell::ff_ = fitArg.getValue();
            }
            else 
                PolyCell::ff_ = 5;
            InitMatrix();
            GENE_COUNT = LoadPrimordialGenes(geneListFile,genesPath);
            std::cout << "Gene count: " << GENE_COUNT << std::endl;
            // if matrix is given
            if (matrixArg.isSet()){
                matrixVec = matrixArg.getValue();
                assert(matrixVec.size()==1);
                ExtractDMSMatrix(matrixVec.front().c_str());
            }
            else{
                PolyCell::useDist_ = true;
                if (gammaArg.isSet()){
                    double shape = alphaArg.getValue();
                    double scale = betaArg.getValue();
                    Gene::initGamma(shape, scale);
                }
                else if (normalArg.isSet()){
                    double mean = alphaArg.getValue();
                    double stddev = betaArg.getValue();
                    Gene::initNormal(mean, stddev);
                }
            }
        }
        else if (inputType == "stability"){
            InitMatrix();
            GENE_COUNT = LoadPrimordialGenes(geneListFile,genesPath);
            PolyCell::ff_ = fitArg.getValue();
            // if DDG matrix is given
            if (matrixArg.isSet()){
                matrixVec = matrixArg.getValue();
                int nMat = matrixVec.size();
                switch (nMat){
                    case 2:
                        bind_DG = ExtractDDGMatrix(matrixVec.front().c_str(),true);
                        std::cout << "Average ∆∆G_binding is " << bind_DG << " ..." << std::endl;
                    case 1:
                        fold_DG = ExtractDDGMatrix(matrixVec.front().c_str(),false);
                        std::cout << "Average ∆∆G_folding is " << fold_DG << " ..." << std::endl;
                        break;
                }
                
            }
            else{
                PolyCell::useDist_ = true;
            }
        }

        enableAnalysis = analysisArg.getValue();
        trackMutations = eventsArg.getValue();
        ENCODING = seqArg.getValue();
        createPop = initArg.getValue();

    }catch (TCLAP::ArgException &e){
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;}

    std::cout << "Opening starting population snapshot ..." << std::endl;
    std::fstream startsnap (startSnapFile.c_str(),std::ios::in|std::ios::binary);
    if (!startsnap.is_open()){
        std::cerr << "File could not be open: "<< startSnapFile << std::endl;
        exit(2);
    }

    if (PolyCell::ff_ == 7) {
        noMut = true;
        std::cout << "Mutations are not active." << std::endl;
    }
    
    // header
    int Total_Cell_Count, dummy;
    double frame_time;
    //read frame time
    startsnap.read((char*)(&frame_time),sizeof(double));
    //read number of cells in file
    startsnap.read((char*)(&Total_Cell_Count),sizeof(int));
    //read file ENCODING
    startsnap.read((char*)(&dummy),sizeof(int));

    sprintf(buffer,"out/%s/snapshots",outDir.c_str());
    std::string outPath = buffer;
    std::cout << "Creating directory " << outPath << " ... " << (makePath(outPath) ? "OK" : "failed") << std::endl;

    sprintf(buffer,"out/%s/command.log",outDir.c_str());
    std::ofstream cmdlog;
    cmdlog.open(buffer, std::ios::out | std::ios::trunc);
    if (!cmdlog.is_open()){
        std::cerr << "Command log file could not be opened" << std::endl;
        exit(1);
    }

    std::string args;
    std::for_each( argv + 1, argv + argc , [&]( const char* c_str ){ args += std::string ( c_str ) + " "; } );
    cmdlog << "sodapop " << args << std::endl;
    cmdlog << std::endl;

    std::vector <PolyCell> Cell_arr;
    double w_sum = 0;

    // IF POPULATION IS INITIALLY MONOCLONAL
    // CREATE VECTOR WITH N IDENTICAL CELLS
    if (createPop){
        std::cout << "Creating a population of " << POP_SIZE << " cells ..." << std::endl;
        cmdlog << "Creating a population of " << POP_SIZE << " cells ..." << std::endl;
        PolyCell A(startsnap, genesPath);
        Cell_arr = std::vector <PolyCell>(POP_SIZE,A);
        //for(auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end(); ++cell_it){
        for (auto& cell : Cell_arr) {
            cell.ch_barcode(getBarcode());
        }
        if (PolyCell::ff_ == 5){
            for (auto& cell : Cell_arr) {
                cell.UpdateRates();
            }
        }
    }
    else{
        // ELSE IT MUST BE POPULATED CELL BY CELL FROM SNAP FILE
        Cell_arr.reserve(POP_SIZE);
        int count = 0;
        std::cout << "Constructing population from source " << startSnapFile.c_str() << " ..." << std::endl;
        cmdlog << "Constructing population from source " << startSnapFile.c_str() << " ..." << std::endl;
        //auto cell_it = Cell_arr.begin();
        while (count <Total_Cell_Count && !startsnap.eof()){
            Cell_arr.emplace_back(startsnap, genesPath);
            count++;  
        }
        if (PolyCell::ff_ == 5){
            for (auto& cell : Cell_arr) {
                cell.UpdateRates();
            }
        }
    }
    startsnap.close();


    std::cout << "Saving initial population snapshot ... " << std::endl;
    cmdlog << "Saving initial population snapshot ... " << std::endl;
    sprintf(buffer,"%s/%s.gen%010d.snap",outPath.c_str(),outDir.c_str(), CURRENT_GEN); 

    // Open snapshot file
    std::fstream OUT2(buffer, std::ios::out | std::ios::binary);
    if (!OUT2.is_open()){
         std::cerr << "Snapshot file could not be opened";
         exit(1);
    }

    Total_Cell_Count = Cell_arr.size();
    OUT2.write((char*)(&frame_time),sizeof(double));
    OUT2.write((char*)(&Total_Cell_Count),sizeof(int));
    OUT2.write((char*)(&ENCODING),sizeof(int));

    int idx;
    switch (ENCODING){
        case 0: //"normal" output format
        case 2: idx=1;
                //for(auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end(); ++cell_it){
                for (const auto& cell : Cell_arr) {
                    w_sum += cell.fitness();
                    cell.dump(OUT2,idx++);
                } 
            break;
        case 1: //"short" output format
                for (const auto& cell : Cell_arr) {
                    w_sum += cell.fitness();
                    cell.dumpShort(OUT2);
                }
            break;
        case 3: //dump with parent data, to be implemented
            break;
    }   

    OUT2.close();
    std::string command = "gzip -f ";
    command += buffer;
    const char *cmd = command.c_str();
    system(cmd);

    std::ofstream MUTATIONLOG;
    if (trackMutations && !noMut){
        // Open MUTATION LOG
        sprintf(buffer, "out/%s/MUTATION_LOG",outDir.c_str());
        MUTATIONLOG.open(buffer);
        if ( !MUTATIONLOG.is_open() ) {
            std::cerr << "Mutation log file could not be opened";
            exit(1);
        }
    }
    
    std::cout << "Starting evolution ..." << std::endl;
    cmdlog << "Starting evolution ..." << std::endl;

    // PSEUDO WRIGHT-FISHER PROCESS
    while (CURRENT_GEN < MAX_GEN){
        printProgress(CURRENT_GEN*1.0/MAX_GEN);
        std::vector<PolyCell> Cell_temp;
        // reserve 2N to allow overflow and prevent segfault
        if (POP_SIZE<10000){
            Cell_temp.reserve(POP_SIZE*5);
        }
        else{
            Cell_temp.reserve(POP_SIZE*2);
        }
        // for each cell in the population
        for (auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end(); ++cell_it){
            // fitness of cell j with respect to sum of population fitness
            double relative_fitness = cell_it->fitness()/w_sum;
            // probability parameter of binomial distribution
            std::binomial_distribution<> binCell(POP_SIZE, relative_fitness);
            // number of progeny k is drawn from binomial distribution with N trials and mean w=relative_fitness
            int n_progeny = binCell(g_rng);
            
            // if nil, the cell will be wiped from the population
            if(n_progeny == 0) continue; 

            // iterator to current available position
            auto it = std::end(Cell_temp);

            // iterator to end position of fill
            auto last = it + n_progeny;

            cell_it->setParent(cell_it - Cell_arr.begin());
            // fill vector with k times the current cell
            std::fill_n(std::back_inserter(Cell_temp),n_progeny,(*cell_it));

            auto link = it;

            do{
                link->linkGenes();
                link++;
            }while(link < last);

            if (!noMut){
            // after filling with children, go through each one for mutation
                do{
                	std::binomial_distribution<> binMut(it->genome_size(), it->mrate());
                	int n_mutations = binMut(g_rng);
                    // attempt n mutations
                    for (int i=0;i<n_mutations;++i){
                        MUTATION_CTR++;
                        if (trackMutations){
                            // mutate and write mutation to file
                            it->ranmut_Gene(MUTATIONLOG,CURRENT_GEN);
                        }
                        else{
                            it->ranmut_Gene();
                        }       
                    }
                    it++;
                }while(it < last);
            }
        }
        // if the population is below N
        // randomly draw from progeny to pad
        while (Cell_temp.size() < POP_SIZE){
            auto cell_it = Cell_temp.begin();
            Cell_temp.emplace_back(*(cell_it + randomNumber()*Cell_temp.size()));
        }

        if (Cell_temp.size() > POP_SIZE){
            std::shuffle(Cell_temp.begin(), Cell_temp.end(), g_rng);
            Cell_temp.resize(POP_SIZE);
        }

        //alternative to shuffling
       /* while(v_size > N){
            int rand_idx = v_size*randomNumber();
            remove_at(Cell_temp,rand_idx);
            v_size--;
        }*/

        Total_Cell_Count = (int)(Cell_temp.size());
        assert(Total_Cell_Count == POP_SIZE);
        
        // swap population with initial vector
        Cell_arr.swap(Cell_temp);

        // reset and update w_sum
        // update Ns and Na for each cell

        w_sum = 0;
        double fittest = 0;
        //for(auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end(); ++cell_it){
        for (auto& cell : Cell_arr) {
            double current = cell.fitness();
            w_sum += current;
            if (current > fittest) 
                fittest = current;
            cell.UpdateNsNa();
        }
        //normalize by fittest individual to prevent overflow
        if (inputType == "s"){
            w_sum = 0;
            for (auto& cell : Cell_arr) {
                w_sum += cell.normalizeFit(fittest);
            }
        }
        
        // update generation counter
        CURRENT_GEN++;
        // save population snapshot every STEP generations
        if( (CURRENT_GEN % STEP) == 0){
            sprintf(buffer,"%s/%s.gen%010d.snap",outPath.c_str(),outDir.c_str(), CURRENT_GEN); 
            //Open snapshot file
            //OUT2 is target output stream
            std::fstream OUT2(buffer, std::ios::out | std::ios::binary);
            if (!OUT2.is_open()){
                 std::cerr << "Snapshot file could not be opened";
                 exit(1);
            }

            double frame_time = CURRENT_GEN;
            OUT2.write((char*)(&frame_time),sizeof(double));
            OUT2.write((char*)(&Total_Cell_Count),sizeof(int));
            OUT2.write((char*)(&ENCODING),sizeof(int));

            int count;
            switch (ENCODING){
                case 0: //"normal" output format 
                case 2: count=1;
                        //for(auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end(); ++cell_it){
                        for (const auto& cell : Cell_arr) {
                            cell.dump(OUT2,count++);
                            //count++;
                        }
                    break;
                case 1: //"short" output format
                        //for(auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end(); ++cell_it){
                        for (const auto& cell : Cell_arr) {
                            cell.dumpShort(OUT2);
                        } 
                    break;
                case 3: //dump with parent data, to be implemented
                    break;
            }   
              
            OUT2.close();
            //compress last written file with gzip
            std::string command = "gzip -f ";
            command += buffer;
            const char *cmd = command.c_str();
            system(cmd);

         }
    }

    printProgress(CURRENT_GEN/MAX_GEN);
    std::cout << std::endl;
    MUTATIONLOG.close();
    std::cout << "Done." << std::endl;
    cmdlog << "Done." << std::endl;
    std::cout << "Total number of mutation events: " << MUTATION_CTR << std::endl;
    cmdlog << "Total number of mutation events: " << MUTATION_CTR << std::endl;
    cmdlog.close();

    // if the user toggled analysis, call shell script
    if (enableAnalysis){
        std::string script = "tools/barcodes.sh";
        std::string command = "/bin/bash "+script+" "+outDir+" "+std::to_string(MAX_GEN)+" "+std::to_string(POP_SIZE)+" "+std::to_string(STEP)+" "+std::to_string(ENCODING)+" "+std::to_string(GENE_COUNT);
        const char *cmd = command.c_str();
        system(cmd);
    }
    return 0;
}