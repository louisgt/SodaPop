#include "Population.h"
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
    int currentGen = 1;
    int maxGen = currentGen + 1;
    int mutationCounter = 0;
    int numberOfGenes = 0;
    Encoding_Type outputEncoding = Encoding_Type::by_default;
    unsigned int targetPopSize = 1;
    int timeStep = 1;

    bool enableAnalysis = false;
    bool trackMutations = false;
    Init_Pop createPop = Init_Pop::from_snapFile;
    bool noMut = false;

    std::string inputType;
    std::string geneListFile;
    std::string genesPath;
    std::string outDir;
    std::string startSnapFile;
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

        TCLAP::ValueArg<std::string> prefixArg("o","prefix","Prefix to be used for snapshot files",false,"sim","filename");
        TCLAP::ValueArg<std::string> geneArg("g","gene-list","Gene list file",true,"null","filename");
        TCLAP::ValueArg<std::string> startArg("p","pop-desc","Population description file",true,"null","filename");
        TCLAP::ValueArg<std::string> libArg("l","gene-lib","Gene library directory",false,"files/genes/","filename");

        TCLAP::MultiArg<std::string> matrixArg("i","input","Input file(s) defining the fitness landscape",false,"filepath(s)");
        TCLAP::ValueArg<int> fitArg("f","fitness","Fitness function",false,5,"integer ID");
        TCLAP::ValueArg<std::string> inputArg("","sim-type","Define simulation type\n<s> (from selection coefficient, DMS or otherwise)\n<stability> (from DDG matrix or distribution)", false,"s","string");
        TCLAP::SwitchArg gammaArg("","gamma","Draw selection coefficients from gamma distribution", cmd, false);
        TCLAP::SwitchArg normalArg("","normal","Draw selection coefficients from normal distribution", cmd, false);
        TCLAP::ValueArg<double> alphaArg("","alpha","Alpha parameter of distribution\nGamma -> shape\nNormal -> mean",false,1,"double");
        TCLAP::ValueArg<double> betaArg("","beta","Beta parameter of distribution\nGamma -> scale\nNormal -> S.D.",false,1,"double");
        TCLAP::SwitchArg initArg("c","create-single","Create initial population on the fly", cmd, false);
        TCLAP::SwitchArg analysisArg("a","analysis","Enable analysis scripts", cmd, false);
        
        TCLAP::SwitchArg eventsArg("e","track-events","Track mutation events", cmd, false);
        TCLAP::ValueArg<int> seqArg("s","seq-output","Sequence output format",false,0,"integer");
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
        maxGen = maxArg.getValue();
        targetPopSize   = popArg.getValue();
        timeStep = dtArg.getValue();

        geneListFile = geneArg.getValue();
        outDir = prefixArg.getValue();
        startSnapFile = startArg.getValue();
        genesPath = libArg.getValue();

        inputType = inputArg.getValue();

        if (seedArg.isSet())
            setRngSeed(seedArg.getValue());

        std::cout << "Begin ... " << std::endl;
        if (inputType == "s"){
            Cell::fromS_ = true;
            if (fitArg.getValue()<5){
                Cell::ff_ = fitArg.getValue();
            }
            else 
                Cell::ff_ = 5;
            InitMatrix();
            numberOfGenes = LoadPrimordialGenes(geneListFile,genesPath);
            std::cout << "Gene count: " << numberOfGenes << std::endl;
            // if matrix is given
            if (matrixArg.isSet()){
                matrixVec = matrixArg.getValue();
                assert(matrixVec.size()==1);
                ExtractDMSMatrix(matrixVec.front().c_str());
            }
            else{
                Cell::useDist_ = true;
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
            numberOfGenes = LoadPrimordialGenes(geneListFile,genesPath);
            Cell::ff_ = fitArg.getValue();
            // if DDG matrix is given
            if (matrixArg.isSet()){
                matrixVec = matrixArg.getValue();
                int nMat = matrixVec.size();
                switch (nMat){
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
        }

        enableAnalysis = analysisArg.getValue();
        trackMutations = eventsArg.getValue();
        outputEncoding = intToEncoding_Type(seqArg.getValue());
        createPop = intToPop_Type(initArg.getValue());

    }catch (TCLAP::ArgException &e){
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;}

    /* general simulation initialization, can be put in a global method
    */
    if (Cell::ff_ == 7) {
        noMut = true;
        std::cout << "Mutations are not active." << std::endl;
    }

    createOutputDir(outDir);

    /* OPEN SIMULATION LOG
    */
    std::ofstream CMDLOG;
    try{
        openCommandLog(CMDLOG, outDir, argv, argc);
    }catch (std::runtime_error &e) {}

    /* OPEN MUTATION LOG
    */
    std::ofstream MUTATIONLOG;
    if (trackMutations && !noMut){
        try{
            openMutationLog(MUTATIONLOG,outDir);
        }catch (std::runtime_error &e) {}
    }

    /*OPEN POPULATION SNAPSHOT
        eventually move block to initializing method*/
    /*READ POPULATION SNAPSHOT:
    */

    std::ifstream startFile;
    try{
        openStartingPop(startSnapFile,startFile);
        readSnapshotHeader(startFile);
    }catch (std::runtime_error &e) {}

    //std::vector <Cell> Cell_arr;
    Population currentPop(startFile, genesPath, targetPopSize, createPop);

    /* should be handled by method initializing the population
    */
    startFile.close();

    Total_Cell_Count = currentPop.getSize();

    std::ofstream OUT;

    try{
        currentPop.saveSnapshot(OUT,outDir,currentGen,outputEncoding);
    }catch (std::runtime_error &e) {}

    OUT.close();
    std::string command = "gzip -f ";
    command += buffer;
    const char *cmd = command.c_str();
    system(cmd);

    /* MAIN SIMULATION BLOCK, can be put in a method outside main
    */
    
    // std::cout << "Starting evolution ..." << std::endl;
    // CMDLOG << "Starting evolution ..." << std::endl;

    // // PSEUDO WRIGHT-FISHER PROCESS
    // while (currentGen < maxGen){
    //     printProgress(currentGen*1.0/maxGen);

    //     // Population newPopulation;
    //     // newPopulation.initialize();

    //     std::vector<Cell> Cell_temp;
    //     // reserve 2N to allow overflow and prevent segfault
    //     if (targetPopSize < 10000){
    //         Cell_temp.reserve(targetPopSize * 5);
    //     }
    //     else{
    //         Cell_temp.reserve(targetPopSize * 2);
    //     }

    //     // Population.divide()
    //     // body of method is just the loop below
    //     // for each cell in the population
    //     for (auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end(); ++cell_it){


    //         // fitness of cell j with respect to sum of population fitness
    //         double relative_fitness = cell_it->fitness()/w_sum;
    //         // probability parameter of binomial distribution
    //         std::binomial_distribution<> binCell(targetPopSize, relative_fitness);
    //         // number of progeny k is drawn from binomial distribution with N trials and mean w=relative_fitness
    //         int n_progeny = binCell(g_rng);
            
    //         // if nil, the cell will be wiped from the population
    //         if(n_progeny == 0) continue; 

    //         // iterator to current available position
    //         auto it = std::end(Cell_temp);

    //         // iterator to end position of fill
    //         auto last = it + n_progeny;

    //         cell_it->setParent(cell_it - Cell_arr.begin());
    //         // fill vector with k times the current cell
    //         std::fill_n(std::back_inserter(Cell_temp),n_progeny,(*cell_it));

    //         auto link = it;

    //         do{
    //             link->linkGenes();
    //             ++link;
    //         }while(link < last);

    //         if (!noMut){
    //         // after filling with children, go through each one for mutation
    //             do{
    //             	std::binomial_distribution<> binMut(it->genome_size(), it->mrate());
    //             	int n_mutations = binMut(g_rng);
    //                 // attempt n mutations
    //                 for (int i=0;i<n_mutations;++i){
    //                     ++mutationCounter;
    //                     if (trackMutations){
    //                         // mutate and write mutation to file
    //                         it->ranmut_Gene(MUTATIONLOG,currentGen);
    //                     }
    //                     else{
    //                         it->ranmut_Gene();
    //                     }       
    //                 }
    //                 ++it;
    //             }while(it < last);
    //         }
    //     }
    //     // if the population is below N
    //     // randomly draw from progeny to pad
    //     while (Cell_temp.size() < targetPopSize ) {
    //         auto cell_it = Cell_temp.begin();
    //         Cell_temp.emplace_back(*(cell_it + randomNumber()*Cell_temp.size()));
    //     }

    //     if (Cell_temp.size() > targetPopSize ) {
    //         std::shuffle(Cell_temp.begin(), Cell_temp.end(), g_rng);
    //         Cell_temp.resize(targetPopSize ) ;
    //     }

    //     //alternative to shuffling
    //    /* while(v_size > N){
    //         int rand_idx = v_size*randomNumber();
    //         remove_at(Cell_temp,rand_idx);
    //         v_size--;
    //     }*/

    //     Total_Cell_Count = static_cast<int>(Cell_temp.size());
    //     assert(Total_Cell_Count == targetPopSize ) ;
        
    //     // swap population with initial vector
    //     Cell_arr.swap(Cell_temp);

    //     // reset and update w_sum
    //     // update Ns and Na for each cell

    //     w_sum = 0;
    //     double fittest = 0;
    //     for (auto& cell : Cell_arr) {
    //         double current = cell.fitness();
    //         w_sum += current;
    //         if (current > fittest) 
    //             fittest = current;
    //         cell.UpdateNsNa();
    //     }
    //     //normalize by fittest individual to prevent overflow
    //     if (inputType == "s"){
    //         w_sum = 0;
    //         for (auto& cell : Cell_arr) {
    //             w_sum += cell.normalizeFit(fittest);
    //         }
    //     }
        
    //     // update generation counter
    //     ++currentGen;
    //     // save population snapshot every timeStep generations
    //     if( (currentGen % timeStep) == 0){
    //         sprintf(buffer,"%s/%s.gen%010d.snap",outPath.c_str(),outDir.c_str(), currentGen); 
    //         //Open snapshot file
    //         //OUT is target output stream
    //         std::ofstream OUT(buffer, std::ios::out | std::ios::binary);
    //         if (!OUT.is_open()){
    //              std::cerr << "Snapshot file could not be opened";
    //              exit(1);
    //         }

    //         double frame_time = currentGen;
    //         OUT.write((char*)(&frame_time),sizeof(double));
    //         OUT.write((char*)(&Total_Cell_Count),sizeof(int));
    //         OUT.write((char*)(&outputEncoding),sizeof(int));

    //         int count;
    //         switch (outputEncoding){
    //             case Encoding_Type::by_default: //"normal" output format 
    //             case Encoding_Type::full: count=1;
    //                     for (const auto& cell : Cell_arr) {
    //                         cell.dump(OUT,count++);
    //                         //count++;
    //                     }
    //                 break;
    //             case Encoding_Type::no_sequence: //"short" output format
    //                     for (const auto& cell : Cell_arr) {
    //                         cell.dumpShort(OUT);
    //                     } 
    //                 break;
    //             case Encoding_Type::other: //dump with parent data, to be implemented
    //                 break;
    //         }   
              
    //         OUT.close();
    //         //compress last written file with gzip
    //         std::string command = "gzip -f ";
    //         command += buffer;
    //         const char *cmd = command.c_str();
    //         system(cmd);

    //      }
    // }

    printProgress(currentGen/maxGen);
    std::cout << std::endl;
    MUTATIONLOG.close();
    std::cout << "Done." << std::endl;
    CMDLOG << "Done." << std::endl;
    std::cout << "Total number of mutation events: " << mutationCounter << std::endl;
    CMDLOG << "Total number of mutation events: " << mutationCounter << std::endl;
    CMDLOG.close();

    // if the user toggled analysis, call shell script
    if (enableAnalysis){
        std::string script = "tools/barcodes.sh";
        std::string command = "/bin/bash "+script+" "+outDir+" "+std::to_string(maxGen)+" "+std::to_string(targetPopSize ) +" "+std::to_string(timeStep)+" "+std::to_string(outputEncoding)+" "+std::to_string(numberOfGenes);
        const char *cmd = command.c_str();
        system(cmd);
    }
    return 0;
}