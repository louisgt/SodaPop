#include "Cell.h"
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

    double lambda_plus = 0;
    double lambda_minus = 0;
    double r_prime = 0;
    double a_for_s_x = 0;
    double b_for_s_x = 0;

    char buffer[200];
    bool enableAnalysis = false;
    bool trackMutations = false;
    Init_Pop createPop = Init_Pop::from_snapFile;
    bool noMut = false;

    bool simul_pangenomes_evolution = false;
    bool track_pangenomes_evolution = false;

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

        // boolean switch to simulate neutral Pangenomes evolution (Gain and loss of genes)
        TCLAP::SwitchArg pangenomes_evo_Arg("V","pangenomes-evolution","simulate neutral Pangenomes evolution (random Gain and loss of genes)", cmd, false);
        
        // boolean switch to track pangenomes evolution events (Gain and loss of genes) and also track the evolution of genome size and loss rate / gain rate ratio
        TCLAP::SwitchArg track_pangenomes_evo_Arg("T","track-PanEv","track pangenomes evolution events (Gain and loss of genes) and also track the evolution of genome size and loss rate / gain rate ratio", cmd, false);

        //parameter a to calculate the selection coefficient s(x)
        TCLAP::ValueArg<double> a_Arg("","aForSx","parameter a to calculate s(x)",false,1,"double");

        //parameter b to calculate the selection coefficient s(x)
        TCLAP::ValueArg<double> b_Arg("","bForSx","parameter b to calculate s(x)",false,1,"double");

        //parameter r_prime to calculate r(x)
        TCLAP::ValueArg<double> r_prime_Arg("","rPrime","parameter rPrime to calculate r(x)",false,1,"double");

        //parameter lambda_plus to calculate alpha(x)
        TCLAP::ValueArg<double> lambda_plus_Arg("","lambdaPlus","parameter lambdaPlus to calculate alpha(x)",false,1,"double");

        //parameter lambda_minus to calculate beta(x)***
        TCLAP::ValueArg<double> lambda_minus_Arg("","lambdaMinus","parameter lambdaMinus to calculate beta(x)",false,1,"double");
        
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

        cmd.add(a_Arg);
        cmd.add(b_Arg);
        cmd.add(r_prime_Arg);
        cmd.add(lambda_plus_Arg);
        cmd.add(lambda_minus_Arg);

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

        simul_pangenomes_evolution = pangenomes_evo_Arg.getValue();
        track_pangenomes_evolution = track_pangenomes_evo_Arg.getValue();
        lambda_plus = lambda_plus_Arg.getValue();
        lambda_minus = lambda_minus_Arg.getValue();
        r_prime = r_prime_Arg.getValue();
        a_for_s_x = a_Arg.getValue();
        b_for_s_x = b_Arg.getValue();

    }catch (TCLAP::ArgException &e){
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;}


    /*OPEN POPULATION SNAPSHOT
        eventually move block to initializing method*/
    /*READ POPULATION SNAPSHOT:
        -put in global method (throws file read error)
        -takes stream to read from
        -wrap in try/catch
    */

    std::ifstream startFile;
    try{
        openStartingPop(startSnapFile,startFile);
        readSnapshotHeader(startFile);

    }catch (std::runtime_error &e) {}


    /* general simulation initialization, can be put in a global method
    */
    if (Cell::ff_ == 7) {
        noMut = true;
        std::cout << "Mutations are not active." << std::endl;
    }

    /* general simulation initialization, can be put in a global method
    */

    sprintf(buffer,"out/%s/snapshots",outDir.c_str());
    std::string outPath = buffer;
    std::cout << "Creating directory " << outPath << " ... " << (makePath(outPath) ? "OK" : "failed") << std::endl;

    /* general simulation initialization, can be put in a global method
    */

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

    /* POPULATION INITIALIZATION
        put in global method
    */

    std::vector <Cell> Cell_arr;
    double w_sum = 0;

    // IF POPULATION IS INITIALLY MONOCLONAL
    // CREATE VECTOR WITH N IDENTICAL CELLS
    if(Init_Pop::from_cellFile){
        std::cout << "Creating a population of " << targetPopSize   << " cells ..." << std::endl;
        cmdlog << "Creating a population of " << targetPopSize   << " cells ..." << std::endl;
        Cell A(startFile, genesPath);
        Cell_arr = std::vector <Cell>(targetPopSize , A);
        for (auto& cell : Cell_arr) {
            cell.ch_barcode(getBarcode());
        }
        if (Cell::ff_ == 5){
            for (auto& cell : Cell_arr) {
                cell.UpdateRates();
            }
        }
    }
    else{
        // ELSE IT MUST BE POPULATED CELL BY CELL FROM SNAP FILE
        Cell_arr.reserve(targetPopSize ) ;
        int count = 0;
        std::cout << "Constructing population from source " << startSnapFile.c_str() << " ..." << std::endl;
        cmdlog << "Constructing population from source " << startSnapFile.c_str() << " ..." << std::endl;
        while (count <Total_Cell_Count && !startFile.eof()){
            Cell_arr.emplace_back(startFile, genesPath);
            ++count;  
        }
        if (Cell::ff_ == 5){
            for (auto& cell : Cell_arr) {
                cell.UpdateRates();
            }
        }
    }

    /* should be handled by method initializing the population
    */
    startFile.close();


    std::cout << "Saving initial population snapshot ... " << std::endl;
    cmdlog << "Saving initial population snapshot ... " << std::endl;
    sprintf(buffer,"%s/%s.gen%010d.snap",outPath.c_str(),outDir.c_str(), currentGen); 

    // Open snapshot file
    std::ofstream OUT2(buffer, std::ios::out | std::ios::binary);
    if (!OUT2.is_open()){
         std::cerr << "Snapshot file could not be opened";
         exit(1);
    }

    Total_Cell_Count = Cell_arr.size();
    OUT2.write((char*)(&frame_time),sizeof(double));
    OUT2.write((char*)(&Total_Cell_Count),sizeof(int));
    OUT2.write((char*)(&outputEncoding),sizeof(int));

    int idx;
    switch (outputEncoding){
        case Encoding_Type::by_default: //"normal" output format
        case Encoding_Type::full: 
                idx=1;
                for (const auto& cell : Cell_arr) {
                    w_sum += cell.fitness();
                    cell.dump(OUT2,idx++);
                } 
            break;
        case Encoding_Type::no_sequence: //"short" output format
                for (const auto& cell : Cell_arr) {
                    w_sum += cell.fitness();
                    cell.dumpShort(OUT2);
                }
            break;
        case Encoding_Type::other: //dump with parent data, to be implemented
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

    /* MAIN SIMULATION BLOCK, can be put in a method outside main
    */
    
    std::cout << "Starting evolution ..." << std::endl;
    cmdlog << "Starting evolution ..." << std::endl;

    //If the pangenomes evolution (gain through HGT and loss) option is activated, confirm that there is more than one species simulated. Otherwhise there can't be HGT and the software will loop indefinitely
    if (simul_pangenomes_evolution){
        std::string barcode_initial = Cell_arr.begin()->barcode();
        bool more_than_one_sp = false;
        for(auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end(); ++cell_it){
            if (cell_it->barcode() != barcode_initial){
                more_than_one_sp = true;
                break;
            }
        }
        if (!more_than_one_sp){
            std::cerr << "Error! You activated the option to simulate gene loss and gene gain through HGT but there is only one species simulated! Add another species to initial snapshot";
            exit(1);
        }
    }

    // PSEUDO WRIGHT-FISHER PROCESS
    while (currentGen < maxGen){
        printProgress(currentGen*1.0/maxGen);
        std::vector<Cell> Cell_temp;
        // reserve 2N to allow overflow and prevent segfault
        if (targetPopSize < 10000){
            Cell_temp.reserve(targetPopSize * 5);
        }
        else{
            Cell_temp.reserve(targetPopSize * 2);
        }

        // for each cell in the population
        for (auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end(); ++cell_it){


            // fitness of cell j with respect to sum of population fitness
            double relative_fitness = cell_it->fitness()/w_sum;
            // probability parameter of binomial distribution
            std::binomial_distribution<> binCell(targetPopSize ,  relative_fitness);
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
                ++link;
            }while(link < last);

            if (!noMut){
            // after filling with children, go through each one for mutation
                do{
                	std::binomial_distribution<> binMut(it->genome_size(), it->mrate());
                	int n_mutations = binMut(g_rng);
                    // attempt n mutations
                    for (int i=0;i<n_mutations;++i){
                        ++mutationCounter;
                        if (trackMutations){
                            // mutate and write mutation to file
                            it->ranmut_Gene(MUTATIONLOG,currentGen);
                        }
                        else{
                            it->ranmut_Gene();
                        }       
                    }
                    ++it;
                }while(it < last);
            }
        }
        // if the population is below N
        // randomly draw from progeny to pad
        while (Cell_temp.size() < targetPopSize ) {
            auto cell_it = Cell_temp.begin();
            Cell_temp.emplace_back(*(cell_it + randomNumber()*Cell_temp.size()));
        }

        if (Cell_temp.size() > targetPopSize ) {
            std::shuffle(Cell_temp.begin(), Cell_temp.end(), g_rng);
            Cell_temp.resize(targetPopSize ) ;
        }

        //alternative to shuffling
       /* while(v_size > N){
            int rand_idx = v_size*randomNumber();
            remove_at(Cell_temp,rand_idx);
            v_size--;
        }*/

        Total_Cell_Count = static_cast<int>(Cell_temp.size());
        assert(Total_Cell_Count == targetPopSize ) ;
        
        // swap population with initial vector
        Cell_arr.swap(Cell_temp);

        // reset and update w_sum
        // update Ns and Na for each cell

        w_sum = 0;
        double fittest = 0;
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
        ++currentGen;
        // save population snapshot every timeStep generations
        if( (currentGen % timeStep) == 0){
            sprintf(buffer,"%s/%s.gen%010d.snap",outPath.c_str(),outDir.c_str(), currentGen); 
            //Open snapshot file
            //OUT2 is target output stream
            std::ofstream OUT2(buffer, std::ios::out | std::ios::binary);
            if (!OUT2.is_open()){
                 std::cerr << "Snapshot file could not be opened";
                 exit(1);
            }

            double frame_time = currentGen;
            OUT2.write((char*)(&frame_time),sizeof(double));
            OUT2.write((char*)(&Total_Cell_Count),sizeof(int));
            OUT2.write((char*)(&outputEncoding),sizeof(int));

            int count;
            switch (outputEncoding){
                case Encoding_Type::by_default: //"normal" output format 
                case Encoding_Type::full: count=1;
                        for (const auto& cell : Cell_arr) {
                            cell.dump(OUT2,count++);
                            //count++;
                        }
                    break;
                case Encoding_Type::no_sequence: //"short" output format
                        for (const auto& cell : Cell_arr) {
                            cell.dumpShort(OUT2);
                        } 
                    break;
                case Encoding_Type::other: //dump with parent data, to be implemented
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

    printProgress(currentGen/maxGen);
    std::cout << std::endl;
    MUTATIONLOG.close();
    std::cout << "Done." << std::endl;
    cmdlog << "Done." << std::endl;
    std::cout << "Total number of mutation events: " << mutationCounter << std::endl;
    cmdlog << "Total number of mutation events: " << mutationCounter << std::endl;
    cmdlog.close();

    // if the user toggled analysis, call shell script
    if (enableAnalysis){
        std::string script = "tools/barcodes.sh";
        std::string command = "/bin/bash "+script+" "+outDir+" "+std::to_string(maxGen)+" "+std::to_string(targetPopSize ) +" "+std::to_string(timeStep)+" "+std::to_string(outputEncoding)+" "+std::to_string(numberOfGenes);
        const char *cmd = command.c_str();
        system(cmd);
    }
    return 0;
}