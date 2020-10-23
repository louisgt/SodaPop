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
    Encoding_Type outputEncoding = Encoding_Type::by_default;
    unsigned int carryingCapacity = 1;
    int timeStep = 1;

    bool enableAnalysis = false;
    bool trackMutations = false;
    Init_Pop createPop = Init_Pop::from_snapFile;
    bool noMut = false;

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
        carryingCapacity = popArg.getValue();
        timeStep = dtArg.getValue();

        geneListFile = geneArg.getValue();
        outDir = prefixArg.getValue();
        startSnapFile = startArg.getValue();
        genesPath = libArg.getValue();

        Population::simType = stringToInput_Type(inputArg.getValue());

        if (seedArg.isSet())
            setRngSeed(seedArg.getValue());

        std::cout << "*** Begin ... " << std::endl;

        if (matrixArg.isSet()){
            matrixVec = matrixArg.getValue();
        }

        // This statement should be a switch
        // make enum for simType
        if (gammaArg.isSet()){
            Gene::initGamma(alphaArg.getValue(), betaArg.getValue());
        }
        else if (normalArg.isSet()){
            Gene::initNormal(alphaArg.getValue(), betaArg.getValue());
        }
        else{
            // default to gamma
            Gene::initGamma(alphaArg.getValue(), betaArg.getValue());
        }

        //call init here
        Population::initLandscape(fitArg.getValue(), matrixVec,geneListFile,genesPath);

        enableAnalysis = analysisArg.getValue();
        trackMutations = eventsArg.getValue();
        outputEncoding = intToEncoding_Type(seqArg.getValue());
        createPop = intToPop_Type(initArg.getValue());

    }catch (TCLAP::ArgException &e){
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;}

    /* general simulation initialization, can be put in a global method
    */
    if (Cell::ff_ == 7) {
        Population::noMut = true;
        std::cout << "*** Mutations are not active." << std::endl;
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

    // create the inoculum population with size = to carrying capacity
    Population inoculumPop(startFile, genesPath, carryingCapacity, createPop);

    // shuffle population to randomize packets
    inoculumPop.shuffle(g_rng);

    // create the microbiota population with size = to first packet
    Population microbiotaPop(carryingCapacity, inoculumPop);

    /* should be handled by method initializing the population
    */
    startFile.close();

    Total_Cell_Count = microbiotaPop.getSize();

    std::ofstream OUT;

    try{
        microbiotaPop.saveSnapshot(OUT,outDir,currentGen,outputEncoding);
    }catch (std::runtime_error &e) {}

    /* MAIN SIMULATION BLOCK, can be put in a method outside main
    */

    int currentSize = microbiotaPop.getSize();

    int targetSize = currentSize < 10000 ? currentSize*5 : currentSize*2;

    bool remaining = true;
    
    std::cout << "*** Starting evolution ..." << std::endl;
    CMDLOG << "*** Starting evolution ..." << std::endl;
    
    // // PSEUDO WRIGHT-FISHER PROCESS
    while (currentGen < maxGen){
            // check if inoculum is empty
            // if TRUE, set flag

            printProgress(currentGen*1.0/maxGen);

            microbiotaPop.divide(targetSize, carryingCapacity, MUTATIONLOG, !remaining);

            // update generation counter
            ++currentGen;

            // save population snapshot every timeStep generations
            if( (currentGen % timeStep) == 0){
                try{
                    microbiotaPop.saveSnapshot(OUT,outDir,currentGen,outputEncoding);
                }catch (std::runtime_error &e) {}
            }

            // check isEmpty flag
            if(remaining){
                // add new packet to microbiotaPop
                remaining = microbiotaPop.addPacket(carryingCapacity,inoculumPop);
            }
            //else{
            //    break;
            //    std::cout << "Stopping WF process and saving microbiota." << std::endl;
            //}

            currentSize = microbiotaPop.getSize();
            targetSize = currentSize < 10000 ? currentSize*5 : currentSize*2;

     }

    //microbiotaPop.reBarcode();
    //microbiotaPop.saveSnapshot(OUT,outDir,0,outputEncoding);

    printProgress(currentGen/maxGen);
    std::cout << std::endl;
    MUTATIONLOG.close();
    std::cout << "*** Done." << std::endl;
    CMDLOG << "Done." << std::endl;
    std::cout << "*** Total number of mutation events: " << microbiotaPop.getMutationCount() << std::endl;
    CMDLOG.close();

    // if the user toggled analysis, call shell script
    if (enableAnalysis){
        std::string script = "tools/barcodes.sh";
        std::string command = "/bin/bash "+script+" "+outDir+" "+std::to_string(maxGen)+" "+std::to_string(carryingCapacity) +" "+std::to_string(timeStep)+" "+std::to_string(outputEncoding)+" "+std::to_string(Population::numberOfGenes);
        const char *cmd = command.c_str();
        system(cmd);
    }
    return 0;
}
