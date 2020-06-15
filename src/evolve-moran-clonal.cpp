#include "ClonalPopulation.h" 
#include <tclap/CmdLine.h>
#include <unistd.h>
#include "rng.h"



int main(int argc, char *argv[])
{
    // these variables will hold the parameters input (or not) by the user
    int currentGen = 1;
    int maxGen = currentGen + 1;
    Encoding_Type outputEncoding = Encoding_Type::by_default;
    unsigned int targetPopSize = 1;
    int timeStep = 1;
    
    double mutationRate = 0; //per genome mutation rate

    std::string outDir;
    std::string startSnapFile;

    // Wrap everything in a try block
    // errors in input are caught and explained to user
    try { 
        // Define the command line object
        TCLAP::CmdLine cmd("SodaPop: a multi-scale model of molecular evolution", ' ', "v1.0");

        // Define value arguments
        TCLAP::ValueArg<int> maxArg("m","maxgen","Number of generations",false,10,"int");
        TCLAP::ValueArg<int> popArg("n","size","Initial population size",false,1,"int");
        TCLAP::ValueArg<int> dtArg("t","dt","Time interval for snapshots",false,1,"int");

        TCLAP::ValueArg<std::string> startArg("i","popdesc","Clonal population description file",false,"null","filename");        
        TCLAP::ValueArg<std::string> prefixArg("o","prefix","Prefix to be used for snapshot files",false,"sim","filename");

        TCLAP::ValueArg<double> mrate("u","mrate","Mutation rate",false,10,"double");

        TCLAP::ValueArg<std::string> inputArg("","sim-type","Define simulation type\n<s> (from selection coefficient, DMS or otherwise)\n<stability> (from DDG matrix or distribution)", false,"s","string");
        
        TCLAP::SwitchArg gammaArg("","gamma","Draw selection coefficients from gamma distribution", cmd, false);
        TCLAP::SwitchArg normalArg("","normal","Draw selection coefficients from normal distribution", cmd, false);
        TCLAP::ValueArg<double> alphaArg("","alpha","Alpha parameter of distribution\nGamma -> shape\nNormal -> mean",false,1,"double");
        TCLAP::ValueArg<double> betaArg("","beta","Beta parameter of distribution\nGamma -> scale\nNormal -> S.D.",false,1,"double");

        TCLAP::ValueArg<unsigned long> seedArg("", "seed", "Seed value for RNG.", false, 0, "unsigned int (64-bit)");
        
        
        // Add the arguments to the CmdLine object.
        cmd.add(maxArg);
        cmd.add(popArg);
        cmd.add(dtArg);
        
        cmd.add(startArg);
        cmd.add(prefixArg);
        cmd.add(mrate);

        cmd.add(alphaArg);
        cmd.add(betaArg);
        cmd.add(inputArg);

        // Parse the argv array.
        cmd.parse(argc, argv);

        // Get values from args. 
        maxGen = maxArg.getValue();
        targetPopSize = popArg.getValue();
        timeStep = dtArg.getValue();

        outDir = prefixArg.getValue();
        startSnapFile = startArg.getValue();
        
        mutationRate = mrate.getValue();

        if (seedArg.isSet())
            setRngSeed(seedArg.getValue());

        std::cout << "Begin ... " << std::endl;        
    /*
        Instantiate random numbers
    */
        
    }catch (TCLAP::ArgException &e){
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }

    /* OPEN SIMULATION LOG
    */

    std::ofstream CMDLOG;
    try{
        openCommandLog(CMDLOG, outDir, argv, argc);
    }catch (std::runtime_error &e) {}

    //std::cout << "Mutation rate:" << mutationRate << std::endl;
    std::cout << "Starting evolution ..." << std::endl;
    CMDLOG << "Starting evolution ..." << std::endl;
    
    
    //open snapshot
    ClonalPopulation CPin(startSnapFile);
    
    
    std::vector<int> C={1,20000,1000000,10000};
    std::vector<double> f={1,1.3,0.5,0.75};
    ClonalPopulation CP(C,f,mutationRate,0);
    
    /* Open clonal population snapshot
    */
    sprintf(buffer,"%s.conal.snap",outDir.c_str());
    std::ofstream outSnap = std::ofstream(buffer, std::ios::out | std::ios::trunc);
    if (outSnap.is_open()){
        // file was opened successfully
        std::cout << "-> Clonal snapshot was opened successfully ..." << std::endl;
    }
    else{
        // error opening file, throw exception
        throw std::runtime_error("Unable to open clonal snapshot.");
    }
    
    std::cout << "Saving population snapshot ..." << std::endl;
    
    CP.sortClonesByFitness();
    CP.dump(outSnap);
    outSnap.close();
    
    sprintf(buffer, "out/%s/MUTATION_LOG",outDir.c_str());

    
    
    //test
       
    
 //----------------
    // Moran process implementation    
    
//-----------------
    
    printProgress(currentGen/maxGen);
    std::cout << std::endl;

    CMDLOG.close();

    return 0;
};