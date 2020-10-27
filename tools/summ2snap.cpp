#include "../src/Cell.h"
#include "../src/Gene.h"

/*
DESCRIPTION: Converts population summary to a snap file.
*/


template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g){
	std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
	std::advance(start, dis(g));
	return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
	static std::random_device rd;
	static std::mt19937 gen(rd());
	return select_randomly(start, end, gen);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "sodasumm <population summary> [ 0-full | 1-single cell | 2-randomize barcodes | 3-introduce variation]\n";
        exit(1);
    }

    int flag = atoi(argv[2]);
    assert((flag == 0) | (flag == 1) | (flag == 2) | (flag == 3));

    char buffer[200];
    std::vector <Cell> Cell_arr;

    // open stream to read population summary
    std::ifstream popf(argv[1]);
    if (!popf.is_open()) {
        std::cerr << "File could not be opened";
        exit(1);
    }

    double s(0);

    int blockSize = 100;

    int blockCounter = 0;

    double limit = 0.90;

    double mu = 5e-5;

    double genomeSize = 4600000;

    int popSize = 10000000;

    int nMu = mu*popSize;

    std::cout << "nMu is " << nMu << std::endl;

    std::default_random_engine generator;
    std::geometric_distribution<int> distribution(0.01);

    if(flag==0){
        Gene::initExponential(100);
        s = Gene::RandomExponential();
	//blockSize = distribution(generator) + 10;
    }

    double nf = 1 + s;
    
    std::string line;
    Cell_arr.reserve(maxPopSize);
    while (!popf.eof()) {
        getline(popf, line);
        std::string word;
        std::istringstream iss(line, std::istringstream:: in );
        iss >> word;
        if (word == "C") {
            iss >> word; //cell count
            int count = atoi(word.c_str());
            iss >> word; //cell files
            std::ifstream temp(word.c_str()); //convert std::string to char
            if (!temp.is_open()) {
                std::cerr << "File could not be open: " << word << std::endl;
                exit(1);
            }
            Cell A(temp);
            switch(flag)
            {
                case 0: A.ch_barcode(getBarcode());
			//s = Gene::RandomExponential();
			//nf = 1 + s;
                        //A.ch_Fitness(nf);   
                        for (int i = 0; i < count; i++) {
                            Cell_arr.push_back(A);
                        }
                break;

                case 1: A.ch_barcode(getBarcode());
                        Cell_arr.push_back(A);
                        temp.close();
                break;

                case 2: for (int i = 0; i < count; i++) {
                            Cell_arr.push_back(A);
                            Cell_arr.back().ch_barcode(getBarcode());
                        }
                break;

                case 3: A.ch_barcode(getBarcode());
                        //change fitness accordingly
                        A.ch_Fitness(nf);
                        for (int k = 0; k < count; ++k) {
                            //draw selection coefficient
                            //std::cout << nf << std::endl;
                            Cell_arr.push_back(A);
			    blockCounter +=1;
			    if (blockCounter == blockSize){
				std::cout << "Changing block... " << std::endl;
				//reset counter
                            	blockCounter = 0;
                            	// new selection coefficient
                            	s = Gene::RandomExponential();
			    	blockSize = distribution(generator) + 10;
				nf = 1 + s;
                        	A.ch_Fitness(nf);
				std::cout << "New fitness for block " << nf << std::endl;
			    }
                        }

                break;
            }
            temp.close();
        }
    }
    popf.close();

    //randomly select cells for beneficial mutations
    for (int k = 0; k < nMu; ++k) {
	s = Gene::RandomExponential();
	nf = 1 + s;
	(*select_randomly(Cell_arr.begin(), Cell_arr.end())).ch_Fitness(nf);
    }

    for (auto& cell : Cell_arr) {
	cell.propagateFitness();
        cell.UpdateRates(); 
    }

    //population snapshot
    int Total_Cell_Count = (int)(Cell_arr.size());

    std::string popname = argv[1];

    std::string rawname = popname.substr(0, popname.find_last_of("."));
    sprintf(buffer, "%s.snap", rawname.c_str());

    // Open stream to write snapshot
    std::ofstream OUT(buffer, std::ios::out);
    if (!OUT.is_open()) {
        std::cerr << "Snapshot file could not be opened";
        exit(1);
    }

    double frame_time = 0;
    int encoding = 0;
    OUT.write((char*)(&frame_time), sizeof(double));
    OUT.write((char*)(&Total_Cell_Count), sizeof(int));
    OUT.write((char*)(&encoding), sizeof(int));

    int count = 0;
    for (auto cell_it = Cell_arr.begin(); cell_it != Cell_arr.end(); ++cell_it) {
        ++count;
        cell_it->dump(OUT, count);
    }
    OUT.close();
}
