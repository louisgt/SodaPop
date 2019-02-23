#include "../src/Cell.h"
#include "../src/Gene.h"

/*
DESCRIPTION: Converts population summary to a snap file.
*/

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "sodasumm <population summary> [ 0-full | 1-single cell | 2-randomize barcodes | 3-introduce variation]\n";
        exit(1);
    }

    int flag = atoi(argv[2]);
    assert((flag == 0) | (flag == 1) | (flag == 2) | (flag == 3));

    char buffer[200];
    std::vector < Cell > Cell_arr;

    // open stream to read population summary
    std::ifstream popf(argv[1]);
    if (!popf.is_open()) {
        std::cerr << "File could not be opened";
        exit(1);
    }

    if(flag==3){
        Gene::initGamma(1.5, 0.0003);
    }

    double limit = 0.90;

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
                        int nClusters = (count + 10 - 1) / 10;
                        int remainder = count % 10;
                        for (int k = 0; k < nClusters - 1; k++) {
                            //edit fitness here (draw selection coefficient for clusters of 10 cells)
                            //fetch selection coefficient
                            double s = Gene::RandomGamma();
                            if (randomNumber()<limit) s *= -1;
                            //change fitness accordingly
                            double nf = A.fitness() + s;
                            //std::cout << nf << std::endl;
                            A.ch_Fitness(nf);
                            for (int i = 0; i < 10; i++) {
                                Cell_arr.push_back(A);
                            }
                        }
                        if(remainder){
                            //edit fitness here (draw selection coefficient for clusters of 10 cells)
                            //fetch selection coefficient
                            double s = Gene::RandomGamma();
                            if (randomNumber()<limit) s *= -1;
                            //change fitness accordingly
                            double nf = A.fitness() + s;
                            //std::cout << nf << std::endl;
                            A.ch_Fitness(nf);
                            for (int i = 0; i < remainder; i++) {
                                Cell_arr.push_back(A);
                            }
                        }
                        
                break;
            }
            temp.close();
        }
    }
    popf.close();

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