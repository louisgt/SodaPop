#include "../src/global.h"

int main(int argc, char * argv[]) {
    if (argc != 3) {
        std::cerr << "sodasnap <snap-binary> <out-ascii>\n";
        exit(1);
    }

    int Total_Cell_Count, encoding;
    double frame_time;

    //open binary file 
    std::fstream IN(argv[1], std::ios:: in | std::ios::binary);
    if (!IN.is_open()) {
        std::cerr << "Binary file could not be opened.\n";
        exit(1);
    }

    std::fstream OUT(argv[2], std::ios::out);
    if (!OUT.is_open()) {
        std::cerr << "Ascii file could not be opened.\n";
        exit(1);
    }

    //Read header
    IN.read((char * )( & frame_time), sizeof(double));
    IN.read((char * )( & Total_Cell_Count), sizeof(int));
    IN.read((char * )( & encoding), sizeof(int));
    OUT << "Generation: " << frame_time << "\n" << "Population size: " << Total_Cell_Count << std::endl;

    //Read Cell array
    for (int i = 0; i < Total_Cell_Count; i++) {
        OUT << std::endl;
        switch(encoding)
            {
                case 0: read_Cell(IN, OUT,false);
                    break;
                case 1: qread_Cell(IN, OUT);
                    break;
                case 2: read_Cell(IN,OUT,true);
                    break;
                case 3: //to be properly implemented
                        //read_Parent(IN,OUT);  
                    break;
            }
    }
    IN.close();
    OUT.close();
}