// Gene.cpp
#include "Cell.h"

Cell::Cell():
    barcode_(getBarcode()),
    ID_(0),
    parent_(0),
    o_mrate_(0),
    c_mrate_(0),
    fitness_(0)
    {
        Gene_L_.reserve(GENECOUNTMAX);
        Gene_arr_.reserve(GENECOUNTMAX);
    }

// Construct from cell file
Cell::Cell(std::fstream & cell_in) {
    char buffer[140];
    Gene_L_.reserve(GENECOUNTMAX);
    Gene_arr_.reserve(GENECOUNTMAX);
    ch_barcode(getBarcode());
    setParent(0);
    std::string line;
    // default genesPath, can be changed by user in command-line
    std::string genesPath = "files/genes/";
    while (!cell_in.eof()) {
        getline(cell_in, line);
        std::string word;
        std::istringstream iss(line, std::istringstream:: in );
        iss >> word;
        if (word == "genes_path") {
            iss >> word;
            genesPath = word.c_str();
        }
        if (word == "org_id") {
            iss >> word;
            ID_ = atoi(word.c_str());
        } else if (word == "mrate") {
            iss >> word;
            o_mrate_ = atof(word.c_str());
            c_mrate_ = atof(word.c_str());
        }
        else if (word == "fitness") {
            iss >> word;
            fitness_ = atof(word.c_str());
        } else if (word == "G") {
            //reading gene files; 
            //concentration and stability from gene file take precedence
            iss >> word;
            //open gene file
            sprintf(buffer, "%s%s.gene", genesPath.c_str(), word.c_str());
            std::fstream gene_data(buffer);
            if (!gene_data.is_open()) {
                std::cerr << "File could not be open: " << buffer << std::endl;
                exit(2);
            }
            Gene A(gene_data,this);
            Gene_arr_.push_back(A);

            //Check if gene is correctly inserted
            std::vector <Gene>::iterator i = Gene_arr_.end();
            i--;
            // if gene fitness is null, assign a randomly fit value
            if((*i).f()==0){
                (*i).ch_f(0.95+randomNumber()*0.02);
            }
            gene_data.close();
        }
    }
}

// Constructs from a unit cell stored in binary 
Cell::Cell(std::fstream & IN,
    const std::string & genesPath) {
    Gene_L_.reserve(GENECOUNTMAX);
    Gene_arr_.reserve(GENECOUNTMAX);

    char buffer[140];
    int cell_id, cell_index, gene_size;
    double m, f;

    IN.read((char*)(&cell_index), sizeof(int));
    IN.read((char*)(&cell_id), sizeof(int));
    ID_ = cell_id;
    parent_ = 0;

    //read barcode
    int l;
    IN.read((char*) & l, sizeof(int));
    //construct vector container with nl elements
    std::vector<char> buf(l);
    IN.read(&buf[0], l);
    ch_barcode(std::string().assign(buf.begin(), buf.end()));

    IN.read((char*)(&f), sizeof(double));
    fitness_ = f;
    IN.read((char*)(&m), sizeof(double));
    c_mrate_ = m;
    IN.read((char*)(&gene_size), sizeof(int));

    //read gene info
    for (int j = 0; j < gene_size; j++) {
        double e, c, dg, f, eff;
        int gene_nid, Ns, Na;
        std::string DNAsequence;

        IN.read((char*)(&gene_nid), sizeof(int));
        IN.read((char*)(&e), sizeof(double));
        IN.read((char*)(&c), sizeof(double));
        IN.read((char*)(&eff), sizeof(double));
        IN.read((char*)(&dg), sizeof(double));
        IN.read((char*)(&f), sizeof(double));
        IN.read((char*)(&Ns), sizeof(int));
        IN.read((char*)(&Na), sizeof(int));

        //read DNA sequence
        int nl;
        IN.read((char*)&nl, sizeof(int));
        //construct vector container with nl elements
        std::vector<char> buf(nl);
        IN.read(& buf[0], nl);
        DNAsequence.assign(buf.begin(), buf.end());

        sprintf(buffer, "%s%d.gene", genesPath.c_str(), gene_nid);
        std::fstream gene_data(buffer);
        if (!gene_data.is_open()) {
            std::cerr << "ERROR: Cannot open gene file " << buffer << std::endl;
            exit(2);
        }
        Gene G(gene_data,this);
        //update gene information
        G.ch_conc(c);
        dg = exp(-dg / kT);
        G.ch_dg(dg);
        G.ch_eff(eff);
        G.ch_f(f);
        G.Update_Sequences(DNAsequence);
        G.ch_Na(Na);
        G.ch_Ns(Ns);
        Gene_arr_.push_back(G);
    }
}

void Cell::linkGenes()
{
    for(auto gene_it = this->Gene_arr_.begin(); gene_it != this->Gene_arr_.end(); gene_it++) {
        gene_it->setCell(this);
    }
}

// // copy constructor
// Cell::Cell(const Cell& C)
// {
//     barcode_ = C.barcode_;
//     ID_ = C.ID_;
//     parent_ = C.parent_;
//     o_mrate_ = C.o_mrate_;
//     c_mrate_ = C.c_mrate_;
//     fitness_ = C.fitness_;
//     Gene_L_ = C.Gene_L_;
//     std::vector < Gene > ::iterator i;
//     for (auto cell_it = C.Gene_arr_.begin(); cell_it != C.Gene_arr_.end(); cell_it++) {
//         Gene A(*cell_it,this);
//         Gene_arr_.push_back(A);
//     }
// }

double Cell::fitness() const
{
    return fitness_;
}

// Initialize the cummulative gene length array
void Cell::FillGene_L() {
    int sum = 0;
    std::vector < Gene > ::iterator i;
    //for (auto cell_it = Gene_arr_.begin(); cell_it != Gene_arr_.end(); ++cell_it) {
    for (const auto& gene : Gene_arr_) {
        sum += gene.length();
        Gene_L_.push_back(sum);
    }
}

// Return total mutation count
// spec:
//  - 0, Ns+Na
//  - 1, Ns
//  - 2. Na
int Cell::total_mutations(const int & spec) {
    assert((spec < 3) && (spec >= 0));

    int sa = 0;
    int s = 0;
    int a = 0;

    //for (auto cell_it = Gene_arr_.begin(); cell_it != Gene_arr_.end(); ++cell_it) {
    for (const auto& gene : Gene_arr_) {
        int Ns = gene.Ns();
        int Na = gene.Na();
        s += Ns;
        a += Na;
        sa += (Ns + Na);
    }

    if (spec == 0) return sa;
    else if (spec == 1) return s;
    else return a;
}