#ifndef CELL_H
#define CELL_H
#include "Gene.h"

/*SodaPop
Copyright (C) 2017 Louis Gauthier

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

class Cell 
{
    protected:
        // organism barcode
        std::string barcode_;
        // organism ID
        int ID_;
        // initial mutation rate
        double o_mrate_;
        // current mutation rate
        double c_mrate_;	
        			
    public:
        //Array of genes
        std::vector<Gene> Gene_arr_;
        //Cummulative sum of gene lengths (i.e. genome size)
        VectInt Gene_L_;
        Cell();
        Cell(std::fstream&);			    
        Cell(std::fstream&, const std::string&);  

        virtual void UpdateRates() = 0;
        virtual void dump(std::fstream&, int) = 0;
        virtual void PrintCell(int) = 0;       
             
        //getters
        const int ID(){return ID_;}
        const double mrate(){return c_mrate_;}
        const int gene_count(){return Gene_arr_.size();}
        const int genome_size(){return Gene_L_.back();}
        const std::string barcode(){return barcode_;}   

        //setters
        void change_ID(int a){ ID_ = a;}
        void ch_barcode(std::string s){barcode_ = s;}
        int total_mutations(const int&);
        void FillGene_L();
};

Cell::Cell():
    barcode_(getBarcode()),
    ID_(0),
    o_mrate_(0),
    c_mrate_(0)
{
  Gene_L_.reserve(GENECOUNTMAX);
  Gene_arr_.reserve(GENECOUNTMAX);
}

// Construct from cell file
Cell::Cell(std::fstream& cell_in)
{
    char buffer[140];
    Gene_L_.reserve(GENECOUNTMAX);
    Gene_arr_.reserve(GENECOUNTMAX);
    ch_barcode(getBarcode());
    std::string line;
    std::string genesPath = "files/genes/";
    while (!cell_in.eof())
    {
        getline(cell_in,line);
        std::string word;
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;
        if ( word == "genes_path" ){
          iss >> word; 
          genesPath = word.c_str();
        }
        if ( word == "org_id" ){
          iss >> word; 
          ID_ = atoi(word.c_str());
        }
        else if ( word == "mrate" ){
          iss >> word; 
          o_mrate_ = atof(word.c_str());
          c_mrate_ = atof(word.c_str());
        }       
        else if ( word == "G" ){
        //reading gene files; 
        //concentration and stability from gene file take precedence
              iss >> word;
              //open gene file
              sprintf(buffer,"%s%s.gene",genesPath.c_str(),word.c_str());
              std::fstream temp (buffer);
              if (!temp.is_open()){
                std::cerr << "File could not be open: " << buffer << std::endl;
                exit(2);
              }
              Gene A(temp);
              Gene_arr_.push_back(A);
              std::cout << "Inserted: "<< word << std::endl;
              
              //Check if gene is correctly inserted
              std::vector<Gene>::iterator i = Gene_arr_.end();
              i--;
              std::cout << (*i).nseq() << std::endl;        
              temp.close();
        }
    }
}

// Constructs from a unit cell stored in binary 
Cell::Cell(std::fstream& IN, const std::string& genesPath)
{
    Gene_L_.reserve(GENECOUNTMAX);
    Gene_arr_.reserve(GENECOUNTMAX);

    char buffer[140];
    int cell_id, cell_index, gene_size;
    double m,m0;
   
    IN.read((char*)(&cell_index),sizeof(int));  
    IN.read((char*)(&cell_id),sizeof(int)); ID_ = cell_id;

    //read barcode
    int l;
    IN.read((char*)&l, sizeof(int));
    //construct vector container with nl elements
    std::vector<char> buf(l);
    IN.read(&buf[0], l);
    ch_barcode(std::string().assign(buf.begin(), buf.end()));

    IN.read((char*)(&m0),sizeof(double)); o_mrate_ = m0;
    IN.read((char*)(&m),sizeof(double));  c_mrate_ = m;
    IN.read((char*)(&gene_size),sizeof(int));
    
    //read gene info
    for(int j=0; j<gene_size; j++){
        double e, c, dg, f;
        int gene_nid, Ns, Na;
        std::string DNAsequence;
     
        IN.read((char*)(&gene_nid),sizeof(int));   
        IN.read((char*)(&e),sizeof(double));
        IN.read((char*)(&c),sizeof(double));
        IN.read((char*)(&dg),sizeof(double));
        IN.read((char*)(&f),sizeof(double));
        IN.read((char*)(&Ns),sizeof(int));
        IN.read((char*)(&Na),sizeof(int));

        //read DNA sequence
        int nl;
        IN.read((char*)&nl, sizeof(int));
        //construct vector container with nl elements
        std::vector<char> buf(nl);
        IN.read(&buf[0], nl);
        DNAsequence.assign(buf.begin(), buf.end());
     
        sprintf(buffer,"%s%d.gene",genesPath.c_str(),gene_nid);
        std::fstream temp (buffer);
        if (!temp.is_open()){
            std::cerr << "ERROR: Cannot open gene file " << buffer << std::endl;
            exit(2);
        }
        Gene G(temp);   
        //update gene information
        G.conc = c;
        dg = exp(-dg/kT);
        G.ch_dg(dg);
        G.ch_f(f);
        G.Update_Sequences(DNAsequence);
        G.ch_Na(Na);
        G.ch_Ns(Ns); 
        Gene_arr_.push_back(G);
    }
}

// Initialize the cummulative gene length array
void Cell::FillGene_L()
{
    int sum = 0;
    std::vector<Gene>::iterator i;
    for(i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        sum+= (*i).length();
        Gene_L_.push_back(sum);
    }
}

// Return total mutation count
// spec:
//  - 0, Ns+Na
//  - 1, Ns
//  - 2. Na
int Cell::total_mutations(const int& spec)
{
    assert( (spec < 3) &&  (spec >= 0));

    int sa = 0;
    int s = 0;
    int a = 0;

    for(std::vector<Gene>::iterator i = Gene_arr_.begin(); i != Gene_arr_.end(); ++i){
        int Ns = i->Ns();
        int Na = i->Na();
        s += Ns;
        a += Na;
        sa += (Ns+Na);
    }

    if(spec == 0) return sa;
    else if (spec == 1) return s;
    else return a;
}

#endif


