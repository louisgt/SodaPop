#include "global.h"

VectStr PrimordialAASeq;
double fold_DG = 0;
double bind_DG = 0;

int Total_Cell_Count = 0;
int dummy = 0;
double frame_time = 0;
std::string outPath = "";
char buffer[200];

double w_sum = 0;

double matrix[gene_number][res_number][20];
double matrix_supp[gene_number][res_number][20];

/******* GENETIC CODE MAPPINGS *******/
// these const mappings are hard-coded and populated at compile-time

struct codon_to_num{
    static std::map<std::string,int> create_map()
        {
          std::map<std::string,int> m;
          m["GCA"] = 1;
          m["GCC"] = 1;
          m["GCG"] = 1;
          m["GCT"] = 1;
          m["TGC"] = 2;
          m["TGT"] = 2;
          m["GAC"] = 3;
          m["GAT"] = 3;
          m["GAA"] = 4;
          m["GAG"] = 4;
          m["TTC"] = 5;
          m["TTT"] = 5;
          m["GGA"] = 6;
          m["GGC"] = 6;
          m["GGG"] = 6;
          m["GGT"] = 6;
          m["CAC"] = 7;
          m["CAT"] = 7;
          m["ATA"] = 8;
          m["ATC"] = 8;
          m["ATT"] = 8;
          m["AAA"] = 9;
          m["AAG"] = 9;
          m["CTA"] = 10;
          m["CTC"] = 10;
          m["CTG"] = 10;
          m["CTT"] = 10;
          m["TTA"] = 10;
          m["TTG"] = 10;
          m["ATG"] = 11;
          m["AAC"] = 12;
          m["AAT"] = 12;
          m["CCA"] = 13;
          m["CCC"] = 13;
          m["CCG"] = 13;
          m["CCT"] = 13;
          m["CAA"] = 14;
          m["CAG"] = 14;
          m["AGA"] = 15;
          m["AGG"] = 15;
          m["CGA"] = 15;
          m["CGC"] = 15;
          m["CGG"] = 15;
          m["CGT"] = 15;
          m["AGC"] = 16;
          m["AGT"] = 16;
          m["TCA"] = 16;
          m["TCC"] = 16;
          m["TCG"] = 16;
          m["TCT"] = 16;
          m["ACA"] = 17;
          m["ACC"] = 17;
          m["ACG"] = 17;
          m["ACT"] = 17;
          m["GTA"] = 18;
          m["GTC"] = 18;
          m["GTG"] = 18;
          m["GTT"] = 18;
          m["TGG"] = 19;
          m["TAC"] = 20;
          m["TAT"] = 20;
          m["TAA"] = 21;
          m["TAG"] = 21;
          m["TGA"] = 21;
          return m;
        }
    static const std::map<std::string,int> cnum;
};

struct codon_to_prot{
    static std::map<std::string,char> create_map()
        {
          std::map<std::string,char> m;
          m["AAA"] = 'K';
          m["AAC"] = 'N';
          m["AAG"] = 'K';
          m["AAT"] = 'N';
          m["ACA"] = 'T';
          m["ACC"] = 'T';
          m["ACG"] = 'T';
          m["ACT"] = 'T';
          m["AGA"] = 'R';
          m["AGC"] = 'S';
          m["AGG"] = 'R';
          m["AGT"] = 'S';
          m["ATA"] = 'I';
          m["ATC"] = 'I';
          m["ATG"] = 'M';
          m["ATT"] = 'I';
          m["CAA"] = 'Q';
          m["CAC"] = 'H';
          m["CAG"] = 'Q';
          m["CAT"] = 'H';
          m["CCA"] = 'P';
          m["CCC"] = 'P';
          m["CCG"] = 'P';
          m["CCT"] = 'P';
          m["CGA"] = 'R';
          m["CGC"] = 'R';
          m["CGG"] = 'R';
          m["CGT"] = 'R';
          m["CTA"] = 'L';
          m["CTC"] = 'L';
          m["CTG"] = 'L';
          m["CTT"] = 'L';
          m["GAA"] = 'E';
          m["GAC"] = 'D';
          m["GAG"] = 'E';
          m["GAT"] = 'D';
          m["GCA"] = 'A';
          m["GCC"] = 'A';
          m["GCG"] = 'A';
          m["GCT"] = 'A';
          m["GGA"] = 'G';
          m["GGC"] = 'G';
          m["GGG"] = 'G';
          m["GGT"] = 'G';
          m["GTA"] = 'V';
          m["GTC"] = 'V';
          m["GTG"] = 'V';
          m["GTT"] = 'V';
          //STOP
          m["TAA"] = 'X';
          m["TAC"] = 'Y';
          //STOP
          m["TAG"] = 'X';
          m["TAT"] = 'Y';
          m["TCA"] = 'S';
          m["TCC"] = 'S';
          m["TCG"] = 'S';
          m["TCT"] = 'S';
          //STOP
          m["TGA"] = 'X';
          m["TGC"] = 'C';
          m["TGG"] = 'W';
          m["TGT"] = 'C';
          m["TTA"] = 'L';
          m["TTC"] = 'F';
          m["TTG"] = 'L';
          m["TTT"] = 'F';
          return m;
        }
    static const std::map<std::string,char> cprot;
};

struct prot_to_num{
    static std::map<char,int> create_map()
        {
          std::map<char,int> m;
          m['A'] = 1;
          m['C'] = 2;
          m['D'] = 3;
          m['E'] = 4;
          m['F'] = 5;
          m['G'] = 6;
          m['H'] = 7;
          m['I'] = 8;
          m['K'] = 9;
          m['L'] = 10;
          m['M'] = 11;
          m['N'] = 12;
          m['P'] = 13;
          m['Q'] = 14;
          m['R'] = 15;
          m['S'] = 16;
          m['T'] = 17;
          m['V'] = 18;
          m['W'] = 19;
          m['Y'] = 20;
          m['X'] = 21;//STOP
          return m;
        }
    static std::map<char,int> const pnum;
};

// populate genetic code mappings
std::map<std::string,int> const codon_to_num::cnum = codon_to_num::create_map();
std::map<char,int> const prot_to_num::pnum = prot_to_num::create_map();
std::map<std::string,char> const codon_to_prot::cprot = codon_to_prot::create_map();

Encoding_Type intToEncoding_Type(int type){
    switch (type){
        case 0:
            return Encoding_Type::by_default;
        case 1:
            return Encoding_Type::no_sequence;
        case 2:
            return Encoding_Type::full;
        case 3:
            return Encoding_Type::other;
        default:
            return Encoding_Type::by_default;
    }
}

Init_Pop intToPop_Type(int type){
    switch (type){
        case 0:
            return Init_Pop::from_snapFile;
        case 1:
            return Init_Pop::from_cellFile;
        default:
            return Init_Pop::from_snapFile;
    }
}

void openCommandLog(std::ofstream& cmdlog, std::string outDir, char *argv[], int argc){
    sprintf(buffer,"out/%s/command.log",outDir.c_str());
    cmdlog = std::ofstream(buffer, std::ios::out | std::ios::trunc);

    if (cmdlog.is_open()){
        // file was opened successfully
        std::cout << "-> Command log was opened successfully ..." << std::endl;
        std::string args;
        std::for_each( argv + 1, argv + argc , [&]( const char* c_str ){ args += std::string ( c_str ) + " "; } );
        cmdlog << "sodapop " << args << std::endl;
        cmdlog << std::endl;
    }
    else{
        // error opening file, throw exception
        throw std::runtime_error("Unable to open command log.");
    }
}

void openMutationLog(std::ofstream& mutlog, std::string outDir){
    sprintf(buffer, "out/%s/MUTATION_LOG",outDir.c_str());
    mutlog = std::ofstream(buffer, std::ios::out | std::ios::trunc);

    if (mutlog.is_open()){
        // file was opened successfully
        std::cout << "-> Mutation log was opened successfully ..." << std::endl;
    }
    else{
        // error opening file, throw exception
        throw std::runtime_error("Unable to open mutation log.");
    }
}

void openStartingPop(std::string filePath, std::ifstream& fileStream){
    std::cout << "Opening starting population snapshot ..." << std::endl;
    fileStream = std::ifstream(filePath.c_str(),std::ios::in|std::ios::binary);
    if (fileStream.is_open()){
        // file was opened successfully
        std::cout << "-> File was opened successfully ..." << std::endl;
    }
    else{
        // error opening file, throw exception
        throw std::runtime_error("Unable to open starting population snapshot.");
    }
}

void readSnapshotHeader(std::ifstream& snapshot)
{
    snapshot.read((char*)(&frame_time),sizeof(double));
    //read number of cells in file
    snapshot.read((char*)(&Total_Cell_Count),sizeof(int));
    //read file outputEncoding
    snapshot.read((char*)(&dummy),sizeof(int));
}

void createOutputDir(std::string dirName){
    sprintf(buffer,"out/%s/snapshots",dirName.c_str());
    outPath = buffer;
    std::cout << "Creating directory " << outPath << " ... " << (makePath(outPath) ? "OK" : "failed") << std::endl;
}

/******* MAPPING FUNCTIONS *******/

// input:  3 nucleotide codon as string
// output: index number of codon
int GetIndexFromCodon(std::string in_codon)
{
    std::map <std::string, int>::const_iterator it;
    it = codon_to_num::cnum.find(in_codon);  
    if (it == codon_to_num::cnum.end()){
        std::cerr << "Invalid codon: "<< in_codon << std::endl;
        std::cerr << "Nucleotide sequence must not contain STOP codons."<< std::endl;
        exit(2);
    }
    return it->second;
}

// input:  3 nucleotide codon sequence as string
// output: amino acid sequence as string
std::string GetProtFromNuc(std::string in_seq)
{
    int ln = in_seq.length();
    if ((ln % 3) != 0){
        std::cerr << "Invalid length for nucleotide sequence: " << ln << std::endl;
        std::cerr << "Nucleotide sequence length must be divisible by 3." << std::endl;
        exit(2);
    }
    int la=ln/3;   
    std::string AA="";
    for (int i=0; i<la;++i){
        std::string temp=in_seq.substr(i*3,3);
        
        //check for valid code
        std::map <std::string, char> :: const_iterator Iter;
        Iter = codon_to_prot::cprot.find(temp);
        if (Iter == codon_to_prot::cprot.end()){
          std::cerr << "Invalid codon: "<< temp << std::endl;
          std::cerr << "Nucleotide sequence must not contain STOP codons."<< std::endl;
          exit(2);
        }   
        AA.push_back(codon_to_prot::cprot.at(temp));
    }
    return AA;
}

// input:  amino acid letter as string
// output: index number of amino acid
int GetIndexFromAA(char aa)
{  
    //check for valid amino acid
    std::map <char, int> :: const_iterator it;
    it = prot_to_num::pnum.find(aa);
    if (it == prot_to_num::pnum.end()){
        std::cerr << "Invalid amino acid: "<< aa << std::endl;
        exit(2);
    }
    return it->second;
}

// returns the next or second to next bp from a given nucleotide
const char AdjacentBP(char a, int j){
 
    if ( j > 2 ){
        std::cerr << "Invalid bp distance. Error in AdjacentBP(). "<< j << std::endl;
        exit(2);
    }

    std::map <char, int> A;
    A['A'] = 0;
    A['T'] = 1;
    A['G'] = 2;
    A['C'] = 3;

    std::map <int, char> B;
    B[0] = 'A';
    B[1] = 'T';
    B[2] = 'G';
    B[3] = 'C';

    int x = (A[a] + j + 1) % 4; 	//get (j+1) nucleotide from a
    return B[x];
}

//Checks if codon b when mutated resulted in a STOP codon a.
//If b is a STOP codon, it is replaced by 
//Input: 	a, new codon
//		    b, old codon
//		    i, mutation site in codon (< 3)
std::string n3_to_n3(std::string a, std::string b, int i){
  //double r = RandomNumber();
  double r = randomNumber();
  double l;

  if ( a == "TAA"){
      switch (i)
      {
          case 0:
              l =  (2*r);

              if (b == "CAA" ){
                if  ( l<1 )	
                  a = "GAA";
                else
                  a = "AAA";
              }else if (b == "GAA" ){
                if  ( l<1 )
                  a = "CAA";
                else
                  a = "AAA";
              }else if (b == "AAA" ){
                if  ( l<1 )
                  a = "GAA";
                else
                  a = "CAA";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

          case 1:
              if 	(b == "TTA")
                a = "TCA";
              else if (b == "TCA")
                a = "TTA";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

            break;

          case 2:
              if 	(b == "TAT")
                a = "TAC";
              else if (b == "TAC")
                a = "TAT";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

            break;

           default:
              std::cerr << "ERROR in translating STOP codon to NON-STOP codon." << std::endl;    
      }
  }
  else if (a == "TAG")
  {
      switch (i)
      {
          case 0:
              l = (2*r);

              if (b == "AAG" ){
                if  ( l<1 )
                  a = "GAG";
                else
                  a = "CAG";
              }else if (b == "GAG" ){
                if  ( l<1 )
                  a = "AAG";
                else
                  a = "CAG";
              }else if (b == "CAG" ){
                if  ( l<1 )
                  a = "GAG";
                else
                  a = "AAG";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

          case 1:
              l = (2*r);

              if (b == "TTG" ){
                if  ( l<1 )	a = "TGG";
                else 		a = "TGC";
              }else if (b == "TGG" ){
                if  ( l<1 )	a = "TTG";
                else 		a = "TGC";
              }else if (b == "TCG" ){
                 if  ( l<1 )	a = "TTG";
                else 		a = "TGG";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

          case 2:
              if 	(b == "TAT") a = "TAC";
              else if (b == "TAC") a = "TAT";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

              break;

             default:
                std::cerr << "ERROR in translating STOP codon to NON-STOP codon." << std::endl;    
       }
   }
   else if (a == "TGA")
   {
      switch (i)
      {
          case 0:
              l =  (2*r);

              if (b == "AGA" ){
                if  ( l<1 )	a = "GGA";
                else 		a = "CGA";
              }else if (b == "GGA" ){
                if  ( l<1 )	a = "AGA";
                else 		a = "CGA";
              }else if (b == "CGA" ){
                 if  ( l<1 )	a = "AGA";
                else 		a = "GGA";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

         case 1:
              if 	(b == "TTA") a = "TCA";
              else if (b == "TCA") a = "TTA";
              else {
                std::cerr << "Invalid starting codon in n3_to_n3()." << std::endl;
                exit(2); 
              }

              break;

          case 2:
              l =  (2*r);

              if (b == "TGT" ){
                if  ( l<1 )	a = "TGG";
                else 		a = "TGC";
              }else if (b == "TGG" ){
                if  ( l<1 )	a = "TGT";
                else 		a = "TGC";
              }else if (b == "TGC" ){
                 if  ( l<1 )	a = "TGT";
                else 		a = "TGG";
              }else{
                std::cerr << "Invalid mutation. n3_to_n3()." << std::endl;
                exit(2);
              }

           break;

           default:
              std::cerr << "ERROR in translating STOP codon to NON-STOP codon." << std::endl;    
      }  
  }
  return a;
}

// generates a random, 15nt barcode
std::string getBarcode()
{
    char seq [16];
    for (int i = 0; i < 15; ++i){
        if (randomNumber()<0.5){
            if (randomNumber()<0.5){
                seq[i] = 'G';
            }
            else seq[i] = 'C';
        }
        else{
            if (randomNumber()<0.5){
                seq[i] = 'A';
            }
            else seq[i] = 'T';
        }
    }
    seq[15] = '\0';
    return std::string(seq);
}

// initializes the 3D matrix for DDG values
void InitMatrix()
{
    std::cout << "Initializing matrix ..." << std::endl;
    for (int i = 0; i != gene_number; ++i)
      for (int j = 0; j != res_number; ++j)
        for (int k = 0; k != 20; ++k){
            matrix[i][j][k] = 1; //i.e., exp(-0/kT) = 1
            matrix_supp[i][j][k] = 1; // maybe it would make more sense here to set it to inf
        }    
}

// extracts values from the DDG file and stores them in the matrix
double ExtractDDGMatrix(std::string filepath, Matrix_Type m)
{
    std::ifstream temp(filepath);
    if (!temp.is_open()){
        std::cerr << "File could not be opened: "<< filepath << std::endl;
        exit(2);
    }
    std::cout << "Extracting DDG matrix ..." << std::endl;
    std::string line;
    int gene_num = 0;
    double sum = 0;
    int idx = 0;
    while(!temp.eof()){
        std::string word;
        getline(temp,line);
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;
        if ( word == "DDG"){
            iss >> word;
            //residue index
            int i = atoi(word.c_str());
            for (int j = 0; iss>>word; ++j){
                //extract DDG values
                double x = atof(word.c_str());
                sum +=x;
                ++idx;
                //if(m){
                    //matrix_supp[gene_num][i-1][j] = exp(-x/kT);
                //}
                //else{
                    matrix[gene_num][i-1][j] = exp(-x/kT);
                //}
            }
        }
        else if ( word == "rCat"){
            iss >> word;
            //residue index
            int i = atoi(word.c_str());
            for (int j = 0; iss>>word; ++j){
                //extract DDG values
                double x = atof(word.c_str());
                sum +=x;
                ++idx;
                matrix_supp[gene_num][i-1][j] = x;
            }
        }
        else if ( word == "Gene_NUM"){
            iss >> word;
            gene_num = atoi(word.c_str());
        }
    }
    temp.close();
    return sum/idx;
}

void ExtractDMSMatrix(std::string filepath)
{
    std::ifstream temp(filepath);
    if (!temp.is_open()){
        std::cerr << "File could not be open: "<< filepath << std::endl;
        exit(2);
    }
    std::cout << "Extracting DMS matrix ..." << std::endl;
    std::string line;
    int gene_num = 0;
    while (!temp.eof()){
        std::string word;
        getline(temp,line);
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;
        if ( word == "DMS"){
            iss >> word;
            //residue index
            int i = atoi(word.c_str());
            for (int j = 0; iss>>word; ++j){
                //extract DMS values
                double x = atof(word.c_str());
                matrix[gene_num][i-1][j] = x;
            }
        }
        else if ( word == "Gene_NUM"){
            iss >> word;
            gene_num = atoi(word.c_str());
        }
    }
    temp.close();
}

double Ran_Gaussian(double const mean, double const sigma)
{
    double x, y, r2;
    do{
        // choose x,y in uniform square [-1,+1]
        x = -1 + 2 * randomNumber();
        y = -1 + 2 * randomNumber();
        // check if it is in the unit circle
        r2 = x * x + y * y;
    }while (r2 > 1.0 || r2 == 0); 
    // Box-Muller transform
    return mean + sigma * y * sqrt (-2.0 * log (r2) / r2);
}

// Loads primordial genes in a VectStr
int LoadPrimordialGenes(const std::string& genelistfile, const std::string& genesPath)
{  
    std::ifstream genelistIN (genelistfile.c_str());
    if (!genelistIN.is_open()){
        std::cerr << "File could not be open: "<< genelistfile <<std::endl;
        exit(2);
    }
    std::cout << "Loading primordial genes file ..." << std::endl;
    int flag_AASeq = 0; 
    int gc = 0;
    while(!genelistIN.eof()){
        std::string word, line;
        getline(genelistIN,line);
        std::istringstream iss(line, std::istringstream::in);
        iss >> word; 
        if ( word=="G" ){
            iss >> word;
            word = genesPath + word;
            std::ifstream genefileIN (word.c_str(), std::ifstream::in | std::ifstream::out);
            if (!genefileIN.is_open()){
                std::cerr << "File could not be open: " << word << std::endl;
                exit(2);
            }
            assert ( gc > 0);

            int gn = -1;
            while(!genefileIN.eof()){
                  std::string w, l;
                  getline(genefileIN,l);
                  std::istringstream iss(l, std::istringstream::in);
                  iss >> w;
                  if (w=="Gene_NUM"){
                      iss >> w;
                      gn = atoi(w.c_str());
                  }
                  else if ( w=="N_Seq" ){
                      assert( gn >= 0);
                      iss >> w;
                      std::string aaseq=GetProtFromNuc(w);
                      //check stop codons in midsequence
                      size_t loc = aaseq.find('X', 0);
                      assert( loc == std::string::npos ); // no match
                      VectStr_iterator iter = PrimordialAASeq.begin();
                      PrimordialAASeq.insert(iter+gn, aaseq); 
                      flag_AASeq += 1;
                  }
            }
        }
        else if ( word=="Gene_Count" ){
            iss >> word;
            gc = atoi(word.c_str());
            PrimordialAASeq.reserve(gc);
        }
    }
    assert( flag_AASeq==gc );
    return gc;
}

// Reads a unit cell stored in binary format using Cell::dumpShort()
void qread_Cell(std::ifstream& IN, std::ofstream& OUT)
{
    char mybuffer[140];
    int na, ns;
    double f;
    std::string barcode;

    int l;
    IN.read((char*)&l, sizeof(int));
    std::vector<char> buf(l);
    IN.read(&buf[0], l);
    barcode.assign(buf.begin(), buf.end());

    OUT << barcode;

    IN.read((char*)(&na),sizeof(int));
    IN.read((char*)(&ns),sizeof(int));
    IN.read((char*)(&f),sizeof(double));
    
    sprintf(mybuffer,"\t%d\t%d\t%e", na, ns, f);
    OUT << mybuffer;
}

// Reads a unit cell stored in binary format using Cell::dumpSeq()
void seqread_Cell(std::ifstream& IN, std::ofstream& OUT)
{
    char mybuffer[140];
    int cell_id, cell_index, gene_size;
    double f,m;
    std::string barcode;

    IN.read((char*)(&cell_index),sizeof(int));
    IN.read((char*)(&cell_id),sizeof(int));

    int l;
    IN.read((char*)&l, sizeof(int));
    std::vector<char> buf(l);
    IN.read(&buf[0], l);
    barcode.assign(buf.begin(), buf.end());

    OUT << barcode;

    IN.read((char*)(&f),sizeof(double));
    IN.read((char*)(&m),sizeof(double));
    IN.read((char*)(&gene_size),sizeof(int));
    
    sprintf(mybuffer,"\t%d\t%e\t%e\t", cell_index, f, m);

    OUT << mybuffer << std::endl;

    for (int j=0; j<gene_size; ++j){
        std::string DNAsequence;   
        int Na, Ns;

        IN.read((char*)(&Na),sizeof(int));
        IN.read((char*)(&Ns),sizeof(int));

        //read DNA sequence
        int nl;
        IN.read((char*)&nl, sizeof(int));
        std::vector<char> buff(nl);
        IN.read(&buff[0], nl);  
        DNAsequence.assign(buff.begin(), buff.end());

        sprintf(mybuffer,"%d\tG\t%d\t%d\t",j,Ns,Na);
        OUT << mybuffer << std::endl;
        OUT << DNAsequence << std::endl;
    }
}

void read_Parent(std::ifstream& IN, std::ofstream& OUT)
{
    uint32_t a;

    IN.read((char*)(&a),sizeof(a));
    
    OUT << a;
}

// Reads a unit cell stored in binary format using Cell::dump()
void read_Cell(std::ifstream& IN, std::ofstream& OUT, bool DNA)
{
    char mybuffer[140];
    int cell_id, cell_index, gene_size;
    double m,f;
    std::string barcode;

    IN.read((char*)(&cell_index),sizeof(int));
    IN.read((char*)(&cell_id),sizeof(int));

    int l;
    IN.read((char*)&l, sizeof(int));
    std::vector<char> buf(l);
    IN.read(&buf[0], l);  
    barcode.assign(buf.begin(), buf.end());

    OUT << barcode;

    IN.read((char*)(&f),sizeof(double));
    IN.read((char*)(&m),sizeof(double));
    IN.read((char*)(&gene_size),sizeof(int));
    
    sprintf(mybuffer,"\t%d\t%.9f\t%e\t", cell_index, f, m);
    OUT << mybuffer << std::endl;

    for (int j=0; j<gene_size; ++j){
        double e, c, dg, f, eff;
        int gene_nid, Ns, Na;
        std::string DNAsequence;

        IN.read((char*)(&gene_nid),sizeof(int));   
        IN.read((char*)(&e),sizeof(double));
        IN.read((char*)(&c),sizeof(double));
        IN.read((char*)(&eff),sizeof(double));
        IN.read((char*)(&dg),sizeof(double));
        IN.read((char*)(&f),sizeof(double));

        IN.read((char*)(&Na),sizeof(int));
        IN.read((char*)(&Ns),sizeof(int));
        
        //read DNA sequence
        int nl;
        IN.read((char*)&nl, sizeof(int));
        std::vector<char> buff(nl);
        IN.read(&buff[0], nl);  
        DNAsequence.assign(buff.begin(), buff.end());

        sprintf(mybuffer,"%d\tG\t%e\t%.8f\t%.9f\t%d\t%d\t",j, c, dg, f, Na, Ns);
        OUT << mybuffer << std::endl;
        if (DNA) OUT << DNAsequence << std::endl;
        else OUT << GetProtFromNuc(DNAsequence) << std::endl;
    } 
}

// Counts the number of dissimilar characters in 2 strings
int StringDiff(const std::string& A, const std::string& B)
{
    unsigned int L = A.length();
    assert(L == B.length());
    int ctr = 0;  
    for (unsigned int i =0; i<L; ++i){
        if ( A.at(i) != B.at(i)) 
          ctr+=1; 
    } 
    return ctr;
}

std::string trim(const std::string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first){
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

// checks if given directory exists
bool isDirExist(const std::string& path)
{
#if defined(_WIN32)
    struct _stat info;
    if (_stat(path.c_str(), &info) != 0){
        return false;
    }
    return (info.st_mode & _S_IFDIR) != 0;
#else 
    struct stat info;
    if (stat(path.c_str(), &info) != 0){
        return false;
    }
    return (info.st_mode & S_IFDIR) != 0;
#endif
}

// creates directory from specified path
bool makePath(const std::string& path)
{
#if defined(_WIN32)
    int ret = _mkdir(path.c_str());
#else
    mode_t mode = 0755;
    int ret = mkdir(path.c_str(), mode);
#endif
    if (ret == 0)
        return true;

    switch (errno){
    case ENOENT:
        // parent didn't exist, try to create it
        {
            std::size_t pos = path.find_last_of('/');
            if (pos == std::string::npos)
#if defined(_WIN32)
                pos = path.find_last_of('\\');
            if (pos == std::string::npos)
#endif
                return false;
            if (!makePath( path.substr(0, pos) ))
                return false;
        }
        // try to create again
#if defined(_WIN32)
        return 0 == _mkdir(path.c_str());
#else 
        return 0 == mkdir(path.c_str(), mode);
#endif
    case EEXIST:
        // done
        return isDirExist(path);

    default:
        return false;
    }
}

void printProgress (double progress)
{
    std::cout << '[';
    int pos = PBWidth * progress;
    for (int i = 0; i < PBWidth; ++i){
        if (i < pos)
          std::cout << '=';
        else if (i == pos) 
          std::cout << '+';
        else 
          std::cout << ' ';
    }

    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}
