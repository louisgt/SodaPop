#include "global.h"

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

class Gene 
{
    private:
        int g_num_;		//numeric ID pointing to primordial gene
        int ln_;		//length nuc seq
        int la_;		//length aa seq

        int Na_;		//number of non-synonymous substitutions
        int Ns_;		//number of sysnonymous substitutions
        
        std::string nucseq_;	//nucleotide sequence
        
        double dg_;		//stability
    public:
        double conc;	//concentration
        double e;		//essentiality: 1-if directly involved in replication, 0-otherwise
   
    public:
        Gene();
        Gene(int, std::string, double);
        Gene(std::fstream&);
        Gene(const Gene&); //copy constructor 
        ~Gene(); 
      
        bool operator==(Gene&);
        Gene& operator=(const Gene&);
      
        double Mutate_BP_Gaussian(int, int);
        std::string Mutate_BP(int, int);

        void Update_Sequences(std::string);
     
        void ch_dg(const double a){ dg_ = a; }
        void ch_Na(const int a){ Na_ = a; }
        void ch_Ns(const int a){ Ns_ = a; }
        void ch_ln(int l){ln_ = l;}
        void ch_la(int l){la_ = l;}
        void ch_gnum(int i){g_num_ = i;}

        const int num(){return g_num_;}
        const int length(){return ln_;}
        const int AAlength(){return la_;}
        const std::string nseq(){return nucseq_;}
        const double dg(){return dg_;}
        const int Ns(){return Ns_;}
        const int Na(){return Na_;}

        double CheckDG();
        double functional();
        double misfolded();
        double Pnat();
};

Gene::Gene(){
      g_num_ = 0;
      ln_ = 0; la_ = 0;
      Na_ = 0; Ns_ = 0;
      nucseq_ = ""; 
      dg_ = 1;
      conc = 1;
      e = 0;
}

//Input: gene number, nuc. sequence, concentration
Gene::Gene(const int i, const std::string a, double c)
{
    if((a.length() % 3) != 0){
        std::cerr << "Invalid length for nucleotide sequence: " << a.length() << std::endl;
        exit(2);
    }
    else{
        g_num_=i;
        nucseq_=a;
        ln_=a.length();
        la_=ln_/3;
        std::string aaseq = GetProtFromNuc(nucseq_);

        //check for stop codons in midsequence
        std::string::size_type loc = aaseq.find("X", 0 );
        if(loc != std::string::npos){
            std::cerr << "ERROR: DNA sequence has STOP codons in the middle"<<std::endl;
            exit(2);
        }           

        dg_= 1;
        conc = c;
        e = 0;
        Na_ = 0;
        Ns_ = 0;
    }
}

//Input: gene file
Gene::Gene(std::fstream& gene_in)
{
    std::string line;
    // Read gene file line by line
    while (!gene_in.eof()){
        getline(gene_in,line);
        std::string word;
        std::istringstream iss(line, std::istringstream::in);
        iss >> word;
        if (word == "Gene_NUM"){ 
            iss>>word; 
            g_num_ = atoi(word.c_str());
        }
        else if (word == "N_Seq"){ 
            iss>>nucseq_;
            ln_=nucseq_.length();
            if((ln_ % 3) != 0){
                  std::cerr << "Invalid length for nucleotide sequence: " << ln_ << std::endl;
                  exit(2);
            }
            else{
                  la_ = ln_/3;
                  std::string aaseq=GetProtFromNuc(nucseq_);

                  //check stop codons in midsequence
                  std::string::size_type loc = aaseq.find("X", 0 );
                  if(loc != std::string::npos){
                        std::cerr << "ERROR: DNA sequence has STOP codons in the middle"<<std::endl;
                        exit(2);
                    }           
            }
        }
        else if ( word=="E" ){
            iss >> word;
            e = atoi(word.c_str());
        }
        else if (word == "CONC"){
            iss >> word;
    	   conc = atof(word.c_str());
        }
        else if (word == "DG")
        { 
            iss>>word; 
      	    dg_ = atof(word.c_str());
            dg_ = exp(-dg_/kT);
        }
        else if (word == "//"){;}//do nothing
    }
    Na_ = 0; //default
    Ns_ = 0;
}

// copy constructor
Gene::Gene(const Gene& G)
{
    g_num_ = G.g_num_;
    ln_ = G.ln_;
    la_ = G.la_;
    nucseq_ = G.nucseq_;
    dg_ = G.dg_; 
    conc = G.conc;
    e = G.e;
    Na_ = G.Na_;
    Ns_ = G.Ns_;
}

Gene::~Gene()
{
}

//Genes are equal if DNA sequence and concentration are equal.
bool Gene::operator== (Gene& G) 
{
    std::string temp = G.nseq();
    if ( (temp.compare(nucseq_) == 0) && (conc == G.conc) ) return true;
    else return false;
}

// assignment overloading
Gene& Gene::operator=(const Gene& A)
{ 
    if (this != &A){
        this->g_num_ = A.g_num_;
        this->ln_ = A.ln_;
        this->la_ = A.la_;
        this->dg_ = A.dg_;
        this->conc = A.conc;
        this->e = A.e;
        this->Na_ = A.Na_;
        this->Ns_ = A.Ns_;
        (this->nucseq_).assign(A.nucseq_);
    }
    return *this;
}

/*
This version of the mutation function draws the DDG value from a gaussian distribution
with a shifting mean to mimic sequence depletion.
*/
double Gene::Mutate_BP_Gaussian(int i, int j)
{ 
    if(i>=ln_){
        std::cerr << "ERROR: Mutation site out of bounds."<< std::endl;
        exit(2);
    }       
    double fNs = 0.775956284; //fraction of non-synonymous substitutions in a typical protein
    double ran = RandomNumber();
       
    if(ran <= fNs){//non-synonymous mutation
        double temp = Ran_Gaussian(1.0, 1.7);
        double x = exp(-temp/kT);

        this->dg_ *= x;
        this->Na_ += 1;

        return x;
    }
    else{
        this->Ns_ += 1;

        return 1;
    }
}

/*
This version of the mutation function gets the DDG value from the DDG matrix
input by the user.
INPUT: 
    i -> site to mutate
    j -> bp to mutate to
*/
std::string Gene::Mutate_BP(int i, int j)
{ 
    // extract codon to be mutated
    int cdn_ndx = (i%3);
    int cdn_start = i - cdn_ndx; 
    int resi = cdn_start/3;

    // fetch current codon
    std::string cdn_curr = nucseq_.substr(cdn_start, 3);
    // fetch current amino acid
    int aa_curr = GetIndexFromCodon(cdn_curr);
    std::string cdn_new = cdn_curr;

    std::string s = PrimordialAASeq.at(g_num_);     
    int aa_primo = GetIndexFromAA(s.at(resi));

    // get mutated bp
    std::string bp = AdjacentBP( cdn_curr.substr(cdn_ndx, 1), j); //new BP
   
    // mutate codon
    cdn_new.replace(cdn_ndx, 1, bp);
    // check for stop codon
    cdn_new = n3_to_n3(cdn_new, cdn_curr, cdn_ndx);
    // get new amino acid
    int aa_new = GetIndexFromCodon(cdn_new);
    
    // get DDG value from matrix
    double x = matrix[g_num_][resi][aa_new-1];

    std::string mutation = std::to_string(g_num_) + '\t' + GetProtFromNuc(cdn_curr) + '\t' + std::to_string(resi) + '\t' + GetProtFromNuc(cdn_new);

    // fetch primordial amino acid

    //Ignore mutations to and from CYSTEINE
    if( (aa_new==2) || (aa_curr==2)){
        return "CYSTEINE\tNA\tNA\tNA";
    }

    //Case unphysical DDG estimate
    if( x>DDG_min || x<DDG_max){
        return "UNPHYSICAL\tNA\tNA\tNA";
    }

    if( aa_curr == aa_new){//SILENT
          nucseq_.replace(cdn_start, 3, cdn_new);
          Ns_ += 1;
          return "SILENT\tNA\tNA\tNA";
    }
    else if(aa_primo == aa_new){//REVERT TO WT

          double x_curr = matrix[g_num_][resi][aa_curr-1];
          assert( x_curr<DDG_min || x_curr>DDG_max); 
          
          dg_ /= x_curr;
          nucseq_.replace(cdn_start, 3, cdn_new);
          Na_ += 1;
          return mutation;
    }
    else{//TYPICAL NON-SYNONYMOUS

          double x_curr = matrix[g_num_][resi][aa_curr-1];
          assert( x_curr<DDG_min || x_curr>DDG_max); 

          // assign new DG value
          dg_ /= x_curr;
          dg_ *= x;
          nucseq_.replace(cdn_start, 3, cdn_new);
          Na_ += 1;
          return mutation;
    }
}


// Updates the current DNA sequence
void Gene::Update_Sequences(const std::string DNAsequence)
{ 
    int l = DNAsequence.length();

    if(l != ln_)
    {
        std::cerr << "ERROR: Replacing DNA sequence with a non-equal length DNA. "<< std::endl;
        exit(2);
    }       

    nucseq_ = DNAsequence;
}

// from Privalov 1979 (see also: Serohijos & Shakhnovich 2013)
// Boltzmann probability of the gene product to be in the native state
double Gene::Pnat()
{
    return this->dg()/(1+this->dg());
}

// Number of functional copies in the cell
double Gene::functional()
{
    return (this->conc)*this->Pnat();
}

// Number of misfolded copies in the cell
double Gene::misfolded()
{
    return (this->conc)*(1-Pnat());
}

//***************************************************
//	ERROR CHECKING ROUTINES
//***************************************************
double Gene::CheckDG()
{
      std::string A0 = PrimordialAASeq.at(g_num_);

      std::string A = GetProtFromNuc(nucseq_);
      std::cout << A0 << std::endl;
      std::cout << A << std::endl;

      double ddg = 1; 
      //calculate DDG
      for(int i=0; i<la_; i++){
        char aa0 = A0.at(i);
        std::cout << i << " " << aa0 << " ";
        char aa = A.at(i);
        std::cout << aa;
        if(aa != aa0){
              int j = GetIndexFromAA(aa);
              double x = matrix[g_num_][i][j-1]; 
              if ( x>DDG_min || x<DDG_max) std::cout << " Ignoring " << -kT*log(x);
              else{
                    ddg*=x;
                    std::cout << " " << -kT*log(x);
              }
        }
        std::cout << std::endl;
      }  
      return ddg; 
}
