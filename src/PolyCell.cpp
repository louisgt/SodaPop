#include "PolyCell.h"

// By default the fitness function is set to neutral
int PolyCell::ff_ = 6;
bool PolyCell::useDist_ = false;
bool PolyCell::fromS_ = false;

PolyCell::PolyCell(){}
PolyCell::PolyCell(std::fstream& f) : Cell(f)
{
    selectFitness();
	// Update current rates
  	this->UpdateRates();  
  	// Fill gene length array
  	this->FillGene_L();
}    
PolyCell::PolyCell(std::fstream& f, const std::string& s) : Cell(f,s)
{
    selectFitness();
	// Update current rates
  	this->UpdateRates();  
  	// Fill gene length array
  	this->FillGene_L();
}

// Initialize the cumulative gene length array
void PolyCell::FillGene_L()
{
    int sum = 0;
    std::vector<Gene>::iterator i;
    //for(i = genomeVec_.begin(); i != genomeVec_.end(); ++i){
    for (auto& gene : genomeVec_){
        sum+= gene.geneLength();
        geneBlocks_.push_back(sum);
    }
}

void PolyCell::selectFitness()
{
    switch (PolyCell::ff_){
        case 1: fit = &PolyCell::fold;
            break;
        case 2: fit = &PolyCell::flux;
            break;
        case 3: fit = &PolyCell::toxicity;
            break;
        case 4: fit = &PolyCell::metabolicOutput;
            break;
        case 5: fit = &PolyCell::multiplicative;
            break;
        case 6: fit = &PolyCell::neutral;
            break;
        case 7: fit = &PolyCell::noMut;
            break;
        case 8: fit = &PolyCell::growthRate;
            break;
        default:;
    }
}

// FOLDING-STABILITY BASED FLUX FITNESS FUNCTION
double PolyCell::fold() const
{
    double sum_func = 0;
    //sum (concentration*Pnat) over all genes
    //for(auto gene_it = genomeVec_.begin(); gene_it != genomeVec_.end(); ++gene_it){
    for (const auto& gene : genomeVec_) {
        sum_func += gene.functional();
    }
    return sum_func;
}

// METABOLIC FLUX FITNESS FUNCTION
double PolyCell::flux() const
{
    double sum_func = 0;
    double a = 0;
    //for(auto gene_it = genomeVec_.begin(); gene_it != genomeVec_.end(); ++gene_it){
    for (const auto& gene : genomeVec_) {
        // efficiency * Pnat * abundance
        sum_func += 1/(gene.eff()*gene.functional());
        a += gene.A_factor();
    }
    return static_cast<double>(a)/sum_func;
}

// MISFOLDING TOXICITY FITNESS FUNCTION
double PolyCell::toxicity() const
{
    double f = 0;
    //sum (concentration*(1-Pnat)) over all genes
    //for(auto gene_it = genomeVec_.begin(); gene_it != genomeVec_.end(); ++gene_it){
    for (const auto& gene : genomeVec_) {
        f += gene.misfolded();
    }
    return exp(-(MISFOLDING_COST*f));
}

// COMBINED METABOLIC OUTPUT FITNESS FUNCTION
double PolyCell::metabolicOutput() const
{
    double flux = 0;
    double toxicity = 0;
    double a = 0;
    //for (auto& it : genomeVec_) {
    for (const auto& gene : genomeVec_) {
        flux += 1 / (gene.eff()*gene.functional());
        toxicity += gene.misfolded();
        a += gene.A_factor();
    }

    flux = a/flux;
    toxicity = MISFOLDING_COST * toxicity;

    double fitness = flux - toxicity;
    // if toxicity > flux, return 0
    return (fitness < 0) ? 0 : fitness;
}

// MULTIPLICATIVE FITNESS FUNCTION
double PolyCell::multiplicative() const
{
    double fitness = 0;
    //for (auto gene_it = genomeVec_.begin(); gene_it != genomeVec_.end(); ++gene_it)
    for (const auto& gene : genomeVec_){
        fitness += gene.f()*gene.e();
    }
    return fitness/gene_count();
}

// NEUTRAL FITNESS FUNCTION
double PolyCell::neutral() const
{
    return 1;
}

// TURN OFF MUTATIONS
double PolyCell::noMut() const
{
    return fitness();
}

// FOLDING-STABILITY BASED FLUX FITNESS FUNCTION
double PolyCell::growthRate() const
{
    double sum = 0;
    //sum (concentration*Pnat) over all genes
    //for(auto gene_it = genomeVec_.begin(); gene_it != genomeVec_.end(); ++gene_it){
    for (const auto& gene : genomeVec_) {
        sum += gene.functional()*gene.eff();
    }
    double fit = PREFACTOR/sum;
    return 1/(fit+1);
}

void PolyCell::UpdateRates()
{
    fitness_ = (this->*fit)();
}

void PolyCell::ranmut_Gene(std::ofstream& log,int ctr)
{
    // get genome size
    int L = geneBlocks_.back();

    // pick random site to mutate

    int site = static_cast<int>( L * randomNumber());

    // find the corresponding gene
    std::vector<Gene>::iterator j = genomeVec_.begin();
    VectInt_citerator k = geneBlocks_.begin();

    if (site >= (*k)){
    // random number generated is greater than
    // the cumulative sum of genes
         for (k = geneBlocks_.begin(); k != geneBlocks_.end(); ++k){
             if( site< (*k) ) 
                break;
             j++; 
         }        
         k--;
         site = site - (*k);        
    }

    std::string mutation = "";

    int bp = static_cast<int>(3 * randomNumber());

    double wi = fitness();
    if (fromS_){
        if (useDist_){
            (*j).Mutate_Select_Dist(site,bp);
        }
        else{
            mutation = (*j).Mutate_Select(site,bp);
        }
    }
    else{
        if (useDist_){
            (*j).Mutate_Stabil_Gaussian(site,bp);
        }
        else{
            mutation = (*j).Mutate_Stabil(site,bp);
        }
    }

    UpdateRates();
    double wf = fitness();
    double s = wf - wi;

    // save beneficial mutations to log
    // we could save all mutations with abs(s) >= some value x
    log << barcode().c_str() << "\t";
    log << std::fixed;
    log << mutation << "\t";
    log << s << "\t";
    log << ctr << std::endl;
}

void PolyCell::ranmut_Gene()
{
    // get genome size
    int L = geneBlocks_.back();
    // pick random site to mutate

    int site = static_cast<int>( L * randomNumber());

    // find the corresponding gene
    std::vector<Gene>::iterator j = genomeVec_.begin();
    VectInt_citerator k = geneBlocks_.begin();


    if (site >= (*k)){
    // random number generated is greater than
    // the cumulative sum of genes
         for (k = geneBlocks_.begin(); k != geneBlocks_.end(); ++k){
             if (site< (*k) )
                break;
             j++; 
         }        
         k--;
         site = site - (*k);        
    }

    int bp = static_cast<int>(3 * randomNumber());
    // what is the input type?
    if (fromS_){
        if (useDist_){
            (*j).Mutate_Select_Dist(site,bp);
        }
        else{
            (*j).Mutate_Select(site,bp);
        }
    }
    else{
        if (useDist_){
            (*j).Mutate_Stabil_Gaussian(site,bp);
        }
        else{
            (*j).Mutate_Stabil(site,bp);
        }
    }
         
    UpdateRates();
}

double PolyCell::normalizeFit(double fittest){
    if (fittest <= 0) {
        std::cerr << "Population collapse, average fitness is null.\n";
        exit(1);
    }
    double newfit = (genomeVec_.begin()->f())/fittest;
    genomeVec_.begin()->ch_f(newfit);
    UpdateRates();
    return newfit;
}

// Dump cell information to binary file
void PolyCell::dump(std::fstream& OUT, int cell_index) const
{
    int x;
    double y;
    
    OUT.write((char*)(&cell_index),sizeof(int));
    
    //cell ID
    OUT.write((char*)(&ID_),sizeof(int));

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    y = fitness();
    OUT.write((char*)(&y),sizeof(double));
    
    y = c_mrate_;		 	 
    OUT.write((char*)(&y),sizeof(double));

    x = static_cast<int>(genomeVec_.size());		 	 
    OUT.write((char*)(&x),sizeof(int));

   //for(auto gene_it = genomeVec_.begin(); gene_it != genomeVec_.end(); ++gene_it){
   for (const auto& gene : genomeVec_) {
        int gene_nid = gene.num();
        double s = gene.e();
        double c = gene.conc();
        double eff = gene.eff();
        double dg = -kT*log(gene.dg());
        double f = gene.f();

        int Ns = gene.Ns();
        int Na = gene.Na();

        OUT.write((char*)(&gene_nid),sizeof(int));
        OUT.write((char*)(&s),sizeof(double));
        OUT.write((char*)(&c),sizeof(double));
        OUT.write((char*)(&eff),sizeof(double));
        OUT.write((char*)(&dg),sizeof(double));
        OUT.write((char*)(&f),sizeof(double));
        OUT.write((char*)(&Na),sizeof(int));
        OUT.write((char*)(&Ns),sizeof(int));

        //Save length of nucleo sequence
        std::string DNAsequence = gene.geneSeq();
        int nl = DNAsequence.length();
        OUT.write((char*)&nl, sizeof(int));
        OUT.write(DNAsequence.data(), nl);
    }
}

// Dump cell summary to binary file
void PolyCell::dumpShort(std::fstream& OUT) const
{
    int x;
    double y;

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    x = Na();
    OUT.write((char*)(&x),sizeof(int));

    x = Ns();
    OUT.write((char*)(&x),sizeof(int));

    y = fitness();         
    OUT.write((char*)(&y),sizeof(double));
}

// Dump cell summary to binary file
void PolyCell::dumpSeq(std::fstream& OUT, int cell_index) const
{
    int x;
    double y;

    OUT.write((char*)(&cell_index),sizeof(int));

    //cell ID
    OUT.write((char*)(&ID_),sizeof(int));

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    y = fitness();         
    OUT.write((char*)(&y),sizeof(double));

    y = c_mrate_;            
    OUT.write((char*)(&y),sizeof(double));

    x = static_cast<int>(genomeVec_.size());             
    OUT.write((char*)(&x),sizeof(int));

    //for(auto gene_it = genomeVec_.begin(); gene_it != genomeVec_.end(); ++gene_it){
    for (const auto& gene : genomeVec_) {

        int Ns = gene.Ns();
        int Na = gene.Na();

        OUT.write((char*)(&Na),sizeof(int));
        OUT.write((char*)(&Ns),sizeof(int));

        //Save length of nucleo sequence
        std::string DNAsequence = gene.geneSeq();
        int nl = DNAsequence.length();
        OUT.write((char*)&nl, sizeof(int));
        OUT.write(DNAsequence.data(), nl);
    }
}

// Dump cell parent to binary file
void PolyCell::dumpParent(std::fstream& OUT) const
{
    uint32_t a;
    a = parent();
    OUT.write((char*)(&a),sizeof(a));
}

void PolyCell::UpdateNsNa()
{
    int new_Na = 0;
    int new_Ns = 0;
    //for(auto gene_it = genomeVec_.begin(); gene_it != genomeVec_.end(); ++gene_it){
    for (auto& gene : genomeVec_){
        new_Na += gene.Na();
        new_Ns += gene.Ns();
    }
    Total_Na_ = new_Na;
    Total_Ns_ = new_Ns;
}

// Print cell information to stdout
void PolyCell::PrintCell(int cell_ndx) const
{
      char buffer[140];
      sprintf(buffer,"C %6d %6d %12e %12e %d", cell_ndx, ID_, o_mrate_, mrate(), static_cast<int>(genomeVec_.size()));  
      std::cout << buffer << std::endl;
      //for(auto gene_it = genomeVec_.begin(); gene_it != genomeVec_.end(); ++gene_it){
      for (const auto& gene : genomeVec_) {
        int gene_nid = gene.num();
        double e = gene.e();
        double c = gene.conc();
        double dg = -kT*log(gene.dg());
        int Ns = gene.Ns();
        int Na = gene.Na();
           
        sprintf(buffer,"G %d% 2.2f %10.8f %10.8f %d %d ", gene_nid, e, c, dg, Ns, Na);
        std::cout << buffer << std::endl;  
      }
      std::cout << std::endl;
}