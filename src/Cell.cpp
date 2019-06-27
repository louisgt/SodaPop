#include "Cell.h"

// By default the fitness function is set to neutral
int Cell::ff_ = 6;
bool Cell::useDist_ = false;
bool Cell::fromS_ = false;

Gene Cell::selected_gene = Gene();

Cell::Cell():
    barcode_(getBarcode()),
    ID_(0),
    parent_(0),
    o_mrate_(0),
    c_mrate_(0),
    fitness_(0)
    {
        geneBlocks_.reserve(maxGeneCount);
        genomeVec_.reserve(maxGeneCount);
    }

// Construct from cell file
Cell::Cell(std::ifstream & cell_in) {
    char mybuffer[140];
    geneBlocks_.reserve(maxGeneCount);
    genomeVec_.reserve(maxGeneCount);
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
            sprintf(mybuffer, "%s%s.gene", genesPath.c_str(), word.c_str());
            std::ifstream gene_data(mybuffer);
            if (!gene_data.is_open()) {
                std::cerr << "File could not be open: " << mybuffer << std::endl;
                exit(2);
            }
            Gene A(gene_data,this);
            genomeVec_.push_back(A);

            //Check if gene is correctly inserted
            auto i = genomeVec_.end();
            i--;
            // if gene fitness is null, assign a randomly fit value
            if ((*i).f()==0){
                (*i).ch_f(0.95+randomNumber()*0.02);
            }
            gene_data.close();
        }
    }
    selectFitness();
    // Update current rates
    this->UpdateRates();  
    // Fill gene length array
    this->FillGene_L();
}

// Constructs from a unit cell stored in binary 
Cell::Cell(std::ifstream & IN,
    const std::string & genesPath) {
    geneBlocks_.reserve(maxGeneCount);
    genomeVec_.reserve(maxGeneCount);

    char mybuffer[140];
    int l(0);
    int cell_id(0);
    int cell_index(0);
    int gene_size(0);
    double m(0);
    double f(0);

    IN.read((char*)(&cell_index), sizeof(int));
    IN.read((char*)(&cell_id), sizeof(int));
    ID_ = cell_id;
    parent_ = 0;

    //read barcode
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
    double e(0);
    double c(0);
    double dg(0);
    double eff(0);

    int gene_nid(0);
    int Ns(0);
    int Na(0);

    std::string DNAsequence;
    for (int j = 0; j < gene_size; ++j) {

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

        sprintf(mybuffer, "%s%d.gene", genesPath.c_str(), gene_nid);
        std::ifstream gene_data(mybuffer);
        if (!gene_data.is_open()) {
            std::cerr << "ERROR: Cannot open gene file " << mybuffer << std::endl;
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
        genomeVec_.push_back(G);
    }
    selectFitness();
    // Update current rates
    this->UpdateRates();  
    // Fill gene length array
    this->FillGene_L();
}

void Cell::linkGenes()
{
    for (auto gene_it = this->genomeVec_.begin(); gene_it != this->genomeVec_.end(); gene_it++) {
        gene_it->setCell(this);
    }
}

void Cell::print_summary_Gene_arr_() {
    for (const auto& gene : genomeVec_) {
         std::cout<<"current gene of cell"<< this->ID()<<" is gene"<< gene.num()<<" and has length "<<gene.geneLength()<<std::endl;
    }
}

void Cell::print_summary_Gene_L_() {
    std::cout<<"cumulative lengths of cell"<< this->ID()<<" is :"<<std::endl;
    for (const auto& cumul : geneBlocks_) {
        std::cout<< cumul <<std::endl;
    }
}

double Cell::fitness() const
{
    return fitness_;
}

// Initialize the cummulative gene length array
void Cell::FillGene_L() {
    int sum = 0;
    for (const auto& gene : genomeVec_) {
        sum += gene.geneLength();
        geneBlocks_.push_back(sum);
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

    for (const auto& gene : genomeVec_) {
        int Ns = gene.Ns();
        int Na = gene.Na();
        s += Ns;
        a += Na;
        sa += (Ns + Na);
    }

    if (spec == 0) 
        return sa;
    else if (spec == 1) 
        return s;
    else 
        return a;
}

void Cell::selectFitness()
{
    switch (Cell::ff_){
        case 1: fit = &Cell::fold;
            break;
        case 2: fit = &Cell::flux;
            break;
        case 3: fit = &Cell::toxicity;
            break;
        case 4: fit = &Cell::metabolicOutput;
            break;
        case 5: fit = &Cell::multiplicative;
            break;
        case 6: fit = &Cell::neutral;
            break;
        case 7: fit = &Cell::noMut;
            break;
        case 8: fit = &Cell::growthRate;
            break;
        default:;
    }
}

// FOLDING-STABILITY BASED FLUX FITNESS FUNCTION
double Cell::fold() const
{
    double f = std::accumulate(begin(genomeVec_), end(genomeVec_), 0.0, [](double i, const Gene& gene)
        { return gene.functional() + i;});
    return f;
}

// METABOLIC FLUX FITNESS FUNCTION
double Cell::flux() const
{
    double sum_func = 0;
    double a = 0;
    for (const auto& gene : genomeVec_) {
        // efficiency * Pnat * abundance
        sum_func += 1/(gene.eff()*gene.functional());
        a += gene.A_factor();
    }
    return static_cast<double>(a)/sum_func;
}

// MISFOLDING TOXICITY FITNESS FUNCTION
double Cell::toxicity() const
{
    // //sum (concentration*(1-Pnat)) over all genes using fold
    double f = std::accumulate(begin(genomeVec_), end(genomeVec_), 0.0, [](double i, const Gene& gene)
        { return gene.misfolded() + i;});
    return exp(-(misfoldingCost*f));
}

// COMBINED METABOLIC OUTPUT FITNESS FUNCTION
double Cell::metabolicOutput() const
{
    double flux = 0;
    double toxicity = 0;
    double a = 0;
    for (const auto& gene : genomeVec_) {
        flux += 1 / (gene.eff()*gene.functional());
        toxicity += gene.misfolded();
        a += gene.A_factor();
    }

    flux = a/flux;
    toxicity = misfoldingCost * toxicity;

    double fitness = flux - toxicity;
    // if toxicity > flux, return 0
    return (fitness < 0) ? 0 : fitness;
}

// MULTIPLICATIVE FITNESS FUNCTION
double Cell::multiplicative() const
{
    double f = std::accumulate(begin(genomeVec_), end(genomeVec_), 0.0, [](double i, const Gene& gene)
        { return gene.f()*gene.e() + i;});
    return f/static_cast<double>(gene_count());
}

// NEUTRAL FITNESS FUNCTION
double Cell::neutral() const
{
    return 1;
}

// TURN OFF MUTATIONS
double Cell::noMut() const
{
    return fitness();
}

double Cell::growthRate() const
{
    double f = std::accumulate(begin(genomeVec_), end(genomeVec_), 0.0, [](double i, const Gene& gene)
        { return gene.functional()*gene.eff() + i;});
    double fit = prefactor/f;
    return 1/(fit+1);
}

void Cell::UpdateRates()
{
    fitness_ = (this->*fit)();
}

void Cell::ranmut_Gene(std::ofstream& log,int ctr)
{
    // get genome size
    int L = geneBlocks_.back();

    // pick random site to mutate

    int site = static_cast<int>( L * randomNumber());

    // find the corresponding gene
    auto j = genomeVec_.begin();
    auto k = geneBlocks_.begin();

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

    //change statement to switch

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

void Cell::ranmut_Gene()
{
    // get genome size
    int L = geneBlocks_.back();

    // pick random site to mutate
    int site = static_cast<int>( L * randomNumber());

    // find the corresponding gene
    auto j = genomeVec_.begin();
    auto k = geneBlocks_.begin();

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

double Cell::normalizeFit(double fittest){
    if (fittest <= 0) {
        std::cerr << "Population collapse, average fitness is null.\n";
        exit(1);
    }
    double newfit = fitness()*1.0/fittest;
    ch_Fitness(newfit);
    UpdateRates();
    return newfit;
}

// Dump cell information to binary file
void Cell::dump(std::ofstream& OUT, int cell_index) const
{   
    OUT.write((char*)(&cell_index),sizeof(int));
    
    //cell ID
    OUT.write((char*)(&ID_),sizeof(int));

    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    int x = 0;
    double y = 0;

    y = fitness();
    OUT.write((char*)(&y),sizeof(double));
    
    y = c_mrate_;            
    OUT.write((char*)(&y),sizeof(double));

    x = static_cast<int>(genomeVec_.size());             
    OUT.write((char*)(&x),sizeof(int));

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
void Cell::dumpShort(std::ofstream& OUT) const
{
    int s = barcode().size();
    OUT.write((char*)&s, sizeof(int));

    OUT.write(barcode().c_str(), s);

    int x = 0;
    double y = 0;

    x = Na();
    OUT.write((char*)(&x),sizeof(int));

    x = Ns();
    OUT.write((char*)(&x),sizeof(int));

    y = fitness();         
    OUT.write((char*)(&y),sizeof(double));
}

// Dump cell parent to binary file
void Cell::dumpParent(std::ofstream& OUT) const
{
    uint32_t a;
    a = parent();
    OUT.write((char*)(&a),sizeof(a));
}

void Cell::UpdateNsNa()
{
    int new_Na = 0;
    int new_Ns = 0;
    for (const auto& gene : genomeVec_){
        new_Na += gene.Na();
        new_Ns += gene.Ns();
    }
    Total_Na_ = new_Na;
    Total_Ns_ = new_Ns;
}

// Print cell information to stdout
void Cell::PrintCell(int cell_ndx) const
{
      char mybuffer[140];
      sprintf(mybuffer,"C %6d %6d %12e %12e %d", cell_ndx, ID_, o_mrate_, mrate(), static_cast<int>(genomeVec_.size()));  
      std::cout << mybuffer << std::endl;
      for (const auto& gene : genomeVec_) {
        int gene_nid = gene.num();
        double e = gene.e();
        double c = gene.conc();
        double dg = -kT*log(gene.dg());
        int Ns = gene.Ns();
        int Na = gene.Na();
           
        sprintf(mybuffer,"G %d% 2.2f %10.8f %10.8f %d %d ", gene_nid, e, c, dg, Ns, Na);
        std::cout << mybuffer << std::endl;  
      }
      std::cout << std::endl;
}