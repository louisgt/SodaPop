#include "include_headers.h"
#include "gamma.h"
#include "deviates.h"

using namespace std;

int main(int argc, char *argv[]){
  if(argc != 5){
    std::cerr <<"command <r> <p> <seed1> <seed2> \n";
    exit(1);
  }

/* We note that
	NB(r, p) = Poisson( Gamma(r, beta))
	beta = (1-p)/p

	NB = Negative Binomial Distribution
	Poisson Distribution
	Gamma Distribution
*/

  	int r = atoi(argv[1]);
	double p = atof(argv[2]);

  	int s1 = atoi(argv[3]);
  	int s2 = atoi(argv[4]);

  	if( (p<=0) | (p>=1) ){
    		std::cerr << "Error in Rnd_NB(): p = "<< p << std::endl;
    		exit(1);
  	}

	//Instantiate gamma and poisson deviate generators
  	Gammadev gd(1,1,s1);
  	Poissondev pd(1, s2);
  	double beta = p/(1-p);
	
	double x;
 	 for(int i=0; i < 50000; ++i)
  	{
		//gamma deviate
		x = gd.dev(r,beta);
    		std::cout << pd.dev(x) << std::endl;
  	}

	return 0;
}
