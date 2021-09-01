#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
using namespace std; 

//#define TEMPERATURE  30

#define BINWIDTH 0.05
#define DX 0.01
#define MAX_X 3.0
#define MAX_R 4.5
#define GAS_CONST    8.31451
#define JOULE_CAL    4.184

#define PRM_D 0.0009
#define PRM_A 10
#define PRM_B 20


int main(int argn, char* argv[]){

  string fname_i_cfg;
  double temperature = atof(argv[1]);
  if(argn>3){
    cerr << "Usage:"<<endl;
    cerr << "./dblwell [temperature]"<<endl;
    cerr << "Many parameters are hard-coded. Check the source."<<endl;
    exit(1);
  }
  //double max_r = sqrt(MAX_X*MAX_X+MAX_X*MAX_X+MAX_X*MAX_X);

  int n_bin=MAX_R/BINWIDTH;
  vector<double> hist(n_bin);
  //vector<double> hist;
  for(int i = 0; i < n_bin; i++) hist[i] = 0.0;
  double temp_fact = GAS_CONST/JOULE_CAL * 1e-3 * temperature;
  double sum_pop = 0.0;
  double max_pop = 0.0;
  for ( double x = 0; x < MAX_X; x += DX){
    for ( double y = 0; y < MAX_X; y += DX){
      for ( double z = 0; z < MAX_X; z += DX){
	double r2=x*x+y*y+z*z;
	double r = sqrt(r2);
	if(r > MAX_R) continue;
	double pot=PRM_D*r2*(r2-PRM_A)*(r2-PRM_B);
	double pop = exp(-pot/temp_fact) * 8 * BINWIDTH;
	sum_pop += pop;
	int i_bin = r/BINWIDTH;
	hist[i_bin] += pop;
	if (hist[i_bin] > max_pop) max_pop = hist[i_bin];
      }
    }
  }
  double min_pmf = 1e10;
  vector<double> pmf(n_bin);
  for(int i = 0; i < n_bin; i++){
    // Normalization with setting the summation to unity.
    hist[i] = hist[i]/sum_pop;
    // Normalization wtih setting the maximum to unity.
    //hist[i] = hist[i]/max_pop;
    pmf[i] = -temp_fact*log(hist[i]);
    if(pmf[i] < min_pmf) min_pmf = pmf[i];
  }
  for(int i = 0; i < n_bin; i++){
    cout << i*BINWIDTH << " " << hist[i] << " " << pmf[i]-min_pmf << endl;
  }
  //cout << "test"<< endl;

  return(0);
}
