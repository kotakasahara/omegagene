#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
using namespace std; 

//#define TEMPERATURE  30

#define BINWIDTH 0.05
#define DX 0.001
//#define MAX_X 4.5
#define GAS_CONST    8.31451
#define JOULE_CAL    4.184

#define PRM_D 0.0005
#define PRM_A 10
#define PRM_B 20


int main(int argn, char* argv[]){

  string fname_i_cfg;

  if(argn!=4){
    cerr << "Usage:"<<endl;
    cerr << "./dblwell_1d [min_x] [max_x] [temperature]"<<endl;
    exit(1);
    }
  //double max_r = sqrt(MAX_X*MAX_X+MAX_X*MAX_X+MAX_X*MAX_X);
  double min_x = atof(argv[1]);
  double max_x = atof(argv[2]);
  double temperature = atof(argv[3]);
  int n_bin=(max_x-min_x)/BINWIDTH;
  vector<double> hist(n_bin);
  //vector<double> hist;
  for(int i = 0; i < n_bin; i++) hist[i] = 0.0;
  double temp_fact = GAS_CONST/JOULE_CAL * 1e-3 * temperature;
  double sum_pop = 0.0;
  for ( double x = min_x; x < max_x; x += DX){
    //for ( double y = 0; y < MAX_X; y += DX){
    //for ( double z = 0; z < MAX_X; z += DX){
    double x2=x*x;
    double pot=PRM_D*x2*(x2-PRM_A)*(x2-PRM_B);
    double pop = exp(-pot/temp_fact)*BINWIDTH;
    sum_pop += pop;
    int i_bin = (x-min_x)/BINWIDTH;
    hist[i_bin] += pop;
  }
  double min_pmf = 1e10;
  vector<double> pmf(n_bin);
  for(int i = 0; i < n_bin; i++){
    hist[i] = hist[i]/sum_pop;
    pmf[i] = -temp_fact*log(hist[i]);
    if(pmf[i] < min_pmf) min_pmf = pmf[i];
  }
  for(int i = 0; i < n_bin; i++){
    cout << i*BINWIDTH+min_x << " " << hist[i] << " " << pmf[i]-min_pmf << endl;
  }

  return(0);
}
