#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
using namespace std; 

//#define TEMPERATURE  30

#define Z_WIDTH 0.5
#define SZ_WIDTH 0.25
#define DX 0.01
#define MAX_X 3.0
#define GAS_CONST    8.31451
#define JOULE_CAL    4.184

#define PRM_D 0.0009
#define PRM_A 10
#define PRM_B 20


int main(int argn, char* argv[]){
  cout << "test"<< endl;
  string fname_i_cfg;
  if(argn!=2){
    cerr << "Usage:"<<endl;
    cerr << "./prob_dblwell_3d [temperature]"<<endl;
    cerr << "Many parameters are hard-coded. Review the source code."<<endl;
    exit(1);
  }
  double temperature = atof(argv[1]);
  

  int n_zones = MAX_X*2/SZ_WIDTH-1;
  int n_sz = MAX_X*2/SZ_WIDTH;
  double temp_fact = GAS_CONST/JOULE_CAL * 1e-3 * temperature;
  double sum_pop = 0.0;
  double max_pop = 0.0;  
  
  for ( int i_x = 0; i_x < n_zones; i_x++){
    double zg_x = -MAX_X + SZ_WIDTH * i_x;
    for ( int i_y = 0; i_y < n_zones; i_y++){
      double zg_y = -MAX_X + SZ_WIDTH * i_y;
      for ( int i_z = 0; i_z < n_zones; i_z++){
	double zg_z = -MAX_X + SZ_WIDTH * i_z;
	//cout << "zg" << zg_x << " " << zg_y << " " << zg_z << endl;
	for ( int sz_x = 0; sz_x < 2; sz_x++){
	  double g_x = zg_x + sz_x * SZ_WIDTH;
	  for ( int sz_y = 0; sz_y < 2; sz_y++){
	    double g_y = zg_y + sz_y * SZ_WIDTH;
	    for ( int sz_z = 0; sz_z < 2; sz_z++){
	      double g_z = zg_z + sz_z * SZ_WIDTH;
	      //cout << "g" << g_x << " " << g_y << " " << g_z << endl;
	      double sz_pop = 0.0;

	      for ( double x = g_x; x < g_x+SZ_WIDTH; x += DX){
		for ( double y = g_y; y < g_y+SZ_WIDTH; y += DX){
		  for ( double z = g_z; z < g_z+SZ_WIDTH; z += DX){
		    
		    double r2=x*x+y*y+z*z;
		    double pot=PRM_D*r2*(r2-PRM_A)*(r2-PRM_B);
		    double pop = exp(-pot/temp_fact) * DX*DX*DX;
		    sz_pop += pop;
		    
		  }
		}
	      }

	      sum_pop += sz_pop;
	      cout << i_x + 1 << " " << i_y+1 << " " << i_z+1 << " ";
	      cout << sz_x + 1 << " " << sz_y+1 << " " << sz_z+1 << " ";
	      cout << sz_pop << endl;
	      
	    }
	  }
	}

      }
    }
  }

  return(0);
}
