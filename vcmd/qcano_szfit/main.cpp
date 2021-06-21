#include <iostream>
#include "Config.h"
#include "Read.h"
#include "VirtualStateCoupled.h"
using namespace std; 

int main(int argn, char* argv[]){
  Config cfg;
  Read reader;
  VirtualStateCoupling vsc;

  string fname_i_cfg;
  
  if(argn<2){
    cerr << "Usage: qcano_szfit [mode]"<<endl;
    cerr << "------------------------------------"<<endl;
    cerr << "mode "<<endl;
    cerr << "  --subzone"<<endl;
    exit(1);
  }
  cout << "conf.setall\n";
  cfg.setAll(argn,argv);
  if(cfg.fname_i_cfg!=string())
    cfg.setAll(reader.read_config(cfg.fname_i_cfg));

  vsc.setup(cfg);

  switch(cfg.mode){
  case M_SUBZONEBASED:
    vsc.subzonebased_mc();
  }

  return(0);
}
