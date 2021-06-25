#include "Config.h"
#include <iostream>
using namespace std;

Config::Config(){
  mode=M_DUMMY;
}

Config::~Config(){

}
void Config::setAll(int argn, char* argv[]){
  vector<string> arg;
  int i;
  for(i=1;i<argn;i++)
    arg.push_back(string(argv[i]));
  setAll(arg);
}

void Config::setAll(vector<string> arg){
  vector<string>::iterator itr;
  string type,val;
  for(itr=arg.begin(); itr!=arg.end(); itr++){
    if(*itr=="--mode"){
      itr++;
      if(*itr=="test")         { mode=M_TEST; }
      else if(*itr=="subzone")         { mode=M_SUBZONEBASED; }
      else if(*itr=="")           { mode=M_DUMMY; }
      else{
        cerr<<"invalid mode ["<<*itr<<"]\n"; exit(1);
      }
    }
    else if(*itr=="--i-cfg"){ itr++; fname_i_cfg=(*itr); }
    else if(*itr=="--i-params"){ itr++; fname_i_params=(*itr); }
    else if(*itr=="--i-qraw-is"){ itr++; fname_i_qraw_is=(*itr); }
    else if(*itr=="--o-qcano"){ itr++; fname_o_qcano=(*itr); }
    else if(*itr=="--o-qweight-opt"){ itr++; fname_o_qweight_opt=(*itr); }

    else if(*itr=="--mc-temp"){ itr++; mc_temp=atof((*itr).c_str()); }
    else if(*itr=="--mc-delta-x"){ itr++; mc_delta_x=atof((*itr).c_str()); }
    else if(*itr=="--mc-steps"){ itr++; mc_steps = atoi((*itr).c_str()); }
    else if(*itr=="--mc-log-interval"){ itr++; mc_log_interval = atoi((*itr).c_str());
      cout << "xxx " << mc_log_interval << endl;
    }
    
    else{
      cerr<<"unknown keyword <"<<*itr<<">"<<endl;
      exit(1);
    }
  }
}
