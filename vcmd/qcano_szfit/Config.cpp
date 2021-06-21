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
    if(*itr=="-mode"){
      itr++;
      if(*itr=="subzone")         { mode=M_SUBZONEBASED; }
      else if(*itr=="")           { mode=M_DUMMY; }
      else{
        cerr<<"invalid mode ["<<*itr<<"]\n"; exit(1);
      }
    }
    else if(*itr=="--i-cfg"){ itr++; fname_i_cfg=(*itr); }
    else if(*itr=="--i-params"){ itr++; fname_i_params=(*itr); }
    else if(*itr=="--o-qcano"){ itr++; fname_o_qcano=(*itr); }
    
    else{
      cerr<<"unknown keyword <"<<*itr<<">"<<endl;
      exit(1);
    }
  }
}
