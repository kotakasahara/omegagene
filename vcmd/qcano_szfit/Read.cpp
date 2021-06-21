#include "Read.h"

Read::Read(){
}

Read::~Read(){
}

int Read::open(string in_filename){
  filename = in_filename;
  ifs.open(filename.c_str());
  if(!ifs){
    cerr <<"Cannot open "<< filename << "." <<endl;
    return 1;
  }
  op=true;
  return 0;
}
int Read::close(){
  ifs.close();
  op=false;
  return 0;
}
vector<string> Read::read_config(string in_filename){
  vector<string> vconf;
  open(in_filename);
  string buf;
  while(getline(ifs,buf)){
    if(buf[0]!='#'){
      stringstream ss(buf);
      string tmp;
      while(ss>>tmp) vconf.push_back(tmp);
    }
  }
  close();
  return vconf;
}
