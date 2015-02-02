#include "Write.h"

Write::Write(string inFn)
  : CelesteObject(){
  
  op=false;
  filename =inFn;
}
Write::Write()
  : CelesteObject(){
  
  op=false;
}
int Write::open(){
  ofs.open(filename.c_str());
  if(!ofs){
    cerr <<"Cannot open "<<filename<<"."<<endl;
    return 1;
  }
  op=true;
  return 0;
}
int Write::openApp(){
  ofs.open(filename.c_str(),ios::app);
  if(!ofs){
    cerr <<"Cannot open "<<filename<<"."<<endl;
    return 1;
  }
  op=true;
  return 0;
}
int Write::close(){
  ofs.close();
  op=false;
  return 0;
}
