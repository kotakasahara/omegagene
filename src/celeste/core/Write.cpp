#include "Write.h"
using namespace std;

Write::Write(string inFn) : CelesteObject() {

    op       = false;
    filename = inFn;
}
Write::Write() : CelesteObject() {
    op = false;
}
int Write::open() {
    ofs.open(filename.c_str());
    if (!ofs) {
        cerr << "Cannot open " << filename << "." << endl;
        return 1;
    }
    op = true;
    return 0;
}
int Write::openApp() {
    ofs.open(filename.c_str(), ios::app);
    if (!ofs) {
        cerr << "Cannot open " << filename << "." << endl;
        return 1;
    }
    op = true;
    return 0;
}
int Write::close() {
    ofs.close();
    op = false;
    return 0;
}

WriteTTPVMcMDLog::WriteTTPVMcMDLog() : Write() {}
WriteTTPVMcMDLog::~WriteTTPVMcMDLog() {}

int WriteTTPVMcMDLog::write_ttpvMcMDLog(int step, int vstate) {
    ofs << step << "\t" << vstate + 1 << endl;
    return 0;
}
int WriteTTPVMcMDLog::write_VcMDLog(int step, std::vector<int> vstate) {
  //ofs << step;
  std::vector<int>::iterator itr;
  for (auto itr: vstate){
    ofs << "\t" << itr+1;
  }
  ofs << std::endl;
  return 0;
}

WriteTableLog::WriteTableLog() : Write() {}
WriteTableLog::~WriteTableLog() {}

int WriteTableLog::write_header() {
    return 0;
}
int WriteTableLog::write_row(int *values) {
    return 0;
}
int WriteTableLog::write_row(real *values) {
    return 0;
}
int WriteTableLog::write_row(std::vector<real> values) {
    return 0;
}

WriteVcMDParam::WriteVcMDParam() : Write() {}
WriteVcMDParam::~WriteVcMDParam() {}
int WriteVcMDParam::write(int interval,
			  std::vector< std::vector<string> > grp_names,
			  std::vector< std::vector<real> > min,
			  std::vector< std::vector<real> > max,
			  std::map< std::vector<int>, real > q_cano){
  ofs << interval << endl;
  ofs << min.size() << endl;
  for ( int i = 0; i < min.size(); i++){
    ofs << min[i].size();
    for ( const auto name : grp_names[i] )
      ofs << " " << name;
    ofs << endl;
    for ( int j = 0; j < min[i].size(); j++){
      ofs << min[i][j] << " " << max[i][j] << endl;
    }
  }
  for ( const auto q : q_cano ){
    for ( const auto vs : q.first ) {
      ofs << vs+1 << " ";
    }
    ofs << q.second << endl;
  }
  ofs<<"END"<<endl;
  return 0;
}

WriteTableLogBinary::WriteTableLogBinary() : WriteTableLog() {}
WriteTableLogBinary::~WriteTableLogBinary() {}
int WriteTableLogBinary::write_header() {
    ofs.write((const char *)&MAGIC_NUMBER, sizeof(int));
    ofs.write((const char *)&REAL_BYTE, sizeof(int));
    ofs.write((const char *)&n_cols, sizeof(int));
    return 0;
}
int WriteTableLogBinary::write_row(int *values) {
    // for(int i=0; i<n_cols; i++){
    // int val = values[i];
    // ofs.write((const char*)&val, sizeof(int));
    //}
    ofs.write((const char *)values, sizeof(int) * n_cols);
    return 0;
}
int WriteTableLogBinary::write_row(real *values) {
    //  for(int i=0; i<n_cols; i++){
    // real val = values[i];
    //    ofs.write((const char*)&val, sizeof(real));
    //  }
    ofs.write((const char *)values, sizeof(real) * n_cols);
    return 0;
}
int WriteTableLogBinary::write_row(std::vector<real> values) {
    //  for(int i=0; i<n_cols; i++){
    // real val = values[i];
    //    ofs.write((const char*)&val, sizeof(real));
    //  }
  std::vector<real>::iterator itr;
  for(auto itr: values){
    ofs.write((const char *)&itr, sizeof(real));
  }
  return 0;
}
WriteTableLogAscii::WriteTableLogAscii() : WriteTableLog() {}
WriteTableLogAscii::~WriteTableLogAscii() {}
int WriteTableLogAscii::write_header() {
    return 0;
}
int WriteTableLogAscii::write_row(int *values) {
    ofs << values[0];
    for (int i = 1; i < n_cols; i++) { ofs << "\t" << values[i]; }
    ofs << endl;
    return 0;
}
int WriteTableLogAscii::write_row(real *values) {
    char buf[1024];
    sprintf(buf, "%14.10e", values[0]);
    ofs << buf;
    for (int i = 1; i < n_cols; i++) {
        sprintf(buf, "%14.10e", values[i]);
        ofs << "\t" << buf;
    }
    ofs << endl;
    return 0;
}
int WriteTableLogAscii::write_row(std::vector<real> values) {
  //cout << "dbg 0303 write_row" << endl;
  //for(const auto v : values){
  //cout << v << endl;
  //}
  char buf[1024];
  sprintf(buf, "%14.10e", values[0]);
  ofs << buf;
  for (int i = 1; i < values.size(); i++) {
    sprintf(buf, "%14.10e", values[i]);
    ofs << "\t" << buf;
  }
  ofs << endl;
  return 0;
}
