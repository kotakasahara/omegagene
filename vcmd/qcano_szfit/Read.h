#ifndef __READ_H__
#define __READ_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

class Read {
 private:
  string filename;
  bool op;
  ifstream ifs;

 public:
  explicit Read();
  ~Read();
  int open(string in_filename);
  int close();
  vector<string> read_config(string in_filename);
};
#endif
