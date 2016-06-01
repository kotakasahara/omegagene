#ifndef __WRITE_H__
#define __WRITE_H__

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "CelesteObject.h"

class Write : public CelesteObject {
  private:
    std::string filename;
    bool        op;

  public:
    std::ofstream ofs;
    Write();
    Write(std::string inFn);
    void set_fn(std::string in_fn) { filename = in_fn; };
    std::string             getFn() { return filename; };
    bool                    is_open() { return op; };
    int                     open();
    int                     openApp();
    int                     close();
};

class WriteTTPVMcMDLog : public Write {
  private:
  public:
    WriteTTPVMcMDLog();
    ~WriteTTPVMcMDLog();
    int write_ttpvMcMDLog(int step, int vstate);
};

class WriteTableLog : public Write {
  private:
  protected:
    int n_cols;

  public:
    WriteTableLog();
    ~WriteTableLog();
    inline void set_ncolumns(int in_n_cols) { n_cols = in_n_cols; };
    virtual int                  write_header();
    virtual int write_row(int *values);
    virtual int write_row(real *values);
};

class WriteTableLogBinary : public WriteTableLog {
  private:
  protected:
  public:
    WriteTableLogBinary();
    ~WriteTableLogBinary();
    virtual int write_header();
    virtual int write_row(int *values);
    virtual int write_row(real *values);
};

class WriteTableLogAscii : public WriteTableLog {
  private:
  protected:
  public:
    WriteTableLogAscii();
    ~WriteTableLogAscii();
    virtual int write_header();
    virtual int write_row(int *values);
    virtual int write_row(real *values);
};

#endif