#ifndef __CELESTE_OBJECT_H__
#define __CELESTE_OBJECT_H__

#define DBG 1

//#include <pair>

// type of real values
//typedef float real;
typedef float real;
// type of real values for pairwise energy calculation in GPU
//typedef float real_pw; 
typedef float real_pw; 
// type of real values for summation of force,energy
typedef double real_fc; 
// type of real values for bonding potentials
typedef double real_bp; 

#ifdef F_MPI
#include <mpi.h>
#define  mpi_real MPI_DOUBLE
// # typedef MPI_FLOAT  mpi_real
// # typedef MPI_DOUBLE mpi_real_pw
#define mpi_real_pw MPI_FLOAT
#endif


//typedef pair<int,int> int_pair;
//typedef pair<int,int> real3d;
#include <string> 
#include <iostream>
#include <iomanip>
using namespace std;

#include "define.h"

class CelesteObject {
 private:
 protected:
 public:
  CelesteObject();
  
  static const int MAX_N_ATOMTYPE;
  
  static const string EXE;
  static const string ABOUT_ME;
  static const real PI;
  static const int MAGIC_NUMBER;
  static const real EPS;
  
  static const real ELEM_CHARGE; // 
  static const real AVOGADRO; // [mol^-1]
  static const real PERMITTIVITY;    // [m^-3 kg^-1 s^4 A^2]
  
  // charge_coeff ELEM_CHARGE**2 * AVOGADRO / (4 * PI * PERMITTIVITY * 1e-10 * 4.184 * 1e+3)
  // [Angestrome cal mol-1]
  static const real CHARGE_COEFF;
  static const real FORCE_VEL;
  static const real GAS_CONST;  
  
  static const real JOULE_CAL;
  static const real KINETIC_COEFF;
  
  template <typename TYPE> inline const TYPE& max(const TYPE& a, const TYPE& b){ return a < b ? b : a; }
  template <typename TYPE> inline const TYPE& min(const TYPE& a, const TYPE& b){ return a > b ? b : a; }
  int cross(const double* a, const double* b, double* ret);
  int cross(const float* a, const float* b, double* ret);
};

#endif
