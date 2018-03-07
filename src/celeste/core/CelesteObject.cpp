#include "CelesteObject.h"

#include <cstdlib>
#include <iostream>
#include <sstream>

#ifndef MACRO_BUILD_VERSION
#error-- MACRO_BUILD_VERSION is not defined; please invoke the compiler with -DMACRO_BUILD_VERSION=<VALUE>!
#endif

#ifndef MACRO_BUILD_TIMESTAMP
#error-- MACRO_BUILD_TIMESTAMP is not defined; please invoke the compiler with -DMACRO_BUILD_TIMESTAMP=<VALUE>!
#endif

using namespace std;

CelesteObject::CelesteObject() {}

// const int CelesteObject::MAX_N_ATOMTYPE = 40;

const string CelesteObject::EXE = "omegagene";
const string CelesteObject::ABOUT_ME =
    CelesteObject::EXE + " " + MACRO_BUILD_VERSION + " (" + MACRO_BUILD_TIMESTAMP + ")";
const string CelesteObject::DESCRIPTION = "";
// const int CelesteObject::REAL_BYTE = 4;
const int    CelesteObject::REAL_BYTE     = sizeof(real);
const real   CelesteObject::PI            = 3.14159265;
const int    CelesteObject::MAGIC_NUMBER  = 66261;
const string CelesteObject::LS_VERSION    = "v.0.39.h";
const real   CelesteObject::EPS           = 1e-10;
const real   CelesteObject::EPS3          = 1e-30;
const real   CelesteObject::ELEM_CHARGE   = 1.60217657e-19;
const real   CelesteObject::AVOGADRO      = 6.0221413e23;
const real   CelesteObject::PERMITTIVITY  = 8.85418782e-12;
const real   CelesteObject::CHARGE_COEFF  = 332.06378;
const real   CelesteObject::FORCE_VEL     = 4.184e-4;
const real   CelesteObject::GAS_CONST     = 8.31451;
const real   CelesteObject::JOULE_CAL     = 4.184;
const real   CelesteObject::KINETIC_COEFF = (1e7 / (JOULE_CAL * 1e3)) * 0.5;
const real   CelesteObject::BOLTZMAN      = 1.380658e-23;
// const int CelesteObject::MAX_N_NB15OFF = 32;

int CelesteObject::cross(const double *a, const double *b, double *ret) {
    ret[0] = a[1] * b[2] - a[2] * b[1];
    ret[1] = a[2] * b[0] - a[0] * b[2];
    ret[2] = a[0] * b[1] - a[1] * b[0];
    return 0;
}
int CelesteObject::cross(const float *a, const float *b, float *ret) {
    ret[0] = a[1] * b[2] - a[2] * b[1];
    ret[1] = a[2] * b[0] - a[0] * b[2];
    ret[2] = a[0] * b[1] - a[1] * b[0];
    return 0;
}
int CelesteObject::error_exit(const string msg, const string error_code) {
    stringstream ss;
    ss << "[ Error : " << error_code << " ]" << endl;
    ss << msg << endl;

    cerr << ss.str();
    exit(1);
    return 0;
}
