/*
Stopwatch.cpp: a basic stopwatch to time runs
*/
#include "Stopwatch.h"

namespace utilities = celeste::utilities;
using namespace utilities;
using namespace std;

void utilities::Stopwatch::reset() {
    start = chrono::steady_clock::now();
}

double utilities::Stopwatch::elapsed_ms() {
    return chrono::duration<double, milli>(chrono::steady_clock::now() - start).count();
}

double utilities::Stopwatch::elapsed_sec() {
    return chrono::duration<double>(chrono::steady_clock::now() - start).count();
}
