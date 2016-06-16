/*
Random.cpp: a basic random number generator that is initialized with a provided seed, and can generate either numbers in
a range or between 0 and 1
*/
#include "Random.h"

namespace randlib = celeste::random;
using namespace randlib;

randlib::Random::Random(int seed) {
    _seed     = seed;
    generator = std::mt19937(seed);
    dist      = std::uniform_real_distribution<double>(0, 1);
}

int randlib::Random::get_seed() {
    return _seed;
}

void randlib::Random::set_seed(int new_seed) {
    generator.seed(new_seed);
}

double randlib::Random::operator()() {
    return dist(generator);
}

double randlib::Random::operator()(double high) {
    return dist(generator) * high;
}

double randlib::Random::operator()(double low, double high) {
    return dist(generator) * (high - low) + low;
}
