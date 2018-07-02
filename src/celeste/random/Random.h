#ifndef __CELESTE_RANDOM_H__
#define __CELESTE_RANDOM_H__

#include <random>

namespace celeste {
    namespace random {
        class Random {
          protected:
            int          _seed = 0;
            std::mt19937 generator;
            // std::uniform_real_distribution<double> dist;
	    
          public:
            Random(int seed);
            Random() = default;
            int  get_seed();
            void set_seed(int new_seed);
            double operator()();
            double operator()(double high);
            double operator()(double low, double high);
	    double normal(double mu, double sigma);
        };
    }
}

#endif
