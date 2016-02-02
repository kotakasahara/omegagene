#ifndef __CELESTE_STOPWATCH__
#define __CELESTE_STOPWATCH__

#include <chrono>

namespace celeste { namespace utilities {
    class Stopwatch {
        std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

      public:
        void reset();
        double elapsed_ms();
        double elapsed_sec();
    };
} }

#endif