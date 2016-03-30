#include "General.h"

using namespace std;

RandomNum::RandomNum() {}

RandomNum::~RandomNum() {}

int RandomNum::set_seed(int in_seed) {
    cout << "Random set_seed " << in_seed << endl;
    mt.seed(in_seed);
    return 0;
}

float RandomNum::get_float_01() {
    unsigned int rand   = mt();
    double       rand01 = rand / float(mt.max());
    // uniform_real_distribution<float> gen(0.0, 1.0);
    // float rand = gen(mt);
    // cout << "DBG01 random : " << rand << " - " << rand01 << " / " << mt.min() <<" " << mt.max()<<endl;
    return rand01;
}
