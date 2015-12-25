#include <iostream>
#include <ctime>
#include "define.h"
#include "Celeste.h"
using namespace std;

int main(int argn, char* argv[]){
    clock_t time_start = clock();
    Celeste celeste;
    celeste.setup(argn,argv);
    celeste.main_stream();
    clock_t time_end = clock();
    cout << "Time: " << (time_end - time_start) / (double)CLOCKS_PER_SEC << endl;
    return 0;
}
