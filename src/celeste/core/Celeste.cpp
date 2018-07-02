#include "Celeste.h"
using namespace std;

int Celeste::setup(int argn, char *argv[]) {
    cout << ABOUT_ME << endl;
    string fn_cfg;
    if (argn < 2) {
        cerr << "Usage: " << CelesteObject::EXE << " [mode]" << endl;
        cerr << "------------------------------------" << endl;
        std::exit(1);
    }
    cout << "conf.set_arguments\n";
    config = Config(argn, argv);
    cout << "/setup\n";
    return 0;
}

int Celeste::main_stream() {
    cout << "mainstream\n";
    switch (config.mode) {
        case M_TEST: test_mode(); break;
        case M_DYNAMICS: dynamics_mode(); break;
        default:
            cout << "Invalid Mode is specified.\n";
            test_mode();
            break;
    }
    cout << "Terminated.\n";
    return 0;
}

int Celeste::test_mode() {
    cout << "test_mode.\n";
    MmSystem mmsys;
    Read(config.fn_inp).load_launch_set(mmsys);
    // mmsys.write_data();
    return 0;
}

int Celeste::dynamics_mode() {
    cout << "dynamics_mode\n";

    // #if defined(F_CUDA) && defined(F_MPI)
    //   cout << "F_CUDA + F_MPI flags = ON" << endl;
    //  MpiGpuDynamicsMode* dynamics = new MpiGpuDynamicsMode;
    // GpuDynamicsMode* dynamics = new GpuDynamicsMode;
    //  MPI_Init(NULL, NULL);
    //#if defined(F_CUDA)
    //  cout << "F_CUDA flag = ON" << endl;
    //  GpuDynamicsMode* dynamics = new GpuDynamicsMode;
    //#else
    DynamicsMode *dynamics;
    if (config.integrator_type == INTGRTR_LEAPFROG_PRESTO) {
        dynamics = new DynamicsModePresto();
    } else if (config.integrator_type == INTGRTR_VELOCITY_VERLET) {
        dynamics = new DynamicsModeVelocityVerlet();
    } else if (config.integrator_type == INTGRTR_ZHANG) {
        dynamics = new DynamicsModeZhang();
    } else if (config.integrator_type == INTGRTR_LANGEVIN) {
        dynamics = new DynamicsModeLangevin();    } else {
        cerr << "Unknown Integrator" << endl;
        std::exit(-1);
    }
    //#endif

    if (DBG >= 1) cout << "DBG1: dynamics->set_config_parameters(cfg)" << endl;

    dynamics->set_config_parameters(&config);
    Read(config.fn_inp).load_launch_set(dynamics->mmsys);

    if (DBG >= 1) cout << "DBG1: dynamics->initial_preprocess()" << endl;

    dynamics->initial_preprocess();

    if (DBG >= 1) cout << "DBG1: dynamics->main_stream()" << endl;
    dynamics->main_stream();

    // dynamics->mmsys.writeData();
    //cout << "dbg20180320 a" << endl;
    dynamics->terminal_process();
    //cout << "dbg20160625 a" << endl;
#if defined(F_MPI)
    MPI_Finalize();
#endif
    //cout << "dbg20160625 b"<< endl;
    delete dynamics;
    cout << "terminate Celeste"<< endl;
    return 0;
}
