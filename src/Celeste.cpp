#include "Celeste.h"

Celeste::Celeste()
  : CelesteObject(){
  
}

int Celeste::setup(int argn, char* argv[]){
  cout << ABOUT_ME << endl;
  string fn_cfg;
  if(argn<2){
    cerr << "Usage: celeste [mode]"<<endl;
    cerr << "------------------------------------"<<endl;
    exit(1);
  }
  cout << "conf.setall\n";
  cfg.set_defaults();
  cfg.setAll(argn,argv);
  if(cfg.fn_cfg!=string())
    cfg.setAll(Read(cfg.fn_cfg).load_config());
  cout <<"/setup\n";
  
  return 0;
}

int Celeste::main_stream(){

  cout<<"mainstream\n";
  switch(cfg.mode){
  case M_TEST: test_mode();           break;
  case M_DYNAMICS: dynamics_mode();   break;
  default:
    cout <<"Invalid Mode is specified.\n";
    test_mode();
    break;
  }
  cout <<"Terminated.\n";
  return 0;
}

int Celeste::test_mode(){
  cout<<"test_mode.\n";
  MmSystem mmsys;
  Read(cfg.fn_inp).load_launch_set(mmsys);
  //mmsys.write_data();
  return 0;
}

int Celeste::dynamics_mode(){
  cout<<"dynamics_mode\n";
  
  // #if defined(F_CUDA) && defined(F_MPI)
  //   cout << "F_CUDA + F_MPI flags = ON" << endl;
  //  MpiGpuDynamicsMode* dynamics = new MpiGpuDynamicsMode;
  //GpuDynamicsMode* dynamics = new GpuDynamicsMode;
  //  MPI_Init(NULL, NULL);
  //#if defined(F_CUDA)
  //  cout << "F_CUDA flag = ON" << endl;
  //  GpuDynamicsMode* dynamics = new GpuDynamicsMode;
  //#else
  DynamicsMode* dynamics;
  if (cfg.integrator_type == INTGRTR_LEAPFROG_PRESTO){
    dynamics = new DynamicsModePresto();
  }else if (cfg.integrator_type == INTGRTR_ZHANG){
    dynamics = new DynamicsModeZhang();
  }else{
    cout << "Unknown Integrator" << endl;
  }
  //#endif
  
  Read(cfg.fn_inp).load_launch_set(dynamics->mmsys);
  
  if(DBG>=1)
    cout << "DBG1: dynamics->set_config_parameters(cfg)" << endl;
  
  dynamics->set_config_parameters(&cfg);
  
  if(DBG>=1)
    cout << "DBG1: dynamics->initial_preprocess()" << endl;
  
  dynamics->initial_preprocess();
  
  if(DBG>=1)
    cout << "DBG1: dynamics->main_stream()" << endl;
  dynamics->main_stream();
  
  //dynamics->mmsys.writeData();
  
  dynamics->terminal_process();  

#if defined(F_MPI)
  MPI_Finalize();
#endif
  delete dynamics;
  return 0;
}


