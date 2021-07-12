#include "Config.h"
#include <iostream>
using namespace std;

Config::Config(){
  mode=M_DUMMY;
}

Config::~Config(){
  mc_target_acc_ratio = 0.0;
  greedy_pivot = 0;
  target_error = 1e-8;
  mc_delta_x_max = 1.0;
  mc_delta_x = 0.1;
  qweight_write_mode = QW_FILE_MODE_RAW;
  process_unsampled_zone = UNSAMPLE_OMIT;

  mc_n_window_trend = 0;
  mc_error_ave_window_size = 0;
  mc_max_temp = 0;
  mc_min_temp = 0;
  mc_delta_temp = 0;

}
void Config::setAll(int argn, char* argv[]){
  vector<string> arg;
  int i;
  for(i=1;i<argn;i++)
    arg.push_back(string(argv[i]));
  setAll(arg);
}

void Config::setAll(vector<string> arg){
  vector<string>::iterator itr;
  string type,val;
  for(itr=arg.begin(); itr!=arg.end(); itr++){
    if(*itr=="--mode"){
      itr++;
      if(*itr=="test")         { mode=M_TEST; }
      else if(*itr=="subzone") { mode=M_SUBZONEBASED; }
      else if(*itr=="")        { mode=M_DUMMY; }
      else{
        cerr<<"invalid mode ["<<*itr<<"]\n"; exit(1);
      }
    }
    else if(*itr=="--i-cfg"){ itr++; fname_i_cfg=(*itr); }
    else if(*itr=="--i-params"){ itr++; fname_i_params=(*itr); }
    else if(*itr=="--i-qraw-is"){ itr++; fname_i_qraw_is=(*itr); }
    else if(*itr=="--o-qcano"){ itr++; fname_o_qcano=(*itr); }
    else if(*itr=="--o-qweight-opt"){ itr++; fname_o_qweight_opt=(*itr); }
    else if(*itr=="--qweight-write-mode"){ itr++;
      if  ((*itr) == "RAW")
	qweight_write_mode = QW_FILE_MODE_RAW;
      else if ((*itr) == "LOG")
	qweight_write_mode = QW_FILE_MODE_LOG;
      else{
	cerr << "invalid option for --qweight-write-mode " << (*itr) << endl;
      }
      //cout << "dbg config qweight_write_mode " << qweight_write_mode << endl;
    }
    else if(*itr=="--process-unsapmled-zone"){ itr++;
      if  ((*itr) == "OMIT")
	process_unsampled_zone = UNSAMPLE_OMIT;
      else if ((*itr) == "MIN")
	process_unsampled_zone = UNSAMPLE_MIN;
      else{
	cerr << "invalid option for --process-unsampled-zone " << (*itr) << endl;
      }
      cout << "dbg config qweight_write_mode " << qweight_write_mode << endl;
    }


    else if(*itr=="--target-error"){ itr++; target_error=atof((*itr).c_str()); }
    else if(*itr=="--mc-temp"){ itr++; mc_temp=atof((*itr).c_str()); }
    else if(*itr=="--mc-delta-x"){ itr++; mc_delta_x=atof((*itr).c_str()); }
    else if(*itr=="--mc-delta-x-max"){ itr++; mc_delta_x_max=atof((*itr).c_str()); }
    else if(*itr=="--mc-steps"){ itr++; mc_steps = atoi((*itr).c_str()); }
    else if(*itr=="--mc-log-interval"){ itr++; mc_log_interval = atoi((*itr).c_str()); }
    else if(*itr=="--mc-target-acc-ratio"){ itr++; mc_target_acc_ratio = atof((*itr).c_str());}
    else if(*itr=="--mc-acc-duration"){ itr++; mc_acc_duration = atoi((*itr).c_str());}
    else if(*itr=="--mc-n-window-trend"){ itr++; mc_n_window_trend = atoi((*itr).c_str());}
    else if(*itr=="--mc-error-ave-window-size"){ itr++; mc_error_ave_window_size = atoi((*itr).c_str());}
    else if(*itr=="--mc-min-temp"){ itr++; mc_min_temp = atof((*itr).c_str());}
    else if(*itr=="--mc-max-temp"){ itr++; mc_max_temp = atof((*itr).c_str());}
    else if(*itr=="--mc-delta-temp"){ itr++; mc_delta_temp = atof((*itr).c_str());}
    else if(*itr=="--greedy-max-steps"){ itr++; greedy_max_steps = atoi((*itr).c_str());}
    else if(*itr=="--greedy-pivot"){ itr++; greedy_pivot = atoi((*itr).c_str());}
    else{
      cerr<<"unknown keyword <"<<*itr<<">"<<endl;
      exit(1);
    }
  }
}
