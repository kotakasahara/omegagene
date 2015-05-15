#include "Expand.h"

Expand::Expand()
  : CelesteObject(){
  write_lambda_interval = 0;
}

Expand::~Expand(){
}

void Expand::set_lambda_interval(int in_lambda_interval){
  write_lambda_interval = in_lambda_interval;
}


