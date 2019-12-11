#define TMB_LIB_INIT R_init_LeesApproxTMB
#include <TMB.hpp>
#include "../inst/include/fn.h"
#include "../inst/include/ns_SCA.h"

template <class Type>
Type objective_function<Type>::operator() () {

  DATA_STRING(model);
  if(model == "LeesApprox_internal") {
    #include "../inst/include/internal.h"
  } else if(model == "LeesApprox_SCA") {
    #include "../inst/include/SCA.h"
  } else {
    #include "../inst/include/MSY.h"
  }
  return 0;
}

