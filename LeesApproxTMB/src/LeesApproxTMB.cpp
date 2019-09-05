#include <TMB.hpp>
#include "../inst/include/LeesApprox_TMB_fn.h"
#include "../inst/include/ns_SCA.h"

template <class Type>
Type objective_function<Type>::operator() () {

  DATA_STRING(model);
  if(model == "LeesApprox_internal") {
    #include "../inst/include/LeesApprox_internal.h"
  } else {
    #include "../inst/include/LeesApprox_SCA.h"
  }
  return 0;
}

