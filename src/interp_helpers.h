#ifndef SURVDISTR_INTERP_HELPER_H
#define SURVDISTR_INTERP_HELPER_H

#include <string>

enum InterpMethod {
  CONST_SURV,
  CONST_DENS,
  CONST_HAZ
};

inline InterpMethod parse_method(const std::string& method) {
  if (method == "const_surv") return CONST_SURV;
  if (method == "const_dens") return CONST_DENS;
  if (method == "const_haz")  return CONST_HAZ;
  Rcpp::stop("Unknown interpolation method.");
}

#endif
