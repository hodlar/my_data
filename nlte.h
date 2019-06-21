#ifndef NLTE_H
#define NLTE_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
#include "hmodel.h"

  double f_hydro(double Y, double R, double Z, double T, double vt, double bx, double by, double bz, double xi);
  void NLTE(Model *model, double error, int hydro, double pz1, double fz1,double dx, int chromospheric_network,int cell, int magnetic);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* FUNCTIONS_H */

