#ifndef HMODEL_H
#define HMODEL_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "atom.h"

typedef struct{
	int id;
	double z;
	double T;
	double P;
	double H;
	double V;
	double vt;
	double ne;
	double ne_lte;
	double bhm; //departure coefficient of H- (TEMPORAL)
	double fz; //reciprocal scale of height eq 93 vernazza 73
	double Bx;
	double By;
	double Bz;
	double xi;
}Atmosphere;


typedef struct{
  int n; //niveles en la atmosfera
  int natom; //numero de atomos en el modelo
  char name[100]; //nombre del modelo.
  Atom *atom;
  Atmosphere atm;
}Model;


  Model newModel(char *name);
  void loadInitValues(Model *modelref, int n);
  void writeModel(Model model, char *name);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* FUNCTIONS_H */

