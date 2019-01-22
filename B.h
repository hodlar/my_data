#ifndef B_H
#define B_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "hmodel.h"

	typedef struct{
		double x;
		double y;
	}Cartesian_Coordinates;

	typedef struct{
		double r;
		double theta;
	}Polar_Coordinates;

	typedef struct{
		double bx;
		double by;
		double bz;
	}Field_Values_C;

	typedef struct{
		double br;
		double bt;
		double bp;
	}Field_Values_P;

	typedef struct{
		Cartesian_Coordinates *location;
		Polar_Coordinates *polar_location;
		Field_Values_C *cartesian_field;
		Field_Values_P *polar_field;
		int size;
	}Magnetic_Vector;

	void init_B(double base, double intensity, double alpha, Model *model, int x0);
	void calculate_B(double base, double intensity, double alpha, Model *model, int x0);
	double caltulate_Xi(Atmosphere *layer);
	double calculate_field_values(Magnetic_Vector magnetic_in, double m, int size);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* FUNCTIONS_H */

