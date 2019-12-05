#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "B.h"

#define N_POINTS 10
#define R_SOL 695700
#define sec(x) while(1) {1/cos(x); break;}

double  calculate_field_values(double height, double magnetic_flux, double base_field, double scale_height){
	double B_total;
	double Bz0, r0z, alpha, d_alpha;

	Bz0 = base_field*exp(-height/scale_height);
//	r0z = sqrt( (3*magnetic_flux)/(M_PI*Bz0) );
//	alpha = magnetic_in.location->x/r0z;

	if(alpha <= 1)
	{
		d_alpha = pow(1 - alpha*alpha, 2);
	}
	else
	{
		d_alpha = 0;
	}

	printf( "bz = %f  \n", Bz0 );


	return Bz0;
}


void init_B(double scale_height, double base_field, Model *model, int x0)
{
	FILE *magnetic_f;
	char fullpath[300];
	int size = 143;
//	double m = 2596400000000;
	double arr[200];
	double magnetic_field_z;
	double magnetic_flux;
	strcpy(fullpath,"data/atmosphere/hydrostatic/C07/magnetic.dat");
	magnetic_f = fopen(fullpath,"w");

	magnetic_flux = 2.8e18;

	magnetic_field_z = calculate_field_values(model->atm.z, magnetic_flux, base_field, scale_height);

	if (magnetic_f==NULL){
		printf("Error 82: File %s could not be created.\n",fullpath);
		exit(0);
	}

	fprintf(magnetic_f,"%le %le\n", (*model).atm.z, magnetic_field_z);

	fclose(magnetic_f);

	(*model).atm.Bx = 0;
	(*model).atm.By = 0;
	(*model).atm.Bz = magnetic_field_z;
}

void calculate_B(double scale_height, double base_field, Model *model, int x0)

{
	//amplitude = distancia del radio solar a donde nace el campo
	//intensity = m
	//alpha = angulo entre la normal y el nacimiento
	
	FILE *magnetic_f;
	char fullpath[300];
	int size = 143;
//	double m = 2596400000000;
	double arr[200];
	double magnetic_field_z;
	double magnetic_flux;

	strcpy(fullpath,"data/atmosphere/hydrostatic/");
	strcat(fullpath,"C07");
	strcat(fullpath,"/magnetic.dat");

	magnetic_flux = 2.8e18;

	magnetic_field_z = calculate_field_values(model->atm.z, magnetic_flux, base_field, scale_height);

	if (magnetic_f==NULL){
		printf("Error 82: File %s could not be created.\n",fullpath);
		exit(0);
	}

	fprintf(magnetic_f,"%le %le\n", (*model).atm.z, magnetic_field_z);

	fclose(magnetic_f);

	(*model).atm.Bx = 0;
	(*model).atm.By = 0;
	(*model).atm.Bz = magnetic_field_z;
}
