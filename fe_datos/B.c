#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "B.h"

#define N_POINTS 10
#define R_SOL 695700
#define sec(x) while(1) {1/cos(x); break;}


double calculate_field_values(double height, double magnetic_flux, double base_field, double scale_height){
	double Bz0;

	height = height + 100;
	if(scale_height == 0)
	{scale_height = 1;}
	Bz0 = base_field * exp(-height/scale_height);

	printf("\nBz0 = %le base_field = %le height = %le\n", Bz0, base_field, height);
	return Bz0;
}


void init_B(double base, double intensity, double alpha, Model *model, int x0)
{
	//scale height = intensity
	//base_field = base
	//alpha magnetic flux
	//h es cercana a 200
	FILE *magnetic_f;
	char fullpath[300];
	Magnetic_Vector my_vec;
	double d;
	int size = 143;
	double m;
	double arr[200];
	double B_total;
	Field_Values_C *result;
	strcpy(fullpath,"data/atmosphere/hydrostatic/C07/magnetic.dat");
	magnetic_f = fopen(fullpath,"w");
	double magnetic_flux = 2.8e18;
	double magnetic_field_vertical;

	magnetic_flux = alpha;

	my_vec.location = (Cartesian_Coordinates*)malloc(size*sizeof(Cartesian_Coordinates));
	my_vec.polar_location = (Polar_Coordinates*)malloc(size*sizeof(Polar_Coordinates));
	my_vec.cartesian_field = (Field_Values_C*)malloc(size*sizeof(Field_Values_C));
	my_vec.polar_field = (Field_Values_P*)malloc(size*sizeof(Field_Values_P));

	magnetic_field_vertical = calculate_field_values(model->atm.z, magnetic_flux, base, intensity);
	model->atm.xi = 1;

	if (magnetic_f==NULL){
		printf("Error 82: File %s could not be created.\n",fullpath);
		exit(0);
	}

	fprintf(magnetic_f,"%le %le\n", (*model).atm.z, magnetic_field_vertical);
	fclose(magnetic_f);


	(*model).atm.Bx = 0;
	(*model).atm.By = 0;
	(*model).atm.Bz = magnetic_field_vertical;
	
	free(my_vec.location);
	free(my_vec.polar_location);
	free(my_vec.cartesian_field);
	free(my_vec.polar_field);
}

void calculate_B(double base, double intensity, double alpha, Model *model, int x0)
{
	//scale height = intensity
	//base_field = base
	//alpha unused
	FILE *magnetic_f;
	char fullpath[300];
	Magnetic_Vector my_vec;
	double d;
	int size = 143;
	double m;
	double arr[200];
	Field_Values_C *result;
	double B_total;
	double magnetic_flux = 2.8e18;
	double magnetic_field_vertical;

	magnetic_flux = alpha;

	strcpy(fullpath,"data/atmosphere/hydrostatic/");
	strcat(fullpath,"C07");
	strcat(fullpath,"/magnetic.dat");

	my_vec.location = (Cartesian_Coordinates*)malloc(size*sizeof(Cartesian_Coordinates));
	my_vec.polar_location = (Polar_Coordinates*)malloc(size*sizeof(Polar_Coordinates));
	my_vec.cartesian_field = (Field_Values_C*)malloc(size*sizeof(Field_Values_C));
	my_vec.polar_field = (Field_Values_P*)malloc(size*sizeof(Field_Values_P));

	magnetic_field_vertical = calculate_field_values(model->atm.z, magnetic_flux, base, intensity);
	model->atm.xi = 1;


	magnetic_f = fopen(fullpath,"a");

	if (magnetic_f==NULL)
	{
		printf("Error 82: File %s not found.\n",fullpath);
		exit(0);
	}
	fprintf(magnetic_f,"%le %le\n", (*model).atm.z, magnetic_field_vertical);
	fclose(magnetic_f);

	(*model).atm.Bx = 0;
	(*model).atm.By = 0;
	(*model).atm.Bz = magnetic_field_vertical;
	
	free(my_vec.location);
	free(my_vec.polar_location);
	free(my_vec.cartesian_field);
	free(my_vec.polar_field);
}
