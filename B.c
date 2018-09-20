#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "B.h"

#define N_POINTS 10
#define R_SOL 695700
#define sec(x) while(1) {1/cos(x); break;}

int calculate_vector_position(double *arr, double alpha, double d, Polar_Coordinates *coordinates, int size){
	int i;
	double h, a, b, w;
	b = d;

	for( i = 0; i < size; i++){
		h = *(arr+i);
		a = h + R_SOL;
		(coordinates + i) ->r = 
			sqrt( pow(a,2) + pow(b,2) - 2*a*b*cos(alpha) );
		w = asin( b*sin(alpha) / (coordinates + i)->r );
		(coordinates + i) ->theta = 
			M_PI/2 - alpha - w;
	}
	return size;
}

int calculate_point_position(double z, double alpha, double d, Polar_Coordinates *coordinates, double x0)
{
	x0 = 5;
	double w,a, b, my_tmp, asqrt, bsqrt, sum, dif82, tmp82;
	
	my_tmp = R_SOL;
	
	b = d;
	a = z + R_SOL;
	//asqrt = pow(a,2.0);
	//bsqrt = pow(b,2.0);
	//tmp82 = -2*a*b*cos(alpha);
	//sum = asqrt + bsqrt;
	//dif82 = sum+tmp82;
	coordinates->r = sqrt( pow(a,2.0) + pow(b,2.0) -2*a*b*cos(alpha) );
	w = b*sin(alpha) / coordinates->r;
	if(w < 1)
	{
		w = asin( w );
	}
	else
	{
		w = 0.001;
	}

	if(coordinates->r < 50)
	{
		coordinates->r = 50;
	}

	coordinates->theta = M_PI/2 - alpha - w;

	printf("z=%lf a=%lf b=%lf asqrt=%le bsqrt=%le dif=%le r=%le\n", z, a, b, asqrt, bsqrt, dif82, coordinates->r);
	if(z != 0){
		coordinates->theta = M_PI/2 - alpha - w;
	}
	else
	{
		coordinates->theta = 1.568;
	}

	/*
	printf("z= %f alpha = %f d = %f\n",z, alpha, d);
	printf("r= %f theta = %f\n",coordinates->r, coordinates->theta);
	printf("a= %f b = %f\n",a, b);
	printf("asqrt = %f bsqrt= %f sum=%f\n",asqrt, bsqrt, sum);
	printf("neg_val = %f dif= %f\n", tmp82, dif82);
	*/

	return 0;
}

int calculate_field_values(Magnetic_Vector magnetic_in, double m, int size){
	int i;
	double a, b, B;
	Polar_Coordinates *polar_field;

	polar_field = magnetic_in.polar_location;

	//Calculo br y bt
	for(i = 0; i < size; i++){
		(magnetic_in.polar_field + i)->br = 
			 (2*m*cos( (polar_field + i)->theta ) ) / pow( (polar_field +i)->r ,3);
		(magnetic_in.polar_field + i)->bt = 
			 (m * sin( (polar_field + i)->theta ) ) / pow( (polar_field +i)->r,3);
		(magnetic_in.polar_field + i)->bp = 0;
		B = pow((magnetic_in.polar_field+i)->br,2) + pow((magnetic_in.polar_field+i)->bt,2);
		printf( "r=%f  t=%f  Bsq=%f\n", (polar_field+i)->r, (polar_field+i)->theta, B);
//		printf( "br = %f  bt = %f  \n", (magnetic_in.polar_field+i)->br, (magnetic_in.polar_field+i)->bt );
	}

	//Calculo bx y by
	for(i = 0; i < size; i++){
		(magnetic_in.cartesian_field + i)->bx = 
			(magnetic_in.polar_field +i)->br * cos( (polar_field + i)->theta ) - 
			(magnetic_in.polar_field +i)->bt * sin( (polar_field+i)->theta );
		(magnetic_in.cartesian_field + i)->by = 
			(magnetic_in.polar_field +i)->br * sin( (polar_field + i)->theta ) + 
			(magnetic_in.polar_field +i)->bt * cos( (polar_field+i)->theta );
		(magnetic_in.cartesian_field + i)->bz = 0;
//		printf( "size =%i bx = %f by = %f \n",size,(magnetic_in.cartesian_field+i)->bx, (magnetic_in.cartesian_field+i)->by);
	}

	//Paso valores de r y theta a x, y
	for(i = 0; i < size; i++){
		(magnetic_in.location + i)->y = 
			 sin( (polar_field + i)->theta ) * (polar_field +i)->r;
		(magnetic_in.location + i)->x = 
			 cos( (polar_field + i)->theta ) * (polar_field +i)->r;
	//	printf( "%f   %f   %i\n", (polar_field+i)->r, (polar_field+i)->theta, i);
	}
	return 0;
}


void init_B(double base, double intensity, double alpha, Model *model, int x0)
{
	FILE *magnetic_f;
	char fullpath[300];
	Magnetic_Vector my_vec;
	double d;
	int size = 143;
//	double m = 2596400000000;
	double m;
	double arr[200];
	double B_total;
	Field_Values_C *result;
	strcpy(fullpath,"data/atmosphere/hydrostatic/C07/magnetic.dat");
	magnetic_f = fopen(fullpath,"w");

	m = intensity;
	d = R_SOL;
	x0 = 10;
	d+=base;

	my_vec.location = (Cartesian_Coordinates*)malloc(size*sizeof(Cartesian_Coordinates));
	my_vec.polar_location = (Polar_Coordinates*)malloc(size*sizeof(Polar_Coordinates));
	my_vec.cartesian_field = (Field_Values_C*)malloc(size*sizeof(Field_Values_C));
	my_vec.polar_field = (Field_Values_P*)malloc(size*sizeof(Field_Values_P));

	//calculate_vector_position(arr,alpha,d,my_vec.polar_location,size);
	//printf("alpha=%le\n",alpha);
	calculate_point_position(model->atm.z, alpha, d, my_vec.polar_location, x0);
	calculate_field_values(my_vec,m,1);
	model->atm.xi = .7*exp( -.00035294 * model->atm.z );
//	B_total = pow((my_vec.polar_field)->br,2) + pow((my_vec.polar_field)->bt,2);
//	printf("B_total = %f\n",B_total);

	if (magnetic_f==NULL){
		printf("Error 82: File %s could not be created.\n",fullpath);
		exit(0);
	}

	//fprintf(magnetic_f,"%s\n%le %le %le %le\n","# Header file to hidrostatic model\n# By: Victor H De la Luz\n# vdelaluz@geofisica.unam.mx\n# 06/09/2008\n# Dummy Pakage to prove te read data\n# id: Identification model\n# z : Height over the photosphere (km)\n# T : Temperature (k)\n# P : Presure\n# H : Hidrogen density (cm-3)\n# [ne]: Calculada a partir de las especies\n# V : Doppler Velocity\n# vt: Turbulent Velocity\n#z		Bx		By		Bz", (*model).atm.z, my_vec.cartesian_field->bx, my_vec.cartesian_field->by, my_vec.cartesian_field->bz);


	fprintf(magnetic_f,"%le %le %le %le\n", (*model).atm.z, my_vec.cartesian_field->bx, my_vec.cartesian_field->by, my_vec.cartesian_field->bz);


	fclose(magnetic_f);


	(*model).atm.Bx = my_vec.cartesian_field->bx;
	(*model).atm.By = my_vec.cartesian_field->by;
	(*model).atm.Bz = my_vec.cartesian_field->bz;
	
	free(my_vec.location);
	free(my_vec.polar_location);
	free(my_vec.cartesian_field);
	free(my_vec.polar_field);
}

void calculate_B(double base, double intensity, double alpha, Model *model, int x0)
{
	//amplitude = distancia del radio solar a donde nace el campo
	//intensity = m
	//alpha = angulo entre la normal y el nacimiento
	FILE *magnetic_f;
	char fullpath[300];
	Magnetic_Vector my_vec;
	double d;
	int size = 143;
//	double m = 2596400000000;
	double m;
	double arr[200];
	Field_Values_C *result;
	double B_total;

	strcpy(fullpath,"data/atmosphere/hydrostatic/");
	strcat(fullpath,"C07");
	strcat(fullpath,"/magnetic.dat");
	
	m = intensity;
	d = R_SOL;
	x0 = 10;
	d+= base;

	x0 = 10;
	my_vec.location = (Cartesian_Coordinates*)malloc(size*sizeof(Cartesian_Coordinates));
	my_vec.polar_location = (Polar_Coordinates*)malloc(size*sizeof(Polar_Coordinates));
	my_vec.cartesian_field = (Field_Values_C*)malloc(size*sizeof(Field_Values_C));
	my_vec.polar_field = (Field_Values_P*)malloc(size*sizeof(Field_Values_P));

	//calculate_vector_position(arr,alpha,d,my_vec.polar_location,size);
	calculate_point_position(model->atm.z, alpha, d, my_vec.polar_location, x0);
	calculate_field_values(my_vec,m,1);

//	B_total = pow((my_vec.polar_field)->br,2) + pow((my_vec.polar_field)->bt,2);
//	printf("B_total = %f\n",B_total);

	//ecuacion de XI
	model->atm.xi = .7*exp( -.00035294 * model->atm.z );
//	printf("Xi = %le\n", model->atm.xi);


	magnetic_f = fopen(fullpath,"a");

	if (magnetic_f==NULL)
	{
		printf("Error 82: File %s not found.\n",fullpath);
		exit(0);
	}

	fprintf(magnetic_f, "%le %le %le %le\n", (*model).atm.z, my_vec.cartesian_field->bx, my_vec.cartesian_field->by, my_vec.cartesian_field->bz);
	fclose(magnetic_f);

	(*model).atm.Bx = my_vec.cartesian_field->bx;
	(*model).atm.By = my_vec.cartesian_field->by;
	(*model).atm.Bz = my_vec.cartesian_field->bz;
	
	free(my_vec.location);
	free(my_vec.polar_location);
	free(my_vec.cartesian_field);
	free(my_vec.polar_field);
/*
	for(i = 0; i < size; i++){
//		fprintf(magnetic_f, "%le %le %le\n",(my_vec.cartesian_field + i)->bx,(my_vec.cartesian_field + i)->by,(my_vec.cartesian_field + i)->bz);
	}
*/
}
