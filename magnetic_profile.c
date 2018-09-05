#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_POINTS 10
#define SIZE 10
#define R_SOL 695700
#define sec(x) while(1) {1/cos(x); break;}

typedef struct vector Vector;
typedef struct cartesian_coordinates Cartesian_Coordinates;
typedef struct field_values_c Field_Values_C;
typedef struct field_values_p Field_Values_P;
typedef struct polar_coordinates Polar_Coordinates;

typedef struct cartesian_coordinates{
	double x;
	double y;
} Cartesian_Coordinates;

typedef struct polar_coordinates{
	double r;
	double theta;
} Polar_Coordinates;

typedef struct field_values_c{
	double bx;
	double by;
	double bz;
} Field_Values_C;

typedef struct field_values_p{
	double br;
	double bt;
	double bp;
} Field_Values_P;

typedef struct magnetic_vector{
	Cartesian_Coordinates *location;
	Polar_Coordinates *polar_location;
	Field_Values_C *cartesian_field;
	Field_Values_P *polar_field;
	int size;
} Magnetic_Vector;

void print_array(double *arr, int size){
	int i;
	for(i = 0; i < size; i++){
		printf("%f  ", *(arr+i));
	}
	printf("\n");
}

void print_happy(Cartesian_Coordinates *arr, Field_Values_C *field, int size, FILE **out_file){
	int i;
	for(i = 0; i < size; i++){
		fprintf(*out_file, "%f %f %f %f\n", (arr+i)->x, (arr+i)->y, (field+i)->bx, (field+i)->by);
	}
}

void print_coord(Cartesian_Coordinates *arr, int size, FILE **out_file){
	int i;
	
	for(i = 0; i < size; i++){
		fprintf(*out_file, "x = %f  y = %f\n", (arr+i)->x, (arr+i)->y);
	}
	fprintf(*out_file,"\n");
}

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
//		printf("r = %f theta = %f h = %f \n", (coordinates+i)->r, (coordinates+i)->theta, h );
	}
	return size;
}


int calculate_field_values(Magnetic_Vector magnetic_in, double m, int size){
	int i;
	double a, b, b_total, B;
	Polar_Coordinates *polar_field;

	polar_field = magnetic_in.polar_location;

	//Calculo br y bt
	for(i = 0; i < size; i++){
		(magnetic_in.polar_field + i)->br = 
			 (2*m*cos( (polar_field + i)->theta ) ) / pow( (polar_field +i)->r ,3);
		(magnetic_in.polar_field + i)->bt = 
			 (m * sin( (polar_field + i)->theta ) ) / pow( (polar_field +i)->r,3);
		(magnetic_in.polar_field + i)->bp = 0;
//		printf( "br = %f  bt = %f  \n", (magnetic_in.polar_field+i)->br, (magnetic_in.polar_field+i)->bt );

		B = pow((magnetic_in.polar_field+i)->br,2) + pow((magnetic_in.polar_field+i)->bt,2);
		printf("B^2 = %f\n", B );
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
//		printf( "%f   \n", (magnetic_in.cartesian_field+i)->by );
	}

	//Paso valores de r y theta a x, y
	for(i = 0; i < size; i++){
		(magnetic_in.location + i)->y = 
			 sin( (polar_field + i)->theta ) * (polar_field +i)->r;
		(magnetic_in.location + i)->x = 
			 cos( (polar_field + i)->theta ) * (polar_field +i)->r;
	//	printf( "%f   %f   %i\n", (polar_field+i)->r, (polar_field+i)->theta, i);
	}
}

int create_semicirle(Polar_Coordinates *coordinates, double r_in, int size){
	int i;
	double h, l, CDp, DDp, BDp;

	for( i = 0; i <= size; i++){
		(coordinates + i) ->r = r_in;
		(coordinates + i) ->theta = i* .3141;
//		printf("%f %f\n", (coordinates+i)->r, (coordinates+i)->theta);
	}
	return size;
}

void generate_profiles(double *arr, double alpha, double no, Polar_Coordinates *coordinates, int size, Magnetic_Vector magnetic_in, double si, int min_d, int max_d, int step_d, long int min_m, long int max_m, long int step_m)
{
	FILE *out_file;
	int i, momentum, j, m, d;
	char filename[50];
	double B_total;
	momentum = (int)m;
	d=R_SOL;
	max_d += R_SOL;

	for(d += min_d; d <= max_d; d+= step_d)
	{
		for(m = min_m; m <= max_m; m += step_m)
		{
			sprintf(filename, "datos/m_%i/d_%i.dat", momentum,d);
			calculate_vector_position(arr, alpha, d, coordinates, size);
			calculate_field_values(magnetic_in, m, size);
			out_file = fopen(filename, "w+");
			print_happy(magnetic_in.location, magnetic_in.cartesian_field, size, &out_file);
			fclose(out_file);
		}
	}

}


int main(){

	//arr es altura en escala distancia/R_SOL 
	Magnetic_Vector my_vec;
	Cartesian_Coordinates x0;
	int i, stages;
	double arr[N_POINTS] = {1000, 1100, 1200, 1300, 1400, 1500, 1600, 1900, 2000, 2100};
	int size = SIZE;
	double alpha = .0000006;
	double m = 500000000;
	FILE *out_file;
	double d, distance32;

	d = R_SOL;

	/*
	printf("alpha? ");
	scanf("%lf", &alpha);
	*/


	scanf("%lf",&distance32);
	scanf("%lf",&m);
	d += distance32;

	printf("base_distance = %lf\n", distance32);
	my_vec.location = (Cartesian_Coordinates*)malloc(2 + size*sizeof(Cartesian_Coordinates));
	my_vec.polar_location = (Polar_Coordinates*)malloc(2+size*sizeof(Polar_Coordinates));
	my_vec.cartesian_field = (Field_Values_C*)malloc(2+size*sizeof(Field_Values_C));
	my_vec.polar_field = (Field_Values_P*)malloc(2+size*sizeof(Field_Values_P));


	calculate_vector_position(arr, alpha, d, my_vec.polar_location, size);
	calculate_field_values(my_vec, m, size);


	out_file = fopen("calculated_values.txt", "w+");
	//Inserta valores del vector
	print_happy(my_vec.location, my_vec.cartesian_field, size, &out_file);
/*
	printf("Quantity of stages to be displayed? ");
	//Esto lo uso para saber cuantos semicirculos quiero
	scanf("%i", &stages);
*/
	stages = 0;
	for(i = 0; i < stages; i++)
	{
		create_semicirle(my_vec.polar_location, 1100 + i*200, size);
		calculate_field_values(my_vec, m, size);
		print_happy(my_vec.location, my_vec.cartesian_field, size, &out_file);
	}

	fclose(out_file);


	return 0;
}
