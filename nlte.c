#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "nlte.h"
#include "constants.h"
#include "LU.h"
#include "utils.h"
#include "depcoef.h"
#include "atom.h"
#include "colisionalexcitation.h"
#include "colisionalionization.h"
#include "energylevels.h"
#include "lines.h"
#include "basicparameters.h"
#include "colisionalexcitation.h"
#include "hmodel.h"
#include "specie.h"

#include <string.h>


//#define GMH 1.11691659e-31
//#define mH  1.673534e-24
//#define kboltz  1.380658e-16

double totalParticles =0.0;
double magnetic_true = 0;
int first_time = 0;

double f_hydro(double Y, double R, double Z, double T, double vt, double bx, double by, double bz, double xi){
 //g * m_H 6.674e-8 cm^3g-1s-2 * 1.673534e-24 g
 // double GMH = 1.11691659e-31;
  //  double mH = 1.673534e-24;
  // double kboltz = 1.380658e-16; //erg K^-1
   double R1,KT,R2,R3,R4,R5,R6,B, magnetic_field;

   char tmp_name[300];
   FILE *tmp_file;
   strcpy(tmp_name,"magnetic_and_pressure.dat");

printf("aqui se rompe el asunto\n\n\n\n");

   magnetic_field = bz;
  
   if(first_time == 0)
   {
	   tmp_file = fopen(tmp_name,"w");
	   magnetic_true = magnetic_field;
	   printf("magnetic_field=%le\n",magnetic_field);
	   fprintf(tmp_file, "%le   %le\n", 0.0, magnetic_field);
	   fclose(tmp_file);
   }

   if(magnetic_true != magnetic_field && first_time)
   {
	   tmp_file = fopen(tmp_name,"a");
	   magnetic_true = magnetic_field;
	   printf("magnetic_field=%le\n",magnetic_field);
	   fprintf(tmp_file, "%le   %le\n", 0.0, xi*magnetic_field/(8*M_PI));
	   fclose(tmp_file);
   }
   first_time = 1;
   B = 1.0+R+Y+Z+( (0.5*pow(vt,2.0)*(1.0+4.0*Y)*mH)/(kboltz*T)) + xi*magnetic_field/(8*M_PI);

//   B = 1.0+R+Y+Z+( (0.5*pow(vt,2.0)*(1.0+4.0*Y)*mH)/(kboltz*T));

   return exp(log(GMH) - log(kboltz*T) + log(1.0+4.0*Y) - log (B));

  /*  R1 = pow(vt,2.0);
R2 = 
  */
  /*  R1 = -30.517684+log(1.0+4.0*Y);
  KT = (1.380658e-16)*T;
  R3 = KT+R*KT+Y*KT+Z*KT;
  R4 = 8.36767e-25*pow(vt,2.0) + (3.347068e-24*pow(vt,2.0))*Y;
  R5 = log(R3+R4);
  R6 = R1-R5;
  */
  // printf("R6=%le\n",pow(10.0,R6));
  //  return pow(10.0,R6);
    //GMH*(1.0+4.0*Y)/((1.0+R+Y+Z)*kboltz*T+ (0.5)*mH*(1.0+4.0*Y)*vt*vt);
}

//esta funcion calcula la cantidad total
//de electrones a partir de los valores
//iniciales, recordemos que en el modelo
//distribuimos los electrones entre cada
//especie 
double ComputeAbsoluteNe(Model *model){
  // Model model = *modelrf;
  int i,j;
  int n = (*model).natom;
  double ne=0.0;

  /***** SYSTEM ********
  printf("Loading electronic density absolute contributions:\n");
  *****************/

  (*model).atm.ne = 0.0;
  for (i=0; i < n;i++){
    ne=0.0;
    //printf("%s\n",(*model).atom[i].name);
    for (j=0; j< (*model).atom[i].ns;j++){
      /*
      (*model).atom[i].specie[j].ne= (*model).atom[i].specie[j].nk*(*model).atom[i].specie[j].nrel*(*model).atom[i].A*(*model).atm.H* ((double)(*model).atom[i].specie[j].ionizationstage);
      */
      (*model).atom[i].specie[j].ne= (*model).atom[i].specie[j].nabs*((double)(*model).atom[i].specie[j].ionizationstage);

      (*model).atom[i].specie[j].ne_lte=(*model).atom[i].specie[j].ne;
      ne += (*model).atom[i].specie[j].ne;

      /********* SYSTEM ***********
      printf("%s:%le\n",(*model).atom[i].specie[j].name,(*model).atom[i].specie[j].ne);
      ***********************/

    //printf("ne:%i\n",(*model).atom[i].specie[j].ionizationstage);
    }
    (*model).atom[i].ne = ne;
    (*model).atom[i].ne_lte = ne;
    (*model).atm.ne += ne;
  }
    (*model).atm.ne_lte = (*model).atm.ne;
    /************ SYSTEM **********
  printf("ne:%le\n",(*model).atm.ne);
    *************/
  return (*model).atm.ne;
}

double computeAbsoluteDensities(Model *model){
  //Model model= *modelref;
  int i,j,l;

  totalParticles =0.0;
  //CALCULANDO DENSIDADES ABSOLUTAS 
  /********** SYSTEM ********
  printf("Loading absolute density contributions:\n");
  **************/
  for(i=0; i < (*model).natom; i++){
    (*model).atom[i].n = (*model).atm.H*(*model).atom[i].A;
    totalParticles += (*model).atom[i].n;
    /***************SYSTEM**********
    printf("%s:%le\n",(*model).atom[i].name,(*model).atom[i].n);
    *****************************/
    for (j=0;j < (*model).atom[i].ns; j++){
      (*model).atom[i].specie[j].nabs = (*model).atom[i].specie[j].nrel*(*model).atm.H*(*model).atom[i].A;
      (*model).atom[i].specie[j].n_lte = (*model).atom[i].specie[j].nabs;
      /************ SYSTEM *************
      printf("%s:%le\n",(*model).atom[i].specie[j].name,(*model).atom[i].specie[j].nabs);
      **********************/

      //Calculando las densidades absolutas de n_l
      for (l=0; l < (*model).atom[i].specie[j].energylevels;l++){
	//TAKE CARE!
	// fixed in hmodel.c line 226
	(*model).atom[i].specie[j].el[l].n = (*model).atom[i].specie[j].el[l].nrel * (*model).atom[i].specie[j].nabs;
	/*********** SYSTEM ********
	printf("%le\t",(*model).atom[i].specie[j].el[l].n);
	******************/
      }
      (*model).atom[i].specie[j].nkabs = (*model).atom[i].specie[j].nk * (*model).atom[i].specie[j].nabs;
      /*********** SYSTEM ********
      printf("%le\n",(*model).atom[i].specie[j].nkabs);
      ********************/
    }
  }
  //******printf("%le\t",totalParticles);
}


//3
double Phi(double T){
  double phi_h = 6.626e-27; //ergs s-1
  double phi_m = 9.109e-28; //g
  double phi_k = 1.38e-16;//erg k-1
  double phi_xh = 13.6; //eV
  double phi_keV = 8.617e-5; //eV K-1
  return pow( (phi_h*phi_h)/(2.0*3.1416*phi_m*phi_k*T), 1.5) * exp(phi_xh/(phi_keV*T));

  //eV
//  pow(10.0,(68606.0/T)-15.38286092-(1.5*log10(T)));
   //pow(10.0,-141.15-1.5*log10(T)+(6859.13/T));
    //return pow((5.555962964e-11/T),1.5)*exp(1.579373098e5/T);
//pow(((4.390479437e-53)/(2.0*PI*me*k*T)),1.5)*exp((2.180562257e-11)/((1.3806505e-16)*T));
}


//CALCULA EN EQUILIBRIO TERMODINAMICO LA ECUACION DE SAHA
//Y REGRESA EL VALOR Z NECESARIO EN EL NLTE

double fZ(Model *model){
  //Necesito encontrar la contribucion real de electrones.
  //Necesito la contribucion de ionizacion y la densidad del nivel de continuo
  //para cada especie
  double **A;
  double *b;
  double *x;
  int N = (*model).natom;
  int i,j,m,n;
  int ai,aj;
  double sum,fsaha,u1,u0,ne,T,Xi;

  T = (*model).atm.T; //La temperatura del sistema
  ne =  (*model).atm.ne_lte;  //La densidad electronica inicial
  //printf("neA=%le\n",ne);

  //Inicio las matrices 
  n=0;
  for(i=0; i< N;i++){
    n+= (*model).atom[i].ns;
  }
  
  A = malloc(sizeof(double *)*n);
  for (i=0; i<n;i++){
    A[i] = malloc(sizeof(double)*n);
  }
  b = malloc(sizeof(double)*n);  
  //printf("n=%i\n",n);

  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      A[i][j] = 0.0;
    }
    b[i] =0.0;
  }

  ai=0;
  aj=0;
  //  printf("n=%i\n",n);
  for (i=0;i <N;i++){
    for (j=0; j < (*model).atom[i].ns;j++){
     A[ai][aj+j] = 1.0;
    }
    b[ai] = (*model).atom[i].n;
    ai++; //un rengloncito abajo

    
    for (j=0; j < (*model).atom[i].ns -1 ;j++){
      //calculando la funcion de particion
      //      printf("gradodeionizacionposible=%i\n",(*model).atom[i].g);
      if ((*model).atom[i].specie[j+1].ionizationstage == (*model).atom[i].g){
	u1=1.0;
      }else{
	//HII
	u1=(*model).atom[i].specie[j+1].basics[0].wl; //basics[m]
      }
      //      printf("u1%s=%le\n",(*model).atom[i].specie[j+1].name,u1);
      u0=(*model).atom[i].specie[j].basics[0].wl; //HI  //basics[m]
      //printf("u0%s=%le\n",(*model).atom[i].specie[j].name,u0);
      Xi=(*model).atom[i].specie[j].basics[0].vkl*1e9*h ; //Ej: HI en su estado base
      fsaha = log10(u1/u0)+15.6826+1.5*log10(T)-Xi*(5039.95/T)-log10(ne);


      //fsaha = -0.1761-log10(ne*k*T)+log10(u1/u0) + 2.5*log10(T)-Xi*(5040.0/T);

      //printf("XI=%le\n",Xi);
      //printf("Pe=%le[dynas] %le[eV]\n",ne*k*T*1.602e-12,ne*k*T);
      //printf("fsaha=%le\n",fsaha);
      A[ai][aj+j] = pow(10.0,fsaha);
      //printf("SAHA_A[%i][%i]=%le\n",ai,aj+j,A[ai][aj+j]);
      A[ai][aj+j+1] = -1.0;
      ai++;
    }

    aj+=j+1;
  }
  x = LU(n,A,b);
  m=0;
  (*model).atm.ne_lte=0.0;

  for (i=0; i<n;i++){
    sum = 0.0;
    for (j=0; j < (*model).atom[i].ns;j++){
      if (x[m] < 0.0){
	x[m] = 0.0;
      }
      (*model).atom[i].specie[j].n_lte = x[m];
      (*model).atom[i].specie[j].ne_lte = x[m]*((double)(*model).atom[i].specie[j].ionizationstage);
      sum += x[m]*((double)(*model).atom[i].specie[j].ionizationstage);
      //printf("X[%i]=%le\n",m,x[m]);  
      m++;
    }
    (*model).atom[i].ne_lte = sum;
    (*model).atm.ne_lte += sum;
  }
for (i=0; i<n;i++){
  free(A[i]);
  }
   free(A);
   free(b);
   free(x);
     return (*model).atm.ne_lte-(*model).atom[0].ne_lte;
   //return (*model).atm.ne_lte;
}

void ComputeNLTE(Model *model){
  double cero =0.0;
  // DUMMY PARA CALCULAR LOS b1 desde b1.dat
  char inter3D_str[512];

  sprintf(inter3D_str,"data/atmosphere/hydrostatic/%s/",(*model).name);

  //NLTE
  inter3D b_input = new_inter3D(inter3D_str,"b1.dat");
  B1 b = giveMeb1(b_input, (*model).atm.T, (*model).atm.H);

  if (b.b1 > cero){
    (*model).atom[0].specie[0].el[0].b = b.b1;
  }else{
    (*model).atom[0].specie[0].el[0].b = 1.0;
  }

  inter3D bhm_input = new_inter3D(inter3D_str,"bhm.dat");
  B1 bhm = giveMeb1(bhm_input, (*model).atm.T, (*model).atm.H);
  if (bhm.b1 > cero){
    (*model).atm.bhm = bhm.b1;
  }else{
    (*model).atm.bhm = 1.0;
  }
  

  //IMPORTANT, for LTE computations:
  // (*model).atom[0].specie[0].el[0].b = 1.0;
  //(*model).atm.bhm = 1.0;

}





void NLTE(Model *model, double error, int hydro, double pz1, double fz1,double dx,int chromospheric_network,int cell, int magnetic){
  // Model model = *modelref;
  //  double GMH = 1.11691659e-31;  
  int l,i,j;
  double nhm;
  double nl,n1,ne,ane,nu,ne_lte,ane_lte;
  double Z,Zt,ZtA;
  double d;
  double nH;
  double b1,bhm;
  double T,dinv;
  double sumA,wl,w1,sumA2;
  FILE *file,*fileatm;
  char comando[500];
  char path_chromosphere[500];
  int N;
  double cuatro = 4.0;
  double left = 1.0 - error;
  double right = 1.0 + error;
  double hleft = 1.0 - error;
  double hright = 1.0 + error;
  double sol;
  double R,P;
  double hsol;
  double vt;
  double fz;
  double Y;
  int N_MAX = 1000;
  int N_I=0;
  int H_MAX = 1000;
  int H_I=0;
  double integral;
  double Zval;
  double B_bx = 0;
  double B_by = 0;
  double B_bz = 0;
  double B_xi = 0;;

  //cuidado con las densidades energeticas, ya que solo se pueden recalcular una vez
  //fixed in this code and hmodel.c line 226 and energylevels.h

  N_MAX=100;
  if (hydro){
	  nH = (*model).atm.H;
	  printf("i\tT\t\tH\t\tP\t\tn_e\n");
	  printf("********************************************************************\n");


	  do{ //hydrogen convergence
		  //printf("hydro computing\n");
		  //    left = 1.0 - error;
		  // right = 1.0 + error;
		  (*model).atm.H = (nH+(*model).atm.H)*0.5 ;
		  computeAbsoluteDensities(&(*model));
		  ane=ne = ComputeAbsoluteNe(&(*model));
		  //  (*model).atm.H = nH;
		  T = (*model).atm.T;
		  ComputeNLTE(&(*model)); //Calculamos bi, ni para cada especie!
		  b1= (*model).atom[0].specie[0].el[0].b;
		  bhm = (*model).atm.bhm;
		  n1 = (*model).atom[0].specie[0].el[0].n; //parece que no se usa!
		  N = (*model).atom[0].specie[0].energylevels;
		  //  b1 =1.0;
		  //  bhm=1.0;
		  d = b1*Phi(T);
		  dinv=pow(d,-1.0)*0.5;
		  N_I=0;

		  if (H_I==0)
			  printf("%i\t%le\t%le\t%le\t%le\n",0, T, nH, (*model).atm.P,ne);

		  H_I++;
		  do{ //electronic density convergence
			  N_I++;
			  (*model).atm.ne = (ane+ne)*0.5;   //Sacamos el promedio entre el anterior y el actual (Punto medio)
			  ne = (*model).atm.ne; //Guardamos el valor para compararlo
			  (*model).atm.ne_lte = ne; //Le pasamos la nueva densidad electronica para que se calcule Z
			  ane= ne;
			  //ne = Z = fZ(&(*model));
			  // Z = 0.0;
			  Z = fZ(&(*model));
			  ne = dinv*sqrt(pow((1.0- (Z*d)) ,2.0)+ (cuatro*d*(nH+Z)))-(1.0-(Z*d));
			  //printf("ne=%.50e\n",cuatro*d*(nH+Z) );
			  //n_H-
			  nhm =  bhm*(1.0345e-16)* ne * (*model).atom[0].specie[0].n_lte * pow(T,-1.5)*exp(8762.0/T);
			  if (nhm < 0.0)
				  nhm=0.0;
			  //OJO
			  ne -= nhm;//electrones - ocupados en H-
			  //printf("NE=%le\n",ne);
			  (*model).atom[0].specie[0].n_lte -= nhm; //Hidrogeno Neutro - ocupdos en H-.
			  sol =  fabs(ne/ane);
			  //    printf("NI=%i\tne=%.25e %.25e ERROR=%e\n",N_I,ne,ane,fabs(ne-ane));
			  //    printf("NE=%.11e\n",(*model).atm.ne_lte);
			  //  }while(absd(ne_lte-ane_lte) > 0.0);
			  // }while(  fabs(ne-ane) > 1.0 );
		  }while( (left > sol || sol > right)  && (N_I <= N_MAX) );

		  //*******  printf("%le\t%le\t",(*model).atm.z,T);
		  // printf("ions\n");
		  sumA=0.0;
		  for (i=0; i< (*model).natom; i++){
			  for (j=0; j< (*model).atom[i].ns;j++){
				  //**********	printf("%le\t",(*model).atom[i].specie[j].n_lte);
				  sumA+=(*model).atom[i].specie[j].n_lte;
			  }
		  }

		  sumA+=nhm;
		  //*******  printf("%le\t%le\t%le\t",nhm,(*model).atm.z,b1,(*model).atm.z,bhm);
		  //ne = (sqrt(pow(1.0 - Z*d,2.0)+4.0*d*(nH+Z))-(1.0-Z*d) )/(2.0*d);
		  //******* printf("%le\t%le\t%le\t",(*model).atm.ne_lte,sumA,Z);
		  //******* printf("%le\t%le\t%le",sol,fabs(ne-ane), totalParticles-sumA );
		  ne_lte = (*model).atm.ne_lte;
		  //******* printf("\n");
		  R = 1.0/(1.0 +d*ne);
		  Y= (*model).atom[1].A;
		  vt = (*model).atm.V*1e5;
		  Zval=Z/(*model).atm.H;
		  B_bx = (*model).atm.Bx;
		  B_by = (*model).atm.By;
		  B_bz = (*model).atm.Bz;
		  B_xi = (*model).atm.xi;
		  fz = f_hydro(Y,R,Zval,T,vt,B_bx,B_by,B_bz,B_xi);
		  //  integral = 0.5*(fz+fz1)*dx*1e5;
		  //   P = pz1*exp( integral  );
		  integral = (fz+fz1)*dx*1e5*0.5;
		  P = exp( log(pz1)- integral);
		  // printf("R=%le\tY=%le\tZ=%le\tvt=%le\n",R,Y,Zval,vt);
		  //    printf("integral=%le\tP=%le\n",integral,P);
		  (*model).atm.P = P;
		  //    nH = (fz*P)/(1.11691659e-31*(1.0+4.0*Y));
		  nH = exp(log(fz)+log(P)-log(GMH)-log(1.0+4.0*Y));
		  //printf("%le\t%le\t%le\t%le\t%le\n",pz1,P,fz1,fz,nH);
		  printf("%i\t%le\t%le\t%le\t%le\n",H_I, T, nH, P,ne);
		  (*model).atm.fz = fz; //for the next computation
		  hsol =  fabs(nH/(*model).atm.H);
		  //******* printf("H=%le \t P=%le\n",nH,P);
	  }while( (hleft > hsol || hsol > hright) && (H_I <= H_MAX) );

	  (*model).atm.H = nH;//+(*model).atm.H)*0.5 ;
	  //printf("H = %le\tP = %le\n",(*model).atm.H,(*model).atm.P);
  }else{
	  //No hydrostatic computations!
	  //   left = 1.0 - error;
	  //right = 1.0 + error;
	  computeAbsoluteDensities(&(*model));
	  ane=ne = ComputeAbsoluteNe(&(*model));
	  nH = (*model).atm.H;
	  T = (*model).atm.T;
	  ComputeNLTE(&(*model)); //Calculamos bi, ni para cada especie!
	  b1= (*model).atom[0].specie[0].el[0].b;
	  bhm = (*model).atm.bhm;
	  n1 = (*model).atom[0].specie[0].el[0].n; //parece que no se usa!
	  N = (*model).atom[0].specie[0].energylevels;
	  //  b1 =1.0;
	  //  bhm=1.0;
	  d = b1*Phi(T);
	  dinv=pow(d,-1.0)*0.5;
	  do{
		  (*model).atm.ne = (ane+ne)*0.5;   //Sacamos el promedio entre el anterior y el actual (Punto medio)
		  ne = (*model).atm.ne; //Guardamos el valor para compararlo
		  (*model).atm.ne_lte = ne; //Le pasamos la nueva densidad electronica para que se calcule Z
		  ane= ne;
		  //ne = Z = fZ(&(*model));
		  // Z = 0.0;
		  Z = fZ(&(*model));
		  ne = dinv*sqrt(pow((1.0- (Z*d)) ,2.0)+ (cuatro*d*(nH+Z)))-(1.0-(Z*d));
		  //printf("ne=%.50e\n",cuatro*d*(nH+Z) );
		  //n_H-
		  nhm =  bhm*(1.0345e-16)* ne * (*model).atom[0].specie[0].n_lte * pow(T,-1.5)*exp(8762.0/T);
		  if (nhm < 0.0)
			  nhm=0.0;
		  //OJO
		  ne -= nhm;//electrones - ocupados en H-
		  //printf("NE=%le\n",ne);
		  (*model).atom[0].specie[0].n_lte -= nhm; //Hidrogeno Neutro - ocupdos en H-.
		  sol =  fabs(ne/ane);
		  //printf("%.25e %.25e ERROR=%e\n",ne,ane,fabs(ne-ane));
		  //    printf("NE=%.11e\n",(*model).atm.ne_lte);
		  //  }while(absd(ne_lte-ane_lte) > 0.0);
		  // }while(  fabs(ne-ane) > 1.0 );
	  }while(left > sol || sol > right );

	  //******* printf("%le\t%le\t",(*model).atm.z,T);
	  // printf("ions\n");
	  sumA=0.0;
	  for (i=0; i< (*model).natom; i++){
		  for (j=0; j< (*model).atom[i].ns;j++){
			  //******* printf("%le\t",(*model).atom[i].specie[j].n_lte);
			  sumA+=(*model).atom[i].specie[j].n_lte;
		  }
	  }
	  sumA+=nhm;
	  //*******  printf("%le\t%le\t%le\t",nhm,(*model).atm.z,b1,(*model).atm.z,bhm);
	  //ne = (sqrt(pow(1.0 - Z*d,2.0)+4.0*d*(nH+Z))-(1.0-Z*d) )/(2.0*d);
	  //*******  printf("%le\t%le\t%le\t",(*model).atm.ne_lte,sumA,Z);
	  //******* printf("%le\t%le\t%le",sol,fabs(ne-ane), totalParticles-sumA );
	  ne_lte = (*model).atm.ne_lte;
	  R = 1.0/(1.0 +d*ne_lte);
	  Y= (*model).atom[1].A;
	  vt = (*model).atm.V*1e5;
	  Zval = Z/(*model).atm.H;
	  B_bx = (*model).atm.Bx;
	  B_by = (*model).atm.By;
	  B_bz = (*model).atm.Bz;
	  B_xi = (*model).atm.xi;
	  fz = f_hydro(Y,R,Zval ,T,vt,B_bx,B_by,B_bx,B_xi);
	  (*model).atm.fz = fz; //for the next computation!
    //******* printf("\n");
    //printf("H = %le\tP = %le\n",(*model).atm.H,(*model).atm.P);
  } //finish computations

  //printf("Z=%le\tT=%le\tH=%le\tP=%le\n",(*model).atm.z,(*model).atm.T,(*model).atm.H,(*model).atm.P);

  //WRITE DATA

  //Building paths!

  if (chromospheric_network){
    if (cell){
        sprintf(path_chromosphere,"data/atmosphere/chromosphere/chromosnet/cell");
    }else{
      sprintf(path_chromosphere,"data/atmosphere/chromosphere/chromosnet/net");
    }
  }else{
    sprintf(path_chromosphere,"data/atmosphere/chromosphere/average");
  }


  
  sprintf(comando,"%s/H.dat",path_chromosphere);
  file = fopen(comando,"a");
  fprintf(file, "%le\t%le\n",(*model).atm.z,(*model).atm.H);
  fflush(file);
  fclose(file);

  sprintf(comando,"%s/P.dat",path_chromosphere);
  file = fopen(comando,"a");
  fprintf(file, "%le\t%le\n",(*model).atm.z,(*model).atm.P);
  fflush(file);
  fclose(file);

  sprintf(comando,"%s/ne.dat",path_chromosphere);
  file = fopen(comando,"a");
  fprintf(file, "%le\t%le\n",(*model).atm.z,(*model).atm.ne_lte);
  fflush(file);
  fclose(file);

  sprintf(comando,"%s/b1.dat",path_chromosphere);
  file = fopen(comando,"a");
  fprintf(file, "%le\t%le\n",(*model).atm.z,b1);
  fflush(file);
  fclose(file);

  sprintf(comando,"%s/bhm.dat",path_chromosphere);
  file = fopen(comando,"a");
  fprintf(file, "%le\t%le\n",(*model).atm.z,bhm);
  fflush(file);
  fclose(file);

  sprintf(comando,"%s/T.dat",path_chromosphere);
  file = fopen(comando,"a");
  fprintf(file, "%le\t%le\n",(*model).atm.z,(*model).atm.T);
  fflush(file);
  fclose(file);
  sprintf(comando,"%s/species.dat",path_chromosphere);
  fileatm = fopen(comando,"w");
  for (i=0; i< (*model).natom; i++){
	  for (j=0; j< (*model).atom[i].ns;j++){
		  fprintf(fileatm,"%s\t%i\n",(*model).atom[i].specie[j].name,(*model).atom[i].specie[j].ionizationstage );
		  sprintf(comando,"%s/%s.dat",path_chromosphere,(*model).atom[i].specie[j].name);
		  file = fopen(comando,"a");
		  fprintf(file, "%le\t%le\n",(*model).atm.z,(*model).atom[i].specie[j].n_lte);
		  fflush(file);
		  fclose(file);
      }
  }
  sprintf(comando,"%s/H-.dat",path_chromosphere);
  file = fopen(comando,"a");
  fprintf(file, "%le\t%le\n",(*model).atm.z,nhm);
  fflush(file);
  fclose(file);
  fflush(fileatm);
  fclose(fileatm);
  /*
  printf("\n%le\t",T);
    sumA=0.0;
    for (i=0; i< (*model).natom; i++){
      for (j=0; j< (*model).atom[i].ns;j++){
	printf("%le\t",(*model).atom[i].specie[j].n_lte);
	sumA+=(*model).atom[i].specie[j].nabs;
      }
    }
      printf("%le\t%le\n",ne,sumA);
  */
}
