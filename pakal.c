/* Pakal 20171009
 * Copyright (C) 2006 by Free Software Foundation, Inc.
 * Autores: De la Luz Rodriguez Victor Hugo
 * Francisco Tapia
 * Carlos Miranda
 *                                 <vdelaluz@geofisica.unam.mx>
 *
 *
 *                          COPYING
 *
 * This  program  is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License  as
 * published by the Free Software Foundation; either version 2, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty  of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See  the
 * GNU General Public License for more details.
 *
 * You  should  have  received  a  copy  of the GNU General Public
 * License along with this program; see the file  GPL.txt. if not,
 * write  to the Free Software Foundation, Inc., 59 Temple Place -
 * Suite 330, Boston, MA 02111-1307, USA.
 *
 * ********************** Revision History ***********************
 *
 * $Log$
 * 18/09/2017
 * Integrate Vector to nu frequencies.
 * 08/01/2016 -cn activated: compute the Chromospheric Network
 * Exaple:
 * mpiexec -n 1 ./pakal -model C07 -compute-ion-profile -cn
 * 06/01/2015 Major improvements in the file outputs.
 * Example 1:
 * mpiexec -n 8 ./pakal -model C07 -nu 43 -min 1e-40 -r 4 -dz 1
 * -ia-step 1 -Rt 7.3e5 -z-begin 0 -z-end 764084 -v 1
 * -pedanting 100
 * Example 2 compute ion profile:
 * mpiexec -n 1 ./pakal -model C07 -compute-ion-profile
 * Example 3 compute hydrostatic model taking into acount the
 * initial model C07:
 * mpiexec -n 1 ./pakal -model C07 -hydro
 * 05/04/10 Vamos a cambiar el flujo de calculo para poder usarlo
 * con los clusters. Primero, leeremos las configuraciones desde
 * un archivo seris.conf, despues mandaremos llamar a jaguar para
 * realizar los calculos necesarios y finalmente lanzaremos a
 * pakal para que calcule la ecuacion de transferencia.
 * 21/V/09 free and malloc revised
 * 17-IV-07 43GHz Customize
 * 16-IV-07 Recuperado!
 * 10-X-06 Inicio
 * ********************** Preview ********************************
 *
 * High Atmosferic Solar Model
 * ************************ Notes ********************************

 * The sun's high is 0 at R_sun.

 * ***************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "chromosnet.h"
#include "functions.h"
#include "geometry.h"
#include "help.h"
#include "lineal.h"
#include "modelonumerico.h"
#include "MPI_pakal.h"
#include "particion.h"
#include "physics.h"
#include "Timming.h"
#include "tulum.h"

#define VERSION "PAKAL Rev 20170918"

/****************************
 *****************************
 ** Ordenamiento Merge Sort **
 *****************************
 *****************************/

void merging(int low, int mid, int high,double *a,double *b, double *c,double *d) {
	int l1, l2, i;

	for(l1 = low, l2 = mid + 1, i = low; l1 <= mid && l2 <= high; i++) {
		if(a[l1] <= a[l2]) {
			b[i] = a[l1++];
			d[i] = c[l1-1];
		}else {
			b[i] = a[l2++];
			d[i] = c[l2-1];
		}
	}

	while(l1 <= mid) {
		b[i++] = a[l1++];
		d[i-1] = c[l1-1];
	}
	while(l2 <= high) {
		b[i++] = a[l2++];
		d[i-1] = c[l2-1];
	}

	for(i = low; i <= high; i++) {
		a[i] = b[i];
		c[i] = d[i];
	}
}

void sort(int low, int high,double *a, double *b, double *c,double *d) {
	int mid;

	if(low < high) {
		mid = (low + high) / 2;
		sort(low, mid,a,b,c,d);
		sort(mid+1, high,a,b,c,d);
		merging(low, mid, high,a,b,c,d);
	}else {
		return;
	}
}


int main(int argc, char **argv){
	int intervalos;
	int n_proc;
	int procesador;
	int miproc, numproc;
	MPI_Status status;
	double t3, t4;
	int pori;
	int porj=10;
	int flag=0;

	int reflejo=0;
	int line=0;
	int line0=0;
	int line1=0;

	//Variables para Ordenamiento Rev 20170918
	int iter;
	double *Tb; //Temperatura de brillo
	double *fc; //Frecuencia
	double *Tba; //Temperatura de brillo auxiliar
	double *fca; //Frecuancia-Auixiliar
	double temp1, temp2;
	int num, max;
	int cont=0;
	FILE * fp;
	char * line_array = NULL;
	size_t len_array=0;
	ssize_t read_array;
	int n_nus=0;
	int i_array=0;
	char *frequencies_file;
	//

	// Variables to restricted Rev 20170918
	double restricted;
	int tau_is_negative = 0;
	int is_restricted = 0;
	char error_mpi[256];
	time_t t_random;
	double random_value;
	FILE *file;

	//Rev 20170926
	double Trmin;


	int isline=0;
	int verbose = 0;
	int puntoenx = 0;
	int puntoeny = 0;
	int solounpunto = 0;
	int xini,xfin,yini,yfin;
	int x,y,m;
	int parai=0;
	double z,z_0,dzDetail;
	double alpha,beta,dz,dzBig;
	int n = 5; //ok
	double min = 1e-25;
	double minx,miny;
	double detail = 0.5;
	double F = -2.0*Rsun; //OK
	double H = 2.0*Rsun; //OK
	double dl = 20.0; //km
	double deltaTau = dl;
	double nu = 43e9; //OK Hz
	double wl = 0.0; //OK
	double r,theta,phi,I,rt,n0,tau,mt;
	double Tb_K; // Rev 20170918

	double Rt = 7.5e5; //km
	double localintensity,localintensityT,localemision;
	double t_temp;
	int cuadrante=0;
	char *linea;
	char *limpia;
	char outputfile[100];
	char outputpath[1024];
	char temperaturefile[100];
	char densityfile[100];
	char hydrogenfile[100];
	char presionfile[100];
	char sahafile[100];
	char comando[200];
	char modelCell[200];
	char modelNet[200];
	char atm_model[200];
	int chromosnet=0; //0 No usar, 1 Usar
	int pedanting=0;
	int nped=0;
	int ismin=1;
	int compute_ion_profile = 0;
	lineal temperature;
	lineal density;
	lineal hydrogen;
	lineal presion;
	FILE *fout;
	int i,j;
	int hydro=0;
	//  double hydro_step=1.0;
	int surface = 0; //no estoy en la superficie
	MPI_pakal_message mp;
	particion pa,ra;
	int nmodel=2;

	//campo magnetico
	double B_amplitude = 0, B_intensity = 0, B_alpha = 0;
	int magnetic = 0;

	/****NUEVAS VARIABLES****/
	int nStep = 100; //Big = nStep*detail
	double utime0, stime0, wtime0,
	       utime1, stime1, wtime1,
	       utime2, stime2, wtime2;
	double nu1,nu2;
	int is_spectrum = 0;
	// int n_nus =1;   //Rev 20170918 defined in top

	MPI_Init (&argc, &argv); /* Inicializar MPI */
	t3 = MPI_Wtime(); /* tiempo de inicio */

	MPI_Comm_rank(MPI_COMM_WORLD,&miproc); /* Determinar el rango del proceso invocado*/
	MPI_Comm_size(MPI_COMM_WORLD,&numproc); /* Determinar el numero de procesos */
	//  printf("\nProceso: %i Total: %i\n",miproc+1,numproc);

	MPI_Barrier (MPI_COMM_WORLD);

	/***************** COMIENZA MAESTRO !!! ******************/

	if (miproc == 0) { /* Si es el maestro lee todo lo que tengas que leer*/
		// printf("\nHola Soy el Proceso maestro (0)\n");
		uswtime(&utime0, &stime0, &wtime0); //tomando el tiempo
		sprintf(outputfile,"sun.dat");


		sprintf(temperaturefile,"data/atmosphere/chromosphere/average/T.dat");
		sprintf(densityfile,"data/average/Ne.dat");
		sprintf(hydrogenfile,"data/average/H.dat");
		sprintf(presionfile,"data/average/P.dat");
		sprintf(atm_model,"NO_COMPUTE");
		sprintf(outputpath,"output");

		/*************************Get Parameters*************************/
		/*
		   printf("Pakal ");
		   printf(VERSION);
		   printf(" GNU/GPL license\n");
		   printf(" By Victor H De la Luz (itztli@gmail.com) [UNAM]\n");
		 */
		if (argc == 1) {
			FILE *fpi;
			int icf =0;
			char tempicf[50];
			char lineaicf[501];
			//      int erroricf =0;
			printf("Modo no interactivo!\n");
			//Leer archivo pakal.input
			fpi = fopen("pakal.input","r");
			if (fpi != NULL) {
				while (fgets(lineaicf, 500, fpi )!= NULL) {
					// printf("%s\n",lineaicf);
					limpia =Clean(lineaicf);
					if ( strlen(limpia) > 0) {

						//printf("%s\n",limpia);

						switch(icf) {
						case 0:
							if (sscanf(limpia,"%lf",&nu) != 1) {
								printf("Error in File pakal.input: Frequency");
								return 0;
							}
							//printf("0\n");
							break;
						case 1:
							if (sscanf(limpia,"%lf",&H) != 1) {
								printf("Error in File pakal.input: Height");
								return 0;
							}
							break;
						case 2:
							if (sscanf(limpia,"%lf",&F) != 1) {
								printf("Error in File pakal.input: Floor");
								return 0;
							}
							break;
						case 3:
							if (sscanf(limpia,"%lf",&Rt) != 1) {
								printf("Error in File pakal.input: Total Radii");
								return 0;
							}
							break;
						case 4:
							if (sscanf(limpia,"%i",&n) != 1) {
								printf("Error in File pakal.input: Resolution");
								return 0;
							}
							break;
						case 5:
							if (sscanf(limpia,"%lf",&detail) != 1) {
								printf("Error in File pakal.input: nDetail");
								return 0;
							}
							break;
						case 6:
							if (sscanf(limpia,"%i",&nStep) != 1) {
								printf("Error in File pakal.input: nStep");
								return 0;
							}
							break;
						case 7:
							if (sprintf(outputfile,"%s",limpia) < 0) {
								printf("Error in File pakal.input: outputdir");
								return 0;
							}
							break;

						case 8:
							if (sprintf(tempicf,"%s",limpia) > 0) {

								strncpy(tempicf,limpia,strlen(limpia)-1);
								tempicf[strlen(limpia)-1] = '\0';


								if (strcmp(tempicf,"single_point")==0) {
									solounpunto=1;
									//printf("Solo un Punto\n");
								}else if (strcmp(tempicf,"line")==0) {
									isline=1;
								}else if (strcmp(tempicf,"image")==0) {
									//nothing
								}else{
									printf("Error in File pakal.input: single_point, line or image");
									return 0;
								}
							}else{
								printf("Error in File pakal.input: single_point, line or image");
								return 0;
							}
							break;

						case 9:
							if (solounpunto) {
								if (sscanf(limpia,"%i %i",&puntoenx,&puntoeny) != 2) {
									printf("Error in File pakal.input: x y");
									return 0;
								}
							}
							if (isline) {
								if (sscanf(limpia,"%i",&line) != 1) {
									printf("Error in File pakal.input: line point");
									return 0;
								}
							}
							break;

						case 10:
							if (sscanf(limpia,"%lf",&min) != 1) {
								printf("Error in File pakal.input: min parameter");
								return 0;
							}
							if (min > 0) {
								ismin=0;
							}
							break;

						case 11:
							if (sscanf(limpia,"%i",&porj) != 1) {
								printf("Error in File pakal.input: verbose (nsave)");
								return 0;
							}
							if (porj > 0) {
								verbose=1;
							}
							break;
						case 12:
							if (sscanf(limpia,"%i",&nped) != 1) {
								printf("Error in File pakal.input: pedanting (npedanting)");
								return 0;
							}
							if (nped > 0) {
								pedanting=1;
							}
							break;
						case 13:
							if (sscanf(limpia,"%i",&nmodel) != 1) {
								printf("Error in File pakal.input: model opacity");
								return 0;
							}
							break;
						case 14:
							/*if (sprintf(atm_model,"%s",limpia) < 0){
							   printf("Error in File pakal.input: chromospheric modelr");
							   return 0;
							   }*/
							strncpy(atm_model,limpia,strlen(limpia)-1);
							atm_model[strlen(limpia)-1] = '\0';
							break;
						case 15:
							if (sscanf(limpia,"%i",&chromosnet) != 1) {
								printf("Error in File pakal.input: chromospheric network");
								return 0;
							}
							break;
						case 16:
							if (sscanf(limpia,"%i",&is_spectrum) != 1) {
								printf("Error in File pakal.input: compute spectrum");
								return 0;
							}
							break;

						case 17:
							if (sscanf(limpia,"%lf %lf %i",&nu1,&nu2,&n_nus) != 3) {
								printf("Error in File pakal.input: spectrum parameters");
								nu1 = nu1*1e9;
								nu2 = nu2*1e9;
								return 0;
							}
							break;
						}//switch
						icf++;
					}//if
				}//while


			}else{//file
				printf("Error leyendo archivo de entrada pakal.input\n");
				return 0;
			}
			fclose(fpi);
			//printf("OK\n");
		}else{//modo iterativo

			for (i=1; i<argc; i++) {
				sprintf(comando,"%s",argv[i]);

				if (strcmp(comando,"-opacities") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%i\n",&nmodel) > 0) {
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				/* Rev 20170918

				   if (strcmp(comando,"-spectrum") == 0){
				   sprintf(comando,"%s",argv[++i]);
				   if (sscanf(comando,"%le\n",&nu1) > 0){
				    //printf(".");
				    sprintf(comando,"%s",argv[++i]);
				    if (sscanf(comando,"%le\n",&nu2) > 0){
				      //printf(".");
				      sprintf(comando,"%s",argv[++i]);
				      if (sscanf(comando,"%i\n",&n_nus) > 0){
				        //printf(".");
				        is_spectrum = 1;
				        nu1 = nu1*1e9;
				        nu2 = nu2*1e9;
				      }else{
				        imprimeInstrucciones();
				        return 0;
				      }
				    }else{
				      imprimeInstrucciones();
				      return 0;
				    }
				   }else{
				    imprimeInstrucciones();
				    return 0;
				   }
				   }

				 */

				if (strcmp(comando,"-spectrum") == 0) {
					//sprintf(comando,"%s",argv[++i]);
					frequencies_file = malloc(sizeof(char)*500);
					if (sprintf(frequencies_file,"%s",argv[++i]) > 0) {
						//printf(".");
						// sprintf(mp.atm_model,"%s",atm_model);

						fp = fopen(frequencies_file, "r");
						if (fp == NULL)
							exit(EXIT_FAILURE);
						n_nus=0;
						while ((read_array = getline(&line_array, &len_array, fp)) != -1) {
							n_nus++;
						}
						fclose(fp);
						//mp.nu_vector_full= malloc(sizeof(float)*n_nus);
						//fp = fopen("frequencies.dat", "r");
						fp = fopen(frequencies_file, "r");
						if (fp == NULL)
							exit(EXIT_FAILURE);
						while ((read_array = getline(&line_array, &len_array, fp)) != -1) {
							//		nu_vector_full[i_array]=atoi(line_array) ;
							sscanf(line_array,"%lf",&mp.nu_vector_full[i_array]);
							//printf("Elemento[%d] = %f\n", i_array, nu_vector_full[i_array] );
							i_array++;
						}
						fclose(fp);
						//Useless
						nu1 = mp.nu_vector_full[0]*1e9;
						nu2 = mp.nu_vector_full[n_nus]*1e9;

						//printf("Reading %s\n",frequencies_file);

						free(frequencies_file);
						for(j=0; j<n_nus; j++) {
							//printf("%lf\n",mp.nu_vector_full[j]);
							mp.nu_vector_full[j]*=1e9;
						}

						is_spectrum=1;
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				if (strcmp(comando,"-compute-ion-profile") == 0) {
					compute_ion_profile = 1;
				}


				if (strcmp(comando,"-model") == 0) {
					if (sprintf(atm_model,"%s",argv[++i]) > 0) {
						//printf(".");
						sprintf(mp.atm_model,"%s",atm_model);
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				if (strcmp(comando,"-pedanting") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%i\n",&nped) > 0) {
						pedanting=1;
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}


				if (strcmp(comando,"-v") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%i\n",&porj) > 0) {
						verbose=1;
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				if (strcmp(comando,"-xy") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%i\n",&puntoenx) > 0) {
						sprintf(comando,"%s",argv[++i]);
						if (sscanf(comando,"%i\n",&puntoeny) > 0) {
							solounpunto=1;
						}else{
							imprimeInstrucciones();
							return 0;
						}
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}


				if (strcmp(comando,"-Rt") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%le\n",&Rt) > 0) {
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				if (strcmp(comando,"-min") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%le\n",&min) > 0) {
						printf("min=%le\n",min);
						ismin=0;
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				if (strcmp(comando,"-l") == 0) {
					sprintf(comando,"%s",argv[++i]);
					//if (sscanf(comando,"%i\n",&line) > 0){
					//  isline=1;
					//printf(".");
					//}

					if (sscanf(comando,"%i\n",&line) > 0) {
						//printf(".");
						sprintf(comando,"%s",argv[++i]);
						if (sscanf(comando,"%i\n",&line0) > 0) {
							//printf(".");
							sprintf(comando,"%s",argv[++i]);
							if (sscanf(comando,"%i\n",&line1) > 0) {
								isline=1;
							}else{
								imprimeInstrucciones();
								return 0;
							}
						}else{
							imprimeInstrucciones();
							return 0;
						}
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}
				if (strcmp(comando,"-help") == 0) {
					imprimeAyuda();
					return 0;
				}
				if (strcmp(comando,"-cn") == 0) {
					printf("Chromosnet activated\n");
					chromosnet = 1;
					//return 0;
				}

				if (strcmp(comando,"-wl") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%lf\n",&wl) > 0) {
						//printf(".");
						nu = C_light/wl;
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}
				if (strcmp(comando,"-nu") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%lf\n",&nu) > 0) {
						nu = nu*1e9;
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}
				if (strcmp(comando,"-z-end") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%lf\n",&H) > 0) {
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}
				if (strcmp(comando,"-z-begin") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%lf\n",&F) > 0) {
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				if (strcmp(comando,"-dz") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%le\n",&detail) > 0) {
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				if (strcmp(comando,"-ia-step") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%i\n",&nStep) > 0) {
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				if (strcmp(comando,"-r") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%i\n",&n) > 0) {
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				if (strcmp(comando,"-hydro") == 0) {
					//	sprintf(comando,"%s",argv[++i]);
					//	if (sscanf(comando,"%lf\n",&hydro_step) > 0)
					hydro = 1;
					//printf(".");
					//	else{
					//  imprimeInstrucciones();
					// return 0;
					//}
				}

				if (strcmp(comando,"-o") == 0) {
					if (sprintf(outputfile,"%s",argv[++i])) {
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				if (strcmp(comando,"-restricted") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%le\n",&restricted) > 0) {
						is_restricted=1;
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}


				if (strcmp(comando,"-path") == 0) {
					if (sprintf(outputpath,"%s",argv[++i])) {
						//printf(".");
					}else{
						imprimeInstrucciones();
						return 0;
					}
				}

				//componente de campo magnetico
				if (strcmp(comando,"-B") == 0) {
					sprintf(comando,"%s",argv[++i]);
					if (sscanf(comando,"%lf\n",&B_amplitude) > 0) {
						sprintf(comando,"%s",argv[++i]);
						if (sscanf(comando,"%lf\n",&B_intensity) > 0) {
							sprintf(comando,"%s",argv[++i]);
							if(sscanf(comando,"%lf\n",&B_alpha) > 0){
							magnetic=1;
							}
							else{
								imprimeInstrucciones();
								return 0;
							}
						}
						else{
							imprimeInstrucciones();
							return 0;
						}
					}
					else{
						imprimeInstrucciones();
						return 0;
					}
				}
			}

		}//modo iterativo

		/*
		   sprintf(perfilesOut,"profiles");    //z vs T,n_e,H,HI,HII,HeI,HeII,HeIII
		   sprintf(emisionOut,"emission"); //z vs kappa, tau, I_l, I_t
		 */

		dl=((double)nStep)*detail;

		printf("\n===== PakalMPI rev 20170918 =====\n");
		printf("By: vdelaluz@geofisica.unam.mx\n");
		printf("GNU/GPL License\n");
		printf("==============================\n");

		if ( (n % 2) == 0) n++;


		/******************************** (END)Get Parameters***********************/

		/******************************** INPUTS FILES ***************************/
		//fout = fopen(outputfile, "w+");

		/*
		   temperature = ilineal(temperaturefile);
		   density = ilineal(densityfile); //electronic density
		   hydrogen = ilineal(hydrogenfile);
		   HI = ilineal("data/HI.dat");
		   HII =  ilineal("data/HII.dat");
		   HeI = ilineal("data/HeI.dat");
		   HeII = ilineal("data/HeII.dat");
		   HeIII = ilineal("data/HeIII.dat");
		 */

		//Saha("saha.dat", 0, 50 , 1000.0, 1.0, abundances, electro);


		/******************************** (END) OUTPUTS FILES *********************/

		if (solounpunto==1) { //Vamos a integrar un punto
			printf("Single Point  (%i, %i)\n",puntoenx,puntoeny);
			xini = xfin = puntoenx;
			yini = yfin =puntoeny;
		}else if (isline==1) { // Vamos a integrar una linea
			if (chromosnet) {
				//need review!!!
				yini= -(n-1)/2;
				yfin = (n-1)/2;
				xini = line-5;
				xfin = line+5;
			}else{
				//yini= 0;
				//yfin = (n-1)/2;
				yini= line0;
				yfin = line1;
				xini = line;
				xfin = line;
			}
			printf("Image Line (%i,%i:%i)\n",line,yini,yfin);
		}else{ //vamos a integrar toda la imagen
			xini= -(n-1)/2;
			xfin = (n-1)/2;
			yini=-(n-1)/2;
			yfin=(n-1)/2;
		}
		mp.nmodel = nmodel;
		mp.pedanting = pedanting;
		mp.nped = nped;
		mp.chromosnet= chromosnet;
		mp.n =n;
		mp.intervalo=n_proc+1;
		mp.nu = nu;
		mp.Rt= Rt;
		mp.F= F;
		mp.H = H;
		mp.dl =dl;
		mp.detail=detail;
		mp.nStep = nStep;
		mp.min= min;
		mp.verbose = verbose;
		mp.porj = porj;
		mp.xini = xini;
		mp.xfin = xfin;
		mp.yini = yini;
		mp.yfin = yfin;
		mp.is_spectrum = is_spectrum;
		mp.restricted = restricted; // Rev 20170918
		sprintf(mp.outputpath,"%s", outputpath);

		/* COMIENZA EL CALCULO DE ESPECIES */
		//printf("ok1\n");
		//if (strcmp(atm_model, "NO_COMPUTE") == 0){	

		if (!(compute_ion_profile)) {
			if (hydro) {

				sprintf(comando,"rm data/atmosphere/chromosphere/average/*.dat");
				system(comando);
				printf("Computing hydrostatic atmosphere.\n");
				//sprintf(comando,"jaguar/jaguar %s 1",atm_model);
				
				if(magnetic){
				sprintf(comando,"jaguar/jaguar %s -hydro -B %le %le %le",atm_model, B_amplitude, B_intensity, B_alpha);}
				else{
				sprintf(comando,"jaguar/jaguar %s -hydro",atm_model);}

				system(comando);		

				printf("Ready\n");
				return 0;
			}else{
				printf("Using previous ion profiles.\n");
			}
		}else{
			//ADD for CN
			if (chromosnet) {
				printf("Computing ion profiles CN activated.\n");
				sprintf(comando,"rm data/atmosphere/chromosphere/chromosnet/cell/*.dat");
				system(comando);
				sprintf(modelCell,"%s-CELL",atm_model);
				printf("Computing %s \n", modelCell);
				sprintf(comando,"jaguar/jaguar %s -cn 1",modelCell); //1
				printf("%s\n",comando);
				system(comando);
				printf("Ready\n");
				printf("Computing ion profiles.\n");
				sprintf(comando,"rm data/atmosphere/chromosphere/chromosnet/net/*.dat");
				system(comando);
				sprintf(modelNet,"%s-NET",atm_model);
				printf("Computing %s \n", modelNet);
				sprintf(comando,"jaguar/jaguar %s -cn 0",modelNet); //0
				printf("%s\n",comando);
				system(comando);
				printf("Ready\n");
				return 0;
			}else{
				printf("Computing ion profiles.\n");
				sprintf(comando,"rm data/atmosphere/chromosphere/average/*.dat");
				system(comando);
				sprintf(comando,"jaguar/jaguar %s",atm_model);
				printf("%s\n",comando);
				system(comando);
				printf("Ready\n");
				return 0;
			}

		}



		//Compute the min value!...

		if (ismin) {
			linea = malloc(sizeof(char)*500);

			fout = fopen(temperaturefile, "r");
			if (fout==NULL) {
				printf("Error 29: File %s not found.\n",temperaturefile);
				exit(0);
			}

			mp.Tmin =Trmin = 1e10;
			while (fgets(linea, 500, fout ) != NULL) {
				limpia = Clean(linea);
				if (strlen(limpia)>0) {
					//printf("i=%i,%s\n",i,limpia);
					if (sscanf(limpia,"%le %le",&minx,&miny) == 2) {
						if (miny < Trmin) {
							Trmin = miny;
							mp.Tmin = miny;
						}
					}
				}
			}

			min = (1e-4) *Trmin* 2.0*K*pow(nu,2.0)/pow(C_light,2.0);
			fclose(fout);
			free(linea);

			//      printf("min=%le\tTmin=%le\n ",min,mp.Tmin);
		}else{
			printf("Min manual\n");
			mp.Tmin = (1e4)*pow(C_light,2.0)*min/(2.0*K*pow(nu,2.0));
			// printf("min=%le\tTmin=%le\n ",min,mp.Tmin);

		}

		printf("Model             %s\n",atm_model);
		if (is_spectrum) {
			printf("Spectrum Mode activated\n");
			printf("Frecuencies       %lf %lf %i [GHz]\n",nu1/1e9, nu2/1e9, n_nus);
		}else{
			if (wl != 0.0)
				printf("Wavelength  %lf\n",wl);
			else
				printf("Frecuency         %lf [GHz]\n",nu/1e9);
		}
		printf("z_begin           %lf [Km]\n",F);
		printf("z_end             %lf [Km]\n",H);
		printf("Rt                %le [Km]\n",Rt);
		printf("dz                %le [Km]\n",detail);
		printf("ia-step           %le [Km]\n",dl);
		printf("Resolution        %ix%i [px]\n",n,n);
		// printf("Minimal parmeter  %le [I]\n",min);
		//printf("(OUT) File Image %s\n",outputfile);

		if (chromosnet)
			printf("Chromospheric Network (chromosnet)\n");
		else
			printf("No Chromospheric Network (average)\n");
		printf("Minimal parmeter  %le [I]\n",min);
		printf("                  %lf [K]\n",mp.Tmin);



		printf("****************************************************\n");
		printf("Index\tProc\tx[px]\ty[px]\tFrequency[GHz]\tTb[K]\n");

/*    if (is_spectrum){
        printf ("Frequency[GHz]\tTb[K]\n");
    }else{

    }
 */

		/*********************************************      **
		 *********************************************      ***
		 ***********  COMIENZA LA INTEGRACION ********      ***
		 *********************************************  ***********
		 *********************************************    *******
		 *********************************************      ***
		 *********************************************       *
		 */

		// Aqui lo chido seria pasar la estructura mp directamente...

	}else{ /***************** TERMINA MAESTRO !!! ******************/
		// printf("%i Esperando...\n",miproc+1);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	//PASANDO VALORES
	if (miproc==0) {
		if (is_spectrum==1) {
			//mp.delta_nu = (nu2-nu1)/(n_nus-1); Revision 20170918
			mp.n_spectrums = n_nus/numproc;
			//printf("Delta_nu=%le\n", mp.delta_nu);
			//printf("Espectro por cada proc=%i\n", mp.n_spectrums);
		}
		//printf("numproc %i\n",numproc);
		mp.nu_begin=mp.n_spectrums; //Rev 20170918
		mp.nu_end=mp.n_spectrums; //Rev 20170918

		for (n_proc=1; n_proc < numproc; n_proc++ ) {
			//      printf("0: Enviando MP a %i\n",n_proc);
			mp.intervalo=n_proc+1;
			mp.nu_end+=mp.n_spectrums; //Rev 20170918

			//mp.nu1 = nu1 + mp.delta_nu*(n_proc*mp.n_spectrums); Rev 20170918
			//mp.nu2 = nu1 + mp.delta_nu*((n_proc+1)*mp.n_spectrums); Rev 20170918



			MPI_Send(&mp, sizeof(mp), MPI_CHARACTER, n_proc, 98, MPI_COMM_WORLD);

			mp.nu_begin=mp.nu_end; // Rev 20170918

			if (is_spectrum==1) {
				//printf("%i: (%le %le)\n",n_proc,mp.nu1,mp.nu2 );
			}

		}
		mp.intervalo=1; //el maestro trabaja.

		if (mp.is_spectrum==1) {
			//n_proc = 0;

			mp.intervalo=1;
			mp.nu_begin=0;
			mp.nu_end=mp.n_spectrums;

			//mp.nu1 = nu1;  Rev 20170918
			//mp.nu2 = mp.nu1 + mp.delta_nu*mp.n_spectrums; Rev 20170918
			//printf("%i: (%le %le)\n",0,mp.nu1,mp.nu2 );

		}

	}else{
		// printf("\n%i: Recibiendo MP!\n",miproc+1);
		MPI_Recv(&mp, sizeof(mp), MPI_CHARACTER, 0, 98, MPI_COMM_WORLD, &status);
		//    printf("\n%i: MP Recibido\n",miproc+1);
	}



	if (miproc==0) {//deleting spectrum intermediate files
		if (is_spectrum==1) {
			for (i = 0; i<numproc; i++) {
				sprintf(comando,"spectrum-%i.dat",i);
				file = fopen(comando,"w");
				// fprintf(file, "%lf\t%lf\n",mp.nu/1e9,Tb_K);
				//fprintf(file, "%lf\n",C_light*C_light*I/(2.0*K*mp.nu*mp.nu));
				fflush(file);
				fclose(file);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD); /* para sincronizar la distribucion de valores iniciales */

	//Aqui ya tenemos los intervalos, ahora necesitamos calcularlos!


	if (mp.is_spectrum == 1) {
		//printf("\nHola, soy %i %i\n",miproc, mp.n_spectrums);

		for (i = 0; i < mp.n_spectrums; i++) {
			// pa.xini = mp.xini;
			//pa.xfin = mp.xfin;
			//pa.yini = mp.yini;
			//pa.yfin = mp.yfin;
			//printf("Calculando Limites\n");
			// ra = calculaLimites(pa,1,mp.intervalo,numproc,1,numproc);
			//printf("ok\n");

			// printf("Particion %i: %i %i %i %i\n",miproc+1,ra.xini,ra.xfin,ra.yini,ra.yfin);

			//printf("\n%i: Parametros Leidos Rt= %le!\n",mp.intervalo, mp.Rt);

			//printf("%i: n=%i\n",mp.intervalo,mp.n);

			x = mp.xini;
			y = mp.yini;
			// for (x = ra.xini; x <= ra.xfin;x++ ){

			alpha = Alpha(mp.Rt,(double)x, mp.n); //ok
			//   for (y= ra.yini; y <= ra.yfin;y++ ){
			beta = Beta(mp.Rt, (double)x,(double)y,mp.n); //ok
			cuadrante = 0;
			if (alpha >= 0.0 && beta >= 0.0) cuadrante=1;
			if (alpha <= 0.0 && beta >= 0.0) cuadrante=2;
			if (alpha <= 0.0 && beta <= 0.0) cuadrante=3;
			if (alpha >= 0.0 && beta <= 0.0) cuadrante=4;
			mp.x =x;
			mp.y = y;
			mp.alpha = alpha;
			mp.beta = beta;
			mp.cuadrante=cuadrante;
			//mp.nu = mp.nu1+i*mp.delta_nu; //Rev 20170918
			mp.nu=mp.nu_vector_full[mp.nu_begin+i]; // Rev 20170918
			//printf("%i:Computing %le\n",miproc,mp.nu);


			mp.min = (1e-4) *mp.Tmin* 2.0*K*pow(nu,2.0)/pow(C_light,2.0); //Rev 20170926

			//printf ("mp.Tmin=%lf\n",mp.Tmin);



			//mp.min = (1e-4) *min* 2.0*K*pow(nu,2.0)/pow(C_light,2.0); //Rev 20170926


			I = depth_integrate(mp);
			/* Edite esta linea (ftapia) 28/03/2016 */
			//original//
			//printf("%i\t%i\t%i\t%lf\t%lf\n",miproc,x,y,mp.nu/1e9,C_light*C_light*I/(2.0*K*mp.nu*mp.nu)); Rev 20170918

			//ftapia
			//printf("%lf\t%lf\n",mp.nu/1e9,C_light*C_light*I/(2.0*K*mp.nu*mp.nu));
			Tb_K =C_light*C_light*I/(2.0*K*mp.nu*mp.nu); //Rev 20170918
			printf("[%lf,%lf]",mp.nu/1e9,Tb_K); // Rev 20170918
			//printf("%le\t%lf\n",mp.nu,C_light*C_light*I/(2.0*K*mp.nu*mp.nu));
			//fprintf(fout,"%i %i %le\n",x,y,I); //revisar que pasa con esto!
			//fflush(fout);
			// }
			// }


			// Rev 20170918 /////////////////////////////////////////////////////////////
			/* Crear un archivo sprectum_%i e imprimir */
			//printf("%i\n",i);
			//hacer esto ---> sprintf(comando,"spectrum_%i.dat",numproc);
			sprintf(comando,"spectrum-%i.dat",miproc);
			file = fopen(comando,"a");
			fprintf(file, "%lf\t%lf\n",mp.nu/1e9,Tb_K);
			//fprintf(file, "%lf\n",C_light*C_light*I/(2.0*K*mp.nu*mp.nu));
			fflush(file);
			fclose(file);
			/////////////////////////////////////////////////////////////////////
		} //for de espectros


		MPI_Barrier(MPI_COMM_WORLD); /* para sincronizar el calculo */

		// printf("%i:pase la barrera\n",miproc);
		//   FILE *fr;
		// fr = fopen("spectrum.dat","w");

		if (miproc == 0) {
			// printf("%i %le\n",0,I);
			for (i = 1; i < numproc; i++) {
				MPI_Recv (&I, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status);
				// fprintf(fr, "%le\t%le\n",mp.nu1+i*mp.delta_nu*);
				// printf("%i %le\n",i,I);
				if ((I < 0.0) && (is_restricted==1)) {
					//MPI_Abort(MPI_COMM_WORLD,-1);
					//MPI_Bcast(&error_mpi, 1, MPI_INT, 0, MPI_COMM_WORLD);
					tau_is_negative = 1;
					//printf("BEFORE:%i\n",tau_is_negative);
				}

			}//end for
			//close(fr);

			//MPI_Barrier(MPI_COMM_WORLD); Rev 20170918
			t4 = MPI_Wtime();
			printf("****************************************************\n");
			// printf("Valor de I: %lf \n",I);
			printf("\nExecution Time: %.3lf seg\n",t4-t3);

			uswtime(&utime1, &stime1, &wtime1);
			printf("\nBenchmarks (sec):\n");
			printf("real  %.3f\n",  wtime1 - wtime0);
			printf("user  %.3f\n",  utime1 - utime0);
			printf("sys   %.3f\n",  stime1 - stime0);
			printf("\n");
			printf("CPU/Wall   %.3f %% \n",
			       100.0 * (utime1 - utime0 + stime1 - stime0) / (wtime1 - wtime0));
			printf("\n");
			MPI_Barrier(MPI_COMM_WORLD); //Rev 20170918

		}else
		{
			MPI_Send(&I, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
			MPI_Barrier (MPI_COMM_WORLD);
		}










		// BEGIN Rev 20170918



		if (miproc==0) {//sorting spectrum intermediate

			if (tau_is_negative == 1) {
				printf("WARNING: TAU < 0.0\n");


				//      mp.nu_vector_full[n_nus]*1e9
				srand((unsigned) time(&t_random));

				printf("Spectrum:\n");
				sprintf(comando,"synthetic-spectrum.dat");    //<---------------------------------------------------
				file = fopen(comando,"w");
				for(i=0; i<numproc; i++) {
					random_value=1.1 + ((double)rand()/(double)(RAND_MAX)) * restricted;
					//random_value= (double)(rand()%100)*0.1*restricted;
					fprintf(file,"%lf\t%lf\n",mp.nu_vector_full[i]/1e9,random_value);
					printf("%lf\t%lf\n",mp.nu_vector_full[i]/1e9,random_value);
				}

				fflush(file);
				fclose(file);


				sprintf(comando,"jaguar/julia-iteration.dat");
				file = fopen(comando,"r");
				//  iter = fgetc(file);
				fscanf(file,"%i",&iter);

				fflush(file);
				fclose(file);

				//      printf("Spectrum:\n");
				sprintf(comando,"synthetic-spectrum-%i.dat",iter);
				file = fopen(comando,"w");
				for(i=0; i<numproc; i++) {
					fprintf(file,"%lf\t%lf\n",mp.nu_vector_full[i]/1e9,restricted);
				}

				fflush(file);
				fclose(file);






			}else{


				//      printf("DEBUG1::Computing synthetic-spectrum.dat\n");
				sprintf(comando,"spectrum-0.dat");
				file = fopen(comando,"r");
				cont=0;
				while(fscanf(file, "%lf\t%lf", &temp1,&temp2) > 0) {
					//a[i] = num;
					//printf("%f\t%f\n",temp1,temp2);
					cont++;
				}
				fclose(file);



				fc=malloc(cont*numproc*sizeof(double));
				fca=malloc(cont*numproc*sizeof(double));
				Tb=malloc(cont*numproc*sizeof(double));
				Tba=malloc(cont*numproc*sizeof(double));


				cont=0;

				for (i = 0; i<numproc; i++) {
					sprintf(comando,"spectrum-%i.dat",i);
					file = fopen(comando,"r");
					fscanf(file, "%lf\t%lf", &fc[cont],&Tb[cont]);
					cont++;

					//while(fscanf(file, "%lf\t%lf", &fc[cont],&Tb[cont]) > 0) {
					//printf("%i\t%lf\t%lf\n",cont,fc[cont], Tb[cont]);
					//  cont++;
					//}
					fclose(file);
				}
				max=cont;
				// fprintf(file, "%lf\t%lf\n",mp.nu/1e9,Tb_K);
				//fprintf(file, "%lf\n",C_light*C_light*I/(2.0*K*mp.nu*mp.nu));

				sort(0,max-1,fc,fca,Tb,Tba);

				printf("Spectrum:\n");
				sprintf(comando,"synthetic-spectrum.dat");    //<---------------------------------------------------
				file = fopen(comando,"w");
				for(i=0; i<max; i++) {
					fprintf(file,"%lf\t%lf\n",fc[i],Tb[i]);
					printf("%lf\t%lf\n",fc[i],Tb[i]);
				}

				fflush(file);
				fclose(file);


				sprintf(comando,"jaguar/julia-iteration.dat");
				file = fopen(comando,"r");
				//  iter = fgetc(file);
				fscanf(file,"%i",&iter);


				fflush(file);
				fclose(file);

				//      printf("Spectrum:\n");
				sprintf(comando,"synthetic-spectrum-%i.dat",iter);
				file = fopen(comando,"w");
				for(i=0; i<max; i++) {
					fprintf(file,"%lf\t%lf\n",fc[i],Tb[i]);
				}

				fflush(file);
				fclose(file);



				free(fc);
				free(fca);
				free(Tb);
				free(Tba);

			}//else tau_is_negative


		}




		MPI_Barrier (MPI_COMM_WORLD);

		// END Rev 20170918

	}else{ //No soy espectro IMPORTANT
		int rank, size;
		MPI_Status status;
		rank = miproc;
		size = numproc;
		if (rank != 0) { // Slaves
			int buf;
			//int Tb = 0;
			buf = rank;
			int flag_start = 1;
			//MPI_Send(&mp, sizeof(mp), MPI_CHARACTER, n_proc, 98, MPI_COMM_WORLD);

			mp.flag_start=flag_start; //first time flag_star==1
			while(1) {
				MPI_Send(&mp, sizeof(mp), MPI_CHARACTER, 0, 0, MPI_COMM_WORLD);
				MPI_Recv(&mp, sizeof(mp), MPI_CHARACTER, 0, 0, MPI_COMM_WORLD, &status);
				//printf("%i\t%i\t%i\t%lf\tComputing...\n",rank,mp.x,mp.y,mp.nu/1e9);
				mp.I_nu = depth_integrate(mp);
				mp.miproc=rank;
				// buf contiene toda la informacion para el calculo del pixel x,y
				// COMPUTE I
				//Tb = compute_in_depth(buf);
			}

		}
		else { // Master
			int sum;
			int flag = -1, res;
			MPI_Request request;
			MPI_Status status;
			int buf;
			int total = res*res;
			int nready=0;
			int apixel,bpixel;
			int x0,x1,y1;
			x0=x=mp.xini;
			x1 = mp.xfin;
			y = mp.yini;
			y1= mp.yfin;
			sum = 2-size;

			while (1) {
				if(flag != 0)
				{
					MPI_Irecv(&mp, sizeof(mp), MPI_CHARACTER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
					//if res != 0 then buf = compute_new_xy();
					//nready++;
					flag = 0;

				}
				MPI_Test(&request, &flag, &status);

				if (flag != 0) {
					//printf("recv : %d, slave : %d\n", res, status.MPI_SOURCE);
					if (mp.flag_start == 1) {
						mp.flag_start = 0;
					}else{
						printf("%i\t%i\t%i\t%i\t%lf\t%lf\n",sum,mp.miproc,mp.x,mp.y,mp.nu/1e9,C_light*C_light*mp.I_nu/(2.0*K*mp.nu*mp.nu));
					}


					if (status.MPI_SOURCE != -1) {
						sum++;
						//buf = res*res;
						alpha = Alpha(mp.Rt,(double)x, mp.n); //ok
						beta = Beta(mp.Rt, (double)x,(double)y,mp.n); //ok
						cuadrante = 0;
						if (alpha >= 0.0 && beta >= 0.0) cuadrante=1;
						if (alpha <= 0.0 && beta >= 0.0) cuadrante=2;
						if (alpha <= 0.0 && beta <= 0.0) cuadrante=3;
						if (alpha >= 0.0 && beta <= 0.0) cuadrante=4;
						mp.alpha = alpha;
						mp.beta = beta;
						mp.cuadrante=cuadrante;
						mp.x= x;
						mp.y= y;

						if (y <= (mp.n-1)/2) {
							MPI_Send(&mp, sizeof(mp), MPI_CHARACTER, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
						}

						x++;
						if (x == (x1+1) ) {
							x=x0;
							y++;
						}
						if (y > (mp.n-1)/2) {
							y = y1+1;
						}
					}
					flag = -1;
				}

				if (sum == (n*n+1) ) {
					//if (nready == total){
					break;
				}
			}
			printf("Pixel computed: %d\n", sum-1);
		}
		printf("Done! \n");


		t4 = MPI_Wtime();
		// printf("Valor de I: %lf \n",I);
		printf("Tiempo de ejecucion: %.3lf seg\n",t4-t3);

		uswtime(&utime1, &stime1, &wtime1);
		printf("\nBenchmarks (sec):\n");
		printf("real  %.3f\n",  wtime1 - wtime0);
		printf("user  %.3f\n",  utime1 - utime0);
		printf("sys   %.3f\n",  stime1 - stime0);
		printf("\n");
		printf("CPU/Wall   %.3f %% \n",
		       100.0 * (utime1 - utime0 + stime1 - stime0) / (wtime1 - wtime0));
		printf("\n");

		MPI_Abort(MPI_COMM_WORLD,MPI_SUCCESS);


	}//NO SOY ESPECTRO






//
//     // initial values are in mp
//
//
//
//
//
//   /********************************************************/
// /**************** EJECUTO EN PARALELO (TODOS) *************/
// /**********************************************************/
//
// //Falta calcular xini, xfin, yini, yfin.
//
//
//
//
//  pa.xini = mp.xini;
//  pa.xfin = mp.xfin;
//  pa.yini = mp.yini;
//  pa.yfin = mp.yfin;
//
//  ra = calculaLimites(pa,1,mp.intervalo,numproc,1,numproc);
//
//  // printf("Particion %i: %i %i %i %i\n",miproc+1,ra.xini,ra.xfin,ra.yini,ra.yfin);
//
//  //printf("\n%i: Parametros Leidos Rt= %le!\n",mp.intervalo, mp.Rt);
//
//  // printf("%i: n=%i\n",mp.intervalo,mp.n);
//   for (x = ra.xini; x <= ra.xfin;x++ ){
//     alpha = Alpha(mp.Rt,(double)x, mp.n); //ok
//     for (y= ra.yini; y <= ra.yfin;y++ ){
//       beta = Beta(mp.Rt, (double)x,(double)y,mp.n); //ok
//       cuadrante = 0;
//       if (alpha >= 0.0 && beta >= 0.0) cuadrante=1;
//       if (alpha <= 0.0 && beta >= 0.0) cuadrante=2;
//       if (alpha <= 0.0 && beta <= 0.0) cuadrante=3;
//       if (alpha >= 0.0 && beta <= 0.0) cuadrante=4;
//       mp.x =x;
//       mp.y = y;
//       mp.alpha = alpha;
//       mp.beta = beta;
//       mp.cuadrante=cuadrante;
//       //      if (miproc==0)
//       //	printf("\nComputing... %i\n",miproc);
//
//       I = integra_Profundidad(mp);
//
//       /* Edite esta linea (ftapia) 28/03/16*/
//       //original//
//       printf("%i\t%i\t%i\t%lf\t%lf\n",miproc,x,y,mp.nu/1e9,C_light*C_light*I/(2.0*K*mp.nu*mp.nu));
//
//       //ftapia
//       //printf("%lf\t%lf\n",mp.nu/1e9,C_light*C_light*I/(2.0*K*mp.nu*mp.nu));
//
//       //fprintf(fout,"%i %i %le\n",x,y,I); //revisar que pasa con esto!
//       //fflush(fout);
//     }
//   }
//
// /**************************************************/
// /**************** EJECUTO EN PARALELO (FIN)********/
// /**************************************************/
//
//
//   MPI_Barrier (MPI_COMM_WORLD);
//
//
//   if (miproc == 0) {
//     // printf("%i %le\n",0,I);
//     for (i = 1; i < numproc; i++){
//       MPI_Recv (&I, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status);
//       //      printf("%i %le\n",i,I);
//     }
//     MPI_Barrier(MPI_COMM_WORLD);
//
//     t4 = MPI_Wtime();
//    // printf("Valor de I: %lf \n",I);
//     printf("Tiempo de ejecucion: %.3lf seg\n",t4-t3);
//
//     uswtime(&utime1, &stime1, &wtime1);
//     printf("\nBenchmarks (sec):\n");
//     printf("real  %.3f\n",  wtime1 - wtime0);
//     printf("user  %.3f\n",  utime1 - utime0);
//     printf("sys   %.3f\n",  stime1 - stime0);
//     printf("\n");
//     printf("CPU/Wall   %.3f %% \n",
//         100.0 * (utime1 - utime0 + stime1 - stime0) / (wtime1 - wtime0));
//     printf("\n");
//
// }else
//     {
//       MPI_Send(&I, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
//       MPI_Barrier (MPI_COMM_WORLD);
//     }
//
//
//   } //No soy espectro
//
	MPI_Finalize ();






	/*********************************************       *
	 *********************************************      ***
	 ***********  TERMINA  LA INTEGRACION ********     *****
	 *********************************************   *********
	 *********************************************      ***
	 *********************************************      ***
	 *********************************************      ***
	 */


	//  fclose(fout);
	return 0;
}


