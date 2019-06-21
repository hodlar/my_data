#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hmodel.h"
#include "B.h"

int main(int argc, char **argv){
  int i,j;
  Model model;
  char env[500];
  char comando[500];
  double hydro_step;
  double pz1;
  double Y;
  double fz1;
  double dx;
  int 	chromospheric_network;// ADD
  int cell;
  int hydro;
  double amplitude, intensity, alpha;
  int magnetic = 0;
  int B_x0;
  
  printf("Loading Atomic Model:\n");


  model = newModel(argv[1]);
  chromospheric_network=0;
  hydro = 0;

   for (i=1; i<argc;i++){
	   sprintf(comando,"%s",argv[i]);
	   if (strcmp(comando,"-hydro") == 0){
		   hydro= 1;
	   }

	   if (strcmp(comando,"-cn") == 0){
		   sprintf(comando,"%s",argv[++i]);
		   if (sscanf(comando,"%i\n",&cell) > 0){
			   chromospheric_network=1;
			   printf("Error: Network or cell is required.\n");
			   return 0;
		   }
	   }
	   if (strcmp(comando, "-B") == 0){
		   sprintf(comando, "%s",argv[++i]);
		   if(sscanf(comando,"%lf\n",&amplitude) > 0){
			   sprintf(comando, "%s", argv[++i]);
			   if(sscanf(comando,"%lf\n",&intensity) > 0){
				   sprintf(comando, "%s", argv[++i]);
				   if(sscanf(comando,"%lf\n",&alpha) > 0){
					   sprintf(comando, "%s", argv[++i]);
					   magnetic = 1;
				   }
			   }
		   }
//		   calculate_B(amplitude, intensity, alpha);
	   }
   }

  //if (argc > 2){ // ADD
//      if (strcmp(argv[2],"-cn") == 0){ // ADD
//	      if (strcmp(argv[3],"1") == 0){ // ADD
//		cell=1;
//		printf("jaguar cell=1\n");
//	      }else{
//		cell=0;
//		printf("jaguar cell = 0\n");
//	      }
//	      
//	chromospheric_network = 1;// ADD
  //    }// ADD
  //  }else{// ADD
 //   	chromospheric_network = 0;// ADD
//  }// ADD

  
  B_x0 = 1; 
  loadInitValues(&model, 0); //we use nH and P as initial values in our model.
  if(magnetic == 1){
	  init_B(amplitude, intensity, alpha, &model, B_x0);
  }
  printf("Layer 0\n");
  NLTE(&model,1e-14,0,0.0,0.0,0.0,chromospheric_network,cell,0);
  writeModel(model,"dummy/");

  for (i=1;i<= model.n;i++){
    pz1 = model.atm.P;
    fz1 = model.atm.fz;
    dx = model.atm.z;
    loadInitValues(&model, i);
	if(magnetic == 1){
//(double amplitude, double intensity, double alpha, Model *model, int x0)
		calculate_B(amplitude, intensity, alpha, &model, B_x0);
	}
    if (!(hydro)){
      printf("Layer %i (ion)\n",i);
      NLTE(&model,1e-14,0,0.0,0.0,0.0,chromospheric_network,cell,0);
    }else{
      printf("Layer %i (hydro)\n",i);
      dx = model.atm.z - dx;
      NLTE(&model,1e-14,1,pz1,fz1,dx,chromospheric_network,cell,magnetic);
    }
    writeModel(model,"dummy/");
    //   for (j=0;j<model.natom;j++)
    //  printf("A:%le\n",model.atom[j].A);
    //printf("X:%i\n",model.atm[11].id);
  }
  
  //printf("%i argc\n",argc);

  return 0;
}
