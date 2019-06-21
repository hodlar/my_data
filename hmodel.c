#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include "atom.h"
#include "utils.h"
#include "hmodel.h"


Model newModel(char *name){
  Model model;
  Atom *atom;
  Atmosphere atm;
  Atom atomprev;
  FILE *f;
  char *linea;
  char *limpia;
  char *fullpath;
  int i=0;
  int j;
  int a1;
  double a2,a3,a4,a5,a6,a7;


  linea = malloc(sizeof(char)*500);
  fullpath = malloc(sizeof(char)*500);

  strcpy(fullpath,"data/atmosphere/hydrostatic/");
  strcat(fullpath,name);
  strcat(fullpath,"/atoms.dat");
  f = fopen(fullpath,"r");
  if (f==NULL){
    printf("Error 16: File %s not found.\n",fullpath);
    exit(0);
  }
  while (fgets(linea, 500, f ) != NULL){
    limpia = Clean(linea);
    if (strlen(limpia)>0)
      i++;
  }
  //printf("%i\n",i);
  rewind(f);
  atom = malloc(sizeof(Atom)*i);
  model.natom=i;
  i=0;
  while (fgets(linea, 500, f ) != NULL){
    limpia = Clean(linea);
    if (strlen(limpia)>0){
      //printf("i=%i,%s\n",i,limpia);
      atomprev = newAtom("data/atoms/",limpia);
      atom[i] = atomprev;//initvalues(atomprev, "models/dummy",0);
      //printf("Name:%s\n",atom[i].name);
      //printf("M=%le A=%le IS=%i\n",atom[i].M,atom[i].A,atom[i].g );
      i++;
    }
  }
  //printf("DEBUG:%le\n",atom[0].specie[0].ce.pce[1].omegalm[2]);

  //free(atom);
  //free(linea);
  model.atom =atom;
  fclose(f);

  //Cargamos el modelo atmosferico
  strcpy(fullpath,"data/atmosphere/hydrostatic/");
  strcat(fullpath,name);
  strcat(fullpath,"/hydrostatic.dat");
  f = fopen(fullpath,"r");
  if (f==NULL){
    printf("Error 17: File %s not found.\n",fullpath);
    exit(0);
  }
  i=0;
  while (fgets(linea, 500, f ) != NULL){
    limpia = Clean(linea);
    if (strlen(limpia)>0)
      i++;
  }
  //printf("%i\n",i);
  rewind(f);
  //atm = malloc(sizeof(Atmosphere)*i);
  model.n= i-1;
  strcpy(model.name,name);

  model.atom =atom;
  model.atm = atm;
  fclose(f);
  return model;
}


void loadInitValues(Model *model, int n){
  FILE *f, *magnetic_f;
  char *linea, *linea_m;
  char *limpia, *limpia_m;
  char fullpath[300], magnetic_path[300];
  int i,j,k,cont;
  double a1,a2;
  //  Model *model = *modelref;
  int b1;
  double b2,b3,b4,b5,b6,b7,bm1,bm2,bm3,bm4;

  linea = (char*)malloc(sizeof(char)*500);

  //LEYENDO LOS DATOS ATMOSFERICOS
  strcpy(fullpath,"data/atmosphere/hydrostatic/");
  strcat(fullpath,(*model).name);
  strcat(fullpath,"/hydrostatic.dat");

  f = fopen(fullpath,"r");
//  magnetic_f = fopen(magnetic_path,"r");
  //lee Campo Magnetico
	if (f==NULL){
		printf("Error 29: File %s not found.\n",fullpath);
		exit(0);
	}

  j=0;
  while (fgets(linea, 500, f ) != NULL){
    limpia = Clean(linea);
    if (strlen(limpia)>0){
      //printf("i=%i,%s\n",i,limpia);
      if (sscanf(limpia,"%i %le %le %le %le %le %le",&b1,&b2,&b3,&b4,&b5,&b6,&b7) == 7){
	if (j==n){
	  (*model).atm.id = b1;
	  (*model).atm.z = b2;
	  (*model).atm.T = b3;
	  (*model).atm.P = b4;
	  (*model).atm.H = b5;
	  (*model).atm.V = b6;
	  (*model).atm.vt = b7;
	}
	j++;
      }
    }
  }
  if (j==0){
    printf("Error 18: in the %s file.",fullpath);
    exit(0);
  }

  //printf("Ruta:%s natom=%i\n",fullpath,model.natom);
  //Leemos el valor inicial de cada atomo para el nivel que nos piden
  if (n==0)
    printf("N_tot\t\tz\t\tT\t\t");

  for (i=0; i<(*model).natom; i++){
    sprintf(fullpath,"data/atmosphere/hydrostatic/%s/%i/%sabundance.dat",(*model).name, (*model).atm.id, (*model).atom[i].name);
    f=fopen(fullpath,"r");
    if (f==NULL){
    printf("Error 20:The file %s not found.\n",fullpath);
    exit(0);
    }
    j=0;
    while (fgets(linea, 500, f ) != NULL){
      limpia = Clean(linea);
      if (strlen(limpia)>0){
	//printf("i=%i,%s\n",i,limpia);
	if (sscanf(limpia,"%le",&a1) == 1){
	  j++;
	  (*model).atom[i].A = a1;  
	}
      }
    }

    fclose(f);

    if (j < 1 || j  > 1){
      printf("Error 21:In the file %s.\n",fullpath);
      exit(0);
    }


	/***** SYSTEM *****

    printf(">%s[ok]\n",(*model).atom[i].name);
    
    ************/

    //ahora vamos a leer las especies!
    //printf("g:%i\n",(*model).atom[i].ns);
    for(j=0;j<(*model).atom[i].ns; j++){
      //AQUI IMPRIMIMOS LOS NOMBRES DE LAS ESPECIES!
  if (n==0)     
    printf("%s%s\t\t", (*model).atom[i].name,ArabicToRoman(j+1));
      sprintf(fullpath,"data/atmosphere/hydrostatic/%s/%i/%s%sabundance.dat",(*model).name, (*model).atm.id, (*model).atom[i].name,ArabicToRoman(j+1));
      f=fopen(fullpath,"r");
      //printf("%s\n",fullpath);
      if (f==NULL){
	printf("Error 21:The file %s not found.\n",fullpath);
	exit(0);
      }
      cont=0;
      while (fgets(linea, 500, f ) != NULL){
	limpia = Clean(linea);
	if (strlen(limpia)>0){
	  //printf("i=%i,%s\n",i,limpia);
	  if (sscanf(limpia,"%le",&a1) == 1){
	    (*model).atom[i].specie[j].nrel = a1;  
	    cont++;
	  }
	}
      }
      
      if(cont != 1){
	printf("Error 22: In the file %s.\n",fullpath);
	exit(0);
      }
      fclose(f);

      //LEYENDO LOS PARAMETROS NONLTE DENTRO DE LAS ESPECIES

      sprintf(fullpath,"data/atmosphere/hydrostatic/%s/%i/%s%snonlte.dat",(*model).name, (*model).atm.id, (*model).atom[i].name,ArabicToRoman(j+1));
      f=fopen(fullpath,"r");
      //printf("%s\n",fullpath);
      if (f==NULL){
	printf("Error 23:The file %s not found.\n",fullpath);
	exit(0);
      }
      cont=0;
      while (fgets(linea, 500, f ) != NULL){
	limpia = Clean(linea);
	if (strlen(limpia)>0){
	  //printf("i=%i,%s\n",i,limpia);
	  if (sscanf(limpia,"%le %le",&a1,&a2) == 2){
	    //Revisamos si esta totalmente ionizado
	    if (cont < (*model).atom[i].specie[j].energylevels){
	      (*model).atom[i].specie[j].el[cont].b = a1;  
	      //hydrostatic version!
	      (*model).atom[i].specie[j].el[cont].nrel = a2;  
	      cont++;
	    }else if  (cont == (*model).atom[i].specie[j].energylevels) {
	      //totalmente ionizado	      
	      //printf("CHIN\n");
	      (*model).atom[i].specie[j].nk = a2;  
	      cont++;
	    }
	  }
	}
      }
      //printf("cont:%i,%i\n",cont,model.atom[i].specie[j].energylevels+1);
      if(cont != ((*model).atom[i].specie[j].energylevels+1)){
	printf("Error 22: In the file %s.\n",fullpath);
      printf("Energy Levels in the model: %i\nEnergy Levels in the file: %i\n",(*model).atom[i].specie[j].energylevels,cont);
	exit(0);
      }
      fclose(f);

      /******* SYSTEM ********
      printf("-->%s[%iL]\n",(*model).atom[i].specie[j].name,cont);
      ******************/

    }//for de las Especies


  }//for de Atomos


  if (n==0)
    printf("H-\t\tb1\t\tb_hm\t\tn_e\t\tA\t\tZ\t\tE1\t\tE2\t\tE3\n");
}


/*
 * Save the model to disk
 *
 */
void writeModel(Model model, char *name){

}




