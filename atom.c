#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "atom.h"
#include "utils.h"
#include "specie.h"

Atom newAtom(char *path, char *name){
  Atom atom;
  FILE *f;
  char fullpath[500] = "";
  char fullpathsearch[500] = "";
  char *roman;
  char *speciestr;
  char *speciename;
  int len = strlen(name);
  char *linea;
  char *limpia;
  int j,i=0;
  int cont = 0;
  Specie *specie; 

  linea = malloc(sizeof(char)*500);



  //Cargando la cabezera del atomo
  strcpy(atom.name,name);
  atom.name[len-1] = '\0';
  strcat(fullpath,path);
  strcat(fullpath,atom.name);
  strcat(fullpath,"/header.dat");
  printf("> %s\n",atom.name);
  f = fopen(fullpath,"r");
  if (f == NULL){
    printf("Error 01: The file %s not found!\n",fullpath);
    printf("You can fix this error in the jaguar.conf file.\n");
    exit(0);
  }
  cont = 0;
  while ((fgets(linea, 500, f )!= NULL) || (cont==0) ){
    //printf("%i:%s\n",i,linea);
    limpia = Clean(linea);
    if (strlen(limpia)>0){
      //printf("%i:%s",i,limpia);
      //La primera linea limpia debe tener la info
      //M		A		IS
      if (sscanf(limpia,"%le %i",&atom.M, &atom.g) == 2){
	atom.n = 0.0;
	atom.A=0.0;
	atom.ne=0.0;
	cont = 1;
      }
    }
    i++;
  }
  free(linea);
  fclose(f);
  //Encontramos la cabecera del atomo?
  if (cont==1){
    //Continuamos (Buscando Especies)
    i = 0;
    j=0;
    cont = 0;
    roman = malloc(sizeof(char)*20);
    speciestr = malloc(sizeof(char)*30);
    speciename = malloc(sizeof(char)*30);

    //Primero revisamos cuantas especies existen
   while (i<=atom.g){
      roman = ArabicToRoman(i+1);
      sprintf(speciestr,"/%s%s/",atom.name,roman);
      i++;
      strcpy(fullpath,"");
      strcat(fullpath,path);
      strcat(fullpath,atom.name);
      strcat(fullpath,speciestr);
      //      strcat(fullpath,"energylevels.dat");
	if (access(fullpath, R_OK) == 0)
	  j++;
	/*
      f = fopen(fullpath,"r");
      if (f != NULL){
	fclose(f);
	j++;
      }
	*/
   }

   //Ahora, reservamos espacio en la memoria
   //y creamos la lista de especies.
   //   printf("Reservando:%i\n",j);
   specie= malloc(sizeof(Specie)*j);
   atom.ns =j; //verificar
   //printf("ns=%i\n",j);
   i=0;
   j=0;
   while (i<=atom.g){
      roman = ArabicToRoman(i+1);
      sprintf(speciestr,"/%s%s/",atom.name,roman);
      sprintf(speciename,"%s%s",atom.name,roman);
      strcpy(fullpath,"");
      strcat(fullpath,path);
      strcat(fullpath,atom.name);
      strcat(fullpath,speciestr);
      strcpy(fullpathsearch,fullpath);
      strcat(fullpathsearch,"energylevels.dat");
      if (i < atom.g){
	f = fopen(fullpathsearch,"r");
	if (f == NULL){
	  printf("-->%s[not found]\n",speciename);
	}else{
	  fclose(f);
	  //printf("Reading %s \n",fullpath);
	  specie[j] = newSpecie(fullpath,speciename,i);
	  j++;
	}
      }else{ // es el ultimo estado (totalmente ionizado)

	if (access(fullpath, R_OK) == 0){
        strcpy(specie[j].name,atom.name);
        strcat(specie[j].name,roman);
	printf("-->%s[1L]\n",speciename);
        specie[j].ionizationstage = i;
	specie[j].energylevels=0;
	}else{
	  printf("-->%s[not found]\n",speciename);
	}
      }
      i++;
   }
   atom.specie =specie;

   //free(speciename);
   //free(specie);
   //free(speciestr);
   //free(roman);
   return atom;
  }else{
    //No encontramos la cabecera del atomo
    printf("Error 02: The file %s has no valid lines.\n", fullpath);
    exit(0);
  }
}
