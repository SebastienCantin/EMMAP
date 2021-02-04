#include "uclib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <dirent.h>

/* -------------------------------------------------------------------------  */
/*              Compute the rate of condensation of each particle             */
/* -------------------------------------------------------------------------  */

void USERFUNCTION_EXPORT
  RateCondensation(CoordReal *result,int size,CoordReal *index_Parcel,CoordReal *residence_Time_Parcel,CoordReal *Time,CoordReal *timeLevel,CoordReal *particleRadius,CoordReal *Volume,CoordReal *Temperature,CoordReal *Fract_Mol_H20,CoordReal *Abs_Tot_Press)
  {

	/************************ GLOBAL VARIABLES ******************************/
	int debug_output =0; 									/* 1 == Mode deboggage */
	extern int saturationRatioThresh;
	extern double timeStep;
	extern int bufferSize;
  extern double meanRadius;
	extern double nbParcel;
  extern char path[];

	int q=0;
	int i=0;
	int Count_rows = 0;
  /* 100 characters per line and 20 characters per data point limit, adjust if necessary */
  char line[100];
	char data[2][20]; /* 2 lines for the Parcel Index and User ComputeRadius */
	int idata=0;
	int read_header = 0;
	int index_Parcel_Rayon[2]={0,0};
	double old_Particle_Radius_CSV[bufferSize]; /* Old particle radius from the CSV file */
	int ID_Parcel_CSV_Int[bufferSize]; /* ID particle from the CSV file */
	int ID_Parcel_Int=0;
	double old_Particle_Radius=0;

	double VaporizationRate,pressVapSat=0;
	double pressPartial,pressVapSatLiq,saturationRatio,saturationRatioLiq,density=0;

  /******* When running in parallel, it might happen that the user ********/
  /******* code does not need to be evaluated for any cells ***************/
  /******* of a particular partition. We can exit immediately in **********/
  /******* this case. *****************************************************/
  if (size==0)
  {
    if (debug_output==1)
      {
        printf("No work to do on this partition. Exiting... \n");
        fflush(stdout);
      }
    return;
  }


  if(timeLevel[0]>=2)
    {
  /* ************ Recover old particle radius from a CSV file ************* */
  /* Header order for check "Parcel Index","User ComputeRadius","X (m)","Y (m)","Z (m)" */
  char filename[100]; /* Name of the file */
  FILE* myfile;
  sprintf(filename,"%s/Particules_CSV/Particule_%e.csv",path,Time[0]-timeStep);
  if (debug_output==1)
    printf("File name is %s \n",filename);
  /* -----Find out how much memory needs to be allocated ------ */
  /* Check number of lines inside the CSV file                  */
  /* Please note: This is only a maximum check as there may be  */
  /* empty lines at the end of the file                         */
  myfile = fopen(filename, "r");
  int ch, numberOfLines = 0;

  /* Comparison of the name file in "Filename" and the path */
  if (strncmp(filename,"%s/Particules_CSV/Particule_%e.csv",path,Time[0]-timeStep)==0)
  printf("Les chaines de characteres sont les memes \n");
  if (debug_output==1)
    {
         printf("Comparaison :\n");
         printf("1:%s\n",filename);
    }
  if (debug_output==1)
      printf("2:%s/Particules_CSV/Particule_%e.csv\n",path,Time[0]-timeStep);

  /* Return of the file does not exist */
    if (myfile==NULL)
      {
        printf("%s: No file found. Check the location of your file or specify absolute path in user code Particule_radius.c \n",filename);
        fflush(stdout);
        return;
      }
  /* Count the number of line in the file */
  do
    {
      ch = fgetc(myfile); /* Read a charactere in the file myfile */
      if(ch == '\n')
      numberOfLines++;
    }
  while (ch != EOF);
  if(ch != '\n' && numberOfLines != 0) /* If the charactere read is different of a line break and if the number of line is different of 0, numberOfLines is incremented */
    numberOfLines++;
  if (debug_output==1)
    {
      printf("%s: %i lines of data were found (this includes the header and the last line might be empty) \n",filename, numberOfLines);
      fflush(stdout);
    }
  fclose(myfile);

  /* **************** Read data set from CSV file ************************* */
  myfile = fopen(filename, "r");
  int consistency_check = 0;
  while (fgets(line,100,myfile) != NULL) /* fgets , read  "line" , in the file myfile, 100 char max per line */
    {
         /* Read the header line if not already done */
         if (read_header==0)
              {
                  read_header = 1;
                  consistency_check=sscanf(line,"%[^,],%[^,]", data[0],data[1]); /* sscanf, read hearder line and store into data[0] and data [1] */
                  if (debug_output==1)
                    {
                      printf("The Following header was found: %s, %s \n", data[0],data[1]);
                      fflush(stdout);
                    }
              /* Store the order of the data for later use */
              /* strncmp compare the hearder line data */
              for (idata = 0; idata < 2; idata++)
                {
                  if (strncmp(data[idata],"\"Parcel Index\"",14)==0) /* data[0]= Parcel Index */
                    {
                      index_Parcel_Rayon[idata] = 0;
                          if (debug_output==1)
                            printf("A \n");
                    }
                  else if(strncmp(data[idata],"\"User ComputeRadius\"",20)==0) /* data[1]= ComputeRadius */
                    {
                      index_Parcel_Rayon[idata] = 1;
                          if (debug_output==1)
                            printf("B \n");
                    }
                }
              if (debug_output==1)
                  {
                    printf("The following positions were detected for the data according to the header file \n");
                    printf("Parcel Index=0, User ComputeRadius=1: %i, %i \n",index_Parcel_Rayon[0],index_Parcel_Rayon[1]);
                    fflush(stdout);
                  }
                } /* end boucle if (read_header==0) */
                /* Read data line */
              else
                {
                  consistency_check=sscanf(line,"%[^,],%[^,],%*[^,],%*[^,],%*s", data[0],data[1]);
                  if (debug_output==1)
                    {
                      printf("Check in line \n %s \n",line);
                      printf("Check in data 0 \n %s \n",data[0]);
                      printf("Check in data 1 \n %s \n",data[1]);
                    }
                  /* Fill the variable (Parcel Index and User ComputeRadius) from CSV file in a static table */
                  for (idata = 0; idata < 2; idata++)
                    {
                      if (index_Parcel_Rayon[idata]==0)
                        {
                          ID_Parcel_CSV_Int[Count_rows] = atoi(data[idata]);/* atof converti un argument string en un nombre int */
                          if (debug_output==1)
                            printf("Check in ID_Parcel_CSV : %d\n",ID_Parcel_CSV_Int[Count_rows]);
                        }
                      else
                        {
                          old_Particle_Radius_CSV[Count_rows] = atof(data[idata]);
                            if (debug_output==1)
                              printf("Check in old_Particle_Radius_CSV : %e\n",old_Particle_Radius_CSV[Count_rows]);
                        }
                     }
                  Count_rows++;


                } /* End loop else */
                  /* Consistency check */
                  if (consistency_check!=2)
                    {
                      printf("Less than 2 items were read in at least one line of %s. Check file format \n",filename);
                      printf("The user function is exiting now ... Please restart STAR-CCM+ to make sure is in a proper state. \n");
                      fflush(stdout);
                      fclose(myfile);
                      return;
                    }
          } /* End Read CSV file, end while loop */
        fclose(myfile);


        /* ****** Sort in ascending order according to the Parcel ID ****** */
        /* Quick Sort function */
        int begin=0;
        int end=Count_rows;

        void recursiveSorting(double ID_Parcel_CSV_Int[],int begin, int end)
        {
          const int pivot_ID=ID_Parcel_CSV_Int[begin]; /* Pivot is the first case of the table */
          const double pivot_rayon=old_Particle_Radius_CSV[begin];
          int pos=begin;
          int g;
          if(begin>=end)
          return;

          for (g=begin;g<end;g++)
          {

            if(ID_Parcel_CSV_Int[g]<pivot_ID)
              {
                ID_Parcel_CSV_Int[pos]=ID_Parcel_CSV_Int[g];
                old_Particle_Radius_CSV[pos]=old_Particle_Radius_CSV[g];
                pos++;
                ID_Parcel_CSV_Int[g]=ID_Parcel_CSV_Int[pos]; /* deplace ID_Parcel_CSV_Int[g] vers la gauche dune case */
                old_Particle_Radius_CSV[g]=old_Particle_Radius_CSV[pos];
                ID_Parcel_CSV_Int[pos]=pivot_ID; /* le pivot est deplace à la place de pos juste */
                old_Particle_Radius_CSV[pos]=pivot_rayon;
              }
          }
          recursiveSorting(ID_Parcel_CSV_Int,begin,pos);
          recursiveSorting(ID_Parcel_CSV_Int,pos+1,end);
      }
   if (debug_output==1)
   {
         int s=0;
         for(s=0;s<Count_rows;s++)
         {
         printf("ID_Parcel_CSV row %d ID %d\n",s,ID_Parcel_CSV_Int[s]);
         printf("old_Particle_Radius_CSV row %d radius %e\n",s,old_Particle_Radius_CSV[s]);
         }
    }
  }
  /* ********************* Compute the rate of vaporization ***************** */

	for (i=0; i != size; i++)
		{
      /* Compute thermodynamics properties  */
      /* Ice */
			if(Temperature[i]<=(-5+273.15))
			   {
				       density=(0.91676-(1.75*0.0001*(Temperature[i]-273.15))-(5*0.0000001*pow((Temperature[i]-273.15),2)))*1000;
				 }
      /* Liquid water */
			else
				 {
								density=(0.9998+0.860*0.0001*(Temperature[i]-273.15)-0.108*0.0001*pow((Temperature[i]-273.15),2))*1000;
				 }


			if(residence_Time_Parcel[i]<=1.5*timeStep)

				{
					VaporizationRate=(density*4*3.14159*(particleRadius[i]-meanRadius)*(pow(particleRadius[i],2))/timeStep)*(nbParcel/Volume[i]);;
      	}
				else if (residence_Time_Parcel[i]>timeStep)
				{
          /* Recover ID Parcel */
          int Stockage_ID=0;
          ID_Parcel_Int=(int)index_Parcel[i];
          Stockage_ID=ID_Parcel_Int; /* Parcel STARCCM+ */
          int k=0;
          if (debug_output==1)
            {
              printf("Evaluation boucle Residence %f vs Pdt %f\n",residence_Time_Parcel[i],timeStep);
              printf("\tNumero de Stockage: %d\n",Stockage_ID);
              printf("\t\tCherche Parcel Starccm+ vs csv\n");
              printf("\t\t\tID_Parcel_CSV : %d\n",ID_Parcel_CSV_Int[k]);
            }
          if (debug_output==1)
            printf("\tNumero de Stockage: %d\n",Stockage_ID);
          while(Stockage_ID!=ID_Parcel_CSV_Int[k]) /* Compare the ID Parcel in StarCCM+ and the ID Parcel from the table */
            {
              k++;
            if (k > bufferSize)
              {
              printf("DEPASSEMENT CAPACITE BUFFER \n");
              return;
              }
            }

          if (debug_output==1)
            printf("Comptage ligne buffer %d \n",k);

          int particule_trouvee=0;
          if(Stockage_ID==ID_Parcel_CSV_Int[k])
          {
            if (debug_output==1)
            printf("ID Parcel : %d \n",ID_Parcel_CSV_Int[k]);
            printf("\t\t\tRayon associe: %e\n",old_Particle_Radius_CSV[k]);
            printf("\t\t\tRecherche Parcel Starccm+ vs csv OK\n");
            particule_trouvee=1;
          }
          else
            {


            if (debug_output==1)
              printf("Particule non trouvee\n");
              return;
            }
          old_Particle_Radius=old_Particle_Radius_CSV[k]; /* Recover the old particle radius */

					if (debug_output==1)
			  			{
				  printf("Numero de ID Parcel CSV: %d\n",index_Parcel[i]);
					printf("Rayon actuel: %f\n",particleRadius[i]);
					printf("Ancien rayon: %f\n",old_Particle_Radius);
							}
          /* Compute vaporization rate */
					VaporizationRate=(density*4*3.14159*(particleRadius[i]-old_Particle_Radius)*(pow(particleRadius[i],2))/timeStep)*(nbParcel/Volume[i]);
					}
					if (debug_output==1)

						if(timeLevel[0]>=2)
								result[i]=VaporizationRate;
						else
								result[i]=0;

		} /* End for loop */

	if (debug_output==1)
	    printf("Fin de la fonction taux de condensation\n");

  }
