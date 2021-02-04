  #include "uclib.h"
  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <sys/types.h>
  #include <dirent.h>
  #include <stdbool.h>
  #include <time.h>

	/* -----------------------------------------------------------------------  */
	/*                 Fourth-order Runge-Kutta user function                   */
	/*                   Compute the growth of ice particle                     */
	/* -----------------------------------------------------------------------  */

  void USERFUNCTION_EXPORT
	ComputeRadius(CoordReal *result,int size,CoordReal *index_Parcel,CoordReal *residence_Time_Parcel, CoordReal *Temperature,CoordReal *Fract_Mol_H20,CoordReal *abs_Tot_Press, CoordReal *Time,CoordReal *timeLevel)
    {

      /************************ GLOBAL VARIABLES ******************************/
			int debug_output =0; 							/* 1 == Mode deboggage */
	    extern int saturationRatioThresh;
	    extern double timeStep;
	    extern int bufferSize;
	    extern int logNormale;
	    extern int nbPartInj;
	    extern double sdtDeviation;
	    extern double maxRadius;
	    extern double minRadius;
	    extern int nbClasses;
  		extern double facteur;
			extern double meanRadius;
	    extern double alpha;
			extern double beta;
      extern char path[];
	    int i=0;
	    int u=0;   /* Variable incrementation distribution lognormale function */

	    int Count_rows = 0;  /*  Count rows in a CSV file */
	    /* 100 characters per line and 20 characters per data point limit, adjust if necessary */
	    char line[100];
	    char data[2][20];    /* 2 lines for the Parcel Index and User ComputeRadius */
	    int idata=0;
	    int read_header = 0;
	    int index_Parcel_Rayon[2]={0,0};

	    double new_Particle_Radius=0;                /* New particle radius computed */
	    double old_Particle_Radius_CSV[bufferSize];  /* Old particle radius from the CSV file */
	    int ID_Parcel_CSV_Int[bufferSize];           /* ID particle from the CSV file */
	    int ID_Parcel_Int=0;
	    double old_Particle_Radius=0;
	    double rayonParticule=0;
	    double rayonLogNormale[nbPartInj];
	    /* Variables Runge Kutta */
	    double F1,F2,k1,k2,k3,k4,A,B,C,D,E=0;
		  /* Thermodynamics variables */
	    double density,enthalpy,pressVapSat,sigma,condThermIce,condThermAir=0;
	    double diffVapor,pressPartial,pressVapSatLiq,pressVapSatIce,saturationRatio,saturationRatioIce,saturationRatioLiq=0;
      double molMassWater=0.028964;
		  double molMassAir=0.018015;
	    double gasConst=8.314;
		  double C_p_Air=(1.4*287.06)/(1.4-1);
      /************************************************************************/

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


  	  /* Generation of a table contaning the log-normale distribution of particle radius where the size is equal at the number of particle injected at each time step */
  	  if(logNormale==1)
			{

	      double deltaRadius=(maxRadius-minRadius)/(nbClasses-1);
	      double classRadius[nbClasses];
	      int radiusFunction[nbClasses];

					    for(i=0;i<nbClasses;i++)
              {
					     /* rayons de chaque classe */
					        classRadius[i]=minRadius+(i*deltaRadius);

					     /* fonction rayon est de la taille du nombre de classe, il contient le nombre de particules associe a chaque classe */
					        radiusFunction[i]=(facteur*exp(-(pow(log(classRadius[i]/meanRadius),2))/(2*(pow(sdtDeviation,2)))));

					     }
					    if (debug_output==1)
					     {
					      for(i=0;i<nbClasses;i++)
					       printf("Value of classRadius %d of table %e\n\n",i, classRadius[i]);

					      for(i=0;i<nbClasses;i++)
					       printf("Value radiusFunction %d of table %d\n\n",i,radiusFunction[i]);
					   		}

                int nbPartTot=0;
			          int g=0;

					      /* Generation of particle radius of each class with the associated radius */
					      for(i=0;i<nbClasses;i++)
					           {

					           for(g=nbPartTot;g<(nbPartTot+radiusFunction[i]);g++)
					                {
					                       rayonLogNormale[g]=classRadius[i];
					                }
					           nbPartTot=nbPartTot+radiusFunction[i];
                     }
      } /* end loop if log normale */
     /* Random soring for particles in the log-normale table */
     if(logNormale==1)
            {
					  int numberChose=0;
					  double temp=0;
					  int a=0;
					  int b=(nbPartInj+1); /* arrondir à la borne superieure */

					  int rand_a_b(int a, int b)
            {
					 		return rand()%(b-a)+a;
            }
					  for(i=0;i<nbPartInj;i++)
            {

					    numberChose=rand_a_b(0,nbPartInj);
					   /* on echange les contenus des cases i et nombre_tire */
					    temp=rayonLogNormale[i];
					    rayonLogNormale[i]=rayonLogNormale[numberChose];
					    rayonLogNormale[numberChose]=temp;
					  }
					  /* verification des rayons */
					  for(i=0;i<nbPartInj;i++)
					  {
					  	if(rayonLogNormale[i]<minRadius)
					   		rayonLogNormale[i]=meanRadius;
					  }
					  if (debug_output==1)
					   {
					     for(i=0;i<nbPartInj;i++)
					      printf("Value sorted %d of table %e\n",i,rayonLogNormale[i]);
					   }
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

	  /* ********************* Compute the particle radius ******************** */

    for (i=0; i != size; i++)
		  {
		  /* Compute the radius of new particle */
		  if(residence_Time_Parcel[i]<=timeStep*1.5)
					  {

							  /* Compute thermodynamics properties  */
							  /* Ice */
							  if(Temperature[i]<=(-5+273.15))
								  {
										density=(0.91676-1.75*0.0001*(Temperature[i]-273.15)-5*0.0000001*pow((Temperature[i]-273.15),2))*1000;
									  sigma=(104.6-0.095*(Temperature[i]-273.15))*0.001;
									  pressVapSat=100*pow(10,((9.5*((Temperature[i]-273.15))/((Temperature[i]-273.15)+265.5))+0.7858));
									  enthalpy=(597.3*(pow((273.15/Temperature[i]),0.167+3.61*0.0001*Temperature[i]))+(79.7+0.485*(Temperature[i]-273.15)-2.5*0.001*pow((Temperature[i]-273.15),2)))*4.18*1000;
                  }
							  /* Liquid water */
							  else
								  {
										density=(0.9998+0.860*0.0001*(Temperature[i]-273.15)-0.108*0.0001*pow((Temperature[i]-273.15),2))*1000;
									  sigma=(76.10-0.155*(Temperature[i]-273.15))*0.001;
									  pressVapSat=100*pow(10,((10.79574*(1-((273.16)/Temperature[i])))-(5.028*log10(Temperature[i]/273.16))+(1.50475*0.0001*(1-pow(10,(-8.2969*((Temperature[i]/273.16)-1)))))+(0.42873*0.001*(pow(10,(4.76955*(1-(273.16/Temperature[i]))))-1))+0.78614));
									  enthalpy=(597.3*4.18*1000*pow((273.15/Temperature[i]),0.167+(3.61*0.0001*Temperature[i])));
									}
							  /* Thermal conductivity of ice */
							  condThermIce=(488.19/Temperature[i])+0.4685;
							  /* Thermal conductivity of dry air */
							  condThermAir=((0.00001*5.69)+0.017*0.00001*(Temperature[i]-273.15))*418;
							  /* Water vapor diffusivity */
							  diffVapor=(0.211*pow((Temperature[i]/273.15),1.94)*1.013*10000/abs_Tot_Press[i])*0.0001;
							  /* Partial pressure of water */
							  pressPartial=Fract_Mol_H20[i]*abs_Tot_Press[i];
							  /* Saturation ratio */
							  saturationRatio=pressPartial/pressVapSat;
							  /* Saturation vapor pressure of liquid water */
							  pressVapSatLiq=100*pow(10,((10.79574*(1-((273.16)/Temperature[i])))-(5.028*log10(Temperature[i]/273.16))+(1.50475*0.0001*(1-pow(10,(-8.2969*((Temperature[i]/273.16)-1)))))+(0.42873*0.001*(pow(10,(4.76955*(1-(273.16/Temperature[i]))))-1))+0.78614));
                /* Saturation vapor pressure of ice */
                pressVapSatIce=100*pow(10,((9.5*((Temperature[i]-273.15))/((Temperature[i]-273.15)+265.5))+0.7858));
                /* Saturation ratio of liquide water */
							  saturationRatioLiq=pressPartial/pressVapSatLiq;
                /* Saturation ratio of ice */
								saturationRatioIce=pressPartial/pressVapSatIce;

                /* Compute Runge Kutta parameters */
							  A=(2*sigma*molMassWater)/(density*gasConst*Temperature[i]);
							  B=(pow(enthalpy,2)*molMassWater*density)/(condThermAir*gasConst*pow(Temperature[i],2));
							  C=(gasConst*Temperature[i]*density)/(pressVapSat*diffVapor*molMassWater);
							  D=(condThermAir*(sqrt(2*3.14159*molMassAir*gasConst*Temperature[i])))/((alpha*abs_Tot_Press[i]*molMassAir*(C_p_Air-(gasConst/2))));
							  E=(diffVapor/beta)*(sqrt(2*3.14159*molMassWater)/(gasConst*Temperature[i]));

             /* Choice of initial particle radius depending of the initial distribution */

					   if(logNormale==0)
					     {
					     	  rayonParticule=meanRadius;
					     }
					   else
					     {
					        rayonParticule=rayonLogNormale[u];
                  u++;
					     }

			          if (debug_output == 1)
                  		{
							  printf("Compare residence Time Step %f vs time step %f\n",residence_Time_Parcel[i],timeStep);
							  printf("Initial radius OK\n");
										  }
					  if (saturationRatioLiq > saturationRatioThresh)
						  {
							  if (debug_output==1)
	        					printf("\tSaturation\n");
                /* Compute Runge Kutta parameters */
							  F1=rayonParticule/(rayonParticule+D);
							  F2=rayonParticule/(rayonParticule+E);
							  k1=timeStep*(1/rayonParticule)*(saturationRatio-exp(A/rayonParticule))/
							  (((B/F1)*exp(A/rayonParticule))+(C/F2));

							  F1=(rayonParticule+(k1/2))/(rayonParticule+(k1/2)+D);
							  F2=(rayonParticule+(k1/2))/(rayonParticule+(k1/2)+E);
							  k2=timeStep*(1/(rayonParticule+(k1/2)))*(saturationRatio-exp(A/(rayonParticule+(k1/2))))/
							  (((B/F1)*exp(A/(rayonParticule+(k1/2))))+(C/F2));

							  F1=(rayonParticule+(k2/2))/(rayonParticule+(k2/2)+D);
							  F2=(rayonParticule+(k2/2))/(rayonParticule+(k2/2)+E);
							  k3=timeStep*(1/(rayonParticule+(k2/2)))*(saturationRatio-exp(A/(rayonParticule+(k2/2))))/
							  (((B/F1)*exp(A/(rayonParticule+(k2/2))))+(C/F2));

							  F1=(rayonParticule+(k3))/(rayonParticule+(k3)+D);
							  F2=(rayonParticule+(k3))/(rayonParticule+(k3)+E);
							  k4=timeStep*(1/(rayonParticule+k3))*(saturationRatio-exp(A/(rayonParticule+k3)))/
							  (((B/F1)*exp(A/(rayonParticule+k3)))+(C/F2));

							  if (debug_output==1)
				    			{
									  printf("Addition plus radius %e\n",(k1+2*k2+2*k3+k4)/6);
									  printf("\t\t%f , %f , %e , %f , %f , %f, %f, %f, %f \n",F1,F2,k1,k2,k3,k4,timeStep,saturationRatio,A);
									}
						  }
					  else
						  {
							  k1=0;
							  k2=0;
							  k3=0;
							  k4=0;
						  }

              new_Particle_Radius=rayonParticule+(k1+2*k2+2*k3+k4)/6;

					  if (debug_output==1)
						  printf("\t\t\tResults loop radius : %d %e\n",i+1,new_Particle_Radius);
				  }
			  /* End first if loop */
			  else if (residence_Time_Parcel[i]>timeStep)
				  {

					  /* Recover ID Parcel */
					  int Stockage_ID=0;
					  ID_Parcel_Int=(int)index_Parcel[i]; /* Parcel STARCCM+ */
					  Stockage_ID=ID_Parcel_Int;
					  int k=0;
					  if (debug_output==1)
			    		{
							  printf("Compare residence Time Step %f vs time step %f\n",residence_Time_Parcel[i],timeStep);
							  printf("\tStore number: %d\n",Stockage_ID);
							  printf("\t\t Search Parcel index in table\n");
							  printf("\t\t\tID_Parcel_CSV : %d\n",ID_Parcel_CSV_Int[k]);
						  }
					  if (debug_output==1)
						  printf("\tStore number: %d\n",Stockage_ID);
					  while(Stockage_ID!=ID_Parcel_CSV_Int[k]) /* Compare the ID Parcel in StarCCM+ and the ID Parcel from the table */
						  {
							  k++;
						  if (k > bufferSize)
							  {
							  printf("OVERFLOW BUFFER \n");
							  return;
							  }
						  }

					  if (debug_output==1)
						  printf("Line in buffer %d \n",k);

					  int particleFound=0;
					  if(Stockage_ID==ID_Parcel_CSV_Int[k])
					  {
						  if (debug_output==1)
						  printf("ID Parcel : %d \n",ID_Parcel_CSV_Int[k]);
              printf("\t\t\tAssociated radius: %e\n",old_Particle_Radius_CSV[k]);
              printf("\t\t\tSearch Parcel index in table OK\n");
						  particleFound=1;
					  }
					  else
						  {


						  if (debug_output==1)
							  printf("Particle not found\n");
                return;
						  }
		        old_Particle_Radius=old_Particle_Radius_CSV[k]; /* Recover the old particle radius */

            /* Compute thermodynamics properties  */
            /* Ice */
            if(Temperature[i]<=(-5+273.15))
              {
                density=(0.91676-1.75*0.0001*(Temperature[i]-273.15)-5*0.0000001*pow((Temperature[i]-273.15),2))*1000;
                sigma=(104.6-0.095*(Temperature[i]-273.15))*0.001;
                pressVapSat=100*pow(10,((9.5*((Temperature[i]-273.15))/((Temperature[i]-273.15)+265.5))+0.7858));
                enthalpy=(597.3*(pow((273.15/Temperature[i]),0.167+3.61*0.0001*Temperature[i]))+(79.7+0.485*(Temperature[i]-273.15)-2.5*0.001*pow((Temperature[i]-273.15),2)))*4.18*1000;
              }
            /* Liquid water */
            else
              {
                density=(0.9998+0.860*0.0001*(Temperature[i]-273.15)-0.108*0.0001*pow((Temperature[i]-273.15),2))*1000;
                sigma=(76.10-0.155*(Temperature[i]-273.15))*0.001;
                pressVapSat=100*pow(10,((10.79574*(1-((273.16)/Temperature[i])))-(5.028*log10(Temperature[i]/273.16))+(1.50475*0.0001*(1-pow(10,(-8.2969*((Temperature[i]/273.16)-1)))))+(0.42873*0.001*(pow(10,(4.76955*(1-(273.16/Temperature[i]))))-1))+0.78614));
                enthalpy=(597.3*4.18*1000*pow((273.15/Temperature[i]),0.167+(3.61*0.0001*Temperature[i])));
              }
            /* Thermal conductivity of ice */
            condThermIce=(488.19/Temperature[i])+0.4685;
            /* Thermal conductivity of dry air */
            condThermAir=((0.00001*5.69)+0.017*0.00001*(Temperature[i]-273.15))*418;
            /* Water vapor diffusivity */
            diffVapor=(0.211*pow((Temperature[i]/273.15),1.94)*1.013*10000/abs_Tot_Press[i])*0.0001;
            /* Partial pressure of water */
            pressPartial=Fract_Mol_H20[i]*abs_Tot_Press[i];
            /* Saturation ratio */
            saturationRatio=pressPartial/pressVapSat;
            /* Saturation vapor pressure of liquid water */
            pressVapSatLiq=100*pow(10,((10.79574*(1-((273.16)/Temperature[i])))-(5.028*log10(Temperature[i]/273.16))+(1.50475*0.0001*(1-pow(10,(-8.2969*((Temperature[i]/273.16)-1)))))+(0.42873*0.001*(pow(10,(4.76955*(1-(273.16/Temperature[i]))))-1))+0.78614));
            /* Saturation vapor pressure of ice */
            pressVapSatIce=100*pow(10,((9.5*((Temperature[i]-273.15))/((Temperature[i]-273.15)+265.5))+0.7858));
            /* Saturation ratio of liquide water */
            saturationRatioLiq=pressPartial/pressVapSatLiq;
            /* Saturation ratio of ice */
            saturationRatioIce=pressPartial/pressVapSatIce;

            /* Compute Runge Kutta parameter */
            A=(2*sigma*molMassWater)/(density*gasConst*Temperature[i]);
            B=(pow(enthalpy,2)*molMassWater*density)/(condThermAir*gasConst*pow(Temperature[i],2));
            C=(gasConst*Temperature[i]*density)/(pressVapSat*diffVapor*molMassWater);
            D=(condThermAir*(sqrt(2*3.14159*molMassAir*gasConst*Temperature[i])))/((alpha*abs_Tot_Press[i]*molMassAir*(C_p_Air-(gasConst/2))));
            E=(diffVapor/beta)*(sqrt(2*3.14159*molMassWater)/(gasConst*Temperature[i]));

					  if (saturationRatioLiq > saturationRatioThresh)

						  {
							  if (debug_output==1)
			    						printf("\tSaturation\n");
                /* Compute Runge Kutta parameter */
							  F1=old_Particle_Radius/(old_Particle_Radius+D);
							  F2=old_Particle_Radius/(old_Particle_Radius+E);
							  k1=(timeStep*(1/old_Particle_Radius)*(saturationRatio-exp(A/old_Particle_Radius))/
							  (((B/F1)*exp(A/old_Particle_Radius))+(C/F2)));

							  F1=(old_Particle_Radius+(k1/2))/(old_Particle_Radius+(k1/2)+D);
							  F2=(old_Particle_Radius+(k1/2))/(old_Particle_Radius+(k1/2)+E);
							  k2=(timeStep*(1/(old_Particle_Radius+(k1/2)))*(saturationRatio-exp(A/(old_Particle_Radius+(k1/2))))/
							  (((B/F1)*exp(A/(old_Particle_Radius+(k1/2))))+(C/F2)));

							  F1=(old_Particle_Radius+(k2/2))/(old_Particle_Radius+(k2/2)+D);
							  F2=(old_Particle_Radius+(k2/2))/(old_Particle_Radius+(k2/2)+E);
							  k3=(timeStep*(1/(old_Particle_Radius+(k2/2)))*(saturationRatio-exp(A/(old_Particle_Radius+(k2/2))))/
							  (((B/F1)*exp(A/(old_Particle_Radius+(k2/2))))+(C/F2)));

							  F1=(old_Particle_Radius+(k3))/(old_Particle_Radius+(k3)+D);
							  F2=(old_Particle_Radius+(k3))/(old_Particle_Radius+(k3)+E);
							  k4=(timeStep*(1/(old_Particle_Radius+k3))*(saturationRatio-exp(A/(old_Particle_Radius+k3)))/
							  (((B/F1)*exp(A/(old_Particle_Radius+k3)))+(C/F2)));
							  if (debug_output==1)
				    		{
									  printf("%e %e %e %e\n",k1,k2,k3,k4);
									  printf("%e %e %e %e\n",timeStep,old_Particle_Radius,saturationRatio,exp(A/old_Particle_Radius));
								}
					   }
					  else
						 {
							  k1=0;
							  k2=0;
							  k3=0;
							  k4=0;
						 }
					  if (debug_output==1)
			    			{
							  printf("Addition plus radius %e\n",(k1+2*k2+2*k3+k4)/6);
							  printf("%e , %e , %e , %e , %e , %e, %e, %e, %e , %e,%e,%e, %e\n",F1,F2,k1,k2,k3,k4,timeStep,saturationRatio,A,B,C,D,E);
							  }

		          new_Particle_Radius=old_Particle_Radius+(k1+2*k2+2*k3+k4)/6;

					  if (debug_output==1)
						  printf("\tResults loop radius : %d %e\n",i+1,new_Particle_Radius);
				  } 	/* fermeture deuxieme boucle if */

				  if (debug_output==1)
	        			printf("line %d : Etat  %d\n",i+1,isnan(new_Particle_Radius!=0));

						 if(new_Particle_Radius<0)
					        new_Particle_Radius=0;

				  result[i]=new_Particle_Radius;

      } /* End  for loop */
	      if (debug_output==1)
	         		printf("End ComputeRadius function\n");


    }
