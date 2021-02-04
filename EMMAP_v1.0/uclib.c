	#include "uclib.h"
	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	/* GLOBAL VARIABLES */
	int saturationRatioThresh=1;
	int logNormale=0;		/* Choice initialization of radius distribution : monodisperse (0) or polydisperse (1) */
		/* VARIABLES FOR POLYDISPERSE */
		int nbPartInj=142071; 		/* Number of particles injected per time step */
		double sdtDeviation=1.6;	/* Standard deviation parameter for the log normale distribution */
		double maxRadius=30E-9;		/* Max particle radius parameter for the log normale distribution */
		double minRadius=5E-9;		/* Min particle radius parameter for the log normale distribution */
		int nbClasses=200;				/* Number of bin radius for the log normale distribution */
		double facteur=743.209;		/* Factor to adjust zith the Excel File, See (Maglaras,2007) for more details about the log normale distribution */
		/* ----------------------------- */
	double meanRadius=20E-9;		/* Initialization particles mean radius */
	double timeStep=0.001;			/* Simulation time step */
	int bufferSize=100000;			/* Size buffer, sufficiently high to contain all the particle in the domain until the end of the simulation */
	double nbParcel=10^7;			/* Number of soot particles per parcel */
	double beta=0.022;			/* Condensation coefficient (Fukuta & Myers, 2007) */
	double alpha=0.68;			/* Accomodation coefficient (Fukuta & Myers, 2007) */
	/* ------------------------------------- */
  char path[]="/home/scantin/Documents/Simulations/Test/Sim1"; /* Path to get the CSV file */
	// Prototypes function //
	void ComputeRadius(CoordReal*, int, CoordReal*);
	void RateCondensation(CoordReal*, int, CoordReal*);
  void SumCondensation(CoordReal*, int, CoordReal*);

void
    USERFUNCTION_EXPORT uclib()
        {

            /* Register user functions here */
            ucfunc(ComputeRadius, "ScalarFieldFunction", "ComputeRadius");
            ucarg(ComputeRadius, "Parcel", "$ParcelId", sizeof(CoordReal));
            ucarg(ComputeRadius, "Parcel", "$ParticleResidenceTime", sizeof(CoordReal));
            ucarg(ComputeRadius, "Parcel","$Temperature",sizeof(CoordReal));
            ucarg(ComputeRadius, "Parcel","$MoleFractionH2O",sizeof(CoordReal));
            ucarg(ComputeRadius, "Parcel","$AbsoluteTotalPressure",sizeof(CoordReal));
            ucarg(ComputeRadius, "Parcel", "$Time", sizeof(CoordReal));
            ucarg(ComputeRadius, "Parcel", "$TimeLevel", sizeof(CoordReal));

            ucfunc(RateCondensation, "ScalarFieldFunction", "RateCondensation");
            ucarg(RateCondensation, "Parcel", "$ParcelId", sizeof(CoordReal));
            ucarg(RateCondensation, "Parcel", "$ParticleResidenceTime", sizeof(CoordReal));
            ucarg(RateCondensation, "Parcel", "$Time", sizeof(CoordReal));
            ucarg(RateCondensation, "Parcel", "$TimeLevel", sizeof(CoordReal));
            ucarg(RateCondensation, "Parcel", "$UserComputeRadius", sizeof(CoordReal));
            ucarg(RateCondensation, "Parcel", "$Volume", sizeof(CoordReal));
	    			ucarg(RateCondensation, "Parcel","$Temperature",sizeof(CoordReal));
	    			ucarg(RateCondensation, "Parcel","$MoleFractionH2O",sizeof(CoordReal));
            ucarg(RateCondensation, "Parcel","$AbsoluteTotalPressure",sizeof(CoordReal));

            ucfunc(SumCondensation, "ScalarFieldFunction", "SumCondensation");
            ucarg(SumCondensation, "Parcel", "$TimeLevel", sizeof(CoordReal));
            ucarg(SumCondensation, "Parcel", "$UserRateCondensation", sizeof(CoordReal));

        }
