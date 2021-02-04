#include "uclib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <dirent.h>


void USERFUNCTION_EXPORT
  SumCondensation(CoordReal *result,int size,CoordReal *timeLevel,CoordReal *RateCondensation)
  {
    int debug_output =0; 							/* 1 == Mode deboggage */
    int i=0;

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

	for (i=0; i != size; i++)
		{

		if(timeLevel[0]>=2)
			{
			result[i]=-(RateCondensation[i]);
			}
		else
			{
			result[i]=0;
			}


		}
  }
