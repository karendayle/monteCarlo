/* This program is for section 5.2.1 of Jacques' 2011 book chapter, page 113.
   https://omlc.org/software/mc/Jacques2011_MonteCarlo_Welch&VanGemert.pdf   

  It simulates the size of the step a photon takes between a single interaction.

   Dayle Kotturi 					October 16, 2019
 */ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N_PHOTONS (1000) /* the more the better. MC works best with 100K events */
#define N_BINS  (50)      /* number of bins * bin width = range of x axis */
#define BIN_WIDTH (0.275) /* this is "ds" in cm. reduce this for better resolution */

int i;
int hist[N_BINS];
double N;
int bin;
int number_of_negatives = 0;
int number_of_overflows = 0;
double myRand;
double myScaledRand;
double s1;
double minS1 = 0.;
double maxS1 = 0.;

int main () 
{
	for (i = 1; i <= N_PHOTONS; i++)
	{
		myRand = rand() + 1.0; /* the +1 is the magic b/c of the log */
		myScaledRand = myRand / (RAND_MAX + 1.0);  /* the +1 is the magic b/c of the log */
		s1 = -1.0 * log(myScaledRand);
		if (s1 < minS1) minS1 = s1;
		if (s1 > maxS1) maxS1 = s1;
		bin = s1/BIN_WIDTH; /* losses here due to conv to int */
                if (bin >= N_BINS) { 
			printf("error bin=%d exceeds max=%d\n", bin, N_BINS);
			number_of_overflows++;
		}
		else {
                	if (bin < 0) {
				printf("error bin=%d less than zero\n", bin);
				number_of_negatives++;
			}
			else {
				hist[bin] += 1; /* add 1 to the histogram for this bin */ 	
			}
		}
		/*printf("myRand: %f, myScaledRand: %f, s1: %f, bin: %d\n", myRand, myScaledRand, s1, bin);*/
	}	
	
	printf("What is in the bins?\n");
	for (i=0; i<N_BINS; i++)
		printf("%6d %6d\n",i, hist[i]);
	printf("Width of each bin (cm): %f\n", BIN_WIDTH);
	printf("Number of times bin number was negative: %d\n", number_of_negatives);
	printf("Number of times bin number exceeded range : %d\n", number_of_overflows);
	printf("Min and max values of S1 modeled as -ln(RAND): %f %f\n", minS1, maxS1);
	return 0;
}
