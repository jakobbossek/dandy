/*
   reads a point set and computes a lower bound on its star discrepancy
*/

/*
  i_tilde inner and outer loops (effect of alpha is considered to be negligible)
  threshold aus Paaren von Nachbarn berechnet
*/

// This time, the file is written so that first element is x[0] :-p

// compute threshold sequence once, then reuse it
// not using alpha

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include "TA_common.h"

#ifndef MC
#define MC 2
#endif

// not recommended: PRINT_ALL_UPDATES
#define PRINT_ALL_UPDATES
// even more ridiculous!
//#define PRINT_UPDATE_CANDIDATES
//#define DISPLAY_CANDIDATES
//#define PRINT_RANGE_DATA

// use THRESH_REPEAT=1 in "production code" (it only helps to denoise testing)
#define THRESH_REPEAT 1

int k_div=0;// "0" means default "4 or 8" setup, other value (from main(), e.g. -k 16) overrides

#define I_TILDE 316    // thresholds to be calculated (sqrt(iterations)), default value 100k
int mc=MC;             //nbr of coordinates to be changed, default value
int i_tilde=I_TILDE;
#define TRIALS 10      //nbr of runs (mean and max will be calculated), default value
int trials=TRIALS;

#ifndef VERSION
// JAKOB: as in Makefile
#define VERSION 2
#endif

// global variables to store info about pointset


// global variables to store info about worst box (also "private")
double real_max_discr=0;
int real_when=0, when=0;
int current_iteration;


//Computes the best of the rounded points -- basic version
//Constant neighbourhood size and mc-values.
//Does not split the search.
//Copies the appropriate "thing" into xc_index (output variable)
double best_of_rounded_delta(int *xn_plus)
{
  double fxc;
  int j, d=n_dimensions;
#if (VERSION == 1) || (VERSION == 2)
  int xn_plus_grow[d];
#endif

  // Growing, shrinking.
#if (VERSION == 3)
  // SNAPUPDATE
  grow_box_randomly(xn_plus);
#elif (VERSION == 1) || (VERSION == 2)
  // Grower, shrinker that copy the point
  for (j=0; j<d; j++)
    xn_plus_grow[j]=xn_plus[j];
  grow_box_randomly(xn_plus_grow);
#endif

  // For version 1: "private update"
#if (VERSION == 1)
  fxc = get_delta(xn_plus_grow);
#ifdef PRINT_UPDATE_CANDIDATES
  fprintf(stderr, "PRIVATE candidate %g (vs %g)\n",
	  fxc, real_max_discr);
#endif
  if (fxc > real_max_discr) {
    real_max_discr = fxc;
    real_when = current_iteration;
#ifdef PRINT_ALL_UPDATES
    fprintf(stderr, "Secret update at %d to %g\n",
	    current_iteration, fxc);
#endif
  }
#endif
  // Ends "private update" block

  // Now, create the official numbers.
#if (VERSION == 2)
  // official update from modified points
  fxc = get_delta(xn_plus_grow);
#else
  // versions 0,3 both compute from the point now given by xn_*
  // version 1 officially reports this as well
  fxc = get_delta(xn_plus);
#endif

  // In delta_only mode: don't copy

  return fxc;
}



double delta_calc(double **pointset, int n, int d, int iter, int max_trials)
{
  trials = max_trials;
  i_tilde = (int)sqrt(iter);

  int k[d], start[d];

  int i, j, p, t;           // loop variables

  double thresh[i_tilde];    //Thresholdsequence
  double T;                            //current Threshold

  double fxc;
  int xc_index[d], xn_plus_index[d];
  int xn_best_index[d];    //Indices of current point, neighbour
  double xglobal[trials+1][d], xbest[d];
  double current, global[trials+1], best, mean;  //current and global best values

  int outerloop=i_tilde, innerloop=i_tilde;

  int anzahl=0;
  int switches[trials+1];
  int global_switches[trials+1];

  //Get pointset from external file
  FILE *datei_ptr=stderr;
  fprintf(datei_ptr,"GLP-Menge %d %d  ",d,n);

  //Sort the grid points, setup global variables
  process_coord_data(pointset, n, d);

  //Algorithm starts here
  for(t=1;t<=trials;t++)
    { //Initialization
      fprintf(stderr, "Trial %d/%d\n", t, trials);

      //Initialize k-value
      for (j=0; j<d; j++) {
	start[j]=(int)((n_coords[j]-1)/2);
      }
      //Initialize mc-value
      mc=2;

      //Initialize iteration count
      current_iteration=0;

      //Generate threshold sequence
      for(i=1;i<=outerloop;i++){

	current_iteration++;
	//Update k-value
	  for (j=0; j<d; j++) {
	    k[j] = start[j]*(((double)outerloop-current_iteration)/(outerloop)) +
	      1*((double)current_iteration/(outerloop));
		  //	    k[j]=start[j] - (int)((3.0/4)*(current_iteration/outerloop)*(start[j]-1));
	  }

        //Update mc-value
	  mc=2+(int)(current_iteration/outerloop*(d-2));


	//generation of random point xc
	generate_xc_delta(xc_index);

	//(Possibly) Snaps the point upwards and computes the fitness
	current = best_of_rounded_delta(xc_index);

	//draw a neighbour of xc
	generate_neighbor_delta(xn_plus_index, xc_index, k, mc);

	//Compute the threshold
	fxc=best_of_rounded_delta(xn_plus_index);
	thresh[i]=0.0-fabs(fxc-current);
      }

      //sort the thresholds in increasing order
      quicksort(1,outerloop,thresh);


      switches[t]=0;
      global_switches[t]=0;
      current=0;
      global[t]=0;
      when=0;
      real_when=0;
      real_max_discr=0;

      //Initialize k-value
      for (j=0; j<d; j++) {
	start[j]=(int)((n_coords[j]-1)/2);
      }
      //Initialize mc-value
      mc=2+(int)(current_iteration/(innerloop*outerloop)*(d-2));


      //draw a random initial point
      generate_xc_delta(xc_index);

      //(Possibly) Snap and compute the best of the rounded points and update current value
      current = best_of_rounded_delta(xc_index);

      global[t] = current;

      current_iteration=0;
      for(i=1;i<=outerloop;i++)
	{
	  T=thresh[i];

	  for(p=1;p<=innerloop;p++)
	    {
	      current_iteration++;

	      //Update k-value
#ifdef PRINT_RANGE_DATA
	      if (p==1)
		fprintf(stderr, "Snapshot: range ");
#endif
	      for (j=0; j<d; j++) {
		k[j] = start[j]*(((double)innerloop*outerloop-current_iteration)/(innerloop*outerloop)) +
		  1*((double)current_iteration/(innerloop*outerloop));
		  //		k[j]=(int)(start[j]-(int)(current_iteration/(innerloop*outerloop)*(start[j]-1)));
#ifdef PRINT_RANGE_DATA
		if (p==1)
		  fprintf(stderr, "%d ", k[j]);
#endif
	      }

	      //Update mc-value
	      mc=2+(int)(current_iteration/(innerloop*outerloop)*(d-2));
#ifdef PRINT_RANGE_DATA
	      if (p==1)
		fprintf(stderr, " threshold %g mc %d\n", T, mc);
#endif
	      //mc=2;

	      //Get random neighbor
	      generate_neighbor_delta(xn_plus_index, xc_index,k,mc);
#ifdef DISPLAY_CANDIDATES
	      fprintf(stderr, "Old: ");
	      for (j=0; j<d; j++)
		fprintf(stderr, "%d ", xc_index[j]);
	      fprintf(stderr, "\nPlus: ");
	      for (j=0; j<d; j++)
		fprintf(stderr, "%d ", xn_plus_index[j]);
	      fprintf(stderr, "\n");
#endif

	      //(Possibly) Snap the points and compute the best of the rounded points
	      fxc = best_of_rounded_delta(xn_plus_index);
#ifdef PRINT_UPDATE_CANDIDATES
	      fprintf(stderr, "Iter. %d candidate %10g (vs %10g best %10g) -- ",
		      current_iteration, fxc, current, global[t]);
#endif
	      //Global update if necessary
	      if(fxc>global[t]){
		global_switches[t]++;
		global[t]=fxc;
		when=current_iteration;
#ifdef PRINT_UPDATE_CANDIDATES
		fprintf(stderr, "global ");
#endif
#ifdef PRINT_ALL_UPDATES
		fprintf(stderr, "%g at %d :", fxc, current_iteration);
		for (j=0; j<d; j++)
		  fprintf(stderr, " %d", xn_plus_index[j]);
		fprintf(stderr, "\n");
#endif
	      }
	      //Update of current best value if necessary
	      if(fxc-current>=T){
#ifdef PRINT_UPDATE_CANDIDATES
		fprintf(stderr, "update\n");
#endif
		switches[t]++;
		current=fxc;
		for(j=0; j<d; j++){
		  xc_index[j]=xn_plus_index[j];
		}
	      }
#ifdef PRINT_UPDATE_CANDIDATES
	      else {
		fprintf(stderr, "skip\n");
	      }
#endif
	    }//innerloop
	}//outerloop
      if (real_max_discr > global[t]) {
	global[t] = real_max_discr;
	when = real_when;
	//	fprintf(stderr, "Max value subsumed\n");
      }
      fprintf(stderr, "Result %g at %d\n", global[t], when);
      fprintf(stdout, "%g\n", global[t]); // To simplify post-execution bookkeeping
    }//trials


  //best calculated value
  best=global[1];
  for(j=0; j<d; j++) xbest[j]=xglobal[1][j];
  for(t=2;t<=trials;t++)
    {
      if(global[t]>best)
	{
	  best=global[t];
	  for(j=0; j<d; j++) xbest[j]=xglobal[t][j];
	}
    }

  for(t=1;t<=trials;t++)
    {
      if(global[t]==best) anzahl++;
    }
  fprintf(datei_ptr,"best %e  ",best);
  // for(j=0; j<d; j++)  fprintf(datei_ptr,"xbest %d coo  %e\n", j,xbest[j]);

  //delta or bar(delta) causing best value?
  //if(best==fabs(delta(xbest,GLP))) fprintf(datei_ptr,"delta\n");
  //else fprintf(datei_ptr,"bar_delta\n");

  //calculation of mean value
  mean=0;
  for(t=1;t<=trials;t++) mean=mean+global[t];
  mean=mean/trials;
  fprintf(datei_ptr,"mean %e  ",mean);
  //fprintf(datei_ptr,"lower_bound %e\n",lower_bound);
  //fprintf(datei_ptr,"upper_bound %e\n",upper_bound);

  //  fprintf(datei_ptr,"Anzahl der Iterationen: %d  ",iteration_count);
  //fprintf(datei_ptr,"Wert von k: %d\n",k);
  // fprintf(datei_ptr,"Wert von Extraminus: %d\n",extraminus);
  fprintf(datei_ptr,"Anzahl best: %d\n",anzahl);
  // for(i=1;i<=outerloop;i++) fprintf(datei_ptr,"Thresh %d = %e\n",i,thresh[i]);

  // for(t=1;t<=trials;t++) {
  //fprintf(datei_ptr,"Anzahl switches in Runde %d: %d\n",t,switches[t]);
  //fprintf(datei_ptr,"Anzahl global_switches in Runde %d: %d\n",t,global_switches[t]);
  //}

  return best;

}


// int main(int argc, char **argv)
// {
//   int dim, npoints,i,j;
//   FILE *pointfile;
//   double **pointset;
//   int pos=1;

//   FILE *random;
//   unsigned int seed;
//   random = fopen("/dev/random", "rb");
//   fread(&seed, 4, 1, random);
//   srand(seed);
//   while (pos < argc) {
//     if (!strcmp(argv[pos], "-kdiv")) {
//       k_div = atoi(argv[++pos]);
//       pos++;
//       fprintf(stderr, "Using k = n/%d\n", k_div);
//     }
//     else if (!strcmp(argv[pos], "-mc")) {
//       mc = atoi(argv[++pos]);
//       pos++;
//       fprintf(stderr, "Using mc = %d\n", mc);
//     }
//     else if (!strcmp(argv[pos], "-iter")) {
//       i_tilde=(int)sqrt(atoi(argv[++pos]));
//       pos++;
//       fprintf(stderr, "Using %d iterations (adj. for sqrt)\n",
// 	      i_tilde*i_tilde);
//     }
//     else if (!strcmp(argv[pos], "-trials")) {
//       trials = atoi(argv[++pos]);
//       pos++;
//       fprintf(stderr, "Doing %d independent trials (currently: times ten thresh. rep.)\n",
// 	      trials);
//       trials*=THRESH_REPEAT;

//     }
//     else
//       break;
//   }
//   switch (argc-pos) {
//   case 0:
//     i=scanf("%d %d reals\n", &dim, &npoints);
//     if (i != 2) {
//       fprintf(stderr, "stdin mode and header line not present\n");
//       exit(EXIT_FAILURE);
//     }
//     pointfile=stdin;
//     break;

//   case 1: // one arg, interpret as file name
//     pointfile = fopen(argv[pos], "r");
//     i=fscanf(pointfile, "%d %d reals\n", &dim, &npoints);
//     if (i != 2) {
//       fprintf(stderr, "stdin mode and header line not present\n");
//       exit(EXIT_FAILURE);
//     }
//     break;

//   case 2: // interpret as dim npoints args
//     dim=atoi(argv[pos++]);
//     npoints=atoi(argv[pos]);
//     pointfile=stdin;
//     break;

//   case 3: // interpret as dim npoints file; file not allowed to have header
//     dim=atoi(argv[pos++]);
//     npoints=atoi(argv[pos++]);
//     pointfile = fopen(argv[pos], "r");
//     break;

//   default:
//     fprintf(stderr, "Usage: calc_discr [dim npoints] [file]\n\nIf file not present, read from stdin. If dim, npoints not present, \nassume header '%%dim %%npoints reals' (e.g. '2 100 reals') in file.\n");
//     exit(EXIT_FAILURE);
//   }

//   fprintf(stderr, "Reading dim %d npoints %d\n", dim, npoints);
//   pointset = malloc(npoints*sizeof(double*));
//   for (i=0; i<npoints; i++) {
//     pointset[i] = malloc(dim*sizeof(double));
//     for (j=0; j<dim; j++) {
//       fscanf(pointfile, "%lg ", &(pointset[i][j]));
//       // newline counts as whitespace
//     }
//   }
//   if (dim<mc)
//     mc=dim;
//   fprintf(stderr, "Calling Carola calculation\n");
//   printf("%g\n", oldmain(pointset, npoints, dim));
//   return EXIT_SUCCESS;

// }

