/*
   reads a point set and computes a lower bound on its star discrepancy
*/

/*
  i_tilde_bardelta inner and outer loops (effect of alpha is considered to be negligible)
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

int k_div_bardelta=0;// "0" means default "4 or 8" setup, other value (from main(), e.g. -k 16) overrides

#define I_TILDE 316    // thresholds to be calculated (sqrt(iterations)), default value 100k
int mc_bardelta=MC;             //nbr of coordinates to be changed, default value
int i_tilde_bardelta=I_TILDE;
#define TRIALS 10      //nbr of runs (mean and max will be calculated), default value
int trials_bardelta=TRIALS;

#ifndef VERSION
// JAKOB: as in Makefile
#define VERSION 2
#endif

// global variables to store info about pointset


// global variables to store info about worst box (also "private")
double real_max_discr_bardelta=0;
int real_when_bardelta=0, when_bardelta=0;
int current_iteration;


//Computes the best of the rounded points -- basic version
//Constant neighbourhood size and mc_bardelta-values.
//Does not split the search.
//Copies the appropriate "thing" into xc_index (output variable)
double best_of_rounded_bardelta(int *xn_minus, int *xn_extraminus, int *xc_index)
{
   double fxn_minus;
  double fxn_extraminus;
  double fxc;
  int j, d=n_dimensions;
  int use_extraminus=0;
#if (VERSION == 1) || (VERSION == 2)
  int xn_minus_snap[d], xn_extraminus_snap[d];
#endif
  for (j=0; j<d; j++)
    if (xn_minus[j] != xn_extraminus[j]) {
      use_extraminus=1;
      break;
    }

  // Growing, shrinking.
#if (VERSION == 3)
  // SNAPUPDATE
  snap_box(xn_minus);
  if (use_extraminus)
    snap_box(xn_extraminus);
#elif (VERSION == 1) || (VERSION == 2)
  // Grower, shrinker that copy the point
  for (j=0; j<d; j++)
    xn_minus_snap[j]=xn_minus[j];
  snap_box(xn_minus_snap);
  if (use_extraminus) {
    for (j=0; j<d; j++)
      xn_extraminus_snap[j]=xn_extraminus[j];
    snap_box(xn_extraminus_snap);
  }
#endif

  // For version 1: "private update"
#if (VERSION == 1)
  fxc = get_bar_delta(xn_minus_snap);
  if (use_extraminus) {
    fxn_extraminus=get_bar_delta(xn_extraminus_snap);
    fxc=max(fxc, fxn_extraminus);
  }
#ifdef PRINT_UPDATE_CANDIDATES
  fprintf(stderr, "PRIVATE candidate %g (vs %g)\n",
	  fxc, real_max_discr_bardelta);
#endif
  if (fxc > real_max_discr_bardelta) {
    real_max_discr_bardelta = fxc;
    real_when_bardelta = current_iteration;
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
  fxc = get_bar_delta(xn_minus_snap);
  if (use_extraminus) {
    fxn_extraminus=get_bar_delta(xn_extraminus_snap);
    fxc=max(fxc, fxn_extraminus);
  }
#else
  // versions 0,3 both compute from the point now given by xn_*
  // version 1 officially reports this as well
  fxc = get_bar_delta(xn_minus);
  if (use_extraminus) {
    fxn_extraminus=get_bar_delta(xn_extraminus);
    fxc=max(fxc, fxn_extraminus);
  }
#endif

  // Remains only to copy the winning point to output variable xc_index.
  if (use_extraminus && (fxn_extraminus >= fxc)) {
    for (j=0; j<d; j++)
      xc_index[j] = xn_extraminus[j];
  }
  else {
    for(j=0; j<d; j++)
      xc_index[j]= xn_minus[j];
  }

  return fxc;
}



double bardelta_calc(double **pointset, int n, int d, int iter, int max_trials)
{
  trials_bardelta = max_trials;
  i_tilde_bardelta = (int)sqrt(iter);

  int k[d], start[d];

  int i, j, p, t;           // loop variables

  double thresh[i_tilde_bardelta];    //Thresholdsequence
  double T;                            //current Threshold

  double fxc;
  int xc_index[d], xn_minus_index[d], xn_extraminus_index[d];
  int xn_best_index[d];    //Indices of current point, neighbour
  double xglobal[trials_bardelta+1][d], xbest[d];
  double current, global[trials_bardelta+1], best, mean;  //current and global best values

  int outerloop=i_tilde_bardelta, innerloop=i_tilde_bardelta;

  int anzahl=0;
  int switches[trials_bardelta+1];
  int global_switches[trials_bardelta+1];

  //Get pointset from external file
  FILE *datei_ptr=stderr;
  fprintf(datei_ptr,"GLP-Menge %d %d  ",d,n);

  //Sort the grid points, setup global variables
  process_coord_data(pointset, n, d);

  //Algorithm starts here
  for(t=1;t<=trials_bardelta;t++)
    { //Initialization
      fprintf(stderr, "Trial %d/%d\n", t, trials_bardelta);

      //Initialize k-value
      for (j=0; j<d; j++) {
	start[j]=(int)((n_coords[j]-1)/2);
      }
      //Initialize mc_bardelta-value
      mc_bardelta=2;

      //Initialize iteration count
      current_iteration=0;

      //Generate threshold sequence   (only once)
      //      fprintf(stderr, "Generating threshold\n");
      for(i=1;i<=outerloop;i++){

	current_iteration++;
	//Update k-value
	  for (j=0; j<d; j++) {
	    k[j] = start[j]*(((double)outerloop-current_iteration)/(outerloop)) +
	      1*((double)current_iteration/(outerloop));
		  //	    k[j]=start[j] - (int)((3.0/4)*(current_iteration/outerloop)*(start[j]-1));
	  }

        //Update mc_bardelta-value
	  mc_bardelta=2+(int)(current_iteration/outerloop*(d-2));


	//generation of random point xc
	generate_xc_bardelta(xn_minus_index, xn_extraminus_index);

	//(Possibly) Snap the points and compute the largest of the rounded values
	current = best_of_rounded_bardelta(xn_minus_index, xn_extraminus_index, xc_index);

	//draw a neighbour of xc
	generate_neighbor_bardelta(xn_minus_index, xn_extraminus_index, xc_index, k, mc_bardelta);

	//Compute the threshold
	fxc=best_of_rounded_bardelta(xn_minus_index, xn_extraminus_index, xc_index);
	thresh[i]=0.0-fabs(fxc-current);
      }

      //sort the thresholds in increasing order
      quicksort(1,outerloop,thresh);


      switches[t]=0;
      global_switches[t]=0;
      current=0;
      global[t]=0;
      when_bardelta=0;
      real_when_bardelta=0;
      real_max_discr_bardelta=0;

      //Initialize k-value
      for (j=0; j<d; j++) {
	start[j]=(int)((n_coords[j]-1)/2);
      }
      //Initialize mc_bardelta-value
      mc_bardelta=2+(int)(current_iteration/(innerloop*outerloop)*(d-2));


      //draw a random initial point
      generate_xc_bardelta(xn_minus_index, xn_extraminus_index);

      //(Possibly) Snap and compute the best of the rounded points and update current value
      current = best_of_rounded_bardelta(xn_minus_index, xn_extraminus_index, xc_index);

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

	      //Update mc_bardelta-value
	      mc_bardelta=2+(int)(current_iteration/(innerloop*outerloop)*(d-2));
#ifdef PRINT_RANGE_DATA
	      if (p==1)
		fprintf(stderr, " threshold %g mc_bardelta %d\n", T, mc_bardelta);
#endif
	      //mc_bardelta=2;

	      //Get random neighbor
	      generate_neighbor_bardelta(xn_minus_index, xn_extraminus_index, xc_index,k,mc_bardelta);
#ifdef DISPLAY_CANDIDATES
	      fprintf(stderr, "Old: ");
	      for (j=0; j<d; j++)
		fprintf(stderr, "%d ", xc_index[j]);
	      fprintf(stderr, "\nMinus: ");
	      for (j=0; j<d; j++)
		fprintf(stderr, "%d ", xn_minus_index[j]);
	      fprintf(stderr, "\nXMinus: ");
	      for (j=0; j<d; j++)
		fprintf(stderr, "%d ", xn_extraminus_index[j]);
	      fprintf(stderr, "\n");
#endif

	      //(Possibly) Snap the points and compute the best of the rounded points
	      fxc = best_of_rounded_bardelta(xn_minus_index, xn_extraminus_index, xn_best_index);
#ifdef PRINT_UPDATE_CANDIDATES
	      fprintf(stderr, "Iter. %d candidate %10g (vs %10g best %10g) -- ",
		      current_iteration, fxc, current, global[t]);
#endif
	      //Global update if necessary
	      if(fxc>global[t]){
		global_switches[t]++;
		global[t]=fxc;
		when_bardelta=current_iteration;
#ifdef PRINT_UPDATE_CANDIDATES
		fprintf(stderr, "global ");
#endif
#ifdef PRINT_ALL_UPDATES
		fprintf(stderr, "%g at %d :", fxc, current_iteration);
		for (j=0; j<d; j++)
		  fprintf(stderr, " %d", xn_best_index[j]);
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
		  xc_index[j]=xn_best_index[j];
		}
	      }
#ifdef PRINT_UPDATE_CANDIDATES
	      else {
		fprintf(stderr, "skip\n");
	      }
#endif
	    }//innerloop
	}//outerloop
      if (real_max_discr_bardelta > global[t]) {
	global[t] = real_max_discr_bardelta;
	when_bardelta = real_when_bardelta;
	//	fprintf(stderr, "Max value subsumed\n");
      }
      fprintf(stderr, "Result %g at %d\n", global[t], when_bardelta);
      fprintf(stdout, "%g\n", global[t]); // To simplify post-execution bookkeeping
    }//trials_bardelta


  //best calculated value
  best=global[1];
  for(j=0; j<d; j++) xbest[j]=xglobal[1][j];
  for(t=2;t<=trials_bardelta;t++)
    {
      if(global[t]>best)
	{
	  best=global[t];
	  for(j=0; j<d; j++) xbest[j]=xglobal[t][j];
	}
    }

  for(t=1;t<=trials_bardelta;t++)
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
  for(t=1;t<=trials_bardelta;t++) mean=mean+global[t];
  mean=mean/trials_bardelta;
  fprintf(datei_ptr,"mean %e  ",mean);
  //fprintf(datei_ptr,"lower_bound %e\n",lower_bound);
  //fprintf(datei_ptr,"upper_bound %e\n",upper_bound);

  //  fprintf(datei_ptr,"Anzahl der Iterationen: %d  ",iteration_count);
  //fprintf(datei_ptr,"Wert von k: %d\n",k);
  // fprintf(datei_ptr,"Wert von Extraminus: %d\n",extraminus);
  fprintf(datei_ptr,"Anzahl best: %d\n",anzahl);
  // for(i=1;i<=outerloop;i++) fprintf(datei_ptr,"Thresh %d = %e\n",i,thresh[i]);

  // for(t=1;t<=trials_bardelta;t++) {
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
//       k_div_bardelta = atoi(argv[++pos]);
//       pos++;
//       fprintf(stderr, "Using k = n/%d\n", k_div_bardelta);
//     }
//     else if (!strcmp(argv[pos], "-mc_bardelta")) {
//       mc_bardelta = atoi(argv[++pos]);
//       pos++;
//       fprintf(stderr, "Using mc_bardelta = %d\n", mc_bardelta);
//     }
//     else if (!strcmp(argv[pos], "-iter")) {
//       i_tilde_bardelta=(int)sqrt(atoi(argv[++pos]));
//       pos++;
//       fprintf(stderr, "Using %d iterations (adj. for sqrt)\n",
// 	      i_tilde_bardelta*i_tilde_bardelta);
//     }
//     else if (!strcmp(argv[pos], "-trials_bardelta")) {
//       trials_bardelta = atoi(argv[++pos]);
//       pos++;
//       fprintf(stderr, "Doing %d independent trials_bardelta (currently: times ten thresh. rep.)\n",
// 	      trials_bardelta);
//       trials_bardelta*=THRESH_REPEAT;

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
//   if (dim<mc_bardelta)
//     mc_bardelta=dim;
//   fprintf(stderr, "Calling Carola calculation\n");
//   printf("%g\n", oldmain(pointset, npoints, dim));
//   return EXIT_SUCCESS;

// }

