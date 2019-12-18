#include <math.h>
#include <stdio.h>

#include "bz_discr.h"
#include "dem_discr.h"

//don't #define SLACK
//don't #define SPAM
#define WORK_OUTPUT

double globallower;

int cmpdbl(double *a, double *b)
{
  if ((*a) < (*b))
    return -1;
  else if ((*a) == (*b))
    return 0;
  return 1;
}

double oydiscr_cell(int npoints, int dim, int rempoints,
		    double **forced, int nforced,
		    double *lowerleft, double *upperright)
{
  double discr, maxdiscr, coordlist[dim][nforced];
  int indexes[dim];
  int i,j,k, status, dimension;
  double biggest[dim][nforced+1], smallest[dim][nforced+1];
  int maxpoints[dim], ntotal = rempoints + nforced;
  // biggest[i][j]: biggest product of coords 0--i for hitting j points
  // smallest[i][j]: smallest product of coords 0--i for hitting j+1 points
  // maxpoints[i]: number of points you get in total from coords 0--i
  double vol1=1.0, vol2=1.0;
  for (i=0; i<dim; i++) {
    vol1 *= lowerleft[i];
    vol2 *= upperright[i];
  }
#ifdef SPAM
  fprintf(stderr, "Parameters: npoints %d, dim %d, rempoints %d, nforced %d\n",
	  npoints, dim, rempoints, nforced);
  fprintf(stderr, "Lower: (%g", lowerleft[0]);
  for (i=1; i<dim; i++)
    fprintf(stderr, ", %g", lowerleft[i]);
  fprintf(stderr, ")\nUpper: (%g", upperright[0]);
  for (i=1; i<dim; i++)
    fprintf(stderr, ", %g", upperright[i]);
  fprintf(stderr, ")\nUncategorized ('forced') points are:\n");
  for (i=0; i<nforced; i++) {
    fprintf(stderr, "(%g", forced[i][0]);
    for (j=1; j<dim; j++)
      fprintf(stderr, ", %g", forced[i][j]);
    fprintf(stderr, ")\n");
  }
#endif

  maxdiscr = vol2 - (double)rempoints/npoints;
  discr = (double)(rempoints+nforced)/npoints - vol1;
  if (discr > maxdiscr)
    maxdiscr = discr;
  if (maxdiscr < globallower)
    return maxdiscr;
  // quicker code for use in some easy cells
  // otherwise, get to work...
  for (i=0; i<dim; i++) {
    indexes[i]=0;
    for (j=0; j<=nforced; j++) {
      smallest[i][j]=1.0;
      biggest[i][j]=0.0;
    }
  }
  for (i=0; i<nforced; i++) {
    status=0;
    for (j=0; j<dim; j++) {
      // order is chosen to handle final box
      if (forced[i][j] <= lowerleft[j])
	continue;
      else if (forced[i][j] >= upperright[j])
	break;
      else { // strictly internal
	if (status) {
	  fprintf(stderr, "PROBLEM: Point occurs as double-internal\n");
	  fflush(stderr);
	  abort();
	}
	status = 1;
	dimension=j;
      }
    }
    if (j==dim) { // else: hit "break", skip
      if (status) {
	coordlist[dimension][indexes[dimension]]=forced[i][dimension];
	indexes[dimension] += 1;
      }
      else { // completely internal
	rempoints++;
      }
    }
  }
  for (i=0; i<dim; i++)
    qsort(&(coordlist[i][0]), indexes[i], sizeof(double), cmpdbl);
  maxpoints[0]=indexes[0];
  for (i=1; i<dim; i++)
    maxpoints[i]=maxpoints[i-1]+indexes[i];

#ifdef SPAM
  fprintf(stderr, "Categorization: %d lower-left, %d internal.\n", rempoints, maxpoints[dim-1]);
  for (i=0; i<dim; i++) {
    if (!indexes[i]) {
      fprintf(stderr, "Direction %d: Nothing.\n", i);
      continue;
    }
    fprintf(stderr, "Direction %d: %g", i, coordlist[i][0]);
    for (j=1; j<indexes[i]; j++)
      fprintf(stderr, ", %g", coordlist[i][j]);
    fprintf(stderr, "\n");
  }
#endif

  // coord 0 first, since there is no recursion for that:
  smallest[0][0]=lowerleft[0];
  for (j=0; j<indexes[0]; j++) {
    smallest[0][j+1]=biggest[0][j]=coordlist[0][j];
  }
  biggest[0][indexes[0]]=upperright[0];
#ifdef SPAM
  fprintf(stderr, "Direction 0 only, biggest: ");
  for (j=0; j<=maxpoints[0]; j++)
    fprintf(stderr, "%g ", biggest[0][j]);
  fprintf(stderr, "\nDirections 0 only, smallest: ");
  for (j=0; j<=maxpoints[0]; j++)
    fprintf(stderr, "%g ", smallest[0][j]);
  fprintf(stderr, "\n");
#endif

  for (i=1; i<dim; i++) {
    // first the special loop for smallest: "nothing contributed by us"
    for (j=0; j<=maxpoints[i-1]; j++)
      smallest[i][j]=smallest[i-1][j]*lowerleft[i];
    // main loop:
    for (j=0; j<indexes[i]; j++) {
      vol1 = coordlist[i][j];
      for (k=0; k<=maxpoints[i-1]; k++) {
	// for biggest: vol1 is coordinate that adds j new points
	vol2=biggest[i-1][k]*vol1;
	if (vol2 > biggest[i][j+k])
	  biggest[i][j+k]=vol2;
	// for smallest: vol1 is coordinate that adds j+1 new points
	vol2=smallest[i-1][k]*vol1;
	if (vol2 < smallest[i][j+k+1])
	  smallest[i][j+k+1]=vol2;
      }
    }
    // last, special loop for biggest: "all of us"
    for (j=0; j<=maxpoints[i-1]; j++) {
      vol1=biggest[i-1][j]*upperright[i];
      if (vol1 > biggest[i][j+indexes[i]])
	biggest[i][j+indexes[i]]=vol1;
    }
#ifdef SPAM
    fprintf(stderr, "Directions 0--%d, biggest: ", i);
    for (j=0; j<=maxpoints[i]; j++)
      fprintf(stderr, "%g ", biggest[i][j]);
    fprintf(stderr, "\nDirections 0--%d, smallest: ", i);
    for (j=0; j<=maxpoints[i]; j++)
      fprintf(stderr, "%g ", smallest[i][j]);
    fprintf(stderr, "\n");
#endif
  }
  // now, use these to find lower, upper limits
  // mode: always contain "rempoints", additionally
  //         pick from smallest[dim-1], biggest[dim-1]
  maxdiscr=0;
  for (i=0; i<=maxpoints[dim-1]; i++) { // i counts internal points
    // small box
    discr = (double)(rempoints+i)/npoints - smallest[dim-1][i];
    if (discr > maxdiscr)
      maxdiscr=discr;
    // big box
    discr = biggest[dim-1][i] - (double)(rempoints+i)/npoints;
    if (discr > maxdiscr)
      maxdiscr = discr;
  }
  if (maxdiscr > globallower) {
#ifdef WORK_OUTPUT
    fprintf(stderr, "Worse bound: %g\n", maxdiscr);
#endif
    globallower=maxdiscr;
  }
#ifdef SPAM
  else
    fprintf(stderr, "Conclusion: %g\n", maxdiscr);
#endif
  return maxdiscr;
}

// "forced" points are points that are strictly between the boundaries in
// one of the cdim earlier dimensions; pointset[0--rempoints-1] are points that
// so far are at most at the lower-left corner in every dimension.
// this includes a final run with lower-left=1.
// being ON a border changes nothing:
//   ON lower-left counts as in (including when lowerleft=1)
//   ON upper-right counts as out (except if previous).
double oydiscr_int(double **pointset, int npoints, int dim, int rempoints,
		   double **forced, int nforced, int cdim,
		   double *lowerleft, double *upperright)
{
  double coord, forcedcoord, lowedge=0.0, highedge;
  int newcount=0, forcedidx, i, j, previdx=0, newforcedidx, resort=0, curridx;
  int newrempoints, wasforced, wasfinal=0;
  double maxdiscr=0.0, discr;
  double **newforced = malloc((nforced+rempoints)*sizeof(double *));
  // internal vars: previdx points at first element excluded from last pass
  //                    (on or above border coordinate)
  //                forcedidx: next unused forced element (a.k.a. counter)
  //                curridx: current pointset point ("i" as loop var)
  //  coords:
  //           coord is value of current pointset-point,
  //           forcedcoord is value of next-up forced boundary,
  //           lowedge is value of last boundary we used
  // newcount counts number of pointset-points since last coord
  if (cdim==dim) {
    free(newforced);
    return oydiscr_cell(npoints, dim, rempoints,
			forced, nforced,
			lowerleft, upperright);
  }
  comparedim=cdim;
  qsort(pointset, rempoints, sizeof(double *), cmpkeyk);
  qsort(forced, nforced, sizeof(double *), cmpkeyk);
  i=0; forcedidx=0;
  while ((i<rempoints) || (forcedidx < nforced)) {
    if (i<rempoints)
      coord=pointset[i][cdim];
    else
      coord=1.0;
    if (forcedidx < nforced)
      forcedcoord=forced[forcedidx][cdim];
    else
      forcedcoord=1.0;
    if ((forcedcoord > coord) && (newcount <= sqrt(npoints))) {
      newcount++;
      i++;
      if ((i<rempoints) || (forcedidx < nforced))
	continue;
      else { // Add one trailing cell
	lowerleft[cdim]=lowedge;
	highedge=upperright[cdim]=1.0;
	wasforced=0;
	wasfinal=1;
      }
    } // below: create new cell
    if (!wasfinal) {
      if (forcedcoord <= coord) {
	lowerleft[cdim]=lowedge;
	highedge=upperright[cdim]=forcedcoord;
	wasforced=1;
      }
      else { // must be count-based border
	lowerleft[cdim]=lowedge;
	highedge=upperright[cdim]=coord;
	wasforced=0;
      }
    } // end "if (!wasfinal)"
    curridx=i; // for better mnemonics
#ifdef WORK_OUTPUT
    if (!cdim)
      fprintf(stderr, "Coord %g\n", highedge);
#endif
    // creating a new cell (or subslab):
    // 1. surviving forced copied
    for (j=0; (j<nforced) && (forced[j][cdim] < highedge); j++)
      newforced[j]=forced[j];
    newforcedidx=j;
    // 2. new (strictly) internal points appended as forced
    j=previdx;
    while ((j<rempoints) && (pointset[j][cdim] <= lowedge))
      j++;
    newrempoints=j;
    for (; (j<rempoints) && (pointset[j][cdim] < highedge); j++)
      newforced[newforcedidx++] = pointset[j];
    if (j>(curridx+1))
      resort=1;
    // 3. make call with properly adjusted boundaries, update variables
    discr = oydiscr_int(pointset, npoints, dim, newrempoints,
			newforced, newforcedidx, cdim+1,
			lowerleft, upperright);
    if (discr > maxdiscr)
      maxdiscr = discr;
    if (resort) {
      comparedim=cdim;
      qsort(pointset, rempoints, sizeof(double *), cmpkeyk);
      resort=0;
    }
    while ((forcedidx < nforced) && (forced[forcedidx][cdim]==highedge))
      forcedidx++;
    while ((i < rempoints) && (pointset[i][cdim]<=highedge))
      i++;
    lowedge=highedge;
    previdx=i;
    newcount=0;
  }
  // one final call to capture the border cases (for boxes containing a point with coord 1.0)
  // 1. new forced == old forced (copy to avoid interfering with previous stages)
  for (j=0; j<nforced; j++)
    newforced[j]=forced[j];
  // 2. per above, we have no new internal/forced points
  // 3. make the call
  lowerleft[cdim]=lowedge;
  upperright[cdim]=1.0;
  discr = oydiscr_int(pointset, npoints, dim, rempoints,
		      newforced, nforced, cdim+1,
		      lowerleft, upperright);
  if (discr>maxdiscr)
    maxdiscr=discr;
  free(newforced);
  return maxdiscr;
}

double oydiscr(double **pointset, int dim, int npoints, double *lower)
{
  // JAKOB: otherwise in subsequent applications of the function
  // in R we get always the same output.
  globallower = 0;

  double lowerleft[dim], upperright[dim];
  double **pre_force = malloc(2*dim*sizeof(double *));
  double discr, *border;
  double maxcoord;
  int maxpos;
  int is_border[npoints];
  double **clone_set = malloc(npoints*sizeof(double *));
  int i,j,k;
  for (i=0; i<dim; i++) {
    border = malloc(dim*sizeof(double));
    for (j=0; j<dim; j++) {
      if (i==j)
	border[j]=1.0;
      else
	border[j]=0.0;
    }
    pre_force[i]=border;
  }
  for (i=0; i<npoints; i++)
    is_border[i]=0;
  for (i=0; i<dim; i++) {
    maxcoord=-1.0;
    maxpos=-1;
    for (j=0; j<npoints; j++)
      if (pointset[j][i] > maxcoord) {
	maxcoord = pointset[j][i];
	maxpos=j;
      }
    is_border[maxpos]=1;
  }
  j=dim; k=0;
  for (i=0; i<npoints; i++)
    if (is_border[i])
      pre_force[j++]=pointset[i];
    else
      clone_set[k++]=pointset[i];
//  discr = oydiscr_int(pointset, npoints, dim, npoints,
//		      pre_force, dim, 0,
//		      lowerleft, upperright);
//  discr = oydiscr_int(clone_set, npoints, dim, k,
//		      pre_force, j, 0,
//		      lowerleft, upperright);
//  discr = oydiscr_int(clone_set, npoints, dim, k,
//		      pre_force[dim], j-dim, 0,
//		      lowerleft, upperright);
  // final version: NOTHING pre-determined.
  discr = oydiscr_int(pointset, npoints, dim, npoints,
		      pre_force, 0, 0,
		      lowerleft, upperright);
  for (i=0; i<dim; i++)
    free(pre_force[i]);
  free(pre_force);
  free(clone_set);
  *lower = globallower;
  return discr;
}
