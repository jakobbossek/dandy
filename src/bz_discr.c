#include "bz_discr.h"
#include <stdio.h>

#define INT_EPS 1e-12
// define SPAM to get trace output

#ifdef SPAM
#define SPIM
#endif
// SPIM is more benign

int comparedim;

double globaldiscr;

double fmax(double a,double b)
{
  return ((a)>(b))?(a):(b);
}

// argument is pointer to pointsetmember, i.e. double **.
int cmpkeyk(double **pt1, double **pt2)
{
  double a=(*pt1)[comparedim], b=(*pt2)[comparedim];
  if (a<b)
    return -1;
  else if (a>b)
    return 1;
  return 0;
}

int intpoints(double **pointset, int dim, int npoints, double *base)
{
  int n=npoints,i,j;
  for (i=0; i<npoints; i++) {
    for (j=0; j<dim; j++) {
      if (pointset[i][j] >= (base[j]-INT_EPS)) {
	n--;
	break;
      }
    }
  }
  return n;
}

// Fixes one dimension at a time, each dimension defined by a point (which
// must be on the boundary of the box).  Maintains list of points that are 
// still possible (first rempoints points) and smallest base compatible with 
// earlier dimension choices (boundary points must not be excluded). 
double int_poly_discr(double **pointset, int dim, int npoints, int rempoints, 
		      double *base, int cdim)
{
  double discr,maxdiscr=0.0, basecopy[dim], *ptcopy[rempoints];
  int i,j, resort=0;

#ifdef SPAM
  fprintf(stderr, "Dim %d points %d/%d base (%g", cdim, rempoints, npoints, base[0]);
  for (i=1; i<dim; i++)
    fprintf(stderr, ", %g", base[i]);
  fprintf(stderr, ")\n");
#endif
  if (cdim==dim) {
    discr = 1.0;
    for (i=0; i<dim; i++)
      discr *= base[i];
    maxdiscr = fmax((double)rempoints/npoints-discr,
		    discr-(double)intpoints(pointset, dim, rempoints, base)/npoints);
#ifdef SPAM
    fprintf(stderr, "Volume %g point-share %g--%g -> discr. %g\n",
	    discr, (double)intpoints(pointset, dim, rempoints, base)/npoints,
	    (double)rempoints/npoints, maxdiscr);
#endif
    if (maxdiscr > globaldiscr) {
      globaldiscr=maxdiscr;
#ifdef SPIM
      fprintf(stderr, "Improved to %g\n", globaldiscr);
#endif
    }
    return maxdiscr;
  }

  for (i=cdim; i<dim; i++)
    basecopy[i]=base[i];
  comparedim = cdim;
  qsort(pointset, rempoints, sizeof(double *), cmpkeyk);    
#ifdef SPAM
  for (i=0; i<rempoints; i++) {
    fprintf(stderr, "Point %04d: (%g", i, pointset[i][0]);
    for (j=1; j<dim; j++) 
      fprintf(stderr, ", %g", pointset[i][j]);
    fprintf(stderr, ")\n");
  }
#endif

  for (i=0; i<rempoints; i++) {
    if (!cdim)
      fprintf(stderr, "%d/%d\n", i, rempoints);
    if (pointset[i][cdim] < base[cdim]-INT_EPS) {
#ifdef SPAM
      fprintf(stderr, "Skip'ng %d: coord too small\n", i);
#endif      
      continue;
    }
    
    base[cdim]=pointset[i][cdim];
    for (j=cdim+1; j<dim; j++) {
      if (pointset[i][j] > base[j])
	base[j]=pointset[i][j];
    }
    if ((i && (pointset[i-1][cdim] == pointset[i][cdim])) ||
	(((i+1)<rempoints) &&
	 (pointset[i+1][cdim] == pointset[i][cdim]))) {
      resort=1;
      for (j=0; j<rempoints; j++)
	ptcopy[j]=pointset[j];
    } 
    j=i+1;
    while ((j < rempoints) && (pointset[j-1][cdim]==pointset[j][cdim]))
      j++;
    
#ifdef SPAM
    fprintf(stderr, "Calling %d\n", i);
#endif
    discr = int_poly_discr(pointset, dim, npoints, j, base, cdim+1);
    if (discr > maxdiscr)
      maxdiscr=discr;
    if (resort) {
      for (j=0; j<rempoints; j++)
	pointset[j]=ptcopy[j];
      resort=0;
    }
    for (j=cdim+1; j<dim; j++)
      base[j]=basecopy[j];
  }
#ifdef SPAM
    fprintf(stderr, "Calling 1.0\n");
#endif
    if (!cdim)
      fprintf(stderr, "Trying 1.0\n");
  base[cdim]=1.0;
  discr = int_poly_discr(pointset, dim, npoints, rempoints, base, cdim+1);
  for (j=cdim; j<dim; j++)
    base[j]=basecopy[j];
  if (discr > maxdiscr)
    maxdiscr=discr;
  
  return maxdiscr;
}  

double poly_discr(double **pointset, int dim, int npoints)
{
  double *base, discr;
  int i;
  globaldiscr=0;
  base = malloc(dim*sizeof(double));
  for (i=0; i<dim; i++)
    base[i]=0.0;
  discr = int_poly_discr(pointset, dim, npoints, npoints, base, 0);
  free(base);
  return discr;
}
