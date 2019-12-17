#include "simple_discr.h"

// change the tuple; return bool: "success"
int step_tuple(int *tuple, int dim, int max)
{
  int i=0;
  while ((i<dim) && (tuple[i] == max-1)) {
    tuple[i++]=0;
  }
  if (i<dim) {
    tuple[i]++;
    return 1;
  }
  else {
    return 0;
  }
}

// pos1 is coord-wise <= pos2
int point_in_box(double *pos1, double *pos2, int dim)
{
  int i;

  for (i=0; i<dim; i++)
    if (pos1[i] > pos2[i])
      return 0;
  return 1;
}

// assuming point in box, is point also on edge?
int point_on_edge(double *pos1, double *pos2, int dim)
{
  int i;

  for (i=0; i<dim; i++)
    if (fabs(pos1[i]-pos2[i])<1e-10)
      return 1;
  return 0;
}


double box_discr(double **pointset, int npoints, int dim,
		 double *point)
{
  double discr=1.0, discr2;
  int i, inbox=0, onedge=0;
  for (i=0; i<dim; i++)
    discr *= point[i];
  for (i=0; i<npoints; i++)
    if (point_in_box(pointset[i], point, dim)) {
      inbox++;
      if (point_on_edge(pointset[i], point, dim))
	onedge++;
    }

  //  fprintf(stderr, "Points %d+%d vol %lg", inbox, onedge, discr);


  discr -= (double)inbox/npoints;
  discr2 = discr + (double)onedge/npoints;
  discr = fabs(discr);
  discr2 = fabs(discr2);

  //  fprintf(stderr, " values %lg,%lg\n", discr, discr2);

  return (discr > discr2) ? discr : discr2;
}

double exact_discr(double **pointset, int npoints, int dim)
{
  double *point, **coords, maxdiscr=0.0, discr;
  int i,j, *tuple;
  point = malloc(dim*sizeof(double));
  tuple = malloc(dim*sizeof(int));
  coords = malloc(dim*sizeof(double *));
  for (i=0; i<dim; i++) {
    tuple[i]=0;
    coords[i] = malloc((npoints+1)*sizeof(double));
    for (j=0; j<npoints; j++)
      coords[i][j] = pointset[j][i];
    coords[i][npoints]=1.0;
  }
  
  tuple[0]=-1;
  while (step_tuple(tuple, dim, npoints+1)) {
    for (i=0; i<dim; i++)
      point[i] = coords[i][tuple[i]];
    discr = box_discr(pointset, npoints, dim, point);
    if (discr > maxdiscr) {
      fprintf(stderr, "Point (%g", point[0]);
      for (i=1; i<dim; i++)
	fprintf(stderr, ", %g", point[i]);
      fprintf(stderr, ") discr %g\n", discr);
      maxdiscr = discr;
    }
  }
  return maxdiscr;
}

double rnd_coord_discr(double **pointset, int npoints, int dim, int trials)
{
  double *point, **coords, maxdiscr=0.0, discr;
  int i,j;
  point = malloc(dim*sizeof(double));
  coords = malloc(dim*sizeof(double *));
  for (i=0; i<dim; i++) {
    coords[i] = malloc((npoints+1)*sizeof(double));
    for (j=0; j<npoints; j++)
      coords[i][j] = pointset[j][i];
    coords[i][npoints]=1.0;
  }

  for (i=0; i<trials; i++) {
    for (j=0; j<dim; j++)
      point[j] = coords[j][random()%(npoints+1)];

    discr = box_discr(pointset, npoints, dim, point);
    if (discr > maxdiscr) {
      fprintf(stderr, "%d/%d: discr %g\n", i, trials, discr);
      maxdiscr = discr;
    }
  }
  return maxdiscr;
}
