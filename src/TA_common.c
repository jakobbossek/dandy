/* Contains common functions (= most) for TA discrepancy search */

#include "TA_common.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <string.h>

// not recommended: PRINT_ALL_UPDATES
#define PRINT_ALL_UPDATES
// even more ridiculous!
//#define PRINT_UPDATE_CANDIDATES
//#define DISPLAY_CANDIDATES
//#define PRINT_RANGE_DATA

// for stupid speedup reasons (alternative: make it static, and take care of initialization somehow)
int *coordinate;

// we want to use C library qsort to sort.  
// I made a replacement stump
int dbl_compare(const void *a, const void *b)
{
  const double *xp=a, *yp=b;
  const double x=*xp, y=*yp;
  if (x<y)
    return -1;
  if (x>y)
    return 1;
  return 0;
} // oops: "x-y" gets cast to integer, rounded wrong

void quicksort(int left, int right, double *arr) 
{
  qsort(&arr[left], right-left+1, sizeof(double), dbl_compare);
}

double get_coord(int d, int i) 
{
  return coord[d][i];
}

int get_index_up(int d, double key)
{
  int i=(int)key*(n_coords[d]-1);
  while (coord[d][i] < key)
    i++;
  while ((i>0) && (coord[d][i-1] >= key))
    i--;
  // termination: at first coordinate which is >= key
  // (if i=0 then key==0)
  return i;
}

int get_index_down(int d, double key)
{
  int bound=n_coords[d]-1;
  int i=key*bound;
  while (coord[d][i] > key)
    i--;
  while ((i < bound) && (coord[d][i+1] <= key))
    i++;
  return i;
}

// The following two functions use an "output variable" provided by the caller
// (anything else would just mess up the C memory management)
void round_point_up(double *point, int *output)
{
  int i;
  for (i=0; i<n_dimensions; i++)
    output[i] = get_index_up(i, point[i]);
}

void round_point_down(double *point, int *output)
{
  int i;
  for (i=0; i<n_dimensions; i++)
    output[i] = get_index_down(i, point[i]);
}

void round_point_extradown(double *point, int *output)
{
  int i;
  for (i=0; i<n_dimensions; i++) {
    output[i] = get_index_down(i, point[i]);
    if(output[i]==0)
      output[i]=n_coords[i]-1;
  }
}

void process_coord_data(double **points, int n, int d)
{
  int i,j;
  double tmp_coords[n+2];
  int n_uniq, idx;
  double prev_coord;
  //initialise n_coords[], coord[][]
  n_dimensions=d;
  n_points=n;
  coordinate=malloc(d*sizeof(int));
  for (i=0; i<d; i++)
    coordinate[i]=i;
  n_coords = malloc(n_dimensions*sizeof(double));
  coord = malloc(n_dimensions*sizeof(double *));
  for (i=0; i<n_dimensions; i++) {
    for (j=0; j<n_points; j++)
      tmp_coords[j+1]=points[j][i];
    tmp_coords[0]=0.0;
    tmp_coords[n+1]=1.0;
    quicksort(1, n, tmp_coords);
    // 1. count
    n_uniq=1;
    prev_coord=0;
    for (j=1; j<=n+1; j++) { // inclusive bound: contains 1.0
      if (prev_coord==tmp_coords[j])
	continue;
      n_uniq++;
      prev_coord=tmp_coords[j];
    }
    // 2. transfer
    coord[i]=malloc(n_uniq*sizeof(double));
    idx=1;
    prev_coord=tmp_coords[0];
    coord[i][0]=prev_coord;
    for (j=1; j<=n+1; j++) {
      if (prev_coord==tmp_coords[j])
	continue;
      prev_coord=tmp_coords[j];
      coord[i][idx++] = prev_coord;
    }
    n_coords[i]=n_uniq;
//    fprintf(stderr, "Coordinates %d, %d: ", i, n_uniq);
//    for (j=0; j<n_uniq; j++)
//      fprintf(stderr, "%g ", coord[i][j]);
//    fprintf(stderr,"\n");
  }
  // finished setup for: n_coords[], coord[][]

  // next: transfer point set to into point_index
  point_index=malloc(n_points*sizeof(int *));
  for (i=0; i<n_points; i++)
    point_index[i]=malloc(n_dimensions*sizeof(int));
  for (i=0; i<n_points; i++)
    for (j=0; j<n_dimensions; j++) {
      idx=get_index_up(j, points[i][j]);
      if (coord[j][idx] != points[i][j]) {
	fprintf(stderr, "ERROR: located incorrect coordinate (%g at %d, wanted %g).\n",
		coord[j][idx], idx, points[i][j]);
	abort();
      }
      point_index[i][j] = idx;
    }
  // setup finished.
}

double volume(int *corner)
{
  double vol=1.0;
  int i;
  for (i=0; i<n_dimensions; i++)
    vol *= get_coord(i, corner[i]);
  return vol;
}


//is x in [0,z[ ?
int open(int *x, int *z)
{
  int i;
  for (i=0; i<n_dimensions; i++)
    if (x[i] >= z[i])
      return 0;
  return 1;
}

//is x in [0,z] ?
int closed(int *x, int *z)
{
  int i;
  for (i=0; i<n_dimensions; i++)
    if (x[i] > z[i])
      return 0;
  return 1;
}

int count_open(int *corner)
{
  int n=n_points;
  int i;
  int res=0;
  for (i=0; i<n; i++)
    res += open(point_index[i], corner);
  return res;
}

int count_closed(int *corner)
{
  int n=n_points;
  int i;
  int res=0;
  for (i=0; i<n; i++)
    res += closed(point_index[i], corner);
  return res;
}

// group of functions for grow/"snap up"
int point_critical_dimension(int *corner, int *point)
{
  int i;
  int crit=-1;
  //  fprintf(stderr, "Point");
  //  for (i=1; i<=d; i++)
  //    fprintf(stderr, " %g", point[i]);
  for (i=0; i<n_dimensions; i++)
    if (point[i] >= corner[i]) {
      if (crit>=0) {
	//	fprintf(stderr, " double out (%d,%d)\n", crit, i);
	return -1;
      }
      else
	crit=i;
    }
  //  fprintf(stderr, " crit %d\n", crit);
  return crit;
}

// "really old" O(nd^2) time alg.
void old_grow_box(int *corner)
{
  int order[n_dimensions];
  int n=n_points, d=n_dimensions;
  int i,j,k, swap;
  int curr_d;
  int max_idx;
  for (i=0; i<d; i++)
    order[i]=i;
  for (i=0; i<d; i++) {
    j = i + random()%(d-i);
    if (i != j) {
      swap=order[i];
      order[i]=order[j];
      order[j]=swap;
    }
  }
  /* 
  fprintf(stderr, "Growing base ");
  for (i=0; i<d; i++)
    fprintf(stderr, "%d ", corner[i]);
  fprintf(stderr, "in order ");
  for (i=0; i<d; i++)
    fprintf(stderr, "%d ", order[i]);
  fprintf(stderr, "\n");
   */

  for (i=0; i<d; i++) {
    curr_d=order[i];
    max_idx=n_coords[curr_d]-1;    //    max_val=1.0;
    /* 
    fprintf(stderr, "Direction %d (base", curr_d);
    for (j=0; j<d; j++)
      fprintf(stderr ," %d", corner[j]);
    fprintf(stderr, "), start %d\n", max_idx);
     */
    for (j=0; j<n; j++) {
      if (point_critical_dimension(corner, point_index[j])==curr_d)
	if (point_index[j][curr_d] < max_idx) {
	  max_idx=point_index[j][curr_d];
	  /* 
	  fprintf(stderr, "Because of point ");
	  for (k=0; k<d; k++)
	    fprintf(stderr,"%d ", point_index[j][k]);
	  fprintf(stderr, "bounded at %d\n", max_idx);
	   */
	}
    }
    corner[curr_d] = max_idx;
  }
}

//void grow_box_newer(int *corner)
void grow_box_randomly(int *corner)
{
  int order[n_dimensions];
  int n=n_points, d=n_dimensions;
  int i,j,k, swap, memo;
  int curr_d;
  int new_box[d];

  for (i=0; i<d; i++)
    order[i]=i;
  for (i=0; i<d; i++) {
    j = i + random()%(d-i);
    if (i != j) {
      swap=order[i];
      order[i]=order[j];
      order[j]=swap;
    }
  }
  for (i=0; i<d; i++)
    new_box[i] = n_coords[i]-1;

  for (i=0; i<n; i++) {
    memo=-1;
    for (j=0; j<d; j++) {
      curr_d=order[j];
      k=point_index[i][curr_d];
      if (k >= corner[curr_d]) {
	if (k < new_box[curr_d]) {
	  if (memo<0)
	    memo=curr_d;
	}
	else {
	  memo=-1;
	  break;
	}
      }
    }
    if (memo >= 0)
      new_box[memo] = point_index[i][memo];
  }

#ifdef DEBUG
  if (count_open(corner) != count_open(new_box)) {
    fprintf(stderr, "ERROR: Went from %d to %d points.\n", count_open(corner), count_open(new_box));
    fprintf(stderr, "Old box: ");
    for (i=0; i<d; i++)
      fprintf(stderr, "%d ", corner[i]);
    fprintf(stderr, "\nNew box: ");
    for (i=0; i<d; i++)
      fprintf(stderr, "%d ", new_box[i]);
    fprintf(stderr, "\n");
    abort();
  }
#endif

  for (i=0; i<d; i++)
    corner[i]=new_box[i];
}

//Older version does not perform "memo" check
void grow_box_older(int *corner)
{
  int order[n_dimensions];
  int n=n_points, d=n_dimensions;
  int i,j,k, swap;
  int curr_d;
  int new_box[d];

  for (i=0; i<d; i++)
    order[i]=i;
  for (i=0; i<d; i++) {
    j = i + random()%(d-i);
    if (i != j) {
      swap=order[i];
      order[i]=order[j];
      order[j]=swap;
    }
  }
  for (i=0; i<d; i++)
    new_box[i] = n_coords[i]-1;

  for (i=0; i<n; i++) {
    for (j=0; j<d; j++) {
      curr_d=order[j];
      k=point_index[i][curr_d];
      if (k >= corner[curr_d]) {
	if (k < new_box[curr_d])
	  new_box[curr_d]=k;
	break;
      }
    }
  }
#ifdef DEBUG
  if (count_open(corner) != count_open(new_box)) {
    fprintf(stderr, "ERROR: Went from %d to %d points.\n", count_open(corner), count_open(new_box));
    fprintf(stderr, "Old box: ");
    for (i=0; i<d; i++)
      fprintf(stderr, "%d ", corner[i]);
    fprintf(stderr, "\nNew box: ");
    for (i=0; i<d; i++)
      fprintf(stderr, "%d ", new_box[i]);
    fprintf(stderr, "\n");
    abort();
  }
#endif

  for (i=0; i<d; i++)
    corner[i]=new_box[i];
}

//#define grow_box_randomly grow_box_older

// Rounds box down to borders given by its contained points.
void snap_box(int *corner)
{
  int d=n_dimensions;
  int n=n_points;
  int max_idx[d];
  int i,j;
  for (i=0; i<d; i++)
    max_idx[i]=-1;
  for (i=0; i<n; i++)
    if (closed(point_index[i], corner))
      for (j=0; j<d; j++)
	if (point_index[i][j] > max_idx[j])
	  max_idx[j]=point_index[i][j];
  for (i=0; i<d; i++)
    corner[i] = (max_idx[i] < 0) ? 0 : max_idx[i];
}



//calculates delta(x)
double get_delta(int *x)
{ 
  int i, op;
  double vol, delta;
  int n=n_points;
  
  //calculates no of points in [0,x[
  op = count_open(x);

  //calcualtes volume of box generated by x
  vol = volume(x);
  
  //calculates delta
  delta = vol - (double)op/n;
  //  fprintf(stderr, "Stand.: vol %g pts %d -> %g error %g\n",
  //	  vol, op, (double)op/n, delta);  
  return delta;
}

//calculates bar(delta)(x)
double get_bar_delta(int *x)
{ 
  int i, cl;
  int n=n_points;
  double vol, bdelta;
  
  //calculates no of points in [0,x]
  cl = count_closed(x);
  
  //calcualtes volume of box generated by x
  vol = volume(x);
  
  //calculates bar(delta)
  bdelta = (double)cl/n - vol;
  //  fprintf(stderr, "Stand.: vol %g pts %d -> %g error %g\n",
  //	  vol, cl, (double)cl/n, bdelta);  
  return bdelta;
}

//Generate random search point xc
//Fills coordinates (indexes, not numbers) into its three arguments
void generate_xc(int *xn_plus, int *xn_minus, int *xn_extraminus)
{
  int j, d=n_dimensions;
  double xn[d];
  double temp;
  for(j=0; j<d; j++)
    {
      temp=(double)((double)rand()/RAND_MAX);
      xn[j]=pow(temp,(double)((double)1/(double)d));
    }
  round_point_up(xn, xn_plus);
  round_point_down(xn, xn_minus);
  round_point_extradown(xn, xn_extraminus);		
}

void generate_xc_delta(int *xn_plus)
{
  int j, d=n_dimensions;
  double xn[d];
  double temp;
  for(j=0; j<d; j++)
    {
      temp=(double)((double)rand()/RAND_MAX);
      xn[j]=pow(temp,(double)((double)1/(double)d));
    }
  round_point_up(xn, xn_plus);
}

void generate_xc_bardelta(int *xn_minus, int *xn_extraminus)
{
  int j, d=n_dimensions;
  double xn[d];
  double temp;
  for(j=0; j<d; j++)
    {
      temp=(double)((double)rand()/RAND_MAX);
      xn[j]=pow(temp,(double)((double)1/(double)d));
    }
  round_point_down(xn, xn_minus);
  round_point_extradown(xn, xn_extraminus);		
}


//Generate a random neighbor of xc
//k[i] == range radius in component i, mc = number of components to change
//the three xn_* variables are filled with indexes
void generate_neighbor (int *xn_plus_index, int *xn_minus_index, int *xn_extraminus_index, 
			int *xc_index, int *k, int mc)
{
  int i, j, q, d=n_dimensions;
  double temp, upper_bound, lower_bound;
  double xn[d];

  //First copy the values of the current search point 
  for(j=0; j<d; j++)
    {
      xn[j] = coord[j][xc_index[j]];
    }

  // find mc different coordinates to be changed
  for (j=0; j<mc; j++) {
    i = j + random()%(d-j);
    if (i != j) {
      q = coordinate[j];
      coordinate[j]=coordinate[i];
      coordinate[i] = q;
    }
  }
  
  //set lower and upper bound to the box from which the random neighbor will be sampled
  for(j=0; j<mc; j++){ 	
    if (xc_index[coordinate[j]]-k[coordinate[j]]>=0) 
      lower_bound = coord[coordinate[j]][xc_index[coordinate[j]]-k[coordinate[j]]];
    else 
      lower_bound=0.0;
    if (xc_index[coordinate[j]]+k[coordinate[j]] < n_coords[coordinate[j]])
      upper_bound = coord[coordinate[j]][xc_index[coordinate[j]]+k[coordinate[j]]];
    else 
      upper_bound=1.0;
    
    //draw a random number in [0,1]
    temp=(double)((double)rand()/RAND_MAX);
    temp=(pow(upper_bound,d)-pow(lower_bound,d))*temp + pow(lower_bound,d);
    xn[coordinate[j]]=pow(temp,(double)((double)1/(double)d));
  }

  round_point_up(xn, xn_plus_index);
  round_point_down(xn, xn_minus_index);
  round_point_extradown(xn, xn_extraminus_index);
}

void generate_neighbor_delta(int *xn_plus_index, 
			int *xc_index, int *k, int mc)
{
  int i, j, q, d=n_dimensions;
  double temp, upper_bound, lower_bound;
  double xn[d];

  //First copy the values of the current search point 
  for(j=0; j<d; j++)
    {
      xn[j] = coord[j][xc_index[j]];
    }

  // find mc different coordinates to be changed
  for (j=0; j<mc; j++) {
    i = j + random()%(d-j);
    if (i != j) {
      q = coordinate[j];
      coordinate[j]=coordinate[i];
      coordinate[i] = q;
    }
  }
  
  //set lower and upper bound to the box from which the random neighbor will be sampled
  for(j=0; j<mc; j++){ 	
    if (xc_index[coordinate[j]]-k[coordinate[j]]>=0) 
      lower_bound = coord[coordinate[j]][xc_index[coordinate[j]]-k[coordinate[j]]];
    else 
      lower_bound=0.0;
    if (xc_index[coordinate[j]]+k[coordinate[j]] < n_coords[coordinate[j]])
      upper_bound = coord[coordinate[j]][xc_index[coordinate[j]]+k[coordinate[j]]];
    else 
      upper_bound=1.0;
    
    //draw a random number in [0,1]
    temp=(double)((double)rand()/RAND_MAX);
    temp=(pow(upper_bound,d)-pow(lower_bound,d))*temp + pow(lower_bound,d);
    xn[coordinate[j]]=pow(temp,(double)((double)1/(double)d));
  }

  round_point_up(xn, xn_plus_index);
}

void generate_neighbor_bardelta(int *xn_minus_index, int *xn_extraminus_index, 
			int *xc_index, int *k, int mc)
{
  int i, j, q, d=n_dimensions;
  double temp, upper_bound, lower_bound;
  double xn[d];

  //First copy the values of the current search point 
  for(j=0; j<d; j++)
    {
      xn[j] = coord[j][xc_index[j]];
    }

  // find mc different coordinates to be changed
  for (j=0; j<mc; j++) {
    i = j + random()%(d-j);
    if (i != j) {
      q = coordinate[j];
      coordinate[j]=coordinate[i];
      coordinate[i] = q;
    }
  }
  
  //set lower and upper bound to the box from which the random neighbor will be sampled
  for(j=0; j<mc; j++){ 	
    if (xc_index[coordinate[j]]-k[coordinate[j]]>=0) 
      lower_bound = coord[coordinate[j]][xc_index[coordinate[j]]-k[coordinate[j]]];
    else 
      lower_bound=0.0;
    if (xc_index[coordinate[j]]+k[coordinate[j]] < n_coords[coordinate[j]])
      upper_bound = coord[coordinate[j]][xc_index[coordinate[j]]+k[coordinate[j]]];
    else 
      upper_bound=1.0;
    
    //draw a random number in [0,1]
    temp=(double)((double)rand()/RAND_MAX);
    temp=(pow(upper_bound,d)-pow(lower_bound,d))*temp + pow(lower_bound,d);
    xn[coordinate[j]]=pow(temp,(double)((double)1/(double)d));
  }

  round_point_down(xn, xn_minus_index);
  round_point_extradown(xn, xn_extraminus_index);
}

