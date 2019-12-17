#include <math.h>
#include "bz_discr.h"

// Name: The data structure is a tree by Overmars and Yap
// (But the algorithm is the Dobkin/Epstein/Mitchell one.)
double oydiscr(double **pointset, int dim, int npoints, double *lower);
