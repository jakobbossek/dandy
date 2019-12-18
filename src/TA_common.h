/* Contains common functions (= most) for TA discrepancy search */


#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

// not recommended: PRINT_ALL_UPDATES
#define PRINT_ALL_UPDATES
// even more ridiculous!
//#define PRINT_UPDATE_CANDIDATES
//#define DISPLAY_CANDIDATES
//#define PRINT_RANGE_DATA

int n_dimensions, n_points;
double *n_coords;
double **coord;
int **point_index;

#define grow_box_randomly grow_box_newer

// we want to use C library qsort to sort.  
// I made a replacement stump
void quicksort(int left, int right, double *arr);


double get_coord(int d, int i);

int get_index_up_lin(int d, double key);

int get_index_down_lin(int d, double key);

// The following functions use an "output variable" provided by the caller
// (anything else would just mess up the C memory management)
void round_point_up(double *point, int *output);
void round_point_down(double *point, int *output);
void round_point_extradown(double *point, int *output);

void process_coord_data(double **points, int n, int d);

double volume(int *corner);

//is x in [0,z[ ?
int open(int *x, int *z);

//is x in [0,z] ?
int closed(int *x, int *z);

int count_open(int *corner);

int count_closed(int *corner);

void grow_box_randomly(int *corner);

void old_grow_box(int *corner);

// Rounds box down to borders given by its contained points.
void snap_box(int *corner);

//calculates delta(x)
double get_delta(int *x);

//calculates bar(delta)(x)
double get_bar_delta(int *x);

//Generate random search point xc
//Fills coordinates (indexes, not numbers) into its three arguments
void generate_xc(int *xn_plus, int *xn_minus, int *xn_extraminus);

void generate_xc_delta(int *xn_plus);

void generate_xc_bardelta(int *xn_minus, int *xn_extraminus);

//Generate a random neighbor of xc
//k[i] == range radius in component i, mc = number of components to change
//the three xn_* variables are filled with indexes
void generate_neighbor (int *xn_plus_index, int *xn_minus_index, int *xn_extraminus_index, 
			int *xc_index, int *k, int mc);

void generate_neighbor_delta(int *xn_plus_index, 
			     int *xc_index, int *k, int mc);

void generate_neighbor_bardelta(int *xn_minus_index, int *xn_extraminus_index, 
			int *xc_index, int *k, int mc);
