/* ==== NEURONS.H ==== */
#ifndef _neurons_h
#define _neurons_h

#include <Python.h> /* for PyMem_ and friends and stdio/stdlib */

#define PY_ARRAY_UNIQUE_SYMBOL    neurons_ARRAY_API 
#include <numpy/arrayobject.h>

#define NPY_DOUBLE                PyArray_DOUBLE
/* MALLOC */
#define EEK(s) { perror((s)); exit(EXIT_FAILURE); }
#define MALLOC(s,t) if(((s) = malloc(t)) == NULL) { EEK("error: malloc() "); }
#define Py_MALLOC(s,t,type) if(((s) = PyMem_New(type,t)) == NULL) { EEK("error: malloc() "); }
#define CALLOC(a, n, s) if(((a) = calloc(n, s)) == NULL) { EEK("error: calloc() "); }
#define SWAP(a,b) {temp = (a); (a) = (b); (b) = temp;}

#define C_ARRAY                   NPY_ALIGNED | NPY_CONTIGUOUS | NPY_FORCECAST

/* ==== GLOBAL DIMENSIONS ==== set in solve_network() */
int NUM_EQNS;      /* ==4*NDIM. We will allot memory for 4 eqns per
		      cell, with excitatory cells padded with zeros
		      for the x variable. this makes indexing
		      easier, as well a matrix mult. */

int NDIM;          /* set to NUM_EX + NUM_IN later */
int NUM_EX;        /* number of excitatory cells */
int NUM_IN;        /* number of inhibitory cells */
int NDIM2;         /* 2*NDIM */
int NDIM3;         /* 3*NDIM */
int TREE;            /* Used in analyzing network in self-similar
			structure: Only average excitatory nodes in
			top-level cluster */

//GSL_RNG_SEED = 123;

/* struct to hold network params */
typedef struct {

  double *eps;
  double bparam;
  double cparam;
  double alpha_In;
  double alpha_x;
  double alpha;
  double beta_In;
  double beta_x;
  double beta;
  double theta_x;
  double theta_In;
  double theta;
  double VbarE;
  double VbarI;
  double *extInput;
  double *randInput;

  // coupling matrices
  double **ItoE;
  double **EtoI;
  double **ItoI;

  double sigma;
  double amp;   // amplitude for oscillatory input. takes value of gamma
  double gamma; // scaling factor for random numbers
  double omega; // frequency of periodic input function 

} Pstruct;


/* ==== Function prototypes for neurons.c ==== */
/* static PyObject *network_fitness(PyObject *self, PyObject *args); */
void initneurons(void);
PyArrayObject *array_swap(PyArrayObject *arr, int *cdims);
double **pymat_to_Cmatptrs(PyArrayObject *arrayin);
double *pyvec_to_Carrayptr(PyArrayObject *arrayin);
long *pyvec_to_Carrayptr_long(PyArrayObject *arrayin);
int *pyvec_to_Carrayptr_int(PyArrayObject *arrayin);
PyArrayObject *pyvector(PyObject *objin);
double **ptrvector(long n);
void free_Carrayptrs(double **vec);
void initialize_Params(PyObject *cell_input, Pstruct *p);
void array_check(PyArrayObject *arr);

/* ==== prototypes for FITZHUGH_NAGUMO.C visible to NEURONS.C ==== */
int dfdt(double t, const double state[], double f[], void *p);
static double sigmoid(double x, double theta);
static void coupling_sum(double a[], 
			 double b[], 
			 double c[], 
			 Pstruct *p, 
			 const double s[]);

#endif
