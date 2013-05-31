/*==================================================================
 * neurons.c
 *
 * Jesse Berwald
 *
 * Opened: Oct. 16, 2008
 *
 * Solves a system of ODE's defining a neuronal network as described
 * in Terman, 2008.
 ===================================================================*/
#include "neurons.h"
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define NPY_DOUBLE                PyArray_DOUBLE

static PyObject *network_fitness(PyObject *self, PyObject *args);

/* ==== Initialize functions to be called by Python ==== */

/* Python methods table */
static PyMethodDef neuronsMethods[] = {
  {"network_fitness", network_fitness, METH_VARARGS},
  {NULL, NULL}
};

/* initialize C function visible to Python */
void initneurons_file_cells(void)
{
  (void) Py_InitModule("neurons_file_cells", neuronsMethods);
  import_array();
}

/* ==== (global) Indexing variables ==== */
int i, j;

/* ==== C function called by Python. ==================================

 *network_fitness():

   This function is called from Python and is used to populate C data
   types with Numpy arrays, and Python scalars, as well as run the GSL
   ode solver (Runge-Kutta 45; i.e. "rkf45")

   Input     -- final time (stopping time for ode solver)
             -- coupling matrices (entries composed of g_syn terms)
	     -- initial values
	     -- (attempted) ode step.
	     -- dictionary of parameter values for FNE's

   Output     
             Items written to disk (/var/tmp/neurons_tmp) to be analyzed Pythonically later).
             -- time vector returned by ode solver.
	     -- vector of time steps taken  
	     -- Matrix of time series values for each step for each variable.
	     These are stored in one matrix that 2d and of size
	     (NUM_EQNS)x(100 . dt . final time)

	     Summary of procedure:
	     
	     Various PyObjects, PyArrayObjects, C array objects, and
	     variables are initialized.

	     Python input *args is parsed (see below). 

	     Dimensions of system of equations is set from Numpy array
	     object (actually from a Numpy array struct member)
	     
	     Length of sampling window
	     and default step size are set. The coupling matrices
	     created in Python module to represent network
	     architecture are pointed to via the C pointer arrays
	     **C_xx. 

	     Next, the initial conditions passed in from the Python
	     module are used to set the gsl_vectors to their initial
	     state.

	     Note: Even though Excitatory cells only have three
	     equations, in order to have a square-matrix system the
	     excitatory portions of solution vectors are zeroed out by
	     simple vector multiplication. Some overhead, but
	     (hopefully) not too bad.

	     Note on args passed in from Python:

	     arg types (passed to solve_network() in this order --
	     this is very important!):

	     tfinal         = double ('d')
	     cpl_ei         = 2D array ('O!')
	     cpl_I_E         = 2D array ('O!')
	     cpl_I_I         = 2D array ('O!')
	     pystate        = 1D array ('O!') (initial values passed in here!)
	     cell_input     = Python dictionary ('O')
	     dt             = Python float ('d')
	     NUM_EX         = Python int ('i')
	     NUM_IN         = Python int ('i')
	     gamma          = Python float ('d')
	     sigma          = Python float ('d')
 ======================================================================= */
static PyObject *network_fitness(PyObject *self, PyObject *args)
{
  PyArrayObject *cpl_E_I, *cpl_I_E, *cpl_I_I, *pystate;
  PyObject *cell_input;
  double tfinal, dt, gamma, sigma, temp;//, current_item;
  int EX_LOOP;

  // File handles and file name
  FILE *fh;
    const char *tmpfile = "/var/tmp/neurons_tmp";
//"/data/jberwald/neurons/tmp_data/neurons_tmp"; //"/data/tmp/neuronsTmp";

  Pstruct Params;
  Pstruct *param_ptr = &Params;
 
  /* ==== parse Python args ==== */
  if (!PyArg_ParseTuple(args, "dO!O!O!O!Odiiddi",
			&tfinal,
			&PyArray_Type, &cpl_E_I,
			&PyArray_Type, &cpl_I_E,
			&PyArray_Type, &cpl_I_I,
			&PyArray_Type, &pystate,
			&cell_input,
			&dt,
			&NUM_EX,
			&NUM_IN,
			&gamma,
			&sigma,
			&TREE)) return NULL;
 
  /* random number generator object */
  // if (sigma != 0) && (gamma != 0) {
      gsl_rng *rand_array;  
      // }

  /* make sure stuff made it in and is of correct type */
  if (NULL == cpl_E_I) return NULL;
  if (NULL == cpl_I_E) return NULL;
  if (NULL == cpl_I_I) return NULL;
  if (!PyDict_Check(cell_input)) {
      // set except context for TypeError
      PyErr_SetString(PyExc_TypeError, 
		      "cell_inputs is not a dictionary!");
      // tell interpreter to raise an exception
      return NULL;
  }
  if (!PyArray_Check(pystate)) {
    printf("pystate object must be a Numpy array\n");
    return NULL;
  }
  if (!PyArray_Check(cpl_E_I)) {
    printf("cpl_E_I object must be a Numpy array\n");
    return NULL;
  }
  if (!PyArray_Check(cpl_I_E)) {
    printf("cpl_I_E object must be a Numpy array\n");
    return NULL;
  }
  if (!PyArray_Check(cpl_I_I)) {
    printf("cpl_I_I object must be a Numpy array\n");
    return NULL;
  }

  /* ===== INITIALIZATION OF VARIABLES AND ARRAYS ===== */

  /* dim of search space -- NDIM defined globally in nmsearch.h
     NUM_EQNS = (4 eqns per cells) * NDIM */
  NDIM = NUM_EX + NUM_IN;
  NDIM2 = 2*NDIM;
  NDIM3 = 3*NDIM;
  NUM_EQNS = 4*NDIM;

  /* /\* indices of state array to store (do we average over E and I */
  /*    cells, or just E cells?). *\/ */
  /* int ihi = NUM_EX; */

  /* Variables for current time and step size within ode stepper. */
  double t = 0.0, h = dt;

  /* ==== Make Numpy arrays contiguous in memory ==== */
  pystate = (PyArrayObject *)
    PyArray_FROM_OTF((PyObject*)pystate, NPY_DOUBLE, C_ARRAY);
  cpl_E_I = (PyArrayObject *)
    PyArray_FROM_OTF((PyObject*)cpl_E_I, NPY_DOUBLE, C_ARRAY);
  cpl_I_E = (PyArrayObject *)
    PyArray_FROM_OTF((PyObject*)cpl_I_E, NPY_DOUBLE, C_ARRAY);
  cpl_I_I = (PyArrayObject *)
    PyArray_FROM_OTF((PyObject*)cpl_I_I, NPY_DOUBLE, C_ARRAY);

  /* fill Params struct members */
  if (sigma > 0) {
    MALLOC(Params.randInput, NDIM * sizeof(double));
  }

  /* Initialize parameters and parameter matrices */
  Params.sigma = sigma;
  initialize_Params(cell_input, param_ptr);
  Params.gamma = Params.amp = gamma;
  Params.EtoI = pymat_to_Cmatptrs(cpl_E_I);
  Params.ItoE = pymat_to_Cmatptrs(cpl_I_E);
  Params.ItoI = pymat_to_Cmatptrs(cpl_I_I);

  /* make space for Cstate array */
  double Cstate[NUM_EQNS];

  /* initialize Cstate  and tmpvec with initial values.
 tmpvec will hold ODE solution at each step. œ*/
  double *tmpvec = (double *)pystate->data;
  for (i=0; i < NUM_EQNS; i++) {
    Cstate[i] = tmpvec[i];
  }


  /* append initial time and avg. exc. state to vstate or write to
     file. */
  fh = fopen( tmpfile, "w" );
  if (fh == NULL) {
    perror("failed to open neurons_tmp!\n");
    return Py_False;
  }

  /* Initial values -- note space after time in first column */
  fprintf( fh, "%g ", t );

 /* Initial excitatory state. If tree structure, only average the
       excitatory cells in the top cluster */
  if (TREE==1) { 
    EX_LOOP = 5;
  }
  else {
    EX_LOOP = NUM_EX;
  }
  for (i=0; i < EX_LOOP; i++){
    fprintf( fh, "%g ", Cstate[i] ); 
  }
  fprintf( fh, "\n" );

  /* =============
     initialize ode solver -- memory allocated in step, control, and
     evolve. Create workspace for state matrix (allocate memory),
     initially NUM_EX x maxdat. state is the array that holds only the
     results that we want to keep (i.e. the voltage of excitatory
     cells). Cstate will hold the vector of current values at each
     time t, but most will not be stored.
     ============== */
  const gsl_odeiv_step_type *ode_solver
    = gsl_odeiv_step_rk4;
    //  = gsl_odeiv_step_rk4imp;
  //    = gsl_odeiv_step_rkf45;

  gsl_odeiv_step *step
    = gsl_odeiv_step_alloc(ode_solver, NUM_EQNS);
  gsl_odeiv_control *control
    =  gsl_odeiv_control_y_new(1e-6, 0.0);
  gsl_odeiv_evolve *evolve
    = gsl_odeiv_evolve_alloc(NUM_EQNS);

  /* ==================
     
     gsl_odeiv_system struct members are 

     int (* function)
     int (* jacobian)
     size_t dimensions
     void * params
     
     Since we don't use a jacobian, pass NULL to that member.

     dfdt is fitshugh-nagumo defined in neurons.h
     ===================  */
  gsl_odeiv_system sys 
    = {dfdt, NULL, NUM_EQNS, param_ptr};

  /* ===== Setup the RNG if sigma > 0 ==== */
  if (sigma > 0) {
    const gsl_rng_type *RNG;
    
    gsl_rng_env_setup();

    RNG = gsl_rng_default;
    rand_array = gsl_rng_alloc (RNG);
  }
  /* ====== ODE solver loop =====

     vstate: each loop append (t, exc_avg, inh_1, inh_2,..., inh_n)

     ============================= */
  int solver_state = 0;
   while (t < tfinal)
    {
      int status = gsl_odeiv_evolve_apply(evolve,
      					  control,
      					  step,
      					  &sys,
      					  &t,
      					  tfinal,
      					  &h,
      					  Cstate);
      /* NOT IMPLEMENTED UNTIL GSL VERSION ~1.15 */
      /* int status = gsl_odeiv_evolve_apply_fixed_step(evolve, */
      /* 					  control, */
      /* 					  step, */
      /* 					  &sys, */
      /* 					  &t, */
      /* 					  tfinal, */
      /* 					  &h, */
      /* 					  Cstate); */

      if (status != GSL_SUCCESS) {
	solver_state = 1;
	printf("Solver failed before tfinal \n");
	break;
      }
    
      /* =====================
	 File format: 
	 
	 time cell1 cell2 ... celln 
	 ======================*/
      /* record current time */
      fprintf( fh, "%g ", t );
    
      /* excitatory cells*/
      for (i=0; i < EX_LOOP; i++) {
	fprintf( fh, "%g ", Cstate[i] );
      }
      fprintf( fh, "\n" );
      
      /* new random vector */
      if (Params.sigma > 0) {
	for (i=0; i < NDIM; i++) {
	  //	  Params.extInput[i] = gsl_ran_gaussian(r, sigma);
	  Params.randInput[i] = gsl_ran_gaussian(rand_array, sigma);
	}
      }
    } // end while loop
   
   /* ===== CLOSE DOWN SHOP. FREE MEMORY AND DECREF PYOBJECTS ===== */
  /* free gsl_odeiv objects */
  gsl_odeiv_evolve_free(evolve);
  gsl_odeiv_control_free(control);
  gsl_odeiv_step_free(step);

  /* decref numpy objects sent in from python */
  free(Params.EtoI);
  free(Params.ItoE);
  free(Params.ItoI);
  
  if (Params.sigma > 0) free(Params.randInput);

  /* decref input arrays */
  Py_XDECREF(pystate);
  Py_XDECREF(cpl_E_I);
  Py_XDECREF(cpl_I_E);
  Py_XDECREF(cpl_I_I);

/* Close down the file */
  fclose( fh );

  if (solver_state == 0) {
    return Py_True;
  }
  else {
    return Py_False;
  }   

}
// END NETWORK_FITNESS()


/* ============================================

   Memory allocation and deallocation functions

   ============================================ */

/* Make a 2D C array from a Numpy array */
double **pymat_to_Cmatptrs(PyArrayObject *arrayin)
{
  double **mat, *a;
  int nx, ny;
  
  nx = arrayin->dimensions[0];
  ny = arrayin->dimensions[1];
  
  mat = ptrvector(nx);  // array of pointers

  a = (double *) arrayin->data; /* ptr to first element of arrayin
				   data (cast as double pointer) */
  for (i=0; i < nx; i++) {
    mat[i] = a + i*ny;
  }
  return mat;
}
  

/* ==== Make a Python Array Obj. from a PyObject, ================
     generates a double vector w/ contiguous memory which may be a new
     allocation if the original was not a double type or contiguous !!
     Must DECREF the object returned from this routine unless it is
     returned to the caller of this routines caller using return
     PyArray_Return(obj) or PyArray_BuildValue with the "N" construct
     !!!
   ================================================
*/
PyArrayObject *pyvector(PyObject *objin)  
{
    return (PyArrayObject *) PyArray_ContiguousFromObject(objin,
        NPY_DOUBLE, 1,1);
}
/* ==== Create 1D Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.             */
double *pyvec_to_Carrayptr(PyArrayObject *arrayin)  
{
    return (double *) arrayin->data;  /* pointer to arrayin data as double */
}

long *pyvec_to_Carrayptr_long(PyArrayObject *arrayin)  
{
  return (long *) arrayin->data;  /* pointer to arrayin data as double */
}

int *pyvec_to_Carrayptr_int(PyArrayObject *arrayin)  
{
  return (int *) arrayin->data;  /* pointer to arrayin data as double */
}

/* ==== Memory allocation and memory freeing functions ==== */

/* allocate memory for vector of pointers */
double **ptrvector(long n)
{
  double **v;
  v = (double **)malloc((size_t) n*sizeof(double *)); 
  if (!v) {
    printf("In **ptrvector, alloc of memory for double array failed!");
    exit(0);
  }
  return v;
}

/* free double *vector (vec of pointers) */
void free_Carrayptrs(double **vec)
{
  free((char*) vec);
}


/* ============================================================
   Set members of paramstruct from cell_input. Params holds values
   that are not changed, except epsilon array (Params.eps), which is
   updated to new vertex of simplex each time through
   solve_network_gs.  
   ============================================================ */
void initialize_Params(PyObject *cell_input, Pstruct *p)
{
  PyArrayObject *ext, *arr;
  PyObject *keystring, *keys, *value;
  double Cvalue; 
  char *Cstr;
  int len_keys;

  /* ===============
     Set cell input parameter values from cell_inputs dict (from
     Python module). PyDict_Keys() return PyList object. 
     =============== */
  keys = PyDict_Keys(cell_input); // new reference

  if (!PyList_Check(keys)){
    printf("cell_input keys not returned as list");
    exit(0);
  }

  len_keys = PyList_Size(keys);

  /* =============== 
     Grab the keys from cell_inputs and put then into a char
     array 
     =============== */
  for (i=0; i < len_keys; i++) {
    keystring = PyList_GetItem(keys, i); 
    Cstr = PyString_AsString(keystring);

    if (strcmp(Cstr, "epsilon") == 0) {
      arr = (PyArrayObject *) PyDict_GetItem(cell_input, keystring);
      (*p).eps = (double *)arr->data;
    }
    else if (strcmp(Cstr, "external_input") == 0) {
      ext = (PyArrayObject *) PyDict_GetItem(cell_input, keystring);
	(*p).extInput = (double *)ext->data;
    }
/*       if (p->sigma == 0) { */
/* 	ext = (PyArrayObject *) PyDict_GetItem(cell_input, keystring); */
/* 	(*p).extInput = (double *)ext->data; */
/*       } */
/*       else continue; */
      
//    }
/*     else if (strcmp(Cstr, "external_I") == 0) { */
/*       // borrowed ref */
/*       extI = (PyArrayObject *)PyDict_GetItem(cell_input, keystring); */
/*       (*p).extInputI = (double *)extI->data; */
/*     } */
    else {
      value = PyDict_GetItem(cell_input, keystring); // borrowed ref
      Cvalue = PyFloat_AsDouble(value);
    }
      
    if (strcmp(Cstr, "alpha_x") == 0) {
      (*p).alpha_x = Cvalue;
    }
    else if (strcmp(Cstr, "beta_x") == 0) {
      (*p).beta_x = Cvalue;
    }
    else if (strcmp(Cstr, "alpha") == 0) {
      (*p).alpha = Cvalue;
    }
    else if (strcmp(Cstr, "beta") == 0) {
      (*p).beta = Cvalue;
    }
    else if (strcmp(Cstr, "alpha_In") == 0) {
      (*p).alpha_In = Cvalue;
    }
    else if (strcmp(Cstr, "beta_In") == 0) {
      (*p).beta_In = Cvalue; 
    }
    else if (strcmp(Cstr, "theta") == 0) {
      (*p).theta = Cvalue;
    }
    else if (strcmp(Cstr, "theta_I") == 0) {
      (*p).theta_In = Cvalue;
    }
    else if (strcmp(Cstr, "theta_x") == 0) {
      (*p).theta_x = Cvalue;
    }
    else if (strcmp(Cstr, "v_In") == 0) {
      (*p).VbarI = Cvalue;
    }
    else if (strcmp(Cstr, "v_Ex") == 0) {
      (*p).VbarE = Cvalue;
    }
    
  } // end cell_input operations

  /* set dw/dt params b and c here */
  (*p).bparam = 0.8;
  (*p).cparam = 0.7;

  /* delete reference to keys list */
  Py_XDECREF(keys); 
}
// END NEURONS.C
