/*==================================================================
 * fitzhugh_nagumo.c
 * 
 * Jesse Berwald
 *
 * Version:    0.1  
 * Repository: .svn_repo/neuralNet  
 *
 * Opened: October 15, 2008
 *
 * This file holds the FitzHugh-Nagumo equations called by
 * solve_network() in sn.c
 * 
 *===================================================================*/
#include "neurons.h"
#include <math.h>
#include <gsl/gsl_errno.h>

/* #ifdef malloc */
/* #undef malloc */
/* #define malloc PyMem_Malloc */
/* #endif */
/* #ifdef free */
/* #undef free */
/* #define free PyMem_Free */
/* #endif */

#define PI 3.141592653589793

/* #define PY_ARRAY_UNIQUE_SYMBOL neurons_ARRAY_API */
/* #define NO_IMPORT_ARRAY */
/* #include <numpy/arrayobject.h> */


/* ==== for indexing ==== */
int i, j;

/* ==== FitzHugh-Nagumo System ==== 

 FitzHugh-Nagumo Equations with coupling between network cells.
     -------------------------------------------------------------
     f[] = dFdt[0,...,4*NUM_CELLS-1], where F==state is composed of
     the four variables that make up our system, as follows:

     dvdt[] = state[0,..., NUM_CELLS-1] 
     dwdt[] = state[NUM_CELLS,..., 2*NUM_CELLS-1]
     dsdt[] = state[2*NUM_CELLS,..., 3*NUM_CELLS-1]
     dxdt[] = state[3*NUM_CELLS,..., 4*NUM_CELLS-1] 

     Variables and their indices in state vector:

     v  -- 0...NUM_CELLS
     w  -- NUM_CELLS...2*NUM_CELLS
     s  -- 2*NUM_CELLS...3*NUM_CELLS
     x  -- 3*NUM_CELLS...4*NUM_CELLS
     
 
 dfdt() takes void *ptr to be cast as struct paramstruct *p 
 ==================================== */
int dfdt(double t, 
	 const double state[], 
	 double f[], 
	 void *params) 
{

  Pstruct *p = (Pstruct *) params;

  /* vectors to hold coupling sum calculation */
  double cei[NUM_EQNS];
  double cie[NUM_EQNS];
  double cii[NUM_EQNS];

  /* initialize coupling sums to zero--probably faster with calloc */
  for (i=0; i < NUM_EQNS; i++) {
    cei[i] = 0.0;
    cie[i] = 0.0;
    cii[i] = 0.0;
  }
  /* calculate the SUM(s_j) terms in FNE's */
  coupling_sum(cei, cie, cii, p, state);

  /* ================= FITZHUGH-NAGUMO EQUATIONS ================

     ==== dvdt ==== 
     vdot = v - (v^3)/3 - w - (g_xy *(v-[VI or VE])*(cxy)
     ============== */
  for (i=0; i < NDIM; i++) {
    f[i] = state[i] - (state[i]*state[i]*state[i])/3 - state[NDIM + i] 
      - (state[i] - (p->VbarI)) * cie[i]
      - (state[i] - (p->VbarE)) * cei[i] 
      - (state[i] - (p->VbarI)) * cii[i]
      + (p->gamma) * (*p).extInput[i]; 
    if ((p->sigma) > 0) f[i] += (p->gamma) * (*p).randInput[i];
  }
  /* ==== dwdt ==== 
     wdot = eps(v - b*w + c)
     ============== */
  for (i=0; i < NDIM; i++) {
    f[NDIM + i] 
      = (*p).eps[i] * (state[i] - (p->bparam) * state[NDIM + i] + (p->cparam));
  }
  /* ==== dsdt ==== */
  for (i=0; i < NUM_EX; i++) {
    f[NDIM2 + i] 
      = (p->alpha) * (1-state[NDIM2 + i]) * sigmoid(state[i], p->theta) 
      - (p->beta) * state[NDIM2 + i];
  }
  /* ==== ds_I/dt ==== */
  for (i=NUM_EX; i < NDIM; i++) {
    f[NDIM2 + i] 
      = (p->alpha_In) * (1-state[NDIM2 + i]) * 
      sigmoid(state[NDIM3 + i], p->theta_x) - (p->beta_In) * state[NDIM2 + i];
  }
  /* ==== dxdt ==== */
  // Exc. terms 
  for (i=0; i < NUM_EX; i++) {
    f[NDIM3 + i] = 0.0;
  }
  // Inhib. terms
  for (i=NUM_EX; i < NDIM; i++) {
    f[NDIM3 + i] 
      = (*p).eps[i] * (p->alpha_x) * (1-state[NDIM3 + i]) * 
      sigmoid(state[i], p->theta_In) - (*p).eps[i] * (p->beta_x) * 
      state[NDIM3 + i];
  }
  return GSL_SUCCESS;
}  // END DFDT

/* =======Functions used by dfdt() ======= */

/* Steep sigmoidal curve */
static double sigmoid(double x, double theta)
{
  return 1.0 / (1.0 + exp(-100*(x-theta)));
}

/* ==============
   Coupling sums
   ============== 
   (These are sparse matrices. We could make this more efficient.) */
static void coupling_sum(double cpl_sum_E_I[], 
			 double cpl_sum_I_E[], 
			 double cpl_sum_I_I[],
			 Pstruct *p, 
			 const double state[])
{
  /* coupling strengths for E->I connections */
/*   for (i=0; i < NUM_EQNS; i++) { */
/*     for (j=0; j < NUM_EQNS; j++) { */
  // inh. cells
  for (i=NUM_EX; i < NDIM; i++) {
    // loop over excitatory connections
    for (j=NDIM2; j < NDIM2 + NUM_EX; j++) { 
      cpl_sum_E_I[i] += (*p).EtoI[i][j] * state[j];
    }
  }
  /* coupling strengths for I->E connections */
/*   for (i=0; i < NUM_EQNS; i++) { */
/*     for (j=0; j < NUM_EQNS; j++) { */
  // exc. cells
  for (i=0; i < NUM_EX; i++) {
    // loop over inhib. connections
    for (j=NDIM2+NUM_EX; j < NDIM3; j++) {
      cpl_sum_I_E[i] += (*p).ItoE[i][j] * state[j];
    }
  }
  /* coupling strengths for I->I connections */
/*   for (i=0; i < NUM_EQNS; i++) { */
/*     for (j=0; j < NUM_EQNS; j++) { */
  // inh. cells
  for (i=NUM_EX; i < NDIM; i++) {
    // inh. connections
    for (j=NDIM2+NUM_EX; j < NDIM3; j++) {
      cpl_sum_I_I[i] += (*p).ItoI[i][j] * state[j];
    }
  }  
}

// EOF
