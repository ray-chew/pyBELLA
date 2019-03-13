
#define TOMMASO

//#define RUPERT

/* auxiliary labels for addressing various solution components */
#define NSPEC 1   /* no of advected scalars */
#define BUOY  0   /* auxiliary pot. temp. perturbation variable */
#define QV    1   /* foreseen for: water vapor */
#define QC    2   /* foreseen for: cloud water */
#define QR    3   /* foreseen for: rain water  */

/* Output parameters */
#define HDFFORMAT  

/* Debugging output options in various places */
#define OUTPUT_SUBSTEPS    0       /* time step after which detailed output is generated - off */
#define OUTPUT_SPLITSTEPS  0       /* on-off */
#define OUTPUT_HYDROSTATES 0       /* on-off */
#define OUTPUT_FLUXES      0       /* on-off */
#define OUTPUT_LAP_CELLS   0       /* on-off */
#define OUTPUT_RHS_CELLS   0       /* on-off */
#define OUTPUT_LAP_NODES   0       /* on-off */
#define OUTPUT_RHS_NODES   0       /* on-off */
#define OUTPUT_ADV_FLUXES  0       /* on-off */
	
/* ============================================= 
 Explicit predictor options
 ============================================= */

/* ============================================= 
 Semi-implicit solver options  
 ============================================= */
/* PRECON:
 0: no preconditioner
 1: diagonal (no gravity) or column-wise (with gravity) preconditioner
 */
#define PRECON 1
#define DIV_CONTROL_LOCAL 0  /* determines norm for convergence test 0: L1, 1: Linfty */

#define P2_FULL_STENCIL 1.0        /* values: 0.0, 1.0;  0.0 = 5/7pt stencil,  1.0=9/27pt stencil */
#define P2_DIAGONAL_FIVE_POINT 1.0 /* 0.0, 1.0; as above but for node-based Poisson op.           */

/*
#define NEUMANN_Y_BOTTOM_BC
//#define IMP_MIDPT_FOR_NODAL_PI
 */
/* TODO: 
 
 1) Make all appearances of "extern ..." disappear except for those of
 User_Data ud;
 double* W0;
 
 2) Whereever the auxiliary double array W0 is used, introduce an external "W0_in_use" flag, so 
    that other routines can query whether or not this aux array is currently occupied.
 DONE
 
 3) Get rid of the full-size dSol arrays
 
 4) Implement   Hydrostatic_Initial_Pressure()  [and lots of other stuff] for 3D.
 
 5) Rewrite  slanted_wall_min() and related routines for cross-wall advective flux  
    rhoYu  instead of for rhou to improve compatibility with the pseudo-incompressible
    model
 DONE
 
 6) Implement Piotr's conjugate residual scheme also for the nodal projection;
    re-invoke a multigrid solver, such as HYPRE
 
 7) Limiter with plateau- instead of extrema-detection
 
 8) Try applying implicit midpoint rule systematically for the nodal pressure, too, 
    while discarding the pressure (correction) computed in the second projection.
    This would allow us to automatically synchronize the nodal and cell-centered
    pi and rhoY = P variables, without having to invoke some additional averaging 
    procedure.
 
9)  Re-implement outer iteration in the nodal Helmholtz solve so as to guarantee
    exact compliance of nodal Exner pressure with the equation of state P = P(pi).
 
10) Implement a version of the nodal projection that always compute pi-increments
    rather than full pi' data. This will allow us to claim working with "full 
    not perturbation variables" only. 
 
11) Try accessing the Arakawa-Konor model. 
 
 */


#include <assert.h> 
