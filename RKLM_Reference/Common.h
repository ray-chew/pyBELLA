#define NSPEC 1   /* no of advected scalars */
#define BUOY  0   /* auxiliary pot. temp. perturbation variable */
#define QV    1   /* foreseen for: water vapor */
#define QC    2   /* foreseen for: cloud water */
#define QR    3   /* foreseen for: rain water  */

/* Output parameters */
#define HDFFORMAT  

/* Output options in main.c for debugging;  1 -> output */
#define OUTPUT_SUBSTEPS 0
#define OUTPUT_SPLITSTEPS 0

/*
#define OUTPUT_FLUXES
 */

#define ONE_POINT_TIME_SERIES

/* ============================================= 
 Initial Data Options
 ============================================= */

/* ============================================= 
 Explicit predictor options
 ============================================= */
/*
 #define EGDE_VELOCITIES_IN_MUSCL_STEP
 #define SYMMETRIC_ADVECTION
 #define FOURTH_ORDER_ADV_FLUXES
 #define SYMMETRIC_ADVECTION
 */ 
#define EGDE_VELOCITIES_IN_MUSCL_STEP
#define SYMMETRIC_ADVECTION
#define ADVECTION

/* ============================================= 
 Semi-implicit solver options  
 ============================================= */
/* 
 #define NONLINEAR_EOS_IN_1st_PROJECTION -- Newton for  P(pi)
 #define NO_PI_SYNC                      -- P doesn't overwrite pi each time step
 */
#define NONLINEAR_EOS_IN_1st_PROJECTION

/* ============================================= 
 Elliptic Solver Options
 ============================================= */

/* solver options ==============================
 #define SOLVER_1_CR2      ->  Piotr's Conjugate Residual
 #define SOLVER_1_BICGSTAB
 
 #define SOLVER_2_BICGSTAB

 Currently, a simple diagonal (no gravity) and a column-wise 
 preconditioner (in the gravity direction) are implemented 
 for both projections/linear implicit steps. They are always
 on and selected on the fly depending on whether or not 
 gravity is on or off.
 
 */

#define SOLVER_1_CR2
#define SOLVER_2_BICGSTAB

/* Hydrostatic solver options 
 #define SURFACE_PRESSURE_CORRECTION
 #define FULL_D_HYDRO_CORRECTION
 */
#define SURFACE_PRESSURE_CORRECTION

/* First projection options */
/* 
 #define P1_ALTERNATIVE_STENCIL_WEIGHT 0.125   value for bilinear p-ansatz fcts
 #define P1_ALTERNATIVE_STENCIL_WEIGHT 0.0     value for standard five-point Laplacian
*/
#define PROJECTION1 1              /* switch for first projection should be on "1" normally       */
#define P1_ALTERNATIVE_STENCIL_WEIGHT 0.125 

/* Second projection options */
#define PROJECTION2 1              /* switch for second projection should be on "1" normally      */
#define DIV_CONTROL_LOCAL          /* if def'd, div is controlled in L_\infty, otherwise in L2    */
#define P2_FULL_STENCIL 1.0        /* values: 0.0, 1.0;  0.0 = 5/7pt stencil,  1.0=9/27pt stencil */
#define P2_DIAGONAL_FIVE_POINT 1.0 /* 0.0, 1.0; as above but for node-based Poisson op.           */


/* TODO: Code cleaning / maintainance
 
 1) Make all appearances of "extern ..." disappear except for those of
 User_Data ud;
 double* W0;
 
 2) Whereever the auxiliary double array W0 is used, introduce an external "W0_in_use" flag, so 
    that other routines can query whether or not this aux array is currently occupied.
 
 3) Get rid of the full-size dSol arrays
 
 4) Can I get away without ever computing the advective theta-perturbation evolution,
    just modifying the momentum balance to include the semi-implicit effects?
 
 5) Implement   Hydrostatic_Initial_Pressure()  for 3D.
 
 6) Get rid of "Level[]"s in the  MPV struct.
 
 7) Make sure, the fourth order computation of the advective fluxes is implemented
    compatibly with rigid wall boundary conditions. 
 
 8) Rewrite  slanted_wall_min() and related routines for cross-wall advective flux  
    rhoYu  instead of for rhou to improve compatibility with the pseudo-incompressible
    model
 
 
 */


#include <assert.h> 
