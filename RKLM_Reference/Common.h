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
#define OUTPUT_SUBSTEPS_PREDICTOR 0

/*
#define OUTPUT_FLUXES
 */

/* ============================================= 
 Initial Data Options
 ============================================= */

/* ============================================= 
 Advection options
 ============================================= */
/**/ 
#define EGDE_VELOCITIES_IN_MUSCL_STEP
 

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

#include <assert.h> 
