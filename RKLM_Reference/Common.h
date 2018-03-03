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


/* ============================================= 
 Initial Data Options
 ============================================= */

/* ============================================= 
 Advection options
 ============================================= */
/* 
#define EGDE_VELOCITIES_IN_MUSCL_STEP
 #define HALF_STEP_FLUX_EXTERNAL
 */

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
 */

#define SOLVER_1_CR2
#define SOLVER_2_BICGSTAB

/* preconditioning options ======================
 #define PRECON
 #define PRECON_DIAGONAL_1ST_PROJ
 #define PRECON_DIAGONAL_2ND_PROJ
 #define PRECON_VERTICAL_COLUMN_1ST_PROJ
 #define PRECON_VERTICAL_COLUMN_2ND_PROJ
*/

#define PRECON
#ifdef PRECON
#define PRECON_VERTICAL_COLUMN_1ST_PROJ
#define PRECON_VERTICAL_COLUMN_2ND_PROJ
#endif


/* First projection options */
/*
 #define TIME_AVERAGED_COEFFS_PROJ1 
*/

#define PROJECTION1 1

/* Second projection options */
#define PROJECTION2 1
#define DIV_CONTROL_LOCAL
#define P2_FULL_STENCIL 1.0 /* 0.0, 1.0 */
#define P2_DIAGONAL_FIVE_POINT 1.0 /* 0.0, 1.0 */

#include <assert.h> 
