#define NSPEC 1
#define QV    0   /* water vapor */
#define QC    1   /* cloud water */
#define QR    2   /* rain water  */

#define NAUX  2
#define PRES  0
#define BUOY  1

/*
 */
#define ADVECT_THETA_PRIME 0.0

#define PERTURBED_WALL

/* Output parameters */
#define HDFFORMAT
#define OUTPUT_SUBSTEPS_PREDICTOR 1
#define OUTPUT_SUBSTEPS 0
#define OUTPUT_SPLITSTEPS 0


/* ============================================= 
 Elliptic Solver Options
 ============================================= */

/* */ 
#define GRAVITY_IMPLICIT
 

/* solver options ==============================
 #define SOLVER_1_CR2
 #define SOLVER_1_BICGSTAB
 
 #define SOLVER_2_CR2
 #define SOLVER_2_BICGSTAB
 */

#define PROJECTION1 1

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
#define PRECON_DIAGONAL_1ST_PROJ
#define PRECON_DIAGONAL_2ND_PROJ
#define PRECON_LEGACY
#endif

/*
 #define NO_UPWIND_PROJ1          IMPORTANT: Upwinding in the correction induces NOISE  
 */
#define NO_UPWIND_PROJ1

/*
 #define THIRD_ORDER_UPWIND_CORRECTION
 #define SECOND_ORDER_CENTRAL_CORRECTION
 #define FIRST_ORDER_UPWIND_CORRECTION
 */
#define SECOND_ORDER_CENTRAL_CORRECTION

#define PROJECTION2 1
#define DIV_CONTROL_LOCAL
#define P2_FULL_STENCIL 1.0 /* 0.0, 1.0 */
#define P2_DIAGONAL_FIVE_POINT 1.0 /* 0.0, 1.0 */

#include <assert.h> 