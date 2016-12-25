#define NSPEC 1
#define BUOY  0
#define QV    1   /* water vapor */
#define QC    2   /* cloud water */
#define QR    3   /* rain water  */

#define YMOM  1

#define NAUX  2
#define PRES  0
#define SOLD  1


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

#define GRAVITY_IMPLICIT

/* 
 #define GRAVITY_IMPLICIT
 #define GRAVITY_IMPLICIT_1   (buoyancy directly from conserved quantities)
 #define GRAVITY_IMPLICIT_2   (buoyancy via auxiliary variable; Piotr's variant)
 */  
#ifdef GRAVITY_IMPLICIT
#define GRAVITY_IMPLICIT_1
#endif

/* solver options ==============================
 #define SOLVER_1_CR2
 #define SOLVER_1_BICGSTAB
 
 #define SOLVER_2_CR2
 #define SOLVER_2_BICGSTAB

 #define CONTROL_PRECONDITIONED_RESIDUAL_PROJ1
 #define CONTROL_PRECONDITIONED_RESIDUAL_PROJ2
 
*/

#define SOLVER_1_BICGSTAB
#define SOLVER_2_BICGSTAB

/* preconditioning options ======================
 #define PRECON
 #define PRECON_DIAGONAL_1ST_PROJ
 #define PRECON_DIAGONAL_2ND_PROJ
 #define PRECON_VERTICAL_COLUMN_1ST_PROJ
 #define PRECON_VERTICAL_COLUMN_2ND_PROJ

 #define UN_AVERAGE_PROJ1
*/

#define PRECON
#ifdef PRECON
#define PRECON_DIAGONAL_1ST_PROJ
#define PRECON_DIAGONAL_2ND_PROJ
#define PRECON_LEGACY
#endif


/* First projection options */
#define PROJECTION1 1
#define CONTROL_PRECONDITIONED_RESIDUAL_PROJ1

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


/* Second projection options */
#define PROJECTION2 1
#define CONTROL_PRECONDITIONED_RESIDUAL_PROJ2
#define DIV_CONTROL_LOCAL
#define P2_FULL_STENCIL 1.0 /* 0.0, 1.0 */
#define P2_DIAGONAL_FIVE_POINT 1.0 /* 0.0, 1.0 */

#include <assert.h> 