#define NSPEC 1
#define BUOY  0
#define QV    1   /* water vapor */
#define QC    2   /* cloud water */
#define QR    3   /* rain water  */

#define YMOM  1

#define NAUX  2
#define PRES  0

#define PERTURBED_WALL

/* Output parameters */
#define HDFFORMAT
#define OUTPUT_SUBSTEPS 0
#define OUTPUT_SPLITSTEPS 0
#define OUTPUT_SUBSTEPS_PREDICTOR 0


/* ============================================= 
 Initial Data Options
 ============================================= */
/*
 */
#define HYDRO_BALANCED_INIT_DATA

/* ============================================= 
 Elliptic Solver Options
 ============================================= */

#define GRAVITY_IMPLICIT

/* solver options ==============================
 #define SOLVER_1_CR2
 #define SOLVER_1_BICGSTAB
 
 #define SOLVER_2_CR2
 #define SOLVER_2_BICGSTAB
 */

#define SOLVER_1_BICGSTAB
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
/* #define PRECON_LEGACY  */
#endif


/* First projection options */
#define PROJECTION1 1
#define TIME_AVERAGED_COEFFS_PROJ1 

/* Second projection options */
#define PROJECTION2 1
#define DIV_CONTROL_LOCAL
#define P2_FULL_STENCIL 1.0 /* 0.0, 1.0 */
#define P2_DIAGONAL_FIVE_POINT 1.0 /* 0.0, 1.0 */

#include <assert.h> 
