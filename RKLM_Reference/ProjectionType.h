/* solver options ==============================
 #define SOLVER_1_CR2
 #define SOLVER_1_BICGSTAB

 #define SOLVER_2_CR2
 #define SOLVER_2_BICGSTAB
 #define SOLVER_2_HYPRE
 #define SOLVER_2_HYPRE_RUPE
 #define SOLVER_2_HYPRE_MICHAEL

 #define SOLVER_HYPRE

 #define MICHAEL_CORRECTION

 #define SHIFTED_COEFFICIENTS_PROJ1

*/

/*         
 #define SECOND_ORDER_ADVECTING_FLUX 1
  */ 

#define PROJECTION1 1
#define SHIFTED_COEFFICIENTS_PROJ1

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
 #define STANDARD_STENCIL_PROJ1
 */


/*
 #define DP_AVERAGE      (used and needed with first GRESHO CFL 0.99 run on 400x400
                          dp2 filtered only before added to Z)
                         -> no longer needed when I run with diagonal 9-pt exact 
                            projection in the second correction step
*/


#define P_AVERAGE  /* Option seems definitely needed when running with gravity; 
                      see Klemp-Skamarock-internal wave test without mean advection */

/*
 #define TEST_LAPLACIAN_PROJ1
 #define TEST_SOLVER
 */

/*
 #define THIRD_ORDER_UPWIND_CORRECTION
 #define SECOND_ORDER_CENTRAL_CORRECTION
 #define FIRST_ORDER_UPWIND_CORRECTION
 */
#define SECOND_ORDER_CENTRAL_CORRECTION


/* cell-centered velocity projection 
 #define SECOND_CORRECTION_IS_EXACT
*/
#define PROJECTION2 1
#define DIV_CONTROL_LOCAL
#define P2_FULL_STENCIL 1.0 /* 0.0, 1.0 */
#define P2_DIAGONAL_FIVE_POINT 1.0 /* 0.0, 1.0 */

/*
 #define HYP_UPDATE_ETA 1
 */