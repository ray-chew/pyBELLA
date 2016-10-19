/*
#define PERTURBATION_POTENTIAL_TEMPERATURE
 */


#define MAC_PROJ_OLDP_WEIGHT  1.0
#define MAC_PROJ_DELP_WEIGHT  1.0
#define SCND_PROJ_OLDP_WEIGHT 1.0
#define SCND_PROJ_DELP_WEIGHT 1.0

/*
 */
#define THETA_EXP_AT_TIME_LEVEL 0.5 /* 0.5 is theoretical value */

#define NSPEC 1
#define QV    0   /* water vapor */
#define QC    1   /* cloud water */
#define QR    2   /* rain water  */

/* #define LIMITER_FROM_MACROS 
 
 Limiter options:
 
 NONE_LIM
 MINMOD_LIM
 SUPERBEE_LIM
 VANLEER_LIM
 SWEBYMUNZ_LIM
 RUPEZ
 NONE_LIM
 
 */
#define LIMITER_FROM_MACROS

#ifdef LIMITER_FROM_MACROS
#define LIMITER_U VANLEER_LIM
#define LIMITER_V VANLEER_LIM
#define LIMITER_W VANLEER_LIM
#define LIMITER_X VANLEER_LIM
#define LIMITER_Y VANLEER_LIM
#define LIMITER_Z VANLEER_LIM
#endif

/* #define HYDROSTATIC_RHOY_INTERPOLATION 
 #define Y_BUOY_FROM_LEFTS_RIGHTS
 */
#define EDGE_FOCUSED_BUOYANCY

/*  
#define BUOYANCY_VIA_OPSPLIT
#define GRAVITY_1
*/
   
/**/
#define THERMCON

/*
 #define NO_UPWIND_PROJ1          IMPORTANT: Upwinding in the correction induces NOISE  
 */
#define NO_UPWIND_PROJ1

#define PERTURBED_WALL

/* Output parameters */
#define HDFFORMAT
#define MATLAB_OUTPUT
#define OUTPUT_SUBSTEPS_PREDICTOR 1
#define OUTPUT_SUBSTEPS 1
#define OUTPUT_SPLITSTEPS 1


#include <assert.h> 