#include "userdata.h"

/* ================================================================================ */

void set_time_integrator_parameters(User_Data *ud)
{
    
    /* Versions prior to Mar 4, 2018 included settings for more integrators 
    switch (ud->time_integrator) {
        case SI_MIDPT:
            ud->tips.dt_frac         = 0.5;
            ud->tips.flux_frac[0][0] = 0.0;
            ud->tips.flux_frac[0][1] = 1.0;
            ud->tips.flux_frac[1][0] = 0.5;
            ud->tips.flux_frac[1][1] = 0.5;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.flux_frac[k][1] = 1e10;
                ud->tips.flux_frac[k][1] = 1e10;
            }
            ud->tips.update_frac[0]  = 0.5;
            ud->tips.update_frac[1]  = 1.0;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.update_frac[k] = 1e10;
            }
            ud->tips.multiD_updt     = WRONG;
            break;
            
        default:
            assert(0);
            break;
    }
     */ 
    switch (ud->advec_time_integrator) {
      case EXPL_MIDPT:
            ud->tips.dt_frac         = 0.5;
            ud->tips.flux_frac[0][0] = 0.0;
            ud->tips.flux_frac[0][1] = 1.0;
            ud->tips.flux_frac[1][0] = -0.5;
            ud->tips.flux_frac[1][1] = 1.0;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.flux_frac[k][1] = 1e10;
                ud->tips.flux_frac[k][1] = 1e10;
            }
            ud->tips.update_frac[0]  = 0.5;
            ud->tips.update_frac[1]  = 1.0;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.update_frac[k] = 1e10;
            }
            ud->tips.multiD_updt     = CORRECT;
            break;
            
        case HEUN:
            ud->tips.dt_frac         = 0.5;
            ud->tips.flux_frac[0][0] = 0.0;
            ud->tips.flux_frac[0][1] = 1.0;
            ud->tips.flux_frac[1][0] = -0.5;
            ud->tips.flux_frac[1][1] = 0.5;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.flux_frac[k][1] = 1e10;
                ud->tips.flux_frac[k][1] = 1e10;
            }
            ud->tips.update_frac[0]  = 1.0;
            ud->tips.update_frac[1]  = 1.0;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.update_frac[k] = 1e10;
            }
            ud->tips.multiD_updt     = CORRECT;
            break;
            
        default: /* STRANG */
            ud->tips.dt_frac         = 0.5;
            ud->tips.flux_frac[0][0] = 0.0;
            ud->tips.flux_frac[0][1] = 1.0;
            ud->tips.flux_frac[1][0] = 0.5;
            ud->tips.flux_frac[1][1] = 0.5;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.flux_frac[k][1] = 1e10;
                ud->tips.flux_frac[k][1] = 1e10;
            }
            ud->tips.update_frac[0]  = 0.5;
            ud->tips.update_frac[1]  = 1.0;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.update_frac[k] = 1e10;
            }
            ud->tips.multiD_updt     = WRONG;
            break;
    }
}

