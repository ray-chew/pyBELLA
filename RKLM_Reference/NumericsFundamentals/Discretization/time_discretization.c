#include "userdata.h"

/* ================================================================================ */

void set_time_integrator_parameters(User_Data *ud)
{
    
    /* set parameters for time integration scheme (so far only for the predictor step) */ 
    switch (ud->time_integrator) {
        case OP_SPLIT:
            ud->tips.dt_frac         = 0.5;
            ud->tips.flux_frac[0][0] = 0.0;
            ud->tips.flux_frac[0][1] = 1.0;
            ud->tips.flux_frac[1][0] = 0.5;
            ud->tips.flux_frac[1][1] = 0.5;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.flux_frac[k][1] = 1e10;
                ud->tips.flux_frac[k][1] = 1e10;
            }
            ud->tips.no_of_stages    = 2;
            ud->tips.update_frac[0]  = 0.5;
            ud->tips.update_frac[1]  = 1.0;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.update_frac[k] = 1e10;
            }
            ud->tips.multiD_updt     = WRONG;
            break;
        case OP_SPLIT_MD_UPDATE:
            ud->tips.dt_frac         = 0.5;
            ud->tips.flux_frac[0][0] = 0.0;
            ud->tips.flux_frac[0][1] = 1.0;
            ud->tips.flux_frac[1][0] = 0.5;
            ud->tips.flux_frac[1][1] = 0.5;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.flux_frac[k][1] = 1e10;
                ud->tips.flux_frac[k][1] = 1e10;
            }
            ud->tips.no_of_stages    = 1;
            ud->tips.update_frac[0]  = 0.5;
            ud->tips.update_frac[1]  = 1.0;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.update_frac[k] = 1e10;
            }
            ud->tips.multiD_updt     = CORRECT;
            break;
        case HEUN:
            ud->tips.dt_frac         = 0.0;
            ud->tips.flux_frac[0][0] = 0.0;
            ud->tips.flux_frac[0][1] = 1.0;
            ud->tips.flux_frac[1][0] = 0.5;
            ud->tips.flux_frac[1][1] = 0.5;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.flux_frac[k][1] = 1e10;
                ud->tips.flux_frac[k][1] = 1e10;
            }
            ud->tips.no_of_stages    = 2;
            ud->tips.update_frac[0]  = 1.0;
            ud->tips.update_frac[1]  = 1.0;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.update_frac[k] = 1e10;
            }
            ud->tips.multiD_updt     = CORRECT;
            break;
        case EXPL_MIDPT:
            ud->tips.dt_frac         = 0.0;
            ud->tips.flux_frac[0][0] = 0.0;
            ud->tips.flux_frac[0][1] = 1.0;
            ud->tips.flux_frac[1][0] = 0.0;
            ud->tips.flux_frac[1][1] = 1.0;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.flux_frac[k][1] = 1e10;
                ud->tips.flux_frac[k][1] = 1e10;
            }
            ud->tips.no_of_stages    = 2;
            ud->tips.update_frac[0]  = 0.5;
            ud->tips.update_frac[1]  = 1.0;
            for (int k=2; k<NO_OF_RK_STAGES; k++) {
                ud->tips.update_frac[k] = 1e10;
            }
            ud->tips.multiD_updt     = CORRECT;
            break;
        case RK3_SKAMA:
            ud->tips.dt_frac         = 0.0;
            ud->tips.flux_frac[0][0] = 0.0;
            ud->tips.flux_frac[0][1] = 1.0;
            ud->tips.flux_frac[1][0] = 0.0;
            ud->tips.flux_frac[1][1] = 1.0;
            ud->tips.flux_frac[2][0] = 0.0;
            ud->tips.flux_frac[2][1] = 1.0;
            for (int k=3; k<NO_OF_RK_STAGES; k++) {
                ud->tips.flux_frac[k][1] = 1e10;
                ud->tips.flux_frac[k][1] = 1e10;
            }
            ud->tips.no_of_stages    = 3;
            ud->tips.update_frac[0]  = 0.3333333333;
            ud->tips.update_frac[1]  = 0.5;
            ud->tips.update_frac[2]  = 1.0;
            for (int k=3; k<NO_OF_RK_STAGES; k++) {
                ud->tips.update_frac[k] = 1e10;
            }
            ud->tips.multiD_updt     = CORRECT;
            break;
        case RK3_TEST:
            ud->tips.dt_frac         = 0.0;
            ud->tips.flux_frac[0][0] = 0.0;
            ud->tips.flux_frac[0][1] = 1.0;
            ud->tips.flux_frac[1][0] = 0.0;
            ud->tips.flux_frac[1][1] = 1.0;
            ud->tips.flux_frac[2][0] = 0.0;
            ud->tips.flux_frac[2][1] = 1.0;
            for (int k=3; k<NO_OF_RK_STAGES; k++) {
                ud->tips.flux_frac[k][1] = 1e10;
                ud->tips.flux_frac[k][1] = 1e10;
            }
            ud->tips.no_of_stages    = 3;
            ud->tips.update_frac[0]  = 0.5;
            ud->tips.update_frac[1]  = 0.5;
            ud->tips.update_frac[2]  = 1.0;
            for (int k=3; k<NO_OF_RK_STAGES; k++) {
                ud->tips.update_frac[k] = 1e10;
            }
            ud->tips.multiD_updt     = CORRECT;
            break;
        default:
            assert(0);
            break;
    }
    
}

