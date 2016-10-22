//
//  Hydrostatics.c
//  RKLM_Reference
//
//  Created by Klein, Rupert on 19/10/16.
//  Copyright © 2016 Klein, Rupert. All rights reserved.
//

#include <math.h>

#include "Hydrostatics.h"
#include "thermodynamic.h"
#include "userdata.h"
#include "variable.h"


/* ================================================================================== */

void Hydrostatics_State(MPV* mpv, double *Yinvbg, const ElemSpaceDiscr* elem) {
    
    extern Thermodynamic th;
    extern User_Data ud;
        
    const double Gamma = th.gm1 / th.gamm;
    const double gamm  = th.gamm;
    const double Gamma_inv = 1.0/Gamma;
        
    const int icy = elem->icy;
    const int igy = elem->igy;
    
    const double rho0 = 1.0;
    
    double y_p, y_m;
    double p0;
    
    double S_integral_m, S_integral_n, S_integral_p, g;
    double S_m, S_p, p_hydro, rhoY_hydro, rhoY_hydro_n;
    
    int j;
    
    g = ud.gravity_strength[1];
    
    /* Hydrostates in bottom dummy cells */
    y_p          = elem->y[igy];
    S_p          = 1.0/stratification(y_p);
    S_integral_p = 0.5 * elem->dy * 0.5*(S_p + 1.0/stratification(0.0));
    p0           = pow(rho0 * stratification(0.0), gamm);
    
    for(j = igy; j >= 0; j--) {
        
        p_hydro     = pow( pow(p0,Gamma) - Gamma*g*S_integral_p ,  Gamma_inv);
        rhoY_hydro   = pow(p_hydro,1.0/th.gamm);
        
        mpv->HydroState->geopot[j] = g * y_p;
        mpv->HydroState->rho0[j]   = rhoY_hydro * S_p;
        mpv->HydroState->p0[j]     = p_hydro;
        mpv->HydroState->p20[j]    = p_hydro/ud.Msq;
        mpv->HydroState->S0[j]     = S_p;
        mpv->HydroState->S10[j]    = 0.0;
        mpv->HydroState->Y0[j]     = 1.0/S_p;
        mpv->HydroState->rhoY0[j]  = rhoY_hydro;
        
        y_m           = y_p;
        y_p           = y_m - elem->dy;
        S_m           = S_p;
        S_p           = 1.0/stratification(y_p);
        S_integral_m  = S_integral_p;
        S_integral_p -= 0.5*elem->dy*(S_m + S_p);
        S_integral_n  = 0.5*(S_integral_m + S_integral_p);
        
        rhoY_hydro_n  = pow( pow(p0,Gamma) - Gamma*g*S_integral_n ,  th.gm1inv);
        mpv->HydroState_n->rhoY0[j] = rhoY_hydro_n;
    }
    
    /* Hydrostates in bulk of domain */
    y_p          = elem->y[igy];
    S_p          = 1.0/stratification(y_p);
    S_integral_p = 0.5 * elem->dy * 0.5*(S_p + 1.0/stratification(0.0));

    for(j = igy; j < icy - igy; j++) {
        
        p_hydro     = pow( pow(p0,Gamma) - Gamma*g*S_integral_p ,  Gamma_inv);
        rhoY_hydro  = pow(p_hydro,1.0/th.gamm);
        
        mpv->HydroState->geopot[j] = g * y_p;
        mpv->HydroState->rho0[j]   = rhoY_hydro * S_p;
        mpv->HydroState->p0[j]     = p_hydro;
        mpv->HydroState->p20[j]    = p_hydro/ud.Msq;
        mpv->HydroState->S0[j]     = S_p;
        mpv->HydroState->S10[j]    = 0.0;
        mpv->HydroState->Y0[j]     = 1.0/S_p;
        mpv->HydroState->rhoY0[j]  = rhoY_hydro;
        
        y_m           = y_p;
        y_p           = y_m + elem->dy;
        S_m           = S_p;
        S_p           = 1.0/stratification(y_p);
        S_integral_m  = S_integral_p;
        S_integral_p += 0.5*elem->dy*(S_m + S_p);
        S_integral_n  = 0.5*(S_integral_m + S_integral_p);
        
        rhoY_hydro_n  = pow( pow(p0,Gamma) - Gamma*g*S_integral_n ,  th.gm1inv);
        mpv->HydroState_n->rhoY0[j+1] = rhoY_hydro_n;
        
    }
    
    /* Hydrostates in top dummy cells */
    for(j = icy-igy; j < icy; j++) {
        
        p_hydro     = pow( pow(p0,Gamma) - Gamma*g*S_integral_p ,  Gamma_inv);
        rhoY_hydro  = pow(p_hydro,1.0/th.gamm);
        
        mpv->HydroState->geopot[j] = g * y_p;
        mpv->HydroState->rho0[j]   = rhoY_hydro * S_p;
        mpv->HydroState->p0[j]     = p_hydro;
        mpv->HydroState->p20[j]    = p_hydro/ud.Msq;
        mpv->HydroState->S0[j]     = S_p;
        mpv->HydroState->S10[j]    = 0.0;
        mpv->HydroState->Y0[j]     = 1.0/S_p;
        mpv->HydroState->rhoY0[j]  = rhoY_hydro;
        
        y_m           = y_p;
        y_p           = y_m + elem->dy;
        S_m           = S_p;
        S_p           = 1.0/stratification(y_p);
        S_integral_m  = S_integral_p;
        S_integral_p += 0.5*elem->dy*(S_m + S_p);
        S_integral_n  = 0.5*(S_integral_m + S_integral_p);
        
        rhoY_hydro_n  = pow( pow(p0,Gamma) - Gamma*g*S_integral_n ,  th.gm1inv);
        mpv->HydroState_n->rhoY0[j+1] = rhoY_hydro_n;
    }
    
    /* HydroStates in an auxiliary field */
    /* */
    for (int k = 0; k < elem->icz; k++) {
        int nk = k*elem->icx*elem->icy;
        
        for (int j = 0; j < elem->icy; j++) {
            int njk = nk + j*elem->icx;
            
            for (int i = 0; i < elem->icx; i++) {
                int nijk = njk + i;
                
                Yinvbg[nijk] = mpv->HydroState->S0[j];
            }
        }
    }

}