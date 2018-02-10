//
//  Hydrostatics.c
//  RKLM_Reference
//
//  Created by Klein, Rupert on 19/10/16.
//  Copyright © 2016 Klein, Rupert. All rights reserved.
//

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Hydrostatics.h"
#include "thermodynamic.h"
#include "userdata.h"
#include "variable.h"
#include "math_own.h"


/* ================================================================================== */

void Hydrostatics_Column(States* HydroState, 
                         States* HydroState_n, 
                         double* Y, 
                         double* Y_n, 
                         const ElemSpaceDiscr* elem, 
                         const NodeSpaceDiscr* node) 
{
    
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
    
    double S_integral_m, S_integral_p, Sn_integral_m, Sn_integral_p, g;
    double S_m, S_p, Sn_m, Sn_p, p_hydro, rhoY_hydro, rhoY_hydro_n;
    
    int j;
    
    g = ud.gravity_strength[1];

    /* reference state near the ground */
    p0                       = pow(rho0 * Y_n[igy], gamm);
    HydroState_n->rho0[igy]  = rho0;
    HydroState_n->rhoY0[igy] = rho0 * Y_n[igy];
    HydroState_n->Y0[igy]    = Y_n[igy];
    HydroState_n->S0[igy]    = 1.0/Y_n[igy];
    HydroState_n->p0[igy]    = p0;
    HydroState_n->p20[igy]   = p0/ud.Msq;
    Sn_p                     = 1.0/Y_n[igy];
    Sn_integral_p            = 0.0;
    
    /* Hydrostates in bottom dummy cells */
    y_p          = elem->y[igy-1];
    S_m          = 1.0/Y_n[igy];
    S_p          = 1.0/Y[igy-1];
    S_integral_p = - 0.5 * elem->dy * 0.5*(S_p + S_m);

    for(j = igy-1; j >= 0; j--) {
        
        p_hydro     = pow( pow(p0,Gamma) - Gamma*g*S_integral_p ,  Gamma_inv);
        rhoY_hydro  = pow(p_hydro,1.0/th.gamm);
        
        HydroState->geopot[j] = g * y_p;
        HydroState->rho0[j]   = rhoY_hydro * S_p;
        HydroState->p0[j]     = p_hydro;
        HydroState->p20[j]    = p_hydro/ud.Msq;
        HydroState->S0[j]     = S_p;
        HydroState->S10[j]    = 0.0;
        HydroState->Y0[j]     = 1.0/S_p;
        HydroState->rhoY0[j]  = rhoY_hydro;
        
        y_m           = y_p;
        y_p           = y_m - elem->dy;
        S_m           = S_p;
        S_p           = 1.0/Y[MAX_own(j-1, 0)];
        S_integral_m  = S_integral_p;
        S_integral_p -= 0.5*elem->dy*(S_m + S_p);
        
        Sn_m           = Sn_p;
        Sn_p           = 1.0/Y_n[j];
        Sn_integral_m  = Sn_integral_p;
        Sn_integral_p -= 0.5*elem->dy*(Sn_m + Sn_p);
        rhoY_hydro_n   = pow( pow(p0,Gamma) - Gamma*g*Sn_integral_p ,  th.gm1inv);
        HydroState_n->rhoY0[j] = rhoY_hydro_n;
        HydroState_n->Y0[j]    = Y_n[j];
        HydroState_n->S0[j]    = 1.0/Y_n[j];
        HydroState_n->p0[j]    = pow(rhoY_hydro_n,th.gamm);
        HydroState_n->p20[j]   = HydroState_n->p0[j]/ud.Msq;
    }
    
    /* Hydrostates in bulk of domain */
    Sn_p           = 1.0/Y_n[igy];
    Sn_integral_p  = 0.0;

    y_p          = elem->y[igy];
    S_p          = 1.0/Y[igy];
    S_m          = 1.0/Y_n[igy];
    S_integral_p = 0.5 * elem->dy * 0.5*(S_p + S_m);
    
    for(j = igy; j < icy; j++) {
        
        p_hydro     = pow( pow(p0,Gamma) - Gamma*g*S_integral_p ,  Gamma_inv);
        rhoY_hydro  = pow(p_hydro,1.0/th.gamm);
        
        HydroState->geopot[j] = g * y_p;
        HydroState->rho0[j]   = rhoY_hydro * S_p;
        HydroState->p0[j]     = p_hydro;
        HydroState->p20[j]    = p_hydro/ud.Msq;
        HydroState->S0[j]     = S_p;
        HydroState->S10[j]    = 0.0;
        HydroState->Y0[j]     = 1.0/S_p;
        HydroState->rhoY0[j]  = rhoY_hydro;
        
        y_m           = y_p;
        y_p           = y_m + elem->dy;
        S_m           = S_p;
        S_p           = 1.0/Y[j+1];
        S_integral_m  = S_integral_p;
        S_integral_p += 0.5*elem->dy*(S_m + S_p);

        Sn_m           = Sn_p;
        Sn_p           = 1.0/Y_n[j+1];
        Sn_integral_m  = Sn_integral_p;
        Sn_integral_p += 0.5*elem->dy*(Sn_m + Sn_p);
        rhoY_hydro_n  = pow( pow(p0,Gamma) - Gamma*g*Sn_integral_p ,  th.gm1inv);
        HydroState_n->rhoY0[j+1] = rhoY_hydro_n;
        HydroState_n->Y0[j+1]    = Y_n[j+1];
        HydroState_n->S0[j+1]    = 1.0/Y_n[j+1];
        HydroState_n->p0[j+1]    = pow(rhoY_hydro_n,th.gamm);
        HydroState_n->p20[j+1]   = HydroState_n->p0[j+1]/ud.Msq;
    }
    
    /* Hydrostates in top dummy cells
    for(j = icy-igy; j < icy; j++) {
        
        p_hydro     = pow( pow(p0,Gamma) - Gamma*g*S_integral_p ,  Gamma_inv);
        rhoY_hydro  = pow(p_hydro,1.0/th.gamm);
        
        HydroState->geopot[j] = g * y_p;
        HydroState->rho0[j]   = rhoY_hydro * S_p;
        HydroState->p0[j]     = p_hydro;
        HydroState->p20[j]    = p_hydro/ud.Msq;
        HydroState->S0[j]     = S_p;
        HydroState->S10[j]    = 0.0;
        HydroState->Y0[j]     = 1.0/S_p;
        HydroState->rhoY0[j]  = rhoY_hydro;
        
        y_m           = y_p;
        y_p           = y_m + elem->dy;
        S_m           = S_p;
        S_p           = 1.0/Y[j+1];
        S_integral_m  = S_integral_p;
        S_integral_p += 0.5*elem->dy*(S_m + S_p);
        S_integral_n  = 0.5*(S_integral_m + S_integral_p);
        
        rhoY_hydro_n  = pow( pow(p0,Gamma) - Gamma*g*S_integral_n ,  th.gm1inv);
        HydroState_n->rhoY0[j+1] = rhoY_hydro_n;
        HydroState_n->Y0[j]      = Y_n[j];
        HydroState_n->S0[j+1]    = 1.0/Y_n[j];
        HydroState_n->p0[j+1]    = pow(rhoY_hydro_n,th.gamm);
    }   
      */
}

/* ================================================================================== */

void Hydrostatics_State(MPV* mpv, 
                        double *Sbg, 
                        const ElemSpaceDiscr* elem, 
                        const NodeSpaceDiscr* node) 
{
    
    extern Thermodynamic th;
    extern User_Data ud;
        
    const double Gamma = th.gm1 / th.gamm;
    const double gamm  = th.gamm;
    const double Gamma_inv = 1.0/Gamma;
        
    const int icy = elem->icy;
    const int igy = elem->igy;
    
    const double rho0 = 1.0;
    
    double p0;
    
    double g;
    double y_p,  y_m,  S_integral_m,  S_integral_p,  S_m,  S_p; 
    double yn_p, yn_m, Sn_integral_m, Sn_integral_p, Sn_m, Sn_p;
    double p_hydro, p_hydro_n, rhoY_hydro, rhoY_hydro_n;
    
    int j;
    
    g = ud.gravity_strength[1];
    
    /* Reference state at the ground */
    p0           = pow(rho0 * stratification(0.0), gamm);
    mpv->HydroState_n->Y0[igy]    = stratification(0.0);
    mpv->HydroState_n->rho0[igy]  = rho0;
    mpv->HydroState_n->rhoY0[igy] = rho0*stratification(0.0);
    mpv->HydroState_n->S0[igy]    = 1.0/mpv->HydroState_n->Y0[igy];
    mpv->HydroState_n->p0[igy]    = p0;
    mpv->HydroState_n->p20[igy]   = mpv->HydroState_n->p0[igy]/ud.Msq;

    /* Hydrostates in bottom dummy cells */
    y_p          = elem->y[igy-1];
    S_p          = 1.0/stratification(y_p);
    S_integral_p = -0.5 * elem->dy * 0.5*(S_p + 1.0/stratification(0.0));

    yn_p          = node->y[igy-1];
    Sn_p          = 1.0/stratification(yn_p);
    Sn_integral_p = -node->dy * 0.5*(Sn_p + 1.0/stratification(0.0));

    for(j = igy-1; j >= 0; j--) {
        
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
        y_p           = y_m - elem->dy;
        S_m           = S_p;
        S_p           = 1.0/stratification(y_p);
        S_integral_m  = S_integral_p;
        S_integral_p -= elem->dy*0.5*(S_m + S_p);
                
        p_hydro_n     = pow( pow(p0,Gamma) - Gamma*g*Sn_integral_p ,  Gamma_inv);
        rhoY_hydro_n  = pow(p_hydro_n,1.0/th.gamm);
        mpv->HydroState_n->rhoY0[j] = rhoY_hydro_n;
        mpv->HydroState_n->Y0[j]    = stratification(0.5*(y_p+y_m));
        mpv->HydroState_n->rho0[j]  = rhoY_hydro_n / mpv->HydroState_n->Y0[j];
        mpv->HydroState_n->S0[j]    = 1.0/mpv->HydroState_n->Y0[j];
        mpv->HydroState_n->p0[j]    = p_hydro_n;
        mpv->HydroState_n->p20[j]   = mpv->HydroState_n->p0[j]/ud.Msq;

        yn_m           = yn_p;
        yn_p           = yn_m - node->dy;
        Sn_m           = Sn_p;
        Sn_p           = 1.0/stratification(yn_p);
        Sn_integral_m  = Sn_integral_p;
        Sn_integral_p -= elem->dy*0.5*(Sn_m + Sn_p);
    }
    
    /* Hydrostates in bulk of domain */
    y_p          = elem->y[igy];
    S_p          = 1.0/stratification(y_p);
    S_integral_p = 0.5 * elem->dy * 0.5*(S_p + 1.0/stratification(0.0));
    
    yn_p          = node->y[igy+1];
    Sn_p          = 1.0/stratification(yn_p);
    Sn_integral_p = node->dy * 0.5*(Sn_p + 1.0/stratification(0.0));

    for(j = igy; j < icy; j++) {
        
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
        
        p_hydro_n     = pow( pow(p0,Gamma) - Gamma*g*Sn_integral_p ,  Gamma_inv);
        rhoY_hydro_n  = pow(p_hydro_n,1.0/th.gamm);
        mpv->HydroState_n->rhoY0[j+1] = rhoY_hydro_n;
        mpv->HydroState_n->Y0[j+1]    = stratification(0.5*(y_p+y_m));
        mpv->HydroState_n->rho0[j+1]  = rhoY_hydro_n / mpv->HydroState_n->Y0[j+1];
        mpv->HydroState_n->S0[j+1]    = 1.0/mpv->HydroState_n->Y0[j+1];
        mpv->HydroState_n->p0[j+1]    = pow(rhoY_hydro_n,th.gamm);
        mpv->HydroState_n->p20[j+1]   = mpv->HydroState_n->p0[j+1]/ud.Msq;

        yn_m           = yn_p;
        yn_p           = yn_m + elem->dy;
        Sn_m           = Sn_p;
        Sn_p           = 1.0/stratification(yn_p);
        Sn_integral_m  = Sn_integral_p;
        Sn_integral_p += 0.5*elem->dy*(Sn_m + Sn_p);
}
    
    /* Hydrostates in top dummy cells 
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
        mpv->HydroState_n->Y0[j+1]    = stratification(0.5*(y_p+y_m));
        mpv->HydroState_n->S0[j+1]    = 1.0/mpv->HydroState_n->Y0[j+1];
        mpv->HydroState_n->p0[j+1]    = pow(rhoY_hydro_n,th.gamm);
        mpv->HydroState_n->p20[j+1]   = mpv->HydroState_n->p0[j+1]/ud.Msq;
    }
     */
    
    /* HydroStates in an auxiliary field */
    /* */
    for (int k = 0; k < elem->icz; k++) {
        int nk = k*elem->icx*elem->icy;
        
        for (int j = 0; j < elem->icy; j++) {
            int njk = nk + j*elem->icx;
            
            for (int i = 0; i < elem->icx; i++) {
                int nijk = njk + i;
                
                Sbg[nijk] = mpv->HydroState->S0[j];
            }
        }
    }
}

/* ================================================================================== */

void Hydrostatic_Exner_pressure(
                                double *pi, 
                                const double pi0, 
                                const double *S, 
                                const double S0,
                                const double dh,
                                const int n, 
                                const int ig) 
{
    
    extern Thermodynamic th;
    extern User_Data ud;
    
    const double Gamma = th.gm1 / th.gamm;
    
    double S_integral_m, S_integral_p, g;
    double S_m, S_p;
    
    int j;
    
    g = ud.gravity_strength[1];
    
    /* Hydrostatic Exner pressure in bottom dummy cells */
    S_p          = S[ig];
    S_integral_p = 0.5 * dh * 0.5*(S_p + S0);
    
    for(j = ig; j >= 0; j--) {
        
        pi[j]         = pi0 - Gamma*g*S_integral_p;
                
        S_m           = S_p;
        S_p           = S[MAX_own(0,j-1)];
        S_integral_m  = S_integral_p;
        S_integral_p -= 0.5*dh*(S_m + S_p);        
    }
    
    /* Hydrostatic Exner pressure in bulk of domain */
    S_p          = S[ig];
    S_integral_p = 0.5 * dh * 0.5*(S_p + S0);
    
    for(j = ig; j < n - ig; j++) {
        
        pi[j] = pi0 - Gamma*g*S_integral_p;
        
        S_m           = S_p;
        S_p           = S[j+1];
        S_integral_m  = S_integral_p;
        S_integral_p += 0.5*dh*(S_m + S_p);
    }
    
    /* Hydrostates in top dummy cells */
    for(j = n-ig; j < n; j++) {
        
        pi[j] = pi0 - Gamma*g*S_integral_p;
        
        S_m           = S_p;
        S_p           = S[MIN_own(n-1, j+1)];
        S_integral_m  = S_integral_p;
        S_integral_p += 0.5*dh*(S_m + S_p);
    }
}

/* ================================================================================== */

void Hydrostatic_Initial_Pressure(ConsVars* Sol, 
                                  MPV* mpv,
                                  const ElemSpaceDiscr *elem,
                                  const NodeSpaceDiscr *node)
{
    /* 
     Computes the pressure field corresponding to the linear hydrostatic
     and pseudo-incompressible approximation for a vertical slice model 
     and x-periodic conditions. 
       The routine assumes that pi = mpv->p2_cells already contains a 
     column-wise hydrostatically balanced Exner pressure field that is 
     horizontally homogeneous at zero height, i.e.,
     
                     pi = int_0^z  Gamma g / theta dz'
     
     The problem can be reduced to an explicit x-integration and the
     choice of one proper integration constant to establish periodicity
     of the surface pressure. 
     */
    extern User_Data ud;
    extern Thermodynamic th;
     
    double *beta, *bdpdx, *pibot, *coeff;
    double dotPU;
    double NoBG = (ud.time_integrator == SI_MIDPT ? 1.0 : 0.0);
    
    assert(elem->ndim == 2);
    
    beta  = (double*)malloc(elem->icx*sizeof(double));
    bdpdx = (double*)malloc(elem->icx*sizeof(double));
    pibot = (double*)malloc(elem->icx*sizeof(double));
    coeff = (double*)malloc(elem->icx*sizeof(double));
    
    int icx = elem->icx;
    int igx = elem->igx;
    int icy = elem->icy;
    int igy = elem->igy;

    int inx = node->icx;

    double dx  = elem->dx;
    double dy  = elem->dy;
    
    double Gammainv = th.Gammainv;
    
    /* vertical averages (see docs; Semi-Implicit-Gravity.tex, section \ref{sec:HydroInit}) */
    memset(beta,0.0,elem->icx*sizeof(double));
    memset(bdpdx,0.0,elem->icx*sizeof(double));
    for (int i=1; i<icx-1; i++) {
        int ni        = i;
        double height = 0.0;        
        for (int j=igy; j<icy-igy; j++) {
            int nij    = ni + j*icx;
            double Pc  = Sol->rhoY[nij];
            double Pm  = Sol->rhoY[nij-1];
            double thc = Sol->rhoY[nij]/Sol->rho[nij];
            double thm = Sol->rhoY[nij-1]/Sol->rho[nij-1];
            beta[i]   += 0.5*(Pm*thm+Pc*thc)*dy;
            bdpdx[i]  += 0.5*(Pm*thm+Pc*thc)*(mpv->p2_cells[nij]-mpv->p2_cells[nij-1])*dy;
            height    += dy;
        }
        beta[i]  *= Gammainv/height;
        bdpdx[i] *= Gammainv/height/dx;
    }

    /* integrate in x */
    memset(pibot,0.0,elem->icx*sizeof(double));
    memset(coeff,0.0,elem->icx*sizeof(double));
    for (int i=igx+1; i<icx-igx+1; i++) {
        coeff[i] = coeff[i-1] + dx/beta[i];
        pibot[i] = pibot[i-1] - dx*bdpdx[i]/beta[i];
    }
    
    /* determine integration constant for periodicity of the surface pressure */
    dotPU = pibot[icx-igx] / coeff[icx-igx];
    
    /*finalize bottom pressure distribution */
    for (int i=igx; i<icx-igx; i++) {
        pibot[i] -= dotPU*coeff[i];
    }
    
    /*reset mpv->p2_cells, mpv->p2_nodes */
    for (int i=igx; i<icx-igx+1; i++) {
        int nci     = i;
        int nni     = i;
        double pic0 = pibot[i];
        double pin0 = 0.5*(pibot[i]+pibot[i-1]);
        for (int j=igy; j<icy-igy+1; j++) {
            int ncij    = nci + j*icx;
            int nnij    = nni + j*inx;
            mpv->p2_cells[ncij] += pic0 - NoBG * mpv->HydroState->p20[j]/mpv->HydroState->rhoY0[j];
            mpv->p2_nodes[nnij] += pin0 - NoBG * mpv->HydroState_n->p20[j]/mpv->HydroState_n->rhoY0[j];
        }
    }
    
    free(beta);
    free(bdpdx);
    free(pibot);
    free(coeff);
}

