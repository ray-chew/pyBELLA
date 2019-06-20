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
#include "boundary.h"


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
    const double gm1   = th.gm1;
    const double Gamma_inv = 1.0/Gamma;
    const double gm1_inv   = 1.0/gm1;
    
    const int icy = elem->icy;
    const int igy = elem->igy;
    
    const double rhoY0 = 1.0;
    
    double y_p, y_m;
    double p0, pi0;
    
    double S_integral_p, Sn_integral_p, g;
    double S_m, S_p, Sn_m, Sn_p, pi_hydro, p_hydro, rhoY_hydro, rhoY_hydro_n, pi_hydro_n;
    
    int j;
    
    g = ud.gravity_strength[1];

    /* reference state near the ground */
    p0                       = pow(rhoY0, gamm);
    pi0                      = pow(rhoY0, gm1);
    HydroState_n->rho0[igy]  = rhoY0/Y_n[igy];
    HydroState_n->rhoY0[igy] = rhoY0;
    HydroState_n->Y0[igy]    = Y_n[igy];
    HydroState_n->S0[igy]    = 1.0/Y_n[igy];
    HydroState_n->p0[igy]    = p0;
    HydroState_n->p20[igy]   = pi0/ud.Msq;
    
    /* Hydrostates in bottom dummy cells */
    y_p          = elem->y[igy-1];
    S_m          = 1.0/Y_n[igy];
    S_p          = 1.0/Y[igy-1];
    S_integral_p = - 0.5 * elem->dy * 0.5*(S_p + S_m);

    /* midpoint rule hydrostatic integration for the nodal pressure */
    Sn_p                     = 1.0/Y[igy];
    Sn_integral_p            = 0.0;
    
    for(j = igy-1; j >= 0; j--) {
        
        pi_hydro    = pi0 - Gamma*g*S_integral_p;
        p_hydro     = pow(pi_hydro,Gamma_inv);
        rhoY_hydro  = pow(pi_hydro,gm1_inv);
        
        HydroState->rho0[j]   = rhoY_hydro * S_p;
        HydroState->p0[j]     = p_hydro;
        HydroState->p20[j]    = pi_hydro/ud.Msq;
        HydroState->S0[j]     = S_p;
        HydroState->S10[j]    = 0.0;
        HydroState->Y0[j]     = 1.0/S_p;
        HydroState->rhoY0[j]  = rhoY_hydro;
        
        y_m           = y_p;
        y_p           = y_m - elem->dy;
        S_m           = S_p;
        S_p           = 1.0/Y[MAX_own(j-1, 0)];
        S_integral_p -= 0.5*elem->dy*(S_m + S_p);

        /* midpoint rule hydrostatic integration for the nodal pressure */
        Sn_m           = Sn_p;
        Sn_p           = 1.0/Y[j];
        Sn_integral_p -= elem->dy * Sn_p;
        pi_hydro_n     = pi0 - Gamma*g*Sn_integral_p;
        rhoY_hydro_n   = pow(pi_hydro_n,gm1_inv);
        HydroState_n->rhoY0[j] = rhoY_hydro_n;
        HydroState_n->Y0[j]    = Y_n[j];
        HydroState_n->S0[j]    = 1.0/Y_n[j];
        HydroState_n->p0[j]    = pow(rhoY_hydro_n,th.gamm);
        HydroState_n->p20[j]   = pi_hydro_n/ud.Msq;
    }
    
    /* Hydrostates in bulk of domain */
    y_p          = elem->y[igy];
    S_p          = 1.0/Y[igy];
    S_m          = 1.0/Y_n[igy];
    S_integral_p = 0.5 * elem->dy * 0.5*(S_p + S_m);
    
    /* midpoint rule hydrostatic integration for the nodal pressure */
    Sn_p           = 1.0/Y[igy];
    Sn_integral_p  = 0.0;
    
    for(j = igy; j < icy; j++) {
        
        pi_hydro    = pi0 - Gamma*g*S_integral_p;
        p_hydro     = pow(pi_hydro,Gamma_inv);
        rhoY_hydro  = pow(pi_hydro,gm1_inv);
        
        HydroState->rho0[j]   = rhoY_hydro * S_p;
        HydroState->p0[j]     = p_hydro;
        HydroState->p20[j]    = pi_hydro/ud.Msq;
        HydroState->S0[j]     = S_p;
        HydroState->S10[j]    = 0.0;
        HydroState->Y0[j]     = 1.0/S_p;
        HydroState->rhoY0[j]  = rhoY_hydro;
        
        y_m           = y_p;
        y_p           = y_m + elem->dy;
        S_m           = S_p;
        S_p           = 1.0/Y[j+1];
        S_integral_p += 0.5*elem->dy*(S_m + S_p);

        /* midpoint rule hydrostatic integration for the nodal pressure */
        Sn_m           = Sn_p;
        Sn_p           = 1.0/Y[j];
        Sn_integral_p += elem->dy * Sn_p;
        
        pi_hydro_n     = pi0 - Gamma*g*Sn_integral_p;
        rhoY_hydro_n   = pow(pi_hydro_n,gm1_inv);
        HydroState_n->rhoY0[j+1] = rhoY_hydro_n;
        HydroState_n->Y0[j+1]    = Y_n[j+1];
        HydroState_n->S0[j+1]    = 1.0/Y_n[j+1];
        HydroState_n->p0[j+1]    = pow(rhoY_hydro_n,th.gamm);
        HydroState_n->p20[j+1]   = pi_hydro_n/ud.Msq;
    }    
}

/* ================================================================================== */

void Hydrostatics_State(MPV* mpv, 
                        const ElemSpaceDiscr* elem, 
                        const NodeSpaceDiscr* node) 
{
    
    extern Thermodynamic th;
    extern User_Data ud;
        
    const double Gamma     = th.gm1 / th.gamm;
    const double gamm      = th.gamm;
    const double gm1       = th.gm1;
    const double Gamma_inv = 1.0/Gamma;
    const double gm1_inv   = 1.0/gm1;
        
    const int icy = elem->icy;
    const int igy = elem->igy;
    
    const double rhoY0 = 1.0;
    
    double p0, pi0;
    
    double g;
    double y_p,  y_m,  S_integral_p,  S_m,  S_p; 
    double yn_p, yn_m, Sn_integral_p, Sn_m, Sn_p;
    double p_hydro, pi_hydro, pi_hydro_n, rhoY_hydro, rhoY_hydro_n;
    
    int j;
    
    g = ud.gravity_strength[1];
    
    /* Reference state at the ground */
    p0           = pow(rhoY0, gamm);
    pi0          = pow(rhoY0, gm1);
    mpv->HydroState_n->Y0[igy]    = stratification(0.0);
    mpv->HydroState_n->rhoY0[igy] = rhoY0;
    mpv->HydroState_n->rho0[igy]  = rhoY0/stratification(0.0);
    mpv->HydroState_n->S0[igy]    = 1.0/mpv->HydroState_n->Y0[igy];
    mpv->HydroState_n->p0[igy]    = p0;
    mpv->HydroState_n->p20[igy]   = pi0/ud.Msq;

    /* Hydrostates in bottom dummy cells */
    y_p          = elem->y[igy-1];
    S_p          = 1.0/stratification(y_p);
    S_integral_p = -0.5 * elem->dy * 0.5*(S_p + 1.0/stratification(0.0));

    /* midpoint rule hydrostatic integration for the nodal pressure */
    yn_p          = node->y[igy-1];
    Sn_p          = 1.0/stratification(elem->y[igy-1]);
    Sn_integral_p = -node->dy * Sn_p;
    
    for(j = igy-1; j >= 0; j--) {
        
        pi_hydro    = pi0 - Gamma*g*S_integral_p;
        p_hydro     = pow(pi_hydro, Gamma_inv);
        rhoY_hydro  = pow(pi_hydro, gm1_inv);
        
        mpv->HydroState->rhoY0[j]  = rhoY_hydro;
        mpv->HydroState->rho0[j]   = rhoY_hydro * S_p;
        mpv->HydroState->p0[j]     = p_hydro;
        mpv->HydroState->p20[j]    = pi_hydro/ud.Msq;
        mpv->HydroState->S0[j]     = S_p;
        mpv->HydroState->S10[j]    = 0.0;
        mpv->HydroState->Y0[j]     = 1.0/S_p;
        
        y_m           = y_p;
        y_p           = y_m - elem->dy;
        S_m           = S_p;
        S_p           = 1.0/stratification(y_p);
        S_integral_p -= elem->dy*0.5*(S_m + S_p);
                
        /* midpoint rule hydrostatic integration for the nodal pressure */
        pi_hydro_n    = pi0 - Gamma*g*Sn_integral_p;
        rhoY_hydro_n  = pow(pi_hydro_n, gm1_inv);
        mpv->HydroState_n->rhoY0[j] = rhoY_hydro_n;
        mpv->HydroState_n->Y0[j]    = stratification(0.5*(y_p+y_m));
        mpv->HydroState_n->rho0[j]  = rhoY_hydro_n / mpv->HydroState_n->Y0[j];
        mpv->HydroState_n->S0[j]    = 1.0/mpv->HydroState_n->Y0[j];
        mpv->HydroState_n->p0[j]    = pow(rhoY_hydro_n,th.gamm);
        mpv->HydroState_n->p20[j]   = pi_hydro_n/ud.Msq;
        yn_m           = yn_p;
        yn_p           = yn_m - node->dy;
        Sn_m           = Sn_p;
        Sn_p           = 1.0/stratification(0.5*(yn_p+yn_m));
        Sn_integral_p -= elem->dy * Sn_p;
    }
    
    /* Hydrostates in bulk of domain */
    y_p          = elem->y[igy];
    S_p          = 1.0/stratification(y_p);
    S_integral_p = 0.5 * elem->dy * 0.5*(S_p + 1.0/stratification(0.0));
    
    /* midpoint rule hydrostatic integration for the nodal pressure */
    yn_p          = node->y[igy+1];
    Sn_p          = 1.0/stratification(y_p);
    Sn_integral_p = node->dy * Sn_p;
    
    for(j = igy; j < icy; j++) {
        
        pi_hydro    = pi0 - Gamma*g*S_integral_p;
        p_hydro     = pow(pi_hydro, Gamma_inv);
        rhoY_hydro  = pow(pi_hydro, gm1_inv);
        
        mpv->HydroState->rho0[j]   = rhoY_hydro * S_p;
        mpv->HydroState->p0[j]     = p_hydro;
        mpv->HydroState->p20[j]    = pi_hydro/ud.Msq;
        mpv->HydroState->S0[j]     = S_p;
        mpv->HydroState->S10[j]    = 0.0;
        mpv->HydroState->Y0[j]     = 1.0/S_p;
        mpv->HydroState->rhoY0[j]  = rhoY_hydro;
        
        y_m           = y_p;
        y_p           = y_m + elem->dy;
        S_m           = S_p;
        S_p           = 1.0/stratification(y_p);
        S_integral_p += 0.5*elem->dy*(S_m + S_p);
        
        /* midpoint rule hydrostatic integration for the nodal pressure */
        pi_hydro_n    = pi0 - Gamma*g*Sn_integral_p;
        rhoY_hydro_n  = pow(pi_hydro_n, gm1_inv);
        mpv->HydroState_n->rhoY0[j+1] = rhoY_hydro_n;
        mpv->HydroState_n->Y0[j+1]    = stratification(0.5*(y_p+y_m));
        mpv->HydroState_n->rho0[j+1]  = rhoY_hydro_n / mpv->HydroState_n->Y0[j+1];
        mpv->HydroState_n->S0[j+1]    = 1.0/mpv->HydroState_n->Y0[j+1];
        mpv->HydroState_n->p0[j+1]    = pow(rhoY_hydro_n,th.gamm);
        mpv->HydroState_n->p20[j+1]   = pi_hydro_n/ud.Msq;
        yn_m           = yn_p;
        yn_p           = yn_m + node->dy;
        Sn_m           = Sn_p;
        Sn_p           = 1.0/stratification(0.5*(yn_p+yn_m));
        Sn_integral_p += elem->dy * Sn_p;
    }
    
    
#if OUTPUT_HYDROSTATES
    {
        extern User_Data ud;
        FILE* phydrofile = NULL;
        char fn2[120];
        sprintf(fn2, "%s/hydrostate/model.hse.igw", ud.file_name);
        phydrofile = fopen(fn2, "w+");
        fprintf(phydrofile, "# npts = %d\n", elem->icy-2*elem->igx);
        fprintf(phydrofile, "# num of variables = 4\n");
        fprintf(phydrofile, "# density\n");
        fprintf(phydrofile, "# temperature\n");
        fprintf(phydrofile, "# pressure\n");
        fprintf(phydrofile, "# X\n");
        
        /* Translation to MAESTRO's cgs-units; my basic unit system being SI */
        double L_ratio   = ud.h_ref   * (1.0 / 0.01);  
        double rho_ratio = ud.rho_ref * (1.0e3 * 1.0e-6);
        double p_ratio   = ud.p_ref   * (1.0e+3 / 1.0e+2);
        double R_ratio   = ud.p_ref/(ud.rho_ref*ud.T_ref) * (1.0);
        
        for(int j=elem->igx; j<elem->icy-elem->igx; j++) {
            fprintf(phydrofile, "%10f %10f %10f %10f %10f\n", \
                    elem->y[j] * L_ratio, \
                    mpv->HydroState->rho0[j] * rho_ratio, \
                    mpv->HydroState->p0[j]/mpv->HydroState->rho0[j] * (p_ratio/(R_ratio*rho_ratio)), \
                    mpv->HydroState->p0[j] * p_ratio, 1.0);
        }
        
        fclose(phydrofile);        
    }
#endif
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
    
    double S_integral_p, g;
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
        S_integral_p -= 0.5*dh*(S_m + S_p);        
    }
    
    /* Hydrostatic Exner pressure in bulk of domain */
    S_p          = S[ig];
    S_integral_p = 0.5 * dh * 0.5*(S_p + S0);
    
    for(j = ig; j < n - ig; j++) {
        
        pi[j] = pi0 - Gamma*g*S_integral_p;
        
        S_m           = S_p;
        S_p           = S[j+1];
        S_integral_p += 0.5*dh*(S_m + S_p);
    }
    
    /* Hydrostates in top dummy cells */
    for(j = n-ig; j < n; j++) {
        
        pi[j] = pi0 - Gamma*g*S_integral_p;
        
        S_m           = S_p;
        S_p           = S[MIN_own(n-1, j+1)];
        S_integral_p += 0.5*dh*(S_m + S_p);
    }
}

/* ================================================================================== */

/* midpoint rule used to obtain hydrostatic nodal pressure */
void Hydrostatic_Initial_Pressure(ConsVars* Sol, 
                                  MPV* mpv,
                                  const ElemSpaceDiscr *elem,
                                  const NodeSpaceDiscr *node)
{
    /* 
     TODO: Currently implemented only for 2D vertical slices;
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
    
    beta  = (double*)malloc(node->icx*sizeof(double));
    bdpdx = (double*)malloc(node->icx*sizeof(double));
    pibot = (double*)malloc(node->icx*sizeof(double));
    coeff = (double*)malloc(node->icx*sizeof(double));
    
    int icx = elem->icx;
    int igx = elem->igx;
    int icy = elem->icy;
    int igy = elem->igy;
    
    int icxn  = node->icx;
    int icyn  = node->icy;
    
    double dx = elem->dx;
    double dy = elem->dy;
    
    double Gammainv = th.Gammainv;
    
    /* At this stage, the nodal pressure field contains not nodal values but
     pressure valudes corresponding to the top/bottom cell faces, i.e., 
     those pressure that directly balance the "weight" of air mass in a cell. 
     Nodal values have to be reconstructed from these. I will do this 
     separately from the cell-centered pressures. 
     */
    
    /* cell centered pressures -------------------------------------------------------------- */
    
    /* vertical averages (see docs; Semi-Implicit-Gravity.tex, section   \ref{sec:HydroInit}) */
    memset(beta,0.0,node->icx*sizeof(double));
    memset(bdpdx,0.0,node->icx*sizeof(double));
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
    
    /*reset mpv->p2_cells */
    for (int i=igx; i<icx-igx+1; i++) {
        int nci     = i;
        double pic0 = pibot[i];
        for (int j=igy; j<icy-igy+1; j++) {
            int ncij    = nci + j*icx;
            mpv->p2_cells[ncij] += pic0 - NoBG * mpv->HydroState->p20[j];
        }
    }
    
    set_ghostcells_p2(mpv->p2_cells, elem, 2);    
    
    /* node centered pressures -------------------------------------------------------------- */
    
    /* vertical averages (see docs; Semi-Implicit-Gravity.tex, section   \ref{sec:HydroInit}) */
    /* mpv->p2_nodes thus far contains pressures at the bottom of primary control volumes, 
     not nodal values yet. Here we extract the hydrostatic part of that pressure and pursue 
     the transfer to the nodes.
     */
    memset(beta,0.0,elem->icx*sizeof(double));
    memset(bdpdx,0.0,elem->icx*sizeof(double));
    for (int i=1; i<icxn-1; i++) {
        int ni        = i;
        double height = 0.0;        
        for (int j=igy; j<icyn-igy; j++) {
            int nij    = ni + j*icx;   
            int nnij   = ni + j*icxn;
            double Pc  = Sol->rhoY[nij];
            double thc = Sol->rhoY[nij]/Sol->rho[nij];
            beta[i]   += Pc*thc*dy;
            bdpdx[i]  += Pc*thc*(mpv->p2_nodes[nnij]-mpv->p2_nodes[nnij-1])*dy;
            height    += dy;
        }
        beta[i]  *= Gammainv/height;
        bdpdx[i] *= Gammainv/height/dx;
    }
    
    /* integrate in x */
    memset(coeff,0.0,node->icx*sizeof(double));
    memset(pibot,0.0,node->icx*sizeof(double));
    for (int i=igx+1; i<icxn-igx; i++) {
        coeff[i] = coeff[i-1] + dx/beta[i];
        pibot[i] = pibot[i-1] - dx*bdpdx[i]/beta[i];
    }
    
    /* determine integration constant for periodicity of the surface pressure */
    dotPU = pibot[icx-igx] / coeff[icx-igx];
    
    /*finalize bottom pressure distribution */
    for (int i=igx; i<icx-igx; i++) {
        pibot[i] -= dotPU*coeff[i];
    }
    
    /*reset mpv->p2_nodes */
    for (int i=igx; i<icx-igx+1; i++) {
        int nni     = i;
        for (int j=igy; j<icy-igy+1; j++) {
            int nnij    = nni + j*icxn;
            mpv->p2_nodes[nnij] += pibot[i] - NoBG * mpv->HydroState_n->p20[j];
        }
    }
    
    /* transfer to the nodes */
    for (int nn=0; nn < node->nc; nn++) {
        mpv->dp2_nodes[nn] = mpv->p2_nodes[nn];
    }

    for (int j=igy; j<icyn-igy; j++) {
        int nnj = j*icxn;
        int sgn;
        double delp2;
        
        /* first guess for left boundary value */
        mpv->p2_nodes[nnj+node->igx] = mpv->dp2_nodes[nnj+node->igx];
        
        /* in unravelling the averaged bottom cell face values to bottom nodal
         values, we start with a recursion from the left 
         */
        for (int i=igx+1; i<icxn-igx; i++) {
            int nnij = nnj + i; 
            mpv->p2_nodes[nnij] = 2.0 * mpv->dp2_nodes[nnij-1] - mpv->p2_nodes[nnij-1];
        }
        
        /* this gives us a new value at the right edge, and now
         we have to enforce periodicity by adding a suitable multiple 
         of the (+1,-1)-mode: This will work only for uneven node numbers
         in the horizontal direction. */
        assert((node->icx+1)%2 == 1);
        delp2 = 0.5*(mpv->p2_nodes[nnj+icxn-igx-1] - mpv->p2_nodes[nnj+igx]);
        sgn   = 1;
        for (int i=igx; i<icxn-igx; i++) {
            int nnij = nnj + i;
            mpv->p2_nodes[nnij] += sgn*delp2;
            sgn *= -1;
        }
    }
    
    set_ghostnodes_p2(mpv->p2_nodes, node, 2);
    
    for (int nn=0; nn<node->nc; nn++) {
        mpv->dp2_nodes[nn] = 0.0;
    }
    
    if (ud.is_compressible) {
        for (int i=igx; i<icx-igx; i++) {
            int nci     = i;
            for (int j=igy; j<icy-igy; j++) {
                int ncij      = nci + j*icx;
                double pi     = ud.Msq * (mpv->p2_cells[ncij] + NoBG * mpv->HydroState->p20[j]);
                double Y      = Sol->rhoY[ncij]/Sol->rho[ncij];
                double rhoold = Sol->rho[ncij];
                Sol->rhoY[ncij] = pow(pi, th.gm1inv);
                Sol->rho[ncij]  = Sol->rhoY[ncij] / Y;
                Sol->rhou[ncij]*= Sol->rho[ncij] / rhoold;
                Sol->rhov[ncij]*= Sol->rho[ncij] / rhoold;
                Sol->rhow[ncij]*= Sol->rho[ncij] / rhoold;
                Sol->rhoe[ncij]*= Sol->rho[ncij] / rhoold;
                for (int nsp=0; nsp<ud.nspec; nsp++) {
                    Sol->rhoX[nsp][ncij]*= Sol->rho[ncij] / rhoold;
                }
            }
        }
    }
    
    free(beta);
    free(bdpdx);
    free(pibot);
    free(coeff);
}

