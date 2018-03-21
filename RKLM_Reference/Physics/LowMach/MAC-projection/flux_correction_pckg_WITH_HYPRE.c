/*******************************************************************************
 File:   flux_correction.c
 Author: Rupert
 Date:   Feb  14  2004
 *******************************************************************************/
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "Common.h"
#include "mpv.h"
#include "memory.h"
#include "flux_correction.h"
#include "variable_coefficient_poisson_cells.h"
#include "BICGSTAB.h"
#include "variable.h"
#include "warning.h"
#include "enumerator.h"
#include "error.h"
#include "io.h"
#include "userdata.h"
#include "Eos.h"
#include "thermodynamic.h"
#include "math_own.h"
#include "recovery.h"
#include "laplacian_cells.h"
#include "boundary.h"

static enum Constraint integral_condition(ConsVars* flux[3],
                                          double* rhs, 
                                          ConsVars* Sol,
                                          const double dt,
                                          const ElemSpaceDiscr* elem,
                                          MPV* mpv);

static void rhs_fix_for_open_boundaries(
                                        double* rhs, 
                                        const ElemSpaceDiscr* elem, 
                                        ConsVars* Sol_new, 
                                        ConsVars* Sol_old, 
                                        ConsVars* flux[3],
                                        double dt,
                                        MPV* mpv);

static void flux_fix_for_open_boundaries(
                                         ConsVars* flux[3],
                                         const ElemSpaceDiscr* elem, 
                                         MPV* mpv);

static void flux_correction_due_to_pressure_gradients(
                                                      ConsVars* flux[3],
                                                      const ElemSpaceDiscr* elem,
                                                      ConsVars* Sol,
                                                      ConsVars* Sol0,
                                                      const MPV* mpv, 
                                                      double* hplus[3],
                                                      double* hS,
                                                      double* dp2,
                                                      const double t,
                                                      const double dt);

void operator_coefficients(double* hplus[3], 
                           double* wcenter, 
                           double* hS, 
                           const ElemSpaceDiscr* elem,
                           const ConsVars* Sol,
                           const ConsVars* Sol0,
                           const MPV* mpv,
                           const double dt); 

#ifdef NONLINEAR_EOS_IN_1st_PROJECTION
double Newton_rhs(double* rhs,
                  const double* dpi,
                  const double* pi,
                  const double* rhoY0,
                  const MPV* mpv,
                  const ElemSpaceDiscr* elem);

#endif

/* ========================================================================== */

#define RHS_OUTPUT 0
#if RHS_OUTPUT
static int rhs_output_count = 0;
#endif

void flux_correction(ConsVars* flux[3],
                     const ElemSpaceDiscr* elem_base_grid,
                     ConsVars* Sol, 
                     ConsVars* Sol0, 
                     const double t,
                     const double dt,
                     const int step) {
    
    extern User_Data ud;
    extern MPV* mpv;
    
    const ElemSpaceDiscr* elem = mpv->Level[0]->elem;
    const NodeSpaceDiscr* node = mpv->Level[0]->node;
    
    double** hplus       = mpv->Level[0]->wplus;
    double*  hcenter     = mpv->Level[0]->wcenter;
    double*  hS          = mpv->Level[0]->wgrav;
    
    double* rhs          = mpv->Level[0]->rhs;
    double* p2		     = mpv->Level[0]->p;
    double* pi           = (double*)malloc(elem->nc*sizeof(double));
    
    double rhsmax;
    
    int n;
    
    for (int nc=0; nc<elem->nc; nc++) pi[nc] = 0.0;
    
    printf("\n\n====================================================");
    printf("\nFirst Projection");
    printf("\n====================================================\n");
    
    operator_coefficients(hplus, hcenter, hS, elem, Sol, Sol0, mpv, dt);
    
    rhsmax = controlled_variable_flux_divergence(rhs, (const ConsVars**)flux, dt, elem);
    printf("\nrhs_max = %e (before projection)\n", rhsmax);
    
    assert(integral_condition(flux, rhs, Sol, dt, elem, mpv) != VIOLATED); 
    if (ud.is_compressible) {
        for (int nc=0; nc<elem->nc; nc++) {
            rhs[nc] += hcenter[nc]*mpv->p2_cells[nc];
        }
    }
    
    rhs_fix_for_open_boundaries(rhs, elem, Sol, Sol0, flux, dt, mpv);
    
#if RHS_OUTPUT
    extern User_Data ud;
    FILE *prhsfile = NULL;
    char fn2[200], fieldname2[90];
    
    sprintf(fn2, "%s/rhs_cells/rhs_cells_00%d.hdf", ud.file_name, rhs_output_count);
    sprintf(fieldname2, "rhs_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             rhs,
             fn2,
             fieldname2);

    rhs_output_count++;
#endif
    
    
    /* Newton iteration for the nonlinear equation of state rhoY = P(\pi) = \pi^{gamma-1} 
     The first sweep computes the half time level Exner pressure based on the linearization
     of P(.). The further iterates yield updates to this quantity until convergence 
     tolerance is reached.
     */
    variable_coefficient_poisson_cells(p2, rhs, (const double **)hplus, hcenter, Sol, elem, node);
    set_ghostcells_p2(p2, elem, elem->igx);
    
#ifdef NONLINEAR_EOS_IN_1st_PROJECTION
    if (ud.is_compressible) {
        /* Nonlinear iteration for the \partial P/\partial t  term in the P-equation */
        double precision_factor = 1.0;
        ud.flux_correction_precision       *= precision_factor;
        ud.flux_correction_local_precision *= precision_factor;
        for (int nc=0; nc<elem->nc; nc++) {
            pi[nc]  = p2[nc];
            p2[nc] -= mpv->p2_cells[nc];        
        }
        rhsmax = Newton_rhs(rhs, (const double*)p2, (const double*)pi, (const double*)Sol->rhoY, (const MPV*)mpv, (const ElemSpaceDiscr*)elem);
        
        while (rhsmax > ud.flux_correction_precision) {
            variable_coefficient_poisson_cells(p2, rhs, (const double **)hplus, hcenter, Sol, elem, node);
            set_ghostcells_p2(p2, elem, elem->igx);
            for (int nc=0; nc<elem->nc; nc++) pi[nc] += p2[nc];
            rhsmax = Newton_rhs(rhs, (const double*)p2, (const double*)pi, (const double*)Sol->rhoY, (const MPV*)mpv, (const ElemSpaceDiscr*)elem);
        }
        
        for (int nc=0; nc<elem->nc; nc++) {
            p2[nc] = pi[nc];        
        }
        ud.flux_correction_precision       /= precision_factor;
        ud.flux_correction_local_precision /= precision_factor;
    }
#endif
    
    
    
    /* Note: flux will contain only the flux-correction after this routine; 
     it is thus overwritten under the SI_MIDPT time integration sequence */
    flux_correction_due_to_pressure_gradients(flux, elem, Sol, Sol0, mpv, hplus, hS, p2, t, dt);
    
    /* test whether divergence is actually controlled */
    rhsmax = controlled_variable_flux_divergence(rhs, (const ConsVars**)flux, dt, elem);
    printf("\nrhs_max = %e (after projection)\n", rhsmax);
    
    flux_fix_for_open_boundaries(flux, elem, mpv);  
    
#if RHS_OUTPUT
    extern double *W0;
    double* lap = W0;
    
    sprintf(fn2, "%s/rhs_cells/rhs_cells_00%d.hdf", ud.file_name, rhs_output_count);
    sprintf(fieldname2, "rhs_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             rhs,
             fn2,
             fieldname2);
    
    EnthalpyWeightedLap_bilinear_p(elem, node,p2, (const double**)hplus, hcenter, Sol, mpv, dt, lap);
    
    sprintf(fn2, "%s/lap_cells/lap_cells_000.hdf", ud.file_name);
    sprintf(fieldname2, "lap_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             lap,
             fn2,
             fieldname2);
    
    rhs_output_count++;

#endif

    /* store results in mpv-fields */
    for(n=0; n<elem->nc; n++) {
        mpv->dp2_cells[n] = p2[n] - mpv->p2_cells[n];
        mpv->p2_cells[n]  = p2[n];
    }
    
    set_ghostcells_p2(mpv->p2_cells, elem, elem->igx);
    Set_Explicit_Boundary_Data(Sol, elem);
    
    free(pi);
}



/* ========================================================================== */

double controlled_variable_flux_divergence(double* rhs, 
                                           const ConsVars* flux[3],
                                           const double dt, 
                                           const ElemSpaceDiscr* elem)
{
    /* right hand side of pressure equation via advective flux divergence */
    memset(rhs, 0.0, elem->nc*sizeof(double));
    
    const int igx = elem->igx;
    const int icx = elem->icx;
    const int ifx = elem->ifx;
    const int igy = elem->igy;
    const int icy = elem->icy;
    const int ify = elem->ify;
    const int igz = elem->igz;
    const int icz = elem->icz;
    const int ifz = elem->ifz;
    
    const double factor = 2.0 / dt;
    const double dx = elem->dx;
    const double dy = elem->dy;
    const double dz = elem->dz;
    
    double rhsmax = 0.0;
    
    for(int i=0; i<elem->nc; i++) rhs[i] = 0.0;
    
    switch (elem->ndim) {
        case 1:
            for(int i = igx; i < icx - igx; i++) {
                int nc  = i;
                int nfx = i;
                rhs[nc] = factor * (flux[0]->rhoY[nfx+1] - flux[0]->rhoY[nfx])/dx;
                rhsmax  = MAX_own(rhsmax, fabs(rhs[nc]));
            }
            break;
            
        case 2:
            for(int j = igy; j < icy - igy; j++) {
                int mc  = j * icx;
                int mfx = j * ifx; 
                int mfy = j;
                for(int i = igx; i < icx - igx; i++) {
                    int nc  = mc  + i;
                    int nfx = mfx + i;
                    int nfy = mfy + i*ify;
                    rhs[nc] = factor * ((flux[0]->rhoY[nfx+1] - flux[0]->rhoY[nfx])/dx + (flux[1]->rhoY[nfy+1] - flux[1]->rhoY[nfy])/dy);
                    rhsmax  = MAX_own(rhsmax, fabs(rhs[nc]));
                }
            }
            break;
            
        case 3:
            for(int k = igz; k < icz - igz; k++) {
                int lc = k*icx*icy;
                int lfx = k*ifx*icy;
                int lfy = k*ify;
                int lfz = k;
                for(int j = igy; j < icy - igy; j++) {
                    int mc  = lc  + j*icx;
                    int mfx = lfx + j*ifx; 
                    int mfy = lfy + j;
                    int mfz = lfz + j*ifz*icx;
                    for(int i = igx; i < icx - igx; i++) {
                        int nc  = mc  + i;
                        int nfx = mfx + i;
                        int nfy = mfy + i*ify*icz;
                        int nfz = mfz + i*ifz;
                        rhs[nc] = factor * (  (flux[0]->rhoY[nfx+1] - flux[0]->rhoY[nfx])/dx \
                                            + (flux[1]->rhoY[nfy+1] - flux[1]->rhoY[nfy])/dy \
                                            + (flux[2]->rhoY[nfz+1] - flux[2]->rhoY[nfz])/dz);
                        rhsmax  = MAX_own(rhsmax, fabs(rhs[nc]));
                    }
                }
            }
            break;
            
        default:
            break;
    }
    
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "rhs_cells.hdf");
    sprintf(fieldname, "rhs-cells");    
    WriteHDF(prhsfile, elem->icx, elem->icy, elem->icz, elem->ndim, rhs, fn, fieldname);
#endif
    
    return rhsmax;
}

/* ========================================================================== */

static enum Constraint integral_condition(
                                          ConsVars* flux[3],
                                          double* rhs,
                                          ConsVars* Sol,
                                          const double dt,
                                          const ElemSpaceDiscr* elem,
                                          MPV* mpv) {
    
    const int ndim = elem->ndim;
    
    switch(ndim) {
            
        case 1: {
            ERROR("function not available");
            break;
        }
            
        case 2: {
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            
            double tmp = 0.0;
            double rhs_max = 0.0;
            double threshold = 1000*DBL_EPSILON;
            int cnt = 0;
            int i,j,m,n;
            for(j = igy; j < icy - igy; j++) {m = j*icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    tmp += rhs[n];
                    rhs_max = MAX_own(fabs(rhs[n]),rhs_max);
                    cnt++;
                }
            }
            tmp /= cnt;
            if(fabs(tmp) < threshold) {
                printf("integral_condition_cells = OK; tmp = %e\n", fabs(tmp));
                return(SATISFIED);
            }
            else {
                printf("integral_condition_cells = VIOLATED; tmp = %e\n", fabs(tmp));
                return(SATISFIED);
            }
            
            break;
        }
            
        case 3: {
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int igz = elem->igz;
            const int icz = elem->icz;
            
            double tmp = 0.0;
            double rhs_max = 0.0;
            double threshold = 1000*DBL_EPSILON;
            int cnt = 0;
            int i,j,k,l,m,n;
            for(k = igz; k < icz - igz; k++) {l = k*icx*icy;
                for(j = igy; j < icy - igy; j++) {m = + + j*icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        tmp += rhs[n];
                        rhs_max = MAX_own(fabs(rhs[n]),rhs_max);
                        cnt++;
                    }
                }
            }
            tmp /= cnt;
            if(fabs(tmp) < threshold) {
                printf("integral_condition_cells = OK; tmp = %e\n", fabs(tmp));
                return(SATISFIED);
            }
            else {
                printf("integral_condition_cells = VIOLATED; tmp = %e\n", fabs(tmp));
                return(VIOLATED);
            }
            
            break;
        }
            
        default: ERROR("ndim not in {1, 2, 3}");
            return(VIOLATED);
            
    }
    
    return(VIOLATED);
}

/* ========================================================================== */

void operator_coefficients(
                           double* hplus[3], 
                           double* wcenter,
                           double* hS,
                           const ElemSpaceDiscr* elem,
                           const ConsVars* Sol,
                           const ConsVars* Sol0,
                           const MPV* mpv, 
                           const double dt) {
    
    extern User_Data ud;
    extern Thermodynamic th;
    
    const double Gammainv = th.Gammainv;
    const int ndim = elem->ndim;
    
    const double ccw = (ud.time_integrator == SI_MIDPT ? 4.0 : 2.0); 
    /* when p2 is perturbation pressure:
     const double ccenter = - ccw * (ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt); 
     const double cexp    = 1.0-th.gamm;
     */
    /* when p2 is Exner pressure: */
    const double ccenter = - ccw * (ud.compressibility*ud.Msq)*th.gm1inv/(mpv->dt*mpv->dt); 
    const double cexp    = 2.0-th.gamm;
    
    switch(ndim) {
        case 1: {    
            ERROR("function not available");
            break;
        }
        case 2: {
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int ifx = elem->ifx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int ify = elem->ify;
            const double dy = elem->dy;
            double* hx = hplus[0];
            double* hy = hplus[1];
            double* hc = wcenter;
            
            double hi, him, hj, hjm, g, gimp, Msq;
            
            int i, j, m, n, ic, icm, jc, jcm;
                        
            for(j = igy-1; j < icy - igy+1; j++) {
                m = j * ifx;
                
                for(i = igx-1; i < ifx - igx+1; i++) {
                    n     = m + i;
                    ic    = n - j;
                    icm   = ic - 1; 
#if 0
                    hi    = Sol0->rhoY[ic] * Sol0->rhoY[ic] / Sol0->rho[ic]    * Gammainv;   
                    him   = Sol0->rhoY[icm] * Sol0->rhoY[icm] / Sol0->rho[icm] * Gammainv;
#else
                    hi    = Sol0->rhoY[ic] * Sol->rhoY[ic] / Sol->rho[ic]    * Gammainv;   
                    him   = Sol0->rhoY[icm] * Sol->rhoY[icm] / Sol->rho[icm] * Gammainv;
#endif
                    hx[n] = 0.5 * (hi + him);
                    
                    assert(hx[n] > 0.0);
                }
            }
            
            Msq = ud.Msq;
            g   = ud.gravity_strength[1];
            
            for(i = 0; i < icx; i++) {
                n = i * ify;
                
                for(j = igy-1; j < ify - igy+1; j++) {
                    m     = n + j;
                    jc    = j * icx + i;
                    jcm   = jc - icx;          
#if 0
                    hj    = Sol0->rhoY[jc] * Sol0->rhoY[jc] / Sol0->rho[jc] * Gammainv;
                    hjm   = Sol0->rhoY[jcm] * Sol0->rhoY[jcm] / Sol0->rho[jcm] * Gammainv;
#else
                    hj    = Sol0->rhoY[jc] * Sol->rhoY[jc] / Sol->rho[jc] * Gammainv;
                    hjm   = Sol0->rhoY[jcm] * Sol->rhoY[jcm] / Sol->rho[jcm] * Gammainv;
#endif
                    double S     = mpv->HydroState->S0[j];
                    double Sm    = mpv->HydroState->S0[j-1];
                    double Y     = 0.5 * (Sol->rhoY[jc]  / Sol->rho[jc]  + Sol0->rhoY[jc]  / Sol0->rho[jc]);
                    double Ym    = 0.5 * (Sol->rhoY[jcm] / Sol->rho[jcm] + Sol0->rhoY[jcm] / Sol0->rho[jcm]);
                    double Nsq   = - (g/Msq) * 0.5*(Y+Ym) * (S-Sm)/dy;
                    double Nsqsc = 0.25*dt*dt*Nsq;
                    
                    gimp  = 1.0 / (1.0 + Nsqsc);
                    
                    hy[m] = 0.5 * (hj + hjm) * gimp;
                    hS[m] = Nsqsc * gimp * Gammainv / (g/Msq) / Y;
                    
                    assert(hy[m] > 0.0);
                }
            }
            
            for(j = igy; j < icy - igy; j++) {m = j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    /*
                     hc[n] = ccenter * pow(Sol0->rhoY[n],cexp);
                     */
                    hc[n] = ccenter * pow(Sol->rhoY[n],cexp);
                }
            }
            
            
            break;
        }
        case 3: {
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int ifx = elem->ifx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int ify = elem->ify;
            const int igz = elem->igz;
            const int icz = elem->icz;
            const int ifz = elem->ifz;
            const double dy = elem->dy;
            double* hx = hplus[0];
            double* hy = hplus[1];
            double* hz = hplus[2];
            double* hc = wcenter;
            
            double hi, him, hj, hjm, hk, hkm, g, gimp, Msq;
            
            int i, j, k, l, m, n, ic, icm, jc, jcm, kc, kcm;
                        
            for(k = igz; k < icz - igz; k++) {l = k * ifx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
                    for(i = igx; i < ifx - igx; i++) {n = m + i;
                        ic  = k*icx*icy + j*icx + i;
                        icm = ic - 1; 
#if 0
                        hi    = (Sol->rhoY[ic] *Sol->rhoY[ic] /Sol->rho[ic] ) * Gammainv;   
                        him   = (Sol->rhoY[icm]*Sol->rhoY[icm]/Sol->rho[icm]) * Gammainv;
#else
                        hi    = (Sol0->rhoY[ic] *Sol->rhoY[ic] /Sol->rho[ic] ) * Gammainv;   
                        him   = (Sol0->rhoY[icm]*Sol->rhoY[icm]/Sol->rho[icm]) * Gammainv;
#endif
                        hx[n] = 0.5 * (hi + him);
                        assert(hx[n] > 0.0);
                    }
                }
            }
            
            Msq = ud.Msq;
            g   = ud.gravity_strength[1];
            
            for(i = igx; i < icx - igx; i++) {l = i * ify*icz;
                for(k = igz; k < icz - igz; k++) {m = l + k * ify;
                    for(j = igy; j < ify - igy; j++) {n = m + j;
                        jc  = k*icx*icy + j*icx + i;
                        jcm = jc - icx;          
#if 0
                        hj       = (Sol->rhoY[jc] *Sol->rhoY[jc] /Sol->rho[jc] ) * Gammainv;
                        hjm      = (Sol->rhoY[jcm]*Sol->rhoY[jcm]/Sol->rho[jcm]) * Gammainv;
#else
                        hj       = (Sol0->rhoY[jc] *Sol->rhoY[jc] /Sol->rho[jc] ) * Gammainv;
                        hjm      = (Sol0->rhoY[jcm]*Sol->rhoY[jcm]/Sol->rho[jcm]) * Gammainv;
#endif
                        double S   = mpv->HydroState->S0[j];
                        double Sm  = mpv->HydroState->S0[j-1];
                        double Y   = 0.5 * (Sol->rhoY[jc]  / Sol->rho[jc]  + Sol0->rhoY[jc]  / Sol0->rho[jc]);
                        double Ym  = 0.5 * (Sol->rhoY[jcm] / Sol->rho[jcm] + Sol0->rhoY[jcm] / Sol0->rho[jcm]);
                        double Nsq = - (g/Msq) * 0.5*(Y+Ym) * (S-Sm)/dy;
                        
                        gimp  = 1.0 / (1.0 + 0.5*dt*dt*Nsq);
                        
                        hy[n] = 0.5 * (hj + hjm) * gimp;
                        
                        assert(hy[n] > 0.0);
                    }
                }
            }
                        
            for(j = igy; j < icy - igy; j++) {l = j * ifz*icx;
                for(i = igx; i < icx - igx; i++) {m = l + i * ifz;
                    for(k = igz; k < ifz - igz; k++) {n = m + k;
                        kc  = k*icx*icy + j*icx + i;
                        kcm = kc - icx*icy;          
#if 0                        
                        hk       = (Sol->rhoY[kc] *Sol->rhoY[kc] /Sol->rho[kc] ) * Gammainv;
                        hkm      = (Sol->rhoY[kcm]*Sol->rhoY[kcm]/Sol->rho[kcm]) * Gammainv;
#else
                        hk       = (Sol0->rhoY[kc] *Sol->rhoY[kc] /Sol->rho[kc] ) * Gammainv;
                        hkm      = (Sol0->rhoY[kcm]*Sol->rhoY[kcm]/Sol->rho[kcm]) * Gammainv;
#endif                        
                        hz[n] = 0.5 * (hk + hkm);
                        assert(hz[n] > 0.0); 
                    }
                }
            }
            
            for(k = igz; k < icz - igz; k++) {l = k*icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j*icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        /*
                        hc[n] = ccenter * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp);
                        hc[n] = ccenter * pow(Sol0->rhoY[n],cexp);
                         */
                        hc[n] = ccenter * pow(Sol->rhoY[n],cexp);
                    }
                }
            }
            
            
            break;
        }
        default: ERROR("ndim not in {1,2,3}");
    }
}

/* ========================================================================== */

static void flux_correction_due_to_pressure_gradients(
                                                      ConsVars* flux[3],
                                                      const ElemSpaceDiscr* elem,
                                                      ConsVars* Sol,
                                                      ConsVars* Sol0,
                                                      const MPV* mpv, 
                                                      double* hplus[3],
                                                      double* hS,
                                                      double* dp2,
                                                      const double t,
                                                      const double dt) {
    
    extern User_Data ud;
    
    const int ndim = elem->ndim;
    
    switch(ndim) {
        case 1: {
            ERROR("function not available");
            break;
        }
            
        case 2: {
            
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int ifx = elem->ifx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int ify = elem->ify;
            
            const double dto2dx = 0.5 * dt / elem->dx;
            const double dto2dy = 0.5 * dt / elem->dy;

            const double b = P1_ALTERNATIVE_STENCIL_WEIGHT; 
            const double a = 1.0-2.0*b;

            ConsVars* f = flux[0];
            ConsVars* g = flux[1];
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            
            for(int j = igy; j < icy - igy; j++) {
                int mc    = j * ifx;                
                for(int i = igx; i < ifx - igx; i++) {
                    int nn   = mc + i + ifx;
                    int nc   = mc + i;
                    int ns   = mc + i - ifx;
                    int ic   = j*icx + i;
                    int icm  = ic - 1;
                    
                    /* It seems this should be the active version, but test against alternative 
                     for IGW before removing this option, looking especially for vertical velocity
                     in the compressible case towards the end of the standard test run. */
                    f->rhoY[nc] -= dto2dx * (  a *   hplusx[nc] * (dp2[ic]     - dp2[icm]    )  
                                             + b * ( hplusx[nc] * (dp2[ic+icx] - dp2[icm+icx])  
                                                   + hplusx[nc] * (dp2[ic-icx] - dp2[icm-icx])  
                                                   ));
                }
            }  
            
            /* fluxes in the y-direction */
            for(int i = igx; i < icx - igx; i++) {
                int nc = i * ify;
                for(int j = igy; j < ify - igy; j++) {
                    int me  = nc + j + ify;
                    int mc  = nc + j;
                    int mw  = nc + j - ify;
                    int jc  = j * icx + i;
                    int jcm = jc - icx;
                    
                    g->rhoY[mc]  -= dto2dy * (  a *   hplusy[mc] * (dp2[jc]   - dp2[jcm]  ) 
                                              + b * ( hplusy[mc] * (dp2[jc+1] - dp2[jcm+1]) 
                                                    + hplusy[mc] * (dp2[jc-1] - dp2[jcm-1]) 
                                                    ));                    
                }   
            }
            
            break;
        }
        case 3: {
            
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int ifx = elem->ifx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int ify = elem->ify;
            const int igz = elem->igz;
            const int icz = elem->icz;
            const int ifz = elem->ifz;
            
            const double dto2dx = 0.5 * dt / elem->dx;
            const double dto2dy = 0.5 * dt / elem->dy;
            const double dto2dz = 0.5 * dt / elem->dz;
            
            const int dix = 1;
            const int diy = icx; 
            const int diz = icx*icy;            
            
            const double a  = 36.0/64.0;
            const double b  =  6.0/64.0;
            const double c  =  1.0/64.0;
            /*
            const double a  =  4.0/16.0;
            const double b  =  2.0/16.0;
            const double c  =  1.0/16.0;
             */
            
            
            ConsVars* fx = flux[0];
            ConsVars* fy = flux[1];
            ConsVars* fz = flux[2];
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            const double* hplusz   = hplus[2];
            
            /* fluxes in the x-direction */
            for(int k = igz; k < icz - igz; k++) {
                int l = k * ifx*icy;
                for(int j = igy; j < icy - igy; j++) {
                    int m = l + j * ifx;
                    for(int i = igx; i < ifx - igx; i++) {
                        int n = m + i;
                        int ic   = k*diz + j*diy + i*dix;
                        int icm  = ic - dix;
                        
                        fx->rhoY[n] -= dto2dx * hplusx[n] * 
                        (     a *  (dp2[ic] - dp2[icm])  
                         +    b * (  (dp2[ic+diy] - dp2[icm+diy]) 
                                   + (dp2[ic-diy] - dp2[icm-diy])  
                                   + (dp2[ic+diz] - dp2[icm+diz])  
                                   + (dp2[ic-diz] - dp2[icm-diz])
                                   )
                         +    c * (  (dp2[ic+diy+diz] - dp2[icm+diy+diz]) 
                                   + (dp2[ic-diy+diz] - dp2[icm-diy+diz])  
                                   + (dp2[ic+diy-diz] - dp2[icm+diy-diz])  
                                   + (dp2[ic-diy-diz] - dp2[icm-diy-diz])
                                   ));   
                    }
                } 
            }
            
            /* fluxes in the y-direction */
            for(int i = igx; i < icx - igx; i++) {
                int l = i * ify*icz;
                for(int k = igz; k < icz - igz; k++) {
                    int m = l + k * ify;
                    for(int j = igy; j < ify - igy; j++) {
                        int n = m + j;
                        int jc   = k*diz + j*diy + i*dix;
                        int jcm  = jc - diy;
                        
                        fy->rhoY[n] -= dto2dy * hplusy[n] * 
                        (     a *    (dp2[jc] - dp2[jcm])  
                         +    b * (  (dp2[jc+diz] - dp2[jcm+diz]) 
                                   + (dp2[jc-diz] - dp2[jcm-diz])  
                                   + (dp2[jc+dix] - dp2[jcm+dix])  
                                   + (dp2[jc-dix] - dp2[jcm-dix])
                                   )
                         +    c * (  (dp2[jc+diz+dix] - dp2[jcm+diz+dix]) 
                                   + (dp2[jc-diz+dix] - dp2[jcm-diz+dix])  
                                   + (dp2[jc+diz-dix] - dp2[jcm+diz-dix])  
                                   + (dp2[jc-diz-dix] - dp2[jcm-diz-dix])
                                   ));
                    }   
                }
            }
            
            /* fluxes in the z-direction */
            for(int j = igy; j < icy - igy; j++) {
                int l = j * ifz*icx;
                for(int i = igx; i < icx - igx; i++) {
                    int m = l + i * ifz;
                    for(int k = igz; k < ifz - igz; k++) {
                        int n   = m + k;
                        int kc  = k*diz + j*diy + i*dix;
                        int kcm = kc - diz;
                        
                        fz->rhoY[n] -= dto2dz * hplusz[n] *  
                        (     a *    (dp2[kc] - dp2[kcm])  
                         +    b * (  (dp2[kc+dix] - dp2[kcm+dix]) 
                                   + (dp2[kc-dix] - dp2[kcm-dix])  
                                   + (dp2[kc+diy] - dp2[kcm+diy])  
                                   + (dp2[kc-diy] - dp2[kcm-diy])
                                   )
                         +    c * (  (dp2[kc+dix+diy] - dp2[kcm+dix+diy]) 
                                   + (dp2[kc-dix+diy] - dp2[kcm-dix+diy])  
                                   + (dp2[kc+dix-diy] - dp2[kcm+dix-diy])  
                                   + (dp2[kc-dix-diy] - dp2[kcm-dix-diy])
                                   ));
                    }
                }
            }			
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
}

/* ========================================================================== */

static void rhs_fix_for_open_boundaries(
                                        double* rhs, 
                                        const ElemSpaceDiscr* elem, 
                                        ConsVars* Sol_new, 
                                        ConsVars* Sol_old, 
                                        ConsVars* flux[3],
                                        double dt,
                                        MPV* mpv) {
    
    extern User_Data ud;
    
    const States* HydroState = mpv->HydroState;
    
    double** hplus = mpv->Level[0]->wplus;
    
    int ndim = elem->ndim;
    
    if (elem->left == OPEN) {
        
        assert(elem->left == elem->right); /* if you get thrown out here: one-sided open domain not yet impltd */
        
        switch(ndim) {
            case 1: {
                ERROR("function not available");
                break;
            }
                
            case 2: {
                
                const double factor = 1.0 / (elem->dx * dt);
                
                const int igx = elem->igx;
                const int icx = elem->icx;
                const int ifx = elem->ifx;
                
                const int igy = elem->igy;
                const int icy = elem->icy;
                
                double* hplusx  = hplus[0];
                
                double rhs_sum, coeff_sum_right, coeff_sum_left, du;
                int j, mc, nc_left, nc_right, mf, nf_left, nf_right, count;
                
                count = 0;
                rhs_sum  = 0.0;
                coeff_sum_left  = 0.0;
                coeff_sum_right = 0.0;
                
                for(j = igy; j < icy - igy; j++) {
                    mf = j * ifx;
                    nf_left  = mf + igx;
                    nf_right = mf + ifx - igx - 1;
                    
                    rhs_sum         += (flux[0]->rhoe[nf_right] - flux[0]->rhoe[nf_left]);
                    coeff_sum_left  += HydroState->rho0[j] * hplusx[nf_left]; 
                    coeff_sum_right += HydroState->rho0[j] * hplusx[nf_right];
                    
                    count++; 
                } 
                
                rhs_sum         /= count;
                coeff_sum_left  /= count;
                coeff_sum_right /= count;
                du = - rhs_sum / (coeff_sum_left + coeff_sum_right);
                
                for(j = igy; j < icy - igy; j++) {
                    mc = j * icx;
                    mf = j * ifx;
                    nc_left  = mc + igx;
                    nc_right = mc + icx - igx - 1;
                    nf_left  = mf + igx;
                    nf_right = mf + ifx - igx - 1;
                    
                    rhs[nc_left]  -= factor * du * HydroState->rho0[j] * hplusx[nf_left];
                    rhs[nc_right] += factor * du * HydroState->rho0[j] * hplusx[nf_right];
                    
                }
                
                mpv->du = du;
                
                break;
            }
            case 3: {
                printf("3D version of rhs_fix_for_open_boundaries() not yet implemented) )");
                break;
            }
            default: ERROR("ndim not in {1, 2, 3}");
        }
    }
}

/* ========================================================================== */

static void flux_fix_for_open_boundaries(
                                         ConsVars* flux[3],
                                         const ElemSpaceDiscr* elem, 
                                         MPV* mpv) {
    
    const States* HydroState = mpv->HydroState;    
    double** hplus           = mpv->Level[0]->wplus;
    double du                = mpv->du;
    double Sdu               = mpv->Sdu;
    
    int ndim = elem->ndim;
    
    if (elem->left == OPEN) {
        
        assert(elem->left == elem->right); /* if you get thrown out here: one-sided open domain not yet impltd */
        
        switch(ndim) {
            case 1: {
                ERROR("function not available");
                break;
            }
                
            case 2: {
                
                const int igx = elem->igx;
                const int ifx = elem->ifx;
                
                const int igy = elem->igy;
                const int icy = elem->icy;
                
                double* hplusx  = hplus[0];
                
                int j, mf, nf_left, nf_right;
                
                for(j = igy; j < icy - igy; j++) {
                    mf = j * ifx;
                    nf_left  = mf + igx;
                    nf_right = mf + ifx - igx - 1;
                    
                    flux[0]->rho[nf_left]   -= Sdu * HydroState->rho0[j];
                    flux[0]->rho[nf_right]  += Sdu * HydroState->rho0[j];
                    flux[0]->rhoe[nf_left]  -= du  * HydroState->rho0[j] * hplusx[nf_left];
                    flux[0]->rhoe[nf_right] += du  * HydroState->rho0[j] * hplusx[nf_right];
                } 
                
                break;
            }
            case 3: {
                printf("3D version of flux_fix_for_open_boundaries() not yet implemented) )");
                break;
            }
            default: ERROR("ndim not in {1, 2, 3}");
        }
    }
}

/* ========================================================================== */

void update_SI_MIDPT_buoyancy(ConsVars* Sol, 
                              const ConsVars* flux[3], 
                              const MPV* mpv,
                              const ElemSpaceDiscr* elem,
                              const double dt)
{
    /* 
     implicit midpoint for gravity implies that the buoyancy variable
     receives, besides its own advection update, another update 
     contribution due to the vertical advection of the background
     stratification. That is what we do here based on the flux
     determined at the half time level, and done for a full time
     step.
     */
    const int ndim = elem->ndim;
    const double lambda = dt/elem->dy;
    
    switch (ndim) {
        case 1:
            printf("Gravity not implemented for 1D\n");
            return;
            break;
            
        case 2: {
            
            const int icx = elem->icx;
            const int igx = elem->igx;
            const int icy = elem->icy;
            const int igy = elem->igy;
            const int ify = elem->ify;
            
            for (int j=igy; j<icy-igy; j++) {
                int ncj = j*icx;
                int nfj = j;
                double S0p = mpv->HydroState_n->S0[j+1];
                double S0m = mpv->HydroState_n->S0[j];
                double S0c = mpv->HydroState->S0[j];
                for (int i=igx; i<icx-igx; i++) {
                    int ncji = ncj+i;
                    int nfji = nfj+i*ify;
                    Sol->rhoX[BUOY][ncji] += -lambda*(flux[1]->rhoY[nfji+1]*(S0p-S0c) + flux[1]->rhoY[nfji]*(S0c-S0m));
                }
            }
            break;
        }
        case 3: {
            
            const int icx = elem->icx;
            const int igx = elem->igx;
            const int icy = elem->icy;
            const int igy = elem->igy;
            const int ify = elem->ify;
            const int icz = elem->icz;
            const int igz = elem->igz;
            
            for (int k=igz; k<icz-igz; k++) {
                int nck = k*icy*icx;
                int nfk = k*ify;
                for (int j=igy; j<icy-igy; j++) {
                    int nckj = nck + j*icx;
                    int nfkj = nfk + j;
                    double S0p = mpv->HydroState_n->S0[j+1];
                    double S0m = mpv->HydroState_n->S0[j];
                    double S0c = mpv->HydroState->S0[j];
                    for (int i=igx; i<icx-igx; i++) {
                        int nckji = nckj+i;
                        int nfkji = nfkj+i*ify*icz;
                        Sol->rhoX[BUOY][nckji] += -lambda*(flux[1]->rhoY[nfkji+1]*(S0p-S0c) + flux[1]->rhoY[nfkji]*(S0c-S0m));
                    }
                }
            }
        }
            break;
            
        default:
            break;
    }
    Set_Explicit_Boundary_Data(Sol, elem);
}

#ifdef NONLINEAR_EOS_IN_1st_PROJECTION
double Newton_rhs(double* rhs,
                  const double* dpi,
                  const double* pi,
                  const double* rhoY0,
                  const MPV* mpv,
                  const ElemSpaceDiscr* elem)
{
    extern User_Data ud;
    extern Thermodynamic th;
    
    double rhsmax = 0.0;
    
    const double cc1  = 4.0/(mpv->dt*mpv->dt);
    const double cc2  = cc1*ud.Msq*th.gm1inv; 
    const double cexp = 2.0-th.gamm;
    
    const int icx = elem->icx;
    const int igx = elem->igx;
    const int icy = elem->icy;
    const int igy = elem->igy;
    const int icz = elem->icz;
    const int igz = elem->igz;
    
    for (int k=igz; k < icz-igz ; k++) {
        int lc = k*icx*icy;
        for (int j=igy; j < icy-igy ; j++) {
            int mc = lc+j*icx;
            for (int i=igx; i < icx-igx ; i++) {
                int nc = mc + i;
                double pexn = ud.Msq*(mpv->HydroState->p20[j]+pi[nc]);
                double pexo = ud.Msq*(mpv->HydroState->p20[j]+pi[nc]-dpi[nc]);
                double rhoYn = pow(pexn,th.gm1inv);
                double rhoYo = pow(pexo,th.gm1inv);
                rhs[nc] = cc1*(rhoYn - rhoYo) - cc2*pow(rhoYn,cexp)*dpi[nc];
                rhsmax  = MAX_own(rhsmax, fabs(rhs[nc]));
            }
        }
    }
    
    return rhsmax;
}
#endif

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
