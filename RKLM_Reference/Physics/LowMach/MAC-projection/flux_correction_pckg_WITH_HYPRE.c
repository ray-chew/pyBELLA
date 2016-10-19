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
#include "ProjectionType.h"
#include "set_ghostcells_p.h"
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

#ifdef SOLVER_1_HYPRE

#include "space_discretization.h"
/* #include "explicit.h" */

#ifdef WALLTRACK
#include "F3DWall-Common.h"
#endif /* WALLTRACK */

#include "/Users/rupert/Documents/Computation/Walls/wallPoissonTest/Version_2012.05.27/rupert_latest_patched/cellPoisson.h"


double grad_bdry[3];

static void recompute_updates(ConsVars* Sol,
                              double *rhs,
                              ConsVars* Sol0,
                              ConsVars* flux[3],
                              const ElemSpaceDiscr* elem,
                              const MPV* mpv,
                              const double dt);

static void fix_solid_wall_boundaries(ConsVars* Sol,
                                      ConsVars* flux[3],
                                      const ElemSpaceDiscr* elem,
                                      const MPV* mpv,
                                      const double dt);

static enum Constraint integral_condition(
                                          ConsVars* flux[3],
                                          const ElemSpaceDiscr* elem);

static enum Constraint integral_condition_rhs(
                                                   double *rhs,
                                                   const ElemSpaceDiscr* elem);


static void controlled_variable_change_explicit(
                                                double* rhs,
                                                const ElemSpaceDiscr* elem,
                                                ConsVars* Sol_new,
                                                ConsVars* Sol_old,
                                                double dt,
                                                const MPV* mpv);


#define CORRECTION_FROM_MG_CODE
#ifdef CORRECTION_FROM_MG_CODE
static void operator_coefficients(
                                  double* hplus[3],
                                  double* wcenter,
                                  double* wgrav,
                                  const ElemSpaceDiscr* elem,
                                  ConsVars* Sol,
                                  ConsVars* Sol0,
                                  const MPV* mpv,
                                  const double dt);

static void flux_correction_due_to_pressure_gradients(
                                                      ConsVars* flux[3],
                                                      VectorField* buoy,
                                                      const ElemSpaceDiscr* elem,
                                                      ConsVars* Sol,
                                                      ConsVars* Sol0,
                                                      const MPV* mpv,
                                                      double* hplus[3],
                                                      double* hgrav,
                                                      double* dp2,
                                                      const double t,
                                                      const double dt,
                                                      const double implicitness);

static void flux_correction_due_to_pressure_values(
                                       ConsVars* flux[3],
                                       VectorField* buoy,
                                       const ElemSpaceDiscr* elem,
                                       ConsVars* Sol,
                                       double* dp2, 
                                       const double dt);
#else /* CORRECTION_FROM_MG_CODE */
static void interface_enthalpy(
                               double* h[3],
                               const ElemSpaceDiscr* elem,
                               ConsVars* flux[3],
                               ConsVars* Sol,
                               ConsVars* Sol0,
                               const MPV* mpv);

static void correction(
                       ConsVars* dSol,
                       ConsVars* Sol,
                       const ConsVars* Sol0,
                       ConsVars* flux[3],
                       double* rhs,
                       const ElemSpaceDiscr* elem,
                       const MPV* mpv,
                       const double t,
                       const double dt);
#endif /* CORRECTION_FROM_MG_CODE */

/* ========================================================================== */

static enum Boolean first_projection_is_initialized = WRONG;

void flux_correction(
                     ConsVars* flux[3],
                     VectorField* buoy,
                     const ElemSpaceDiscr* elem,
                     const NodeSpaceDiscr* node,
                     ConsVars* Sol,
                     ConsVars* Sol0,
                     MPV* mpv,
                     const double t,
                     const double dt,
                     const double implicitness) {
    
    extern User_Data ud;
    extern Thermodynamic th;
    extern ConsVars* dSol;
    extern double *W0, *W1;
    
    const double dx = elem->dx;
    const double dy = elem->dy;
    
    const int icx     = elem->icx;
    const int icy     = elem->icy;
    const int ndim    = elem->ndim;
    const int nfx     = elem->nfx;
    const int nfy     = elem->nfy;
    const int nfz     = elem->nfz;
    const int icx_inn = elem->icx-2*elem->igx;
    const int icy_inn = elem->icy-2*elem->igy;
    const int nfx_inn = (icx_inn+1)*icy_inn;
    const int nfy_inn = (icy_inn+1)*icx_inn;
    const int nc      = elem->nc;
    
    double*  rhs      = mpv->Level[0]->rhs;
    double*  dp2      = mpv->Level[0]->p;
    double** hplus    = mpv->Level[0]->wplus;
    double*  hcenter  = mpv->Level[0]->wcenter;
    double*  hgrav    = mpv->Level[0]->wgrav;

    
    double *h[3];

    h[0] = (double*)malloc(2 * nfx * sizeof(double));
    if(ndim > 1) h[1] = (double*)malloc(2 * nfy * sizeof(double));
    if(ndim > 2) h[2] = (double*)malloc(2 * nfz * sizeof(double));
    
    int n;
    
    if (first_projection_is_initialized == WRONG) {
        const int inx = node->icx;
        const int iny = node->icy;
        const int inx_inn = inx-2*node->igx;
        const int iny_inn = iny-2*node->igy;
        double* fakeG = (double*)malloc(node->nc * sizeof(double));
        
        for (int jn = 0; jn < iny; jn++) {int mn = jn*inx;
            for (int in = 0; in < inx; in++) {int nn = mn + in;
                fakeG[nn] = -1.0 - sqrt((dx*(in-0.5*inx))*(dx*(in-0.5*inx)) + (dy*(jn-0.5*iny))*(dy*(jn-0.5*iny)));
            }
        }
        
        initCellPoisson(node->x[node->igx], node->y[node->igy], node->dx, node->dy, inx_inn, iny_inn, fakeG);
        free(fakeG);
        first_projection_is_initialized = CORRECT;
    }

    printf("\n\n====================================================");
    printf("\n First Projection");
    printf("\n====================================================\n");
    
    for (int k=0; k<elem->nc; k++) {
        buoy->x[k]=buoy->y[k]=0.0;
    }
    

#ifdef CORRECTION_FROM_MG_CODE
    operator_coefficients(hplus, hcenter, hgrav, elem, Sol, Sol0, mpv, dt);
    map2D_to_michaels_memory_opcoeffs_cells(h, hplus, elem);
#else
    interface_enthalpy(h, elem, flux, Sol, Sol0, mpv);
#endif

    controlled_variable_change_explicit(rhs, elem, Sol, Sol0, dt, mpv);
    assert(integral_condition_rhs(rhs, elem));
    
    /* flip stuff into Michael's memory scheme */
    map2D_to_michaels_memory_cells(rhs, elem, W0);
    map2D_to_michaels_memory_faces_x(h[0], elem, W0, W1);
    map2D_to_michaels_memory_faces_y(h[1], elem, W0, W1);
    
    solveCellPoisson(dx, dy, icx_inn, icy_inn, h[0], h[1], rhs, dp2);
    
#if 0
    extern User_Data ud;
    FILE *prhsfile = NULL;
    char fn2[100], fieldname2[90];
    
    sprintf(fn2, "%s/Tests/rhs_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "rhs_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->icy-2*mpv->Level[0]->elem->igy,
             mpv->Level[0]->elem->icx-2*mpv->Level[0]->elem->igx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             rhs,
             fn2,
             fieldname2);
    
    sprintf(fn2, "%s/Tests/p2_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "p2_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->icy-2*mpv->Level[0]->elem->igy,
             mpv->Level[0]->elem->icx-2*mpv->Level[0]->elem->igx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             dp2,
             fn2,
             fieldname2);

    sprintf(fn2, "%s/Tests/hx_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hx_c");
    
    WriteHDF(prhsfile,
             2*(mpv->Level[0]->elem->icy-2*mpv->Level[0]->elem->igy),
             mpv->Level[0]->elem->ifx-2*mpv->Level[0]->elem->igx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             h[0],
             fn2,
             fieldname2);

    sprintf(fn2, "%s/Tests/hy_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hy_c");
    
    WriteHDF(prhsfile,
             2*(mpv->Level[0]->elem->ify-2*mpv->Level[0]->elem->igy),
             mpv->Level[0]->elem->icx-2*mpv->Level[0]->elem->igx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             h[1],
             fn2,
             fieldname2);

    sprintf(fn2, "%s/Tests/hplusx_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hplusx_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->ifx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             hplus[0],
             fn2,
             fieldname2);
    
    sprintf(fn2, "%s/Tests/hplusy_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hplusy_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->ify,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             hplus[1],
             fn2,
             fieldname2);
#endif
    


#ifdef CORRECTION_FROM_MG_CODE
    map2D_to_ruperts_memory_cells(dp2,elem,W0);
    map2D_to_ruperts_memory_interface_coeffs_x(h[0], elem, W0);
    map2D_to_ruperts_memory_interface_coeffs_y(h[1], elem, W0);

#if 0
    extern User_Data ud;
    FILE *pfile = NULL;
    char fn3[100], fieldname3[90];
    
    sprintf(fn3, "%s/Tests/p2_c_rupesmem.hdf", ud.file_name);
    sprintf(fieldname3, "p2_c_rupesmem");
    
    WriteHDF(pfile,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             dp2,
             fn3,
             fieldname3);
    
    sprintf(fn3, "%s/Tests/hx_c_rupesmem.hdf", ud.file_name);
    sprintf(fieldname3, "hx_rupesmem");
    
    WriteHDF(pfile,
             mpv->Level[0]->elem->ifx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             h[0],
             fn3,
             fieldname3);
    
    sprintf(fn3, "%s/Tests/hy_c_rupesmem.hdf", ud.file_name);
    sprintf(fieldname3, "hy_rupesmem");
    
    WriteHDF(pfile,
             mpv->Level[0]->elem->ify,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             h[1],
             fn3,
             fieldname3);
    
#endif
    
    set_ghostcells_p2(dp2, h, hgrav, elem, elem->igx);

    flux_correction_due_to_pressure_gradients(flux, buoy, elem, Sol, Sol0, mpv, hplus, hgrav, dp2, t, dt, implicitness);
    if (ud.p_flux_correction) {
        flux_correction_due_to_pressure_values(flux, buoy, elem, Sol, dp2, dt); 
    }
    
    
    /*
     flux_fix_for_open_boundaries(flux, elem, mpv);
     */
#else /* CORRECTION_FROM_MG_CODE */
    /**/
    double* dp2Grid;
    dp2Grid = (double*)malloc(2 * nc * sizeof(double));
    memset(dp2Grid, 0, nc * sizeof(double));

    getCellIntegralsGrid(dx, dy, icx_inn, icy_inn, hplus[0], hplus[1], dp2, dp2Grid, f[0], f[1]);

    map2D_to_ruperts_memory_cells(dp2,elem,W0);
    set_ghostcells_p2(dp2, h, hgrav, elem, elem->igx);

    for (n = 0; n < nfx_inn; n++) {
        flux[0]->rhoY[n] = f[0][2*n] + f[0][2*n+1];
    }
    for (n = 0; n < nfy_inn; n++) {
        flux[1]->rhoY[n] = f[1][2*n] + f[1][2*n+1];
    }
    
    map2D_to_ruperts_memory_cells(dp2Grid,elem,W0);
    map2D_to_ruperts_memory_cells(rhs,elem,W0);
    
    /* here I am overwriting the explicitly predicted rhoY-flux
     with its pressure gradient correction */
    map2D_to_ruperts_memory_faces_x(flux[0]->rhoY,elem,W0);
    map2D_to_ruperts_memory_faces_y(flux[1]->rhoY,elem,W0);
    
    correction(dSol, Sol, Sol0, flux, rhs, elem, mpv, t, dt);
    free(dp2Grid);
    if(ndim > 2) free(f[2]);
    if(ndim > 1) free(f[1]);
    free(f[0]);
#endif  /* CORRECTION_FROM_MG_CODE */
    
    memcpy(mpv->dp2_cells, dp2, elem->nc*sizeof(double));
    
    /* store results in mpv-fields */
    for(n=0; n<elem->nc; n++) mpv->p2_cells[n] = MAC_PROJ_OLDP_WEIGHT*mpv->p2_cells[n] + MAC_PROJ_DELP_WEIGHT*dp2[n];
    
    set_ghostcells_p2(mpv->p2_cells, h, hgrav, elem, elem->igx);

    if(ndim > 2) free(h[2]);
    if(ndim > 1) free(h[1]);
    free(h[0]);
}

/* ========================================================================== */

static void fix_solid_wall_boundaries(ConsVars* Sol,
                                      ConsVars* flux[3],
                                      const ElemSpaceDiscr* elem,
                                      const MPV* mpv,
                                      const double dt) {
    

    extern User_Data ud;

    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
    
    const int ifx = elem->ifx;
    const int ify = elem->ify;
    
    const int igx = elem->igx;
    const int igy = elem->igy;
    const int igz = elem->igz;
    
    double uu;
    double lambda_x = dt / elem->dx;
    double lambda_y = dt / elem->dy;
    
    int i, j, k, l, m, n, lf, mf, nf;
    
    if (elem->ndim > 0) {
        /* left */
        if (ud.bdrytype_min[0] == WALL) {
            for (k = igz; k < icz-igz; k++) {l = k*icx*icy; lf = k*ifx*icy;
                for (j = elem->igy; j < elem->icy-elem->igy; j++) {m = l + j*icx; mf = lf + j*ifx;
                    i = igx;
                    n = m + i; nf = mf + i;
                    uu = 0.5 * (Sol->rhou[n]/Sol->rho[n] + Sol->rhou[n-1]/Sol->rho[n-1]);
                    Sol->rho[n]    -= lambda_x * flux[0]->rho[nf];
                    Sol->rhou[n]   -= lambda_x * flux[0]->rho[nf] * uu;
                    Sol->rhov[n]   -= lambda_x * flux[0]->rhov[nf];
                    Sol->rhow[n]   -= lambda_x * flux[0]->rhow[nf];
                    Sol->rhoe[n]   -= lambda_x * flux[0]->rhoe[nf];
                    Sol->rhoY[n]   -= lambda_x * flux[0]->rhoY[nf];
                    Sol->rhoZ[n]   -= lambda_x * flux[0]->rhoZ[nf];
                    
                    
                    flux[0]->rho[nf]  = 0.0;
                    flux[0]->rhou[nf] = 0.0; /* pressure-momentum flux not included */
                    flux[0]->rhov[nf] = 0.0;
                    flux[0]->rhow[nf] = 0.0;
                    flux[0]->rhoe[nf] = 0.0;
                    flux[0]->rhoY[nf] = 0.0;
                    flux[0]->rhoZ[nf] = 0.0;
                }
            }
        }
        
        /* right */
        if (ud.bdrytype_max[0] == WALL) {
            for (k = igz; k < icz-igz; k++) {l = k*icx*icy; lf = k*ifx*icy;
                for (j = elem->igy; j < elem->icy-elem->igy; j++) {m = l + j*icx; mf = lf + j*ifx;
                    i = icx-igx-1;
                    n = m + i; nf = mf + i + 1;
                    uu = 0.5 * (Sol->rhou[n]/Sol->rho[n] + Sol->rhou[n+1]/Sol->rho[n+1]);
                    Sol->rho[n]    -= -lambda_x * flux[0]->rho[nf];
                    Sol->rhou[n]   -= -lambda_x * flux[0]->rho[nf] * uu;
                    Sol->rhov[n]   -= -lambda_x * flux[0]->rhov[nf];
                    Sol->rhow[n]   -= -lambda_x * flux[0]->rhow[nf];
                    Sol->rhoe[n]   -= -lambda_x * flux[0]->rhoe[nf];
                    Sol->rhoY[n]   -= -lambda_x * flux[0]->rhoY[nf];
                    Sol->rhoZ[n]   -= -lambda_x * flux[0]->rhoZ[nf];
                    
                    
                    flux[0]->rho[nf]  = 0.0;
                    flux[0]->rhou[nf] = 0.0;
                    flux[0]->rhov[nf] = 0.0;
                    flux[0]->rhow[nf] = 0.0;
                    flux[0]->rhoe[nf] = 0.0;
                    flux[0]->rhoY[nf] = 0.0;
                    flux[0]->rhoZ[nf] = 0.0;
                }
            }
        }
    }
    
    if (elem->ndim > 1) {
        /* bottom */
        if (ud.bdrytype_min[1] == WALL) {
            for (k = igz; k < icz-igz; k++) {l = k*icx*icy; lf = k*icx*ify;
                for (i = elem->igx; i < elem->icx-elem->igx; i++) {m = l + i; mf = lf + i*ify;
                    j = igy;
                    n = m + j*icx; nf = mf + j;
                    uu = 0.5 * (Sol->rhov[n]/Sol->rho[n] + Sol->rhov[n-1]/Sol->rho[n-1]);
                    Sol->rho[n]    -= lambda_y * flux[1]->rho[nf];
                    Sol->rhou[n]   -= lambda_y * flux[1]->rhou[nf];
                    Sol->rhov[n]   -= lambda_y * flux[1]->rho[nf] * uu;
                    Sol->rhow[n]   -= lambda_y * flux[1]->rhow[nf];
                    Sol->rhoe[n]   -= lambda_y * flux[1]->rhoe[nf];
                    Sol->rhoY[n]   -= lambda_y * flux[1]->rhoY[nf];
                    Sol->rhoZ[n]   -= lambda_y * flux[1]->rhoZ[nf];
                    
                    flux[1]->rho[nf]  = 0.0;
                    flux[1]->rhou[nf] = 0.0;
                    flux[1]->rhov[nf] = 0.0;
                    flux[1]->rhow[nf] = 0.0;
                    flux[1]->rhoe[nf] = 0.0;
                    flux[1]->rhoY[nf] = 0.0;
                    flux[1]->rhoZ[nf] = 0.0;
                }
            }
        }
        
        /* top */
        if (ud.bdrytype_max[1] == WALL) {
            for (k = igz; k < icz-igz; k++) {l = k*icx*icy; lf = k*icx*ify;
                for (i = elem->igx; i < elem->icx-elem->igx; i++) {m = l + i; mf = lf + i*ify;
                    j = icy-igy-1;
                    n = m + j*icx; nf = mf + j + 1;
                    uu = 0.5 * (Sol->rhov[n]/Sol->rho[n] + Sol->rhov[n+1]/Sol->rho[n+1]);
                    Sol->rho[n]    -= -lambda_y * flux[1]->rho[nf];
                    Sol->rhou[n]   -= -lambda_y * flux[1]->rhou[nf];
                    Sol->rhov[n]   -= -lambda_y * flux[1]->rho[nf] * uu;
                    Sol->rhow[n]   -= -lambda_y * flux[1]->rhow[nf];
                    Sol->rhoe[n]   -= -lambda_y * flux[1]->rhoe[nf];
                    Sol->rhoY[n]   -= -lambda_y * flux[1]->rhoY[nf];
                    Sol->rhoZ[n]   -= -lambda_y * flux[1]->rhoZ[nf];
                    
                    flux[1]->rho[nf]  = 0.0;
                    flux[1]->rhou[nf] = 0.0;
                    flux[1]->rhov[nf] = 0.0;
                    flux[1]->rhow[nf] = 0.0;
                    flux[1]->rhoe[nf] = 0.0;
                    flux[1]->rhoY[nf] = 0.0;
                    flux[1]->rhoZ[nf] = 0.0;
                }
            }
        }
    }
    
    if (elem->ndim > 2) {
        ERROR("fix_solid_wall_boundaries() not implemented for 3D \n");
    }
}

/* ========================================================================== */

static void recompute_updates(ConsVars* Sol,
                              double *rhs,
                              ConsVars* Sol0,
                              ConsVars* flux[3],
                              const ElemSpaceDiscr* elem,
                              const MPV* mpv,
                              const double dt) {
    
    /* complete flux-based recomputation of updates; test, whether
     this creates a problem solving the projection step.
     */
    extern User_Data ud;
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    
    const int ifx = elem->ifx;
    const int ify = elem->ify;
    
    const int igx = elem->igx;
    const int igy = elem->igy;
    const double dx = elem->dx;
    const double dy = elem->dy;
    
    double lambda_x = dt / elem->dx;
    double lambda_y = dt / elem->dy;
    
    const double factor = 2.0 * dx * dy / (dt * dt);
    
    double delta, deltax, deltay, deltaSol, delmax, ddelmax;
    int i, j, mc, nc, mfx, nfx, nfxp, mfy, nfy, nfyp;
    
    switch (elem->ndim) {
        case 1:
            ERROR("\n\n1D-case of recompute_updates() not implemented");
            break;
        case 2:
            
            delmax  = 0.0;
            ddelmax = 0.0;
            
            for (j = igy; j < icy-igy; j++) {
                mc  = j*icx;
                mfx = j*ifx;
                mfy = j;
                for (i = igx; i < icx-igx; i++) {
                    nc   = mc + i;
                    nfx  = mfx + i;
                    nfxp = nfx + 1;
                    nfy  = mfy + i*ify;
                    nfyp = nfy + 1;
                    
                    deltax        = - lambda_x * (flux[0]->rhoY[nfxp]   - flux[0]->rhoY[nfx]);
                    deltay        = - lambda_y * (flux[1]->rhoY[nfyp]   - flux[1]->rhoY[nfy]);
                    delta         = deltax + deltay;
                    deltaSol      = Sol->rhoY[nc] - Sol0->rhoY[nc];
                    
                    delmax  = MAX_own(delmax,fabs(deltaSol));
                    ddelmax = MAX_own(ddelmax, fabs(delta-deltaSol));
                    
                    Sol->rhoY[nc] = Sol0->rhoY[nc] + delta;
                    rhs[nc]       = - factor * delta;
                    
                }
            }
            break;
        case 3:
            ERROR("\n\n3D-case of recompute_updates() not implemented");
            break;
        default:
            ERROR("\n\nwe are not doing string theory");;
    }
}

/* ========================================================================== */

static enum Constraint integral_condition(
                                          ConsVars* flux[3],
                                          const ElemSpaceDiscr* elem) {
    
    extern User_Data ud;
    extern int ndim;
    
    if(ud.acoustic_timestep == 0) {
#ifdef CHEMSOURCE
        ERROR("function not available");
#endif
        
        switch(ndim) {
                
            case 1: {
                ERROR("function not available");
                break;
            }
                
            case 2: {
                int i, j, jifx, iify;
                const int igx = elem->igx;
                const int icx = elem->icx;
                const int ifx = elem->ifx;
                const int igy = elem->igy;
                const int icy = elem->icy;
                const int ify = elem->ify;
                const ConsVars* f = flux[0];
                const ConsVars* g = flux[1];
                double F, G, FY, GY;
                
                for(j = igy, F = 0.0; j < icy - igy; j++) {
                    jifx = j * ifx;
                    F  += f->rhoe[jifx + icx - igx] - f->rhoe[jifx + igx];
                }
                for(i = igx, G = 0.0; i < icx - igx; i++) {
                    iify =  i * ify;
                    G  += g->rhoe[iify + icy - igy] - g->rhoe[iify + igy];
                }
                
                for(j = igy, FY = 0.0; j < icy - igy; j++) {
                    jifx = j * ifx;
                    FY  += f->rhoY[jifx + icx - igx] - f->rhoY[jifx + igx];
                }
                for(i = igx, GY = 0.0; i < icx - igx; i++) {
                    iify =  i * ify;
                    GY  += g->rhoY[iify + icy - igy] - g->rhoY[iify + igy];
                }
                
                if(fabs(F * elem->dy + G * elem->dx) > sqrt(DBL_EPSILON) || fabs(FY * elem->dy + GY * elem->dx)  > sqrt(DBL_EPSILON))
                    return VIOLATED;
                else
                    return SATISFIED;
                break;
            }
                
            case 3: {
                int i, j, k, l, m, n;
                const int igx = elem->igx;
                const int icx = elem->icx;
                const int ifx = elem->ifx;
                const int igy = elem->igy;
                const int icy = elem->icy;
                const int ify = elem->ify;
                const int igz = elem->igz;
                const int icz = elem->icz;
                const int ifz = elem->ifz;
                const int ifxicy = ifx * icy;
                const int ifyicz = ify * icz;
                const int ifzicx = ifz * icx;
                const ConsVars* f = flux[0];
                const ConsVars* g = flux[1];
                const ConsVars* h = flux[2];
                double F, G, H;
                
                for(k = igz, F = 0.0; k < icz - igz; k++) {l = k * ifxicy;
                    for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
                        F += f->rhoe[m + icx - igx] - f->rhoe[m +  igx];
                    }
                }
                for(i = igx, G = 0.0; i < icx - igx; i++) {m = i * ifyicz;
                    for(k = igz; k < icz - igz; k++) {n = m + k * ify;
                        G += g->rhoe[n + icy - igy] - g->rhoe[n + igy];
                    }
                }
                for(j = igy, H = 0.0; j < icy - igy; j++) {n = j * ifzicx;
                    for(i = igx; i < icx - igx; i++) {l = n + i * ifz;
                        H += h->rhoe[l + icz - igz] - h->rhoe[l + igz];
                    }
                }
                if(fabs(F * elem->dy + G * elem->dx * H * elem->dz) > DBL_EPSILON) {
                    return VIOLATED;
                }
                else
                    return SATISFIED;
                break;
            }
                
            default: ERROR("ndim not in {1, 2, 3}");
                
        }
    }
    else {
        ERROR("function not available");
    }
    return VIOLATED;
}

/* ========================================================================== */

static enum Constraint integral_condition_rhs(
                                                   double *rhs,
                                                   const ElemSpaceDiscr* elem) {
    
    
    const double res_tol = 10000*DBL_EPSILON;
    
    switch(elem->ndim) {
            
        case 1: {
            ERROR("integral_condition_rhs not implemented fro ndim = 1\n");
            break;
        }
            
        case 2: {
            int i, j, m, n;
            const int icx = elem->icx;
            const int icy = elem->icy;
            const int igx = elem->igx;
            const int igy = elem->igy;
            double res, vol;
            
            res = 0.0;
            vol = 0.0;
            
            /* regular cell contributions */
            for(j = igy; j < icy - igy; j++) {m = j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    res += rhs[n];
                    vol += 1.0;
                }
            }
            
#if 0
            {
                double rescorr;
                
                /* improve */
                rescorr = res/vol;
                for(j = igy; j < icy - igy; j++) {m = j * icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        rhs[n] -= (Wall_old->type[n] == 0 ? rescorr : 0.0);
                    }
                }
                for (imc = 0; imc < WallCs->nMCs; imc++) {
                    n = WallCs->MCi[0].ic_from_imc[imc];
                    rhs[n] -= WallCs->VF[imc] * rescorr;
                }
                
                res = 0.0;
                /* check improvement */
                for(j = igy; j < icy - igy; j++) {m = j * icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        res += rhs[n];
                    }
                }
                
                /* printf("sum of rhs [P1] impr = %e\n", res); */
            }
#endif
            
            if(fabs(res) > res_tol) {
                printf("\nintegral condition 1st projection: res = %e", fabs(res));
                return VIOLATED;
            }
            else
                return SATISFIED;
            break;
        }
            
        case 3: ERROR("integral_condition_rhs not implemented fro ndim = 3\n");
            
        default: ERROR("ndim not in {1, 2, 3}");
            
    }
    
    return VIOLATED;
}

/* ========================================================================== */

static void interface_enthalpy(
                               double* h[3],
                               const ElemSpaceDiscr* elem,
                               ConsVars* flux[3],
                               ConsVars* Sol,
                               ConsVars* Sol0,
                               const MPV* mpv) {
    
    const int ndim = elem->ndim;
    const double tol = DBL_MAX;
    
    switch(ndim) {
        case 1: {
            ERROR("function not available");
            break;
        }
        case 2: {
            int i, j, m, n;
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int ifx = elem->ifx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int ify = elem->ify;
            double* hx = h[0];
            double* hy = h[1];
            double h0min = 9999.0;
            double h1min = 9999.0;
            
            for(n=0; n<icy*ifx; n++)  hx[2*n] = hx[2*n+1] = 1.0;
            for(j = igy; j < icy - igy; j++) {m = j * ifx;
                for(i = igx; i < ifx - igx; i++) {n = m + i;
                    const int ic  = MAX_own(0, n  - j);
                    const int icm = MAX_own(0, ic - 1);
                    
                    const double hi  = 0.5 * (Sol->rhoY[ic] /Sol->rho[ic]  + Sol0->rhoY[ic] /Sol0->rho[ic] );
                    const double him = 0.5 * (Sol->rhoY[icm]/Sol->rho[icm] + Sol0->rhoY[icm]/Sol0->rho[icm]);

                    hx[2*n] = hx[2*n+1] = 0.5 * (hi + him);
                    h0min = MIN_own(h0min,0.5 * (hi + him));
                }
            }

            for(m=0; m<icx*ify; m++)  hy[2*m] = hy[2*m+1] = 1.0;
            for(i = igx; i < icx - igx; i++) {n = i * ify;
                for(j = igy; j < ify - igy; j++) {m = n + j;
                    const int jc  = MAX_own(0, j*icx + i);
                    const int jcm = MAX_own(0,  jc - icx);
                    
                    const double hj  = 0.5 * (Sol->rhoY[jc]  / Sol->rho[jc]  + Sol0->rhoY[jc]  / Sol0->rho[jc]);
                    const double hjm = 0.5 * (Sol->rhoY[jcm] / Sol->rho[jcm] + Sol0->rhoY[jcm] / Sol0->rho[jcm]);

                    hy[2*m] = hy[2*m+1] = 0.5 * (hj + hjm);
                    h1min = MIN_own(h1min,0.5 * (hj + hjm));
                }
            }

            break;
        }
        
        default: ERROR("ndim not in {1, 2,3}");
    }
}

/* ========================================================================== */

static void controlled_variable_change_explicit(
                                                double* rhs,
                                                const ElemSpaceDiscr* elem,
                                                ConsVars* Sol_new,
                                                ConsVars* Sol_old,
                                                double dt,
                                                const MPV* mpv) {
    
    int ndim = elem->ndim;
    
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
            const double dx = elem->dx;
            const double dy = elem->dy;
            
            const double factor = 2.0 * dx * dy / (dt * dt);
            
            int i, j, m, n;
            
            for(i=0; i<elem->nc; i++) rhs[i] = 0.0;
            
            for(j = igy; j < icy - igy; j++) {
                m = j * icx;
                
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    rhs[n] = - factor * ( (Sol_new->rhoY[n] - Sol_old->rhoY[n]));
                }
            }
            
            break;
        }
        case 3: {
            ERROR("ndim = 3 - option not implemented (density_change_explicit() )");
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
    
}


#ifdef CORRECTION_FROM_MG_CODE
/* ========================================================================== */

static void operator_coefficients(
                                  double* hplus[3],
                                  double* wcenter,
                                  double* wgrav,
                                  const ElemSpaceDiscr* elem,
                                  ConsVars* Sol,
                                  ConsVars* Sol0,
                                  const MPV* mpv,
                                  const double dt) {
    
    extern User_Data ud;
    extern Thermodynamic th;
    
    const int ndim = elem->ndim;
    
    const int impl_grav_th = ud.implicit_gravity_theta;
    const int impl_grav_pr = ud.implicit_gravity_press;
    
    const double implicitness = ud.implicitness;
    const double ccw = 2.0; /* ccenterweight   4.0 */
    const double ccenter = - ccw * (ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt);
    /* for ccw = 4.0:
     first factor 2.0 from reinterpretation of Helmholtz-solution as
     (p^{n+1}-p^{n}) and SECOND-ORDER update to t^{n+1/2} in the first projection
     */

    const double cexp    = 1.0-th.gamm;
    
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
            const double dx = elem->dx;
            const double dy = elem->dy;
            double* hx = hplus[0];
            double* hy = hplus[1];
            double* hc = wcenter;
            double* hg = wgrav;
            
            double hi, him, hj, hjm, g, gimp, thet, thetm, Msq;
            
            int i, j, m, n, ic, icm, jc, jcm;
            
            Msq = ud.Msq;
            g   = ud.gravity_strength[0];
            
            for(j = igy; j < icy - igy; j++) {m = j * ifx;
                for(i = igx; i < ifx - igx; i++) {n = m + i;
                    ic  = n - j;
                    icm = ic - 1;
                    
                    hi       = 0.5 * ( Sol->rhoY[ic] / Sol->rho[ic]
                                      + Sol0->rhoY[ic] / Sol0->rho[ic] );
                    
                    him      = 0.5 * ( Sol->rhoY[icm] / Sol->rho[icm]
                                      + Sol0->rhoY[icm] / Sol0->rho[icm] );
                    
                    thet  = Sol->rhoY[ic]  / Sol->rho[ic] ;
                    thetm = Sol->rhoY[icm] / Sol->rho[icm];
                    gimp   = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dx*0.5*(thet+thetm)));
                    
                    hx[n] = 0.5 * (hi + him) * implicitness * gimp;
                    assert(hx[n] > 0.0);
                }
            }
            
            
            g = ud.gravity_strength[1];
            
            for(i = igx; i < icx - igx; i++) {n = i * ify;
                for(j = igy; j < ify - igy; j++) {m = n + j;
                    jc = j * icx + i;
                    jcm = jc - icx;
                    
                    hj       = 0.5 * ( Sol->rhoY[jc] / Sol->rho[jc]
                                      + Sol0->rhoY[jc] / Sol0->rho[jc]);
                    
                    hjm      = 0.5 * ( Sol->rhoY[jcm] / Sol->rho[jcm]
                                      + Sol0->rhoY[jcm] / Sol0->rho[jcm]);
                    
                    thet  = Sol->rhoY[jc]  / Sol->rho[jc] ;
                    thetm = Sol->rhoY[jcm] / Sol->rho[jcm];
                    gimp  = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dy*0.5*(thet+thetm)));
                    
                    hy[m] = 0.5 * (hj + hjm) * implicitness * gimp;
                    
                    hg[m] = 0.0;
                    
                    assert(hy[m] > 0.0);
                }
            }
            
            for(j = igy; j < icy - igy; j++) {m = j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    hc[n] = ccenter * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp);
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
            const double dx = elem->dx;
            const double dy = elem->dy;
            const double dz = elem->dz;
            double* hx = hplus[0];
            double* hy = hplus[1];
            double* hz = hplus[2];
            double* hc = wcenter;
            double* hg = wgrav;
            
            double hi, him, hj, hjm, hk, hkm, g, gimp, thet, thetm, Msq;
            
            int i, j, k, l, m, n, ic, icm, jc, jcm, kc, kcm;
            
            Msq = ud.Msq;
            g   = ud.gravity_strength[0];
            assert(g==0.0); /* implicit gravity only for y-direction */
            
            for(k = igz; k < icz - igz; k++) {l = k * ifx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
                    for(i = igx; i < ifx - igx; i++) {n = m + i;
                        ic  = k*icx*icy + j*icx + i;
                        icm = ic - 1;
                        
                        hi       = 0.5 * ( Sol->rhoY[ic] / Sol->rho[ic]
                                          + Sol0->rhoY[ic] / Sol0->rho[ic] );
                        
                        him      = 0.5 * ( Sol->rhoY[icm] / Sol->rho[icm]
                                          + Sol0->rhoY[icm] / Sol0->rho[icm] );
                        
                        thet  = Sol->rhoY[ic]  / Sol->rho[ic] ;
                        thetm = Sol->rhoY[icm] / Sol->rho[icm];
                        gimp   = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dx*0.5*(thet+thetm)));
                        
                        hx[n] = 0.5 * (hi + him) * implicitness * gimp;
                        assert(hx[n] > 0.0);
                    }
                }
            }
            
            
            g = ud.gravity_strength[1];
            
            for(i = igx; i < icx - igx; i++) {l = i * ify*icz;
                for(k = igz; k < icz - igz; k++) {m = l + k * ify;
                    for(j = igy; j < ify - igy; j++) {n = m + j;
                        jc  = k*icx*icy + j*icx + i;
                        jcm = jc - icx;
                        
                        hj       = 0.5 * ( Sol->rhoY[jc] / Sol->rho[jc]
                                          + Sol0->rhoY[jc] / Sol0->rho[jc]);
                        
                        hjm      = 0.5 * ( Sol->rhoY[jcm] / Sol->rho[jcm]
                                          + Sol0->rhoY[jcm] / Sol0->rho[jcm]);
                        
                        thet  = Sol->rhoY[jc]  / Sol->rho[jc] ;
                        thetm = Sol->rhoY[jcm] / Sol->rho[jcm];
                        gimp   = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dy*0.5*(thet+thetm)));
                        
                        hy[n] = 0.5 * (hj + hjm) * implicitness * gimp;
                        
                        hg[m] = th.gamminv * pow(0.5*(Sol->rhoY[jc]+Sol->rhoY[jcm]),cexp) * g * implicitness * gimp;
                        /* hg[m] = th.gamminv * pow(0.5*(Sol->rhoY[jc]+Sol->rhoY[jcm]),cexp) * (g/Msq) * implicitness * gimp; */
                        
                        assert(hy[n] > 0.0);
                    }
                }
            }
            
            g = ud.gravity_strength[2];
            assert(g==0.0); /* implicit gravity only for y-direction */
            
            for(j = igy; j < icy - igy; j++) {l = j * ifz*icx;
                for(i = igx; i < icx - igx; i++) {m = l + i * ifz;
                    for(k = igz; k < ifz - igz; k++) {n = m + k;
                        kc  = k*icx*icy + j*icx + i;
                        kcm = kc - icx*icy;
                        
                        hk       = 0.5 * ( Sol->rhoY[kc] / Sol->rho[kc]
                                          + Sol0->rhoY[kc] / Sol0->rho[kc]);
                        
                        hkm      = 0.5 * ( Sol->rhoY[kcm] / Sol->rho[kcm] 
                                          + Sol0->rhoY[kcm] / Sol0->rho[kcm]);
                        
                        thet  = Sol->rhoY[kc]  / Sol->rho[kc] ;
                        thetm = Sol->rhoY[kcm] / Sol->rho[kcm];
                        gimp   = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dz*0.5*(thet+thetm)));
                        
                        hz[n] = 0.5 * (hk + hkm) * implicitness * gimp;
                        assert(hz[n] > 0.0); 
                    }
                }
            }
            
            for(k = igz; k < icz - igz; k++) {l = k*icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j*icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        hc[n] = ccenter * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp);
                    }
                }
            }
            
            
            break;
        }
        default: ERROR("ndim not in {1,2,3}");
    }
}

/* ========================================================================== */
#ifdef THIRD_ORDER_UPWIND_CORRECTION
#define INTERPOL(l,c,cr) (FIVE_SIXTHS * c + ONE_THIRD * cr - ONE_SIXTH * l)
#endif
#ifdef SECOND_ORDER_CENTRAL_CORRECTION
#define INTERPOL(l,c,cr) (0.5 * (c + cr))
#endif
#ifdef FIRST_ORDER_UPWIND_CORRECTION
#define INTERPOL(l,c,cr) (c)
#endif

static void flux_correction_due_to_pressure_gradients(
                                                      ConsVars* flux[3],
                                                      VectorField* buoy,
                                                      const ElemSpaceDiscr* elem,
                                                      ConsVars* Sol,
                                                      ConsVars* Sol0,
                                                      const MPV* mpv, 
                                                      double* hplus[3],
                                                      double* hgrav,
                                                      double* dp2,
                                                      const double t,
                                                      const double dt,
                                                      const double implicitness) {
    
    extern User_Data ud;
    
    const int ndim = elem->ndim;
    
    int nsp;
    
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
            
            const double dto2dx = implicitness * 0.5 * dt / elem->dx;
            const double dto2dy = implicitness * 0.5 * dt / elem->dy;
            
            ConsVars* f = flux[0];
            ConsVars* g = flux[1];
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            
            double oorhoi, ui, vi, wi, Yi, Zi, Hi, oorhoim, uim, vim, wim, Yim, Zim, Him;
            double Xi[NSPEC], Xim[NSPEC];
            double oorhoj, uj, vj, wj, Yj, Zj, Hj, oorhojm, ujm, vjm, wjm, Yjm, Zjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];
            double us, vs, ws, Hs, Ys, Zs;
            double Xs[NSPEC];
            double tmpx, tmpy, frhoY, grhoY;
            
            double p0_c, p0_s;
            double upwind;
            
            int i, j, m, n, ic, icm, icmm, icp, jc, jcm, jcmm, jcp;
            
            for(j = igy; j < icy - igy; j++) {
                m    = j * ifx;
                p0_c = mpv->HydroState->p0[j];
                
                for(i = igx; i < ifx - igx; i++) {
                    n    = m + i;
                    ic   = j*icx + i;
                    icp  = ic + 1;
                    icm  = ic - 1;
                    icmm = ic - 2;
                    
                    buoy->x[ic] = 0.0;
                    
                    /* compute interface values of u, v, Y, Z, H */
                    oorhoi  = 1.0 / INTERPOL(Sol->rhoY[icp],  Sol->rhoY[ic],  Sol->rhoY[icm]); 
                    ui      = oorhoi * INTERPOL(Sol->rhou[icp], Sol->rhou[ic], Sol->rhou[icm]);
                    vi      = oorhoi * INTERPOL(Sol->rhov[icp], Sol->rhov[ic], Sol->rhov[icm]);
                    wi      = oorhoi * INTERPOL(Sol->rhow[icp], Sol->rhow[ic], Sol->rhow[icm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xi[nsp]      = oorhoi * INTERPOL(Sol->rhoX[nsp][icp], Sol->rhoX[nsp][ic], Sol->rhoX[nsp][icm]);
                    }
                    Yi      = oorhoi * INTERPOL(Sol->rho[icp], Sol->rho[ic], Sol->rho[icm]);
                    Zi      = oorhoi * INTERPOL(Sol->rhoZ[icp], Sol->rhoZ[ic], Sol->rhoZ[icm]);
                    Hi      = oorhoi * (INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]) + p0_c);
                    
                    oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                    uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                    vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                    wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                    }
                    Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
                    Zim     = oorhoim * INTERPOL(Sol->rhoZ[icmm], Sol->rhoZ[icm], Sol->rhoZ[ic]);
                    Him     = oorhoim * (INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]) + p0_c);
                    
#ifdef STANDARD_STENCIL_PROJ1
                    tmpx = - dto2dx * ( hplusx[n] * (dp2[ic] - dp2[icm]) );
#else
                    tmpx = - dto2dx * (  hplusx[n] * (dp2[ic] - dp2[icm] )  
                                       + 0.125 * ( hplusx[n] * ((dp2[ic+icx] - dp2[icm+icx] ) - (dp2[ic] - dp2[icm]))  
                                                  - hplusx[n] * ((dp2[ic] - dp2[icm] ) - (dp2[ic-icx] - dp2[icm-icx]))  
                                                  ) 
                                       );
#endif
                    
                    frhoY = f->rhoY[n] + tmpx;
                    
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                    
                    us = upwind * uim + (1.0 - upwind) * ui;
                    vs = upwind * vim + (1.0 - upwind) * vi;  
                    ws = upwind * wim + (1.0 - upwind) * wi;  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xs[nsp]      = upwind * Xim[nsp] + (1.0 - upwind) * Xi[nsp];
                    }
                    Ys = upwind * Yim + (1.0 - upwind) * Yi;
                    Zs = upwind * Zim + (1.0 - upwind) * Zi;
                    Hs = upwind * Him + (1.0 - upwind) * Hi;
                    
                    f->rho[n]  = Ys * tmpx;
                    f->rhou[n] = us * tmpx + us * tmpx;
                    f->rhov[n] = vs * tmpx;  
                    f->rhow[n] = ws * tmpx;  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        f->rhoX[nsp][n] = Xs[nsp] * tmpx; 
                    }
                    f->rhoe[n] = 0.0 * Hs * tmpx;
                    f->rhoY[n] = tmpx;
                    f->rhoZ[n] = Zs * tmpx;
                }
            }  
            
            
            /* fluxes in the y-direction */
            for(i = igx; i < icx - igx; i++) {n = i * ify;
                for(j = igy; j < ify - igy; j++) {m = n + j;
                    jc   = j * icx + i;
                    jcp  = jc + icx;
                    jcm  = jc - icx;
                    jcmm = jc - 2*icx;
                    
                    p0_c = mpv->HydroState->p0[j];
                    p0_s = mpv->HydroState->p0[j-1];
                    
                    buoy->y[jc] = 0.0;
                    
                    oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
                    uj = oorhoj * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
                    vj = oorhoj * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
                    wj = oorhoj * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                    }
                    Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
                    Zj = oorhoj * INTERPOL(Sol->rhoZ[jcp], Sol->rhoZ[jc], Sol->rhoZ[jcm]);
                    Hj = oorhoj * (INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]) + p0_c);
                    
                    oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                    ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                    vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                    wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                    }
                    Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
                    Zjm = oorhojm * INTERPOL(Sol->rhoZ[jcmm], Sol->rhoZ[jcm], Sol->rhoZ[jc]);
                    Hjm = oorhojm * (INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]) + p0_s);
                    
#ifdef STANDARD_STENCIL_PROJ1
                    tmpy = - dto2dy * ( hplusy[m] * (dp2[jc] - dp2[jcm] ) );
#else
                    tmpy = - dto2dy * (   hplusy[m] * (dp2[jc]   - dp2[jcm] ) 
                                       + 0.125 * ( hplusy[m] * ((dp2[jc+1] - dp2[jcm+1] ) - (dp2[jc  ] - dp2[jcm  ] )) 
                                                  - hplusy[m] * ((dp2[jc  ] - dp2[jcm  ] ) - (dp2[jc-1] - dp2[jcm-1] )) 
                                                  ) 
                                       );
#endif					
                    tmpy += - 0.5 * dt * hgrav[m] * 0.5 * (dp2[jc]+dp2[jcm]); /* implicit gravity contribution */                    
                    grhoY       = g->rhoY[m] + tmpy;
                    
                    
                    
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(grhoY, 0.01) ); 
#endif
                    
                    us = upwind * ujm + (1.0 - upwind) * uj;
                    vs = upwind * vjm + (1.0 - upwind) * vj;
                    ws = upwind * wjm + (1.0 - upwind) * wj;
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xs[nsp]      = upwind * Xjm[nsp] + (1.0 - upwind) * Xj[nsp];
                    }
                    Ys = upwind * Yjm + (1.0 - upwind) * Yj;
                    Zs = upwind * Zjm + (1.0 - upwind) * Zj;
                    Hs = upwind * Hjm + (1.0 - upwind) * Hj;
                    
                    g->rho[m]  = Ys * tmpy;
                    g->rhou[m] = us * tmpy;  
                    g->rhov[m] = vs * tmpy + vs * tmpy; 
                    g->rhow[m] = ws * tmpy; 
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        g->rhoX[nsp][m] = Xs[nsp] * tmpy; 
                    }
                    g->rhoe[m] = 0.0 * Hs * tmpy;
                    g->rhoY[m] = tmpy;
                    g->rhoZ[m] = Zs * tmpy;
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
            
            const double dto2dx = implicitness * 0.5 * dt / elem->dx;
            const double dto2dy = implicitness * 0.5 * dt / elem->dy;
            const double dto2dz = implicitness * 0.5 * dt / elem->dz;
            
            const int dix = 1;
            const int diy = icx; 
            const int diz = icx*icy;            
            
            
#ifdef STANDARD_STENCIL_PROJ1
            const double cstencil  = 1.0;
#else
            const double cstencil  = 1.0/64.0;
#endif
            
            ConsVars* fx = flux[0];
            ConsVars* fy = flux[1];
            ConsVars* fz = flux[2];
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            const double* hplusz   = hplus[2];
            
            double oorhoi, ui, vi, wi, Yi, Zi, Hi; 
            double oorhoim, uim, vim, wim, Yim, Zim, Him;
            double Xi[NSPEC], Xim[NSPEC];
            
            double oorhoj, uj, vj, wj, Yj, Zj, Hj; 
            double oorhojm, ujm, vjm, wjm, Yjm, Zjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];
            
            double oorhok, uk, vk, wk, Yk, Zk, Hk; 
            double oorhokm, ukm, vkm, wkm, Ykm, Zkm, Hkm;
            double Xk[NSPEC], Xkm[NSPEC];
            
            double us, vs, ws, Hs, Ys, Zs;
            double Xs[NSPEC];
            
            double tmpx, tmpy, tmpz, frhoY;
            
            double p0_c, p0_s;
            double upwind;
            
            int i, j, k, l, m, n;
            
            int ic, icm, icmm, icp;
            int jc, jcm, jcmm, jcp;
            int kc, kcm, kcmm, kcp;
            
            for(k = igz; k < icz - igz; k++) {l = k * ifx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
                    p0_c = mpv->HydroState->p0[j];
                    for(i = igx; i < ifx - igx; i++) {n = m + i;
                        ic   = k*diz + j*diy + i*dix;
                        icp  = ic + dix; /* these appear unused as they only are accessed in macro INTERPOL() */
                        icm  = ic - dix;
                        icmm = ic - 2*dix; /* these appear unused as they only are accessed in macro INTERPOL() */
                        
                        /* compute interface values of u, v, Y, Z, H */
                        
                        oorhoi  = 1.0 / INTERPOL(Sol->rhoY[icp],  Sol->rhoY[ic],  Sol->rhoY[icm]); 
                        ui      = oorhoi * INTERPOL(Sol->rhou[icp], Sol->rhou[ic], Sol->rhou[icm]);
                        vi      = oorhoi * INTERPOL(Sol->rhov[icp], Sol->rhov[ic], Sol->rhov[icm]);
                        wi      = oorhoi * INTERPOL(Sol->rhow[icp], Sol->rhow[ic], Sol->rhow[icm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xi[nsp]      = oorhoi * INTERPOL(Sol->rhoX[nsp][icp], Sol->rhoX[nsp][ic], Sol->rhoX[nsp][icm]);
                        }
                        Yi      = oorhoi * INTERPOL(Sol->rho[icp], Sol->rho[ic], Sol->rho[icm]);
                        Zi      = oorhoi * INTERPOL(Sol->rhoZ[icp], Sol->rhoZ[ic], Sol->rhoZ[icm]);
                        Hi      = oorhoi * (INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]) + p0_c);
                        
                        oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                        uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                        vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                        wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                        }
                        Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
                        Zim     = oorhoim * INTERPOL(Sol->rhoZ[icmm], Sol->rhoZ[icm], Sol->rhoZ[ic]);
                        Him     = oorhoim * (INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]) + p0_c);
                        
#ifdef STANDARD_STENCIL_PROJ1
                        tmpx = - dto2dx * hplusx[n] * cstencil * (dp2[ic] - dp2[icm]);  
#else
                        tmpx = - dto2dx * hplusx[n] * cstencil *
                        (  36.0 *  (dp2[ic] - dp2[icm])  
                         +  6.0 * (  (dp2[ic+diy] - dp2[icm+diy]) 
                                   + (dp2[ic-diy] - dp2[icm-diy])  
                                   + (dp2[ic+diz] - dp2[icm+diz])  
                                   + (dp2[ic-diz] - dp2[icm-diz])
                                   )
                         +  1.0 * (  (dp2[ic+diy+diz] - dp2[icm+diy+diz]) 
                                   + (dp2[ic-diy+diz] - dp2[icm-diy+diz])  
                                   + (dp2[ic+diy-diz] - dp2[icm+diy-diz])  
                                   + (dp2[ic-diy-diz] - dp2[icm-diy-diz])
                                   )
                         );
#endif
                        
                        frhoY = fx->rhoY[n] + tmpx;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * uim + (1.0 - upwind) * ui;
                        vs = upwind * vim + (1.0 - upwind) * vi;  
                        ws = upwind * wim + (1.0 - upwind) * wi;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xim[nsp] + (1.0 - upwind) * Xi[nsp];
                        }
                        Ys = upwind * Yim + (1.0 - upwind) * Yi;
                        Zs = upwind * Zim + (1.0 - upwind) * Zi;
                        Hs = upwind * Him + (1.0 - upwind) * Hi;
                        
                        fx->rho[n]  = Ys * tmpx;
                        fx->rhou[n] = us * tmpx + us * tmpx;
                        fx->rhov[n] = vs * tmpx;  
                        fx->rhow[n] = ws * tmpx;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fx->rhoX[nsp][n] = Xs[nsp] * tmpx; 
                        }
                        fx->rhoe[n] = 0.0 * Hs * tmpx;
                        fx->rhoY[n] = tmpx;
                        fx->rhoZ[n] = Zs * tmpx;
                    }
                } 
            }
            
            
            /* fluxes in the y-direction */
            for(i = igx; i < icx - igx; i++) {l = i * ify*icz;
                for(k = igz; k < icz - igz; k++) {m = l + k * ify;
                    for(j = igy; j < ify - igy; j++) {n = m + j;
                        jc   = k*diz + j*diy + i*dix;
                        jcp  = jc + diy;  
                        jcm  = jc - diy;
                        jcmm = jc - 2*diy; 
                        
                        p0_c = mpv->HydroState->p0[j];
                        p0_s = mpv->HydroState->p0[j-1];
                        
                        buoy->y[jc] = 0.0;
                        
                        oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
                        uj = oorhoj  * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
                        vj = oorhoj  * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
                        wj = oorhoj  * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                        }
                        Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
                        Zj = oorhoj * INTERPOL(Sol->rhoZ[jcp], Sol->rhoZ[jc], Sol->rhoZ[jcm]);
                        Hj = oorhoj * (INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]) + p0_c);
                        
                        oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                        ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                        vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                        wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                        }
                        Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
                        Zjm = oorhojm * INTERPOL(Sol->rhoZ[jcmm], Sol->rhoZ[jcm], Sol->rhoZ[jc]);
                        Hjm = oorhojm * (INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]) + p0_s);
                        
#ifdef STANDARD_STENCIL_PROJ1
                        tmpy = - dto2dy * hplusy[n] * cstencil * (dp2[jc] - dp2[jcm]);  
#else
                        tmpy = - dto2dy * hplusy[n] * cstencil *
                        (  36.0 *    (dp2[jc] - dp2[jcm])  
                         +  6.0 * (  (dp2[jc+diz] - dp2[jcm+diz]) 
                                   + (dp2[jc-diz] - dp2[jcm-diz])  
                                   + (dp2[jc+dix] - dp2[jcm+dix])  
                                   + (dp2[jc-dix] - dp2[jcm-dix])
                                   )
                         +  1.0 * (  (dp2[jc+diz+dix] - dp2[jcm+diz+dix]) 
                                   + (dp2[jc-diz+dix] - dp2[jcm-diz+dix])  
                                   + (dp2[jc+diz-dix] - dp2[jcm+diz-dix])  
                                   + (dp2[jc-diz-dix] - dp2[jcm-diz-dix])
                                   )
                         );
#endif                                                
                        frhoY       = fy->rhoY[n] + tmpy;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * ujm + (1.0 - upwind) * uj;
                        vs = upwind * vjm + (1.0 - upwind) * vj;
                        ws = upwind * wjm + (1.0 - upwind) * wj;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xjm[nsp] + (1.0 - upwind) * Xj[nsp];
                        }
                        Ys = upwind * Yjm + (1.0 - upwind) * Yj;
                        Zs = upwind * Zjm + (1.0 - upwind) * Zj;
                        Hs = upwind * Hjm + (1.0 - upwind) * Hj;
                        
                        fy->rho[n]  = Ys * tmpy;
                        fy->rhou[n] = us * tmpy;  
                        fy->rhov[n] = vs * tmpy + vs * tmpy; 
                        fy->rhow[n] = ws * tmpy;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fy->rhoX[nsp][n] = Xs[nsp] * tmpy; 
                        }
                        fy->rhoe[n] = 0.0 * Hs * tmpy;
                        fy->rhoY[n] = tmpy;
                        fy->rhoZ[n] = Zs * tmpy;
                    }   
                }
            }
            
            /* fluxes in the z-direction */
            for(j = igy; j < icy - igy; j++) {l = j * ifz*icx;
                for(i = igx; i < icx - igx; i++) {m = l + i * ifz;
                    for(k = igz; k < ifz - igz; k++) {n = m + k;
                        kc   = k*diz + j*diy + i*dix;
                        kcp  = kc + diz;   
                        kcm  = kc - diz;
                        kcmm = kc - 2*diz;  
                        
                        p0_c = mpv->HydroState->p0[j];
                        
                        oorhok = 1.0 /INTERPOL(Sol->rhoY[kcp], Sol->rhoY[kc], Sol->rhoY[kcm]);
                        uk = oorhok * INTERPOL(Sol->rhou[kcp], Sol->rhou[kc], Sol->rhou[kcm]);
                        vk = oorhok * INTERPOL(Sol->rhov[kcp], Sol->rhov[kc], Sol->rhov[kcm]);
                        wk = oorhok * INTERPOL(Sol->rhow[kcp], Sol->rhow[kc], Sol->rhow[kcm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xk[nsp]      = oorhok * INTERPOL(Sol->rhoX[nsp][kcp], Sol->rhoX[nsp][kc], Sol->rhoX[nsp][kcm]);
                        }
                        Yk = oorhok * INTERPOL(Sol->rho[kcp], Sol->rho[kc], Sol->rho[kcm]);
                        Zk = oorhok * INTERPOL(Sol->rhoZ[kcp], Sol->rhoZ[kc], Sol->rhoZ[kcm]);
                        Hk = oorhok * (INTERPOL(Sol->rhoe[kcp], Sol->rhoe[kc], Sol->rhoe[kcm]) + p0_c);
                        
                        oorhokm = 1.0 / INTERPOL(Sol->rhoY[kcmm], Sol->rhoY[kcm], Sol->rhoY[kc]);
                        ukm = oorhokm * INTERPOL(Sol->rhou[kcmm], Sol->rhou[kcm], Sol->rhou[kc]);
                        vkm = oorhokm * INTERPOL(Sol->rhov[kcmm], Sol->rhov[kcm], Sol->rhov[kc]);
                        wkm = oorhokm * INTERPOL(Sol->rhow[kcmm], Sol->rhow[kcm], Sol->rhow[kc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xkm[nsp]      = oorhokm * INTERPOL(Sol->rhoX[nsp][kcmm], Sol->rhoX[nsp][kcm], Sol->rhoX[nsp][kc]);
                        }
                        Ykm = oorhokm * INTERPOL(Sol->rho[kcmm], Sol->rho[kcm], Sol->rho[kc]);
                        Zkm = oorhokm * INTERPOL(Sol->rhoZ[kcmm], Sol->rhoZ[kcm], Sol->rhoZ[kc]);
                        Hkm = oorhokm * (INTERPOL(Sol->rhoe[kcmm], Sol->rhoe[kcm], Sol->rhoe[kc]) + p0_c);
                        
#ifdef STANDARD_STENCIL_PROJ1
                        tmpz = - dto2dz * hplusz[n] * cstencil * (dp2[kc] - dp2[kcm]);  
#else
                        tmpz = - dto2dz * hplusz[n] * cstencil *
                        (  36.0 *    (dp2[kc] - dp2[kcm])  
                         +  6.0 * (  (dp2[kc+diz] - dp2[kcm+diz]) 
                                   + (dp2[kc-diz] - dp2[kcm-diz])  
                                   + (dp2[kc+dix] - dp2[kcm+dix])  
                                   + (dp2[kc-dix] - dp2[kcm-dix])
                                   )
                         +  1.0 * (  (dp2[kc+diz+dix] - dp2[kcm+diz+dix]) 
                                   + (dp2[kc-diz+dix] - dp2[kcm-diz+dix])  
                                   + (dp2[kc+diz-dix] - dp2[kcm+diz-dix])  
                                   + (dp2[kc-diz-dix] - dp2[kcm-diz-dix])
                                   )
                         );
#endif                        
                        frhoY       = fz->rhoY[n] + tmpz;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * ukm + (1.0 - upwind) * uk;
                        vs = upwind * vkm + (1.0 - upwind) * vk;
                        ws = upwind * wkm + (1.0 - upwind) * wk;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xkm[nsp] + (1.0 - upwind) * Xk[nsp];
                        }
                        Ys = upwind * Ykm + (1.0 - upwind) * Yk;
                        Zs = upwind * Zkm + (1.0 - upwind) * Zk;
                        Hs = upwind * Hkm + (1.0 - upwind) * Hk;
                        
                        fz->rho[n]  = Ys * tmpz;
                        fz->rhou[n] = us * tmpz;  
                        fz->rhov[n] = vs * tmpz; 
                        fz->rhow[n] = ws * tmpz + ws * tmpz;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fz->rhoX[nsp][n] = Xs[nsp] * tmpz; 
                        }
                        fz->rhoe[n] = 0.0 * Hs * tmpz;
                        fz->rhoY[n] = tmpz;
                        fz->rhoZ[n] = Zs * tmpz;
                    }   
                }
            }			
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
}

/* ========================================================================== */

static void flux_correction_due_to_pressure_values(
                                                   ConsVars* flux[3],
                                                   VectorField* buoy,
                                                   const ElemSpaceDiscr* elem,
                                                   ConsVars* Sol,
                                                   double* dp2, 
                                                   const double dt) {
    
    extern User_Data ud;
    
    const int ndim = elem->ndim;
    
    switch(ndim) {
        case 1: {
            ERROR("function not available");
            break;
        }
            
        case 2: {
            int i, j, m, n, ic, icm, jc, jcm;
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int ifx = elem->ifx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int ify = elem->ify;
            
            ConsVars* f = flux[0];
            ConsVars* g = flux[1];
            
            /* flux_post_correction    change from   flux ... +=  corr    to    flux ... = corr */
            
            /* Momentum flux corrections */

            /* fluxes in the x-direction */
            for(j = igy; j < icy - igy; j++) {
                
                m = j * ifx;
                
                for(i = igx; i < ifx - igx; i++) {
                    
                    n   = m + i;
                    ic  = n - j;
                    icm = ic - 1;
                    
                    f->rhou[n] += ud.p_extrapol*0.5*(
                                                     (ud.latw[0]*dp2[ic+icx]  + ud.latw[1]*dp2[ic]  + ud.latw[2]*dp2[ic-icx]) 
                                                     + (ud.latw[0]*dp2[icm+icx] + ud.latw[1]*dp2[icm] + ud.latw[2]*dp2[icm-icx]) 
                                                     );
                }
            }  
            
            
            /* flux_post_correction    change from   flux ... +=  corr    to    flux ... = corr */
            
            /* fluxes in the y-direction */
            for(i = igx; i < icx - igx; i++) {
                
                n = i * ify;
                
                for(j = igy; j < ify - igy; j++) {
                    
                    m   = n + j;
                    jc  = j * icx + i;
                    jcm = jc - icx;
                    
                    g->rhov[m] += ud.p_extrapol*0.5*(  (ud.latw[0]*dp2[jc+1]  + ud.latw[1]*dp2[jc]  + ud.latw[2]*dp2[jc-1])  
                                                     + (ud.latw[0]*dp2[jcm+1] + ud.latw[1]*dp2[jcm] + ud.latw[2]*dp2[jcm-1])
                                                     );
                }   
            }
            
            break;
        }
        case 3: {
            ERROR("function not available");
            break;
        }
        default: ERROR("ndim not in {1, 2,3}");
    }
#if 0
    printf("leaving flux_correction_due_to_pressure_values()\n");
#endif
}


#else /* CORRECTION_FROM_MG_CODE */

/* ========================================================================== */
#if 1
#define CORR factor * \
( (1-acfl)  * \
(  f->rhoY[nf+1] * (upw_ri * q_ce + (1-upw_ri) * q_ri) \
- f->rhoY[nf]   * (upw_le * q_le + (1-upw_le) * q_ce) \
+ g->rhoY[ng+1] * (upw_to * q_ce + (1-upw_to) * q_to) \
- g->rhoY[ng]   * (upw_bo * q_bo + (1-upw_bo) * q_ce) \
) \
+  acfl * q_ce * delrhoY \
)
#endif

static void correction(
                       ConsVars* dSol,
                       ConsVars* Sol,
                       const ConsVars* Sol0,
                       ConsVars* flux[3],
                       double* rhs,
                       const ElemSpaceDiscr* elem,
                       const MPV* mpv,
                       const double t,
                       const double dt) {
    
    /*
     extern Thermodynamic th;
     const double gm1 = th.gm1;
     const double tol = DBL_MAX;
     */
    extern User_Data ud;
    
    const int compressibility = ud.acoustic_timestep;
    const double dx = elem->dx;
    const double dy = elem->dy;
    
    double factor = 0.5*dt*dt / (dx*dy);
    double del_norm = 0.0;
    double acfl, delrhoY;
    int count;
    
    switch(elem->ndim) {
        case 1: {
            ERROR("correction from first projection not implemented for 1D");
            break;
        }
            
        case 2: {
            int i, j, m, n, mf, nf, mg, ng;
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int ifx = elem->ifx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int ify = elem->ify;
            
            const ConsVars* f = flux[0];
            const ConsVars* g = flux[1];
            
            for(j = igy; j < icy - igy; j++) {m = j * icx; mf = j*ifx; mg = j;
                for(i = igx; i < icx - igx; i++) {n = m + i; nf = mf + i; ng = mg + i*ify;
                    double drhoY_r = Sol->rhoY[n] - Sol0->rhoY[n];
                    double drhoY_c = factor * ((f->rhoY[nf+1] - f->rhoY[nf]) + (g->rhoY[ng+1] - g->rhoY[ng]));
                    double del, del2;
                    double q_ce, q_le, q_ri, q_bo, q_to;
                    
                    acfl = 0.0;
                    delrhoY = 0.0;
                    
                    int upw_le = (1 + SIGNint(-f->rhoY[nf]  )) / 2; 
                    int upw_ri = (1 + SIGNint(-f->rhoY[nf+1])) / 2; 
                    int upw_bo = (1 + SIGNint(-g->rhoY[ng]  )) / 2; 
                    int upw_to = (1 + SIGNint(-g->rhoY[ng+1])) / 2; 
                    
                    /* correction of the to-be-controlled variable */
                    del = drhoY_r + drhoY_c;
                    
                    /* advective corrections for all remaining conservative variables */
                    q_ce = 1.0;
                    q_le = 1.0;
                    q_ri = 1.0;
                    q_bo = 1.0;
                    q_to = 1.0;
                    
                    dSol->rhoY[n] = CORR;
                    del2 = (Sol->rhoY[n] + dSol->rhoY[n]) - Sol0->rhoY[n];
                    
                    del_norm += del2*del2;
                    count++;
                    
                    q_ce = Sol->rho[n]     / Sol->rhoY[n];
                    q_le = Sol->rho[n-1]   / Sol->rhoY[n-1];
                    q_ri = Sol->rho[n+1]   / Sol->rhoY[n+1];
                    q_bo = Sol->rho[n-icx] / Sol->rhoY[n-icx];
                    q_to = Sol->rho[n+icx] / Sol->rhoY[n+icx];
                    
                    dSol->rho[n] = CORR;
                    
                    /* momentum correction due to dp2-gradients?? */
                    
                    q_ce = Sol->rhou[n]     / Sol->rhoY[n];
                    q_le = Sol->rhou[n-1]   / Sol->rhoY[n-1];
                    q_ri = Sol->rhou[n+1]   / Sol->rhoY[n+1];
                    q_bo = Sol->rhou[n-icx] / Sol->rhoY[n-icx];
                    q_to = Sol->rhou[n+icx] / Sol->rhoY[n+icx];
                    
                    dSol->rhou[n] = CORR;
                    
                    q_ce = Sol->rhov[n]     / Sol->rhoY[n];
                    q_le = Sol->rhov[n-1]   / Sol->rhoY[n-1];
                    q_ri = Sol->rhov[n+1]   / Sol->rhoY[n+1];
                    q_bo = Sol->rhov[n-icx] / Sol->rhoY[n-icx];
                    q_to = Sol->rhov[n+icx] / Sol->rhoY[n+icx];
                    
                    dSol->rhov[n] = CORR;
                    
                    /* advective corrections for all remaining conservative variables */
                    q_ce = Sol->rhow[n]     / Sol->rhoY[n];
                    q_le = Sol->rhow[n-1]   / Sol->rhoY[n-1];
                    q_ri = Sol->rhow[n+1]   / Sol->rhoY[n+1];
                    q_bo = Sol->rhow[n-icx] / Sol->rhoY[n-icx];
                    q_to = Sol->rhow[n+icx] / Sol->rhoY[n+icx];
                    
                    dSol->rhow[n] = CORR;
                    
                    /* advective corrections for all remaining conservative variables */
                    q_ce = Sol->rhoZ[n]     / Sol->rhoY[n];
                    q_le = Sol->rhoZ[n-1]   / Sol->rhoY[n-1];
                    q_ri = Sol->rhoZ[n+1]   / Sol->rhoY[n+1];
                    q_bo = Sol->rhoZ[n-icx] / Sol->rhoY[n-icx];
                    q_to = Sol->rhoZ[n+icx] / Sol->rhoY[n+icx];
                    
                    dSol->rhoZ[n] = compressibility * CORR + (1-compressibility)*((Sol->rhoZ[n] / Sol->rho[n]) * dSol->rho[n]);
                    
                }
            }  
            
            break;
        }
        case 3: 
            ERROR("correction from first projection not implemented for 3D");
            
        default: ERROR("ndim not in {1, 2, 3}");
    }
    
    Explicit_step_update(Sol, elem->nc);
}
#endif  /* CORRECTION_FROM_MG_CODE */

#else /* SOLVER_1_HYPRE */

static enum Constraint integral_condition(
										  ConsVars* flux[3],
										  double* rhs, 
										  ConsVars* Sol,
										  const double dt,
										  const ElemSpaceDiscr* elem,
										  MPV* mpv);

static void controlled_variable_change_explicit(
												double* rhs, 
												const ElemSpaceDiscr* elem, 
												ConsVars* Sol_new, 
												ConsVars* Sol_old, 
												double dt,
												const MPV* mpv);

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
													  VectorField* buoy,
													  const ElemSpaceDiscr* elem,
													  ConsVars* Sol,
													  ConsVars* Sol0,
													  const MPV* mpv, 
													  double* hplus[3],
													  double* hgrav,
													  double* dp2,
													  const double t,
													  const double dt,
													  const double implicitness);

static void flux_correction_due_to_pressure_values(
												   ConsVars* flux[3],
												   VectorField* buoy,
												   const ElemSpaceDiscr* elem,
												   ConsVars* Sol,
												   double* dp2, 
												   const double dt);

/* ========================================================================== */

#if 1
static int rhs_output_count = 0;
#endif

void flux_correction(
					 ConsVars* flux[3],
					 VectorField* buoy,
					 const ElemSpaceDiscr* elem_base_grid,
					 ConsVars* Sol, 
					 ConsVars* Sol0, 
					 const double t,
					 const double dt,
					 const double implicitness,
                     const int step) {
	
    extern User_Data ud;
	extern MPV* mpv;
	
	const ElemSpaceDiscr* elem = mpv->Level[0]->elem;
    const NodeSpaceDiscr* node = mpv->Level[0]->node;

	double** hplus       = mpv->Level[0]->wplus;
	double*  hcenter     = mpv->Level[0]->wcenter;
	double*  hgrav       = mpv->Level[0]->wgrav;
	
	double* rhs          = mpv->Level[0]->rhs;
	double* dp2		     = mpv->Level[0]->p;
    
    rhs_output_count     = 0;
    
	int n;
    
#ifdef NEWTON_FOR_EQ_OF_STATE
    extern Thermodynamic th;
    extern double *W0;
    double *ddp = W0;
    double Msq  = ud.Msq;
#endif
    
    printf("\n\n====================================================");
    printf("\nFirst Projection");
    printf("\n====================================================\n");
	    
	operator_coefficients(hplus, hcenter, hgrav, elem, Sol, Sol0, mpv, dt);
    
	/* obtain finest grid data for the MG hierarchy */
	controlled_variable_change_explicit(rhs, elem, Sol, Sol0, dt, mpv);
    
    assert(integral_condition(flux, rhs, Sol, dt, elem, mpv) != VIOLATED); 
    rhs_fix_for_open_boundaries(rhs, elem, Sol, Sol0, flux, dt, mpv);
    
#if 1
    extern User_Data ud;
    FILE *prhsfile = NULL;
    char fn2[200], fieldname2[90];
    
    sprintf(fn2, "%s/rhs_cells/rhs_cells_00%d.hdf", ud.file_name, rhs_output_count);
    rhs_output_count++;
    sprintf(fieldname2, "rhs_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             rhs,
             fn2,
             fieldname2);
    
    /*
    sprintf(fn2, "%s/Tests/p2_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "p2_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             dp2,
             fn2,
             fieldname2);

    sprintf(fn2, "%s/Tests/hx_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hx_c");
    
    WriteHDF(prhsfile,
             2*(mpv->Level[0]->elem->icy-2*mpv->Level[0]->elem->igy),
             mpv->Level[0]->elem->ifx-2*mpv->Level[0]->elem->igx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             h[0],
             fn2,
             fieldname2);
    
    sprintf(fn2, "%s/Tests/hy_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hy_c");
    
    WriteHDF(prhsfile,
             2*(mpv->Level[0]->elem->ify-2*mpv->Level[0]->elem->igy),
             mpv->Level[0]->elem->icx-2*mpv->Level[0]->elem->igx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             h[1],
             fn2,
             fieldname2);
    
    sprintf(fn2, "%s/Tests/hplusx_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hplusx_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->ifx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             hplus[0],
             fn2,
             fieldname2);
    
    sprintf(fn2, "%s/Tests/hplusy_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hplusy_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->ify,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             hplus[1],
             fn2,
             fieldname2);
    
 */
#endif

    
    /**/
#ifdef NEWTON_FOR_EQ_OF_STATE
    {
        double del_nwt = 1.0e10;
        double tol_nwt = 1.0e-6;
        
        for (int ic=0; ic<elem->nc; ic++) dp2[ic] = ddp[ic] = 0.0;

        while (del_nwt > tol_nwt) {
            
            /* make sure rhs contains the residue for the linear solve upon exit */
            variable_coefficient_poisson_cells(ddp, rhs, (const double **)hplus, hcenter, hgrav, Sol, elem, node, implicitness);
            set_ghostcells_p2(ddp, (const double **)hplus, hgrav, elem, elem->igx);
            for (int ic=0; ic<elem->nc; ic++) {
                double pold  = pow(Sol0->rhoY[ic],th.gamm);
                double pnext = pold+Msq*(dp2[ic]+ddp[ic]);
                double plast = pold+Msq*dp2[ic];
                double dPdp  = th.gamminv / pow(pold, th.Gamma);
                rhs[ic] += (pow(pnext, th.gamminv)-pow(plast, th.gamminv) - dPdp*Msq*ddp[ic])/dt;
                del_nwt  = MAX_own(del_nwt, rhs[ic]); 
                dp2[ic] += ddp[ic];
                ddp[ic]  = 0.0;
            }
            printf("Newton: del_nwt = %e\n", del_nwt);
        }
    }
#else
    variable_coefficient_poisson_cells(dp2, rhs, (const double **)hplus, hcenter, hgrav, Sol, elem, node, implicitness);
#endif

    set_ghostcells_p2(dp2, (const double **)hplus, hgrav, elem, elem->igx);

    flux_correction_due_to_pressure_gradients(flux, buoy, elem, Sol, Sol0, mpv, hplus, hgrav, dp2, t, dt, implicitness);
    if (ud.p_flux_correction) {
        flux_correction_due_to_pressure_values(flux, buoy, elem, Sol, dp2, dt); 
    }

    flux_fix_for_open_boundaries(flux, elem, mpv);  

    memcpy(mpv->dp2_cells, dp2, elem->nc*sizeof(double));

    /* store results in mpv-fields */        
    for(n=0; n<elem->nc; n++) mpv->p2_cells[n] = MAC_PROJ_OLDP_WEIGHT*mpv->p2_cells[n] + MAC_PROJ_DELP_WEIGHT*dp2[n];
    
	set_ghostcells_p2(mpv->p2_cells, (const double **)hplus, hgrav, elem, elem->igx);
}

/* ========================================================================== */
 
static void controlled_variable_change_explicit(
                                                double* rhs, 
                                                const ElemSpaceDiscr* elem, 
                                                ConsVars* Sol_new, 
                                                ConsVars* Sol_old, 
                                                double dt,
                                                const MPV* mpv) {
        
    memset(rhs, 0.0, elem->nc*sizeof(double));
    
    const int igx = elem->igx;
    const int icx = elem->icx;
    const int igy = elem->igy;
    const int icy = elem->icy;
    const int igz = elem->igz;
    const int icz = elem->icz;
    
    const double factor = 2.0 / (dt * dt);
        
    int i, j, k, l, m, n;
    
    for(i=0; i<elem->nc; i++) rhs[i] = 0.0;
        
    for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                double rhs_value = - factor * (Sol_new->rhoY[n] - Sol_old->rhoY[n]);
                rhs[n]   = rhs_value;
            }
        }
    }
    
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "rhs_cells.hdf");
    sprintf(fieldname, "rhs-cells");    
    WriteHDF(prhsfile, elem->icx, elem->icy, elem->icz, elem->ndim, rhs, fn, fieldname);
#endif
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
			int i, j, jifx, iify;
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int ifx = elem->ifx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int ify = elem->ify;
			const ConsVars* f = flux[0];
			const ConsVars* g = flux[1];
			double F, G;
			double F_rho, G_rho;
			double F_rhoY, G_rhoY;
						
			/* this implementation does not allow for non-zero global flux divergence */
			for(j = igy, F = 0.0, F_rho = 0.0, F_rhoY = 0.0; j < icy - igy; j++) {
				jifx = j * ifx;
				F      += f->rhoe[jifx + icx - igx] - f->rhoe[jifx +  igx];
				F_rho  += f->rho[jifx + icx - igx]  - f->rho[jifx +  igx];
				F_rhoY += f->rhoY[jifx + icx - igx] - f->rhoY[jifx +  igx];
			} 
			for(i = igx, G = 0.0, G_rho = 0.0, G_rhoY = 0.0; i < icx - igx; i++) {
				iify =  i * ify;
				G      += g->rhoe[iify + icy - igy] - g->rhoe[iify + igy];
				G_rho  += g->rho[iify + icy - igy]  - g->rho[iify + igy];
				G_rhoY += g->rhoY[iify + icy - igy] - g->rhoY[iify + igy];
			}  
						
            /*
            if(0){
				printf("F_rhoY = %e, G_rhoY = %e, DBL_EPSILON = %e\n", F_rhoY, G_rhoY, DBL_EPSILON);
				return VIOLATED;
			}
			else
             */
            {
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
					printf("integral_condition_cells = OK; tmp = %e, \n F = %e, F_rho = %e, F_rhoY = %e, \n G = %e, G_rho = %e, G_rhoY = %e\n", fabs(tmp), F, F_rho, F_rhoY, G, G_rho, G_rhoY);
					return(SATISFIED);
				}
				else {
					printf("integral_condition_cells = VIOLATED; tmp = %e, rhs_max = %e\n F = %e, F_rho = %e, F_rhoY = %e \n G = %e, G_rho = %e, G_rhoY = %e\n", fabs(tmp), rhs_max, F, F_rho, F_rhoY, G, G_rho, G_rhoY);
					return(SATISFIED);
			
                }
			}	
			break;
		}
			
		case 3: {
			int i, j, k, l, m, n;
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int ifx = elem->ifx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int ify = elem->ify;
			const int igz = elem->igz;
			const int icz = elem->icz;
			const int ifz = elem->ifz;
			const int ifxicy = ifx * icy;
			const int ifyicz = ify * icz;
			const int ifzicx = ifz * icx;
			const ConsVars* f = flux[0];
			const ConsVars* g = flux[1];
			const ConsVars* h = flux[2];
			double F, G, H;
			
			printf("\nSanity check for rho-flux and rhoe-flux not implemented in 3D \n");
			
			for(k = igz, F = 0.0; k < icz - igz; k++) {l = k * ifxicy;
				for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
					F += f->rhoY[m + icx - igx] - f->rhoY[m +  igx];
				} 
			}
			for(i = igx, G = 0.0; i < icx - igx; i++) {m = i * ifyicz;
				for(k = igz; k < icz - igz; k++) {n = m + k * ify;
					G += g->rhoY[n + icy - igy] - g->rhoY[n + igy];
				} 
			}
			for(j = igy, H = 0.0; j < icy - igy; j++) {n = j * ifzicx;
				for(i = igx; i < icx - igx; i++) {l = n + i * ifz;
					H += h->rhoY[l + icz - igz] - h->rhoY[l + igz];
				}
			}  
			if(fabs(F * elem->dy + G * elem->dx * H * elem->dz) > 100*DBL_EPSILON) 
				return VIOLATED;
			else
				return SATISFIED;
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
                           double* wgrav,
                           const ElemSpaceDiscr* elem,
                           const ConsVars* Sol,
                           const ConsVars* Sol0,
                           const MPV* mpv, 
                           const double dt) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	
	const int ndim = elem->ndim;
	
    const int impl_grav_th = ud.implicit_gravity_theta;
    const int impl_grav_pr = ud.implicit_gravity_press;
    
	const double implicitness = ud.implicitness;
    const double ccw = 2.0; /* ccenterweight   4.0 */
	const double ccenter = - ccw * (ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt); 
    /* for ccw = 4.0: 
       first factor 2.0 from reinterpretation of Helmholtz-solution as 
       (p^{n+1}-p^{n}) and SECOND-ORDER update to t^{n+1/2} in the first projection 
     */
	const double cexp    = 1.0-th.gamm;
	
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
            const double dx = elem->dx;
            const double dy = elem->dy;
			double* hx = hplus[0];
			double* hy = hplus[1];
			double* hc = wcenter;
			double* hg = wgrav;
			
			double hi, him, hj, hjm, g, gimp, thet, thetm, Msq;
			
			int i, j, m, n, ic, icm, jc, jcm;
            
            Msq = ud.Msq;
            g   = ud.gravity_strength[0];
            
			for(j = igy-1; j < icy - igy+1; j++) {
                m = j * ifx;
				
                for(i = igx-1; i < ifx - igx+1; i++) {
                    n     = m + i;
                    ic    = n - j;
					icm   = ic - 1; 
					
#ifdef OLD_TIME_COEFFS
                    hi    = Sol0->rhoY[ic]  / Sol0->rho[ic];   
                    him   = Sol0->rhoY[icm] / Sol0->rho[icm];
                    thet  = Sol0->rhoY[ic]  / Sol0->rho[ic] ;
                    thetm = Sol0->rhoY[icm] / Sol0->rho[icm];
#else
                    hi    = 0.5 * ( Sol->rhoY[ic] / Sol->rho[ic] + Sol0->rhoY[ic] / Sol0->rho[ic] );   
					him   = 0.5 * ( Sol->rhoY[icm] / Sol->rho[icm] + Sol0->rhoY[icm] / Sol0->rho[icm] );
                    thet  = Sol->rhoY[ic]  / Sol->rho[ic] ;
                    thetm = Sol->rhoY[icm] / Sol->rho[icm];
#endif
                    gimp  = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dx*0.5*(thet+thetm)));
                    
					hx[n] = 0.5 * (hi + him) * implicitness * gimp;
                    
					assert(hx[n] > 0.0);
				}
			}
			
            
            g = ud.gravity_strength[1];
            
			for(i = 0; i < icx; i++) {
                n = i * ify;
				
                for(j = igy-1; j < ify - igy+1; j++) {
                    m     = n + j;
                    jc    = j * icx + i;
					jcm   = jc - icx;          
					
#ifdef OLD_TIME_COEFFS
                    hj    = Sol0->rhoY[jc]  / Sol0->rho[jc];
                    hjm   = Sol0->rhoY[jcm] / Sol0->rho[jcm];
                    thet  = Sol0->rhoY[jc]  / Sol0->rho[jc] ;
                    thetm = Sol0->rhoY[jcm] / Sol0->rho[jcm];
#else
					hj    = 0.5 * ( Sol->rhoY[jc] / Sol->rho[jc] + Sol0->rhoY[jc] / Sol0->rho[jc]);
					hjm   = 0.5 * ( Sol->rhoY[jcm] / Sol->rho[jcm] + Sol0->rhoY[jcm] / Sol0->rho[jcm]);
                    thet  = Sol->rhoY[jc]  / Sol->rho[jc] ;
                    thetm = Sol->rhoY[jcm] / Sol->rho[jcm];
                    /*
                    thet  = mpv->HydroState->Y0[j];
                    thetm = mpv->HydroState->Y0[j-1];
                     */
#endif
                     
                    gimp  = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dy*0.5*(thet+thetm)));
                    
					hy[m] = 0.5 * (hj + hjm) * implicitness * gimp;
                    
                    hg[m] = th.gamminv * pow(0.5*(Sol->rhoY[jc]+Sol->rhoY[jcm]),cexp) * g * impl_grav_pr * gimp; 
                    /* hg[m] = th.gamminv * pow(0.5*(Sol0->rhoY[jc]+Sol0->rhoY[jcm]),cexp) * (g/Msq) * impl_grav_pr * gimp; */
                    
					assert(hy[m] > 0.0);
				}
			}
            
			for(j = igy; j < icy - igy; j++) {m = j * icx;
				for(i = igx; i < icx - igx; i++) {n = m + i;
					hc[n] = ccenter * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp);
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
            const double dx = elem->dx;
            const double dy = elem->dy;
            const double dz = elem->dz;
			double* hx = hplus[0];
			double* hy = hplus[1];
			double* hz = hplus[2];
			double* hc = wcenter;
			double* hg = wgrav;
			
			double hi, him, hj, hjm, hk, hkm, g, gimp, thet, thetm, Msq;
			
			int i, j, k, l, m, n, ic, icm, jc, jcm, kc, kcm;
            
            Msq = ud.Msq;
            g   = ud.gravity_strength[0];
            assert(g==0.0); /* implicit gravity only for y-direction */
            
			for(k = igz; k < icz - igz; k++) {l = k * ifx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
                    for(i = igx; i < ifx - igx; i++) {n = m + i;
                        ic  = k*icx*icy + j*icx + i;
                        icm = ic - 1; 
					
                        hi       = 0.5 * ( Sol->rhoY[ic] / Sol->rho[ic] 
                                          + Sol0->rhoY[ic] / Sol0->rho[ic] );   
					
                        him      = 0.5 * ( Sol->rhoY[icm] / Sol->rho[icm] 
                                          + Sol0->rhoY[icm] / Sol0->rho[icm] );
					
                        thet  = Sol->rhoY[ic]  / Sol->rho[ic] ;
                        thetm = Sol->rhoY[icm] / Sol->rho[icm];
                        gimp   = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dx*0.5*(thet+thetm)));
                    
                        hx[n] = 0.5 * (hi + him) * implicitness * gimp;
                        assert(hx[n] > 0.0);
                    }
                }
            }
			
            
            g = ud.gravity_strength[1];
            
			for(i = igx; i < icx - igx; i++) {l = i * ify*icz;
                for(k = igz; k < icz - igz; k++) {m = l + k * ify;
                    for(j = igy; j < ify - igy; j++) {n = m + j;
                        jc  = k*icx*icy + j*icx + i;
                        jcm = jc - icx;          
					
                        hj       = 0.5 * ( Sol->rhoY[jc] / Sol->rho[jc] 
                                          + Sol0->rhoY[jc] / Sol0->rho[jc]);
					
                        hjm      = 0.5 * ( Sol->rhoY[jcm] / Sol->rho[jcm] 
                                          + Sol0->rhoY[jcm] / Sol0->rho[jcm]);
					
                        thet  = Sol->rhoY[jc]  / Sol->rho[jc] ;
                        thetm = Sol->rhoY[jcm] / Sol->rho[jcm];
                        gimp   = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dy*0.5*(thet+thetm)));
                    
                        hy[n] = 0.5 * (hj + hjm) * implicitness * gimp;
                        
                        hg[m] = th.gamminv * pow(0.5*(Sol->rhoY[jc]+Sol->rhoY[jcm]),cexp) * g * implicitness * gimp;
                        /* hg[m] = th.gamminv * pow(0.5*(Sol->rhoY[jc]+Sol->rhoY[jcm]),cexp) * (g/Msq) * implicitness * gimp; */

                        assert(hy[n] > 0.0);
                    }
                }
            }

            g = ud.gravity_strength[2];
            assert(g==0.0); /* implicit gravity only for y-direction */
            
			for(j = igy; j < icy - igy; j++) {l = j * ifz*icx;
				for(i = igx; i < icx - igx; i++) {m = l + i * ifz;
                    for(k = igz; k < ifz - igz; k++) {n = m + k;
                        kc  = k*icx*icy + j*icx + i;
                        kcm = kc - icx*icy;          
					
                        hk       = 0.5 * ( Sol->rhoY[kc] / Sol->rho[kc] 
                                          + Sol0->rhoY[kc] / Sol0->rho[kc]);
					
                        hkm      = 0.5 * ( Sol->rhoY[kcm] / Sol->rho[kcm] 
                                          + Sol0->rhoY[kcm] / Sol0->rho[kcm]);
					
                        thet  = Sol->rhoY[kc]  / Sol->rho[kc] ;
                        thetm = Sol->rhoY[kcm] / Sol->rho[kcm];
                        gimp   = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dz*0.5*(thet+thetm)));
                    
                        hz[n] = 0.5 * (hk + hkm) * implicitness * gimp;
                        assert(hz[n] > 0.0); 
                    }
                }
			}
            
            for(k = igz; k < icz - igz; k++) {l = k*icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j*icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        hc[n] = ccenter * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp);
                    }
                }
            }
			
			
			break;
		}
		default: ERROR("ndim not in {1,2,3}");
	}
}

/* ========================================================================== */
#ifdef THIRD_ORDER_UPWIND_CORRECTION
#define INTERPOL(l,c,cr) (FIVE_SIXTHS * c + ONE_THIRD * cr - ONE_SIXTH * l)
#endif
#ifdef SECOND_ORDER_CENTRAL_CORRECTION
#define INTERPOL(l,c,cr) (0.5 * (c + cr))
#endif
#ifdef FIRST_ORDER_UPWIND_CORRECTION
#define INTERPOL(l,c,cr) (c)
#endif

#ifdef SHIFTED_COEFFICIENTS_PROJ1

static void flux_correction_due_to_pressure_gradients(
                                                      ConsVars* flux[3],
                                                      VectorField* buoy,
                                                      const ElemSpaceDiscr* elem,
                                                      ConsVars* Sol,
                                                      ConsVars* Sol0,
                                                      const MPV* mpv, 
                                                      double* hplus[3],
                                                      double* hgrav,
                                                      double* dp2,
                                                      const double t,
                                                      const double dt,
                                                      const double implicitness) {
    
    extern User_Data ud;
    
    const int ndim = elem->ndim;
    
    int nsp;
    
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
            
            const double dto2dx = implicitness * 0.5 * dt / elem->dx;
            const double dto2dy = implicitness * 0.5 * dt / elem->dy;
            
            ConsVars* f = flux[0];
            ConsVars* g = flux[1];
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            
            double oorhoi, ui, vi, wi, Yi, Zi, Hi, oorhoim, uim, vim, wim, Yim, Zim, Him;
            double Xi[NSPEC], Xim[NSPEC];
            double oorhoj, uj, vj, wj, Yj, Zj, Hj, oorhojm, ujm, vjm, wjm, Yjm, Zjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];
            double us, vs, ws, Hs, Ys, Zs;
            double Xs[NSPEC];
            double tmpx, tmpy, frhoY, grhoY;
            
            double p0_c, p0_s;
            double upwind;
            
            int i, j, mc, me, mw, nc, nn, ns, ic, icm, icmm, icp, jc, jcm, jcmm, jcp;
            
            for(j = igy; j < icy - igy; j++) {
                mc    = j * ifx;
                p0_c = mpv->HydroState->p0[j];
                
                for(i = igx; i < ifx - igx; i++) {
                    nn   = mc + i + ifx;
                    nc   = mc + i;
                    ns   = mc + i - ifx;
                    ic   = j*icx + i;
                    icp  = ic + 1;
                    icm  = ic - 1;
                    icmm = ic - 2;
                    
                    buoy->x[ic] = 0.0;
                    
                    /* compute interface values of u, v, Y, Z, H */
                    oorhoi  = 1.0 / INTERPOL(Sol->rhoY[icp],  Sol->rhoY[ic],  Sol->rhoY[icm]); 
                    ui      = oorhoi * INTERPOL(Sol->rhou[icp], Sol->rhou[ic], Sol->rhou[icm]);
                    vi      = oorhoi * INTERPOL(Sol->rhov[icp], Sol->rhov[ic], Sol->rhov[icm]);
                    wi      = oorhoi * INTERPOL(Sol->rhow[icp], Sol->rhow[ic], Sol->rhow[icm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xi[nsp]      = oorhoi * INTERPOL(Sol->rhoX[nsp][icp], Sol->rhoX[nsp][ic], Sol->rhoX[nsp][icm]);
                    }
                    Yi      = oorhoi * INTERPOL(Sol->rho[icp], Sol->rho[ic], Sol->rho[icm]);
                    Zi      = oorhoi * INTERPOL(Sol->rhoZ[icp], Sol->rhoZ[ic], Sol->rhoZ[icm]);
                    Hi      = oorhoi * (INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]) + p0_c);
                    
                    oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                    uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                    vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                    wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                    }
                    Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
                    Zim     = oorhoim * INTERPOL(Sol->rhoZ[icmm], Sol->rhoZ[icm], Sol->rhoZ[ic]);
                    Him     = oorhoim * (INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]) + p0_c);
                    
#ifdef STANDARD_STENCIL_PROJ1
                    tmpx = - dto2dx * ( hplusx[nc] * (dp2[ic] - dp2[icm]) );
#else
                    tmpx = - dto2dx * (  0.75  *   hplusx[nc] * (dp2[ic]     - dp2[icm]    )  
                                       + 0.125 * ( hplusx[nn] * (dp2[ic+icx] - dp2[icm+icx])  
                                                 + hplusx[ns] * (dp2[ic-icx] - dp2[icm-icx])  
                                                  ) 
                                       );
#endif
                    
                    frhoY = f->rhoY[nc] + tmpx;
                    
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                    
                    us = upwind * uim + (1.0 - upwind) * ui;
                    vs = upwind * vim + (1.0 - upwind) * vi;  
                    ws = upwind * wim + (1.0 - upwind) * wi;  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xs[nsp]      = upwind * Xim[nsp] + (1.0 - upwind) * Xi[nsp];
                    }
                    Ys = upwind * Yim + (1.0 - upwind) * Yi;
                    Zs = upwind * Zim + (1.0 - upwind) * Zi;
                    Hs = upwind * Him + (1.0 - upwind) * Hi;
                    
                    f->rho[nc]  = Ys * tmpx;
                    f->rhou[nc] = us * tmpx + us * tmpx;
                    f->rhov[nc] = vs * tmpx;  
                    f->rhow[nc] = ws * tmpx;  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        f->rhoX[nsp][nc] = Xs[nsp] * tmpx; 
                    }
                    f->rhoe[nc] = 0.0 * Hs * tmpx;
                    f->rhoY[nc] = tmpx;
                    f->rhoZ[nc] = Zs * tmpx;
                }
            }  
            
            
            /* fluxes in the y-direction */
            for(i = igx; i < icx - igx; i++) {
                nc = i * ify;
                
                for(j = igy; j < ify - igy; j++) {
                    me = nc + j + ify;
                    mc = nc + j;
                    mw = nc + j - ify;
                    jc   = j * icx + i;
                    jcp  = jc + icx;
                    jcm  = jc - icx;
                    jcmm = jc - 2*icx;
                    
                    p0_c = mpv->HydroState->p0[j];
                    p0_s = mpv->HydroState->p0[j-1];
                    
                    buoy->y[jc] = 0.0;
                    
                    oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
                    uj = oorhoj * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
                    vj = oorhoj * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
                    wj = oorhoj * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                    }
                    Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
                    Zj = oorhoj * INTERPOL(Sol->rhoZ[jcp], Sol->rhoZ[jc], Sol->rhoZ[jcm]);
                    Hj = oorhoj * (INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]) + p0_c);
                    
                    oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                    ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                    vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                    wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                    }
                    Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
                    Zjm = oorhojm * INTERPOL(Sol->rhoZ[jcmm], Sol->rhoZ[jcm], Sol->rhoZ[jc]);
                    Hjm = oorhojm * (INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]) + p0_s);
                    
#ifdef STANDARD_STENCIL_PROJ1
                    tmpy = - dto2dy * ( hplusy[mc] * (dp2[jc] - dp2[jcm] ) );
#else
                    tmpy = - dto2dy * (  0.75  *   hplusy[mc] * (dp2[jc]   - dp2[jcm]  ) 
                                       + 0.125 * ( hplusy[me] * (dp2[jc+1] - dp2[jcm+1]) 
                                                 + hplusy[mw] * (dp2[jc-1] - dp2[jcm-1]) 
                                                  ) 
                                       );
#endif					
                    
                    tmpy += - 0.5 * dt * hgrav[mc] * 0.5 * (dp2[jc]+dp2[jcm]); /* implicit gravity contribution */                    
                    grhoY       = g->rhoY[mc] + tmpy;
                    
                    
                    
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(grhoY, 0.01) ); 
#endif
                    
                    us = upwind * ujm + (1.0 - upwind) * uj;
                    vs = upwind * vjm + (1.0 - upwind) * vj;
                    ws = upwind * wjm + (1.0 - upwind) * wj;
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xs[nsp]      = upwind * Xjm[nsp] + (1.0 - upwind) * Xj[nsp];
                    }
                    Ys = upwind * Yjm + (1.0 - upwind) * Yj;
                    Zs = upwind * Zjm + (1.0 - upwind) * Zj;
                    Hs = upwind * Hjm + (1.0 - upwind) * Hj;
                    
                    g->rho[mc]  = Ys * tmpy;
                    g->rhou[mc] = us * tmpy;  
                    g->rhov[mc] = vs * tmpy + vs * tmpy; 
                    g->rhow[mc] = ws * tmpy; 
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        g->rhoX[nsp][mc] = Xs[nsp] * tmpy; 
                    }
                    g->rhoe[mc] = 0.0 * Hs * tmpy;
                    g->rhoY[mc] = tmpy;
                    g->rhoZ[mc] = Zs * tmpy;
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
            
            const double dto2dx = implicitness * 0.5 * dt / elem->dx;
            const double dto2dy = implicitness * 0.5 * dt / elem->dy;
            const double dto2dz = implicitness * 0.5 * dt / elem->dz;
            
            const int dix = 1;
            const int diy = icx; 
            const int diz = icx*icy;            
            
            
#ifdef STANDARD_STENCIL_PROJ1
            const double cstencil  = 1.0;
#else
            const double cstencil  = 1.0/64.0;
#endif
            
            ConsVars* fx = flux[0];
            ConsVars* fy = flux[1];
            ConsVars* fz = flux[2];
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            const double* hplusz   = hplus[2];
            
            double oorhoi, ui, vi, wi, Yi, Zi, Hi; 
            double oorhoim, uim, vim, wim, Yim, Zim, Him;
            double Xi[NSPEC], Xim[NSPEC];
            
            double oorhoj, uj, vj, wj, Yj, Zj, Hj; 
            double oorhojm, ujm, vjm, wjm, Yjm, Zjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];
            
            double oorhok, uk, vk, wk, Yk, Zk, Hk; 
            double oorhokm, ukm, vkm, wkm, Ykm, Zkm, Hkm;
            double Xk[NSPEC], Xkm[NSPEC];
            
            double us, vs, ws, Hs, Ys, Zs;
            double Xs[NSPEC];
            
            double tmpx, tmpy, tmpz, frhoY;
            
            double p0_c, p0_s;
            double upwind;
            
            int i, j, k, l, m, n;
            
            int ic, icm, icmm, icp;
            int jc, jcm, jcmm, jcp;
            int kc, kcm, kcmm, kcp;
            
            for(k = igz; k < icz - igz; k++) {l = k * ifx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
                    p0_c = mpv->HydroState->p0[j];
                    for(i = igx; i < ifx - igx; i++) {n = m + i;
                        ic   = k*diz + j*diy + i*dix;
                        icp  = ic + dix; /* these appear unused as they only are accessed in macro INTERPOL() */
                        icm  = ic - dix;
                        icmm = ic - 2*dix; /* these appear unused as they only are accessed in macro INTERPOL() */
                        
                        /* compute interface values of u, v, Y, Z, H */
                        
                        oorhoi  = 1.0 / INTERPOL(Sol->rhoY[icp],  Sol->rhoY[ic],  Sol->rhoY[icm]); 
                        ui      = oorhoi * INTERPOL(Sol->rhou[icp], Sol->rhou[ic], Sol->rhou[icm]);
                        vi      = oorhoi * INTERPOL(Sol->rhov[icp], Sol->rhov[ic], Sol->rhov[icm]);
                        wi      = oorhoi * INTERPOL(Sol->rhow[icp], Sol->rhow[ic], Sol->rhow[icm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xi[nsp]      = oorhoi * INTERPOL(Sol->rhoX[nsp][icp], Sol->rhoX[nsp][ic], Sol->rhoX[nsp][icm]);
                        }
                        Yi      = oorhoi * INTERPOL(Sol->rho[icp], Sol->rho[ic], Sol->rho[icm]);
                        Zi      = oorhoi * INTERPOL(Sol->rhoZ[icp], Sol->rhoZ[ic], Sol->rhoZ[icm]);
                        Hi      = oorhoi * (INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]) + p0_c);
                        
                        oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                        uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                        vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                        wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                        }
                        Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
                        Zim     = oorhoim * INTERPOL(Sol->rhoZ[icmm], Sol->rhoZ[icm], Sol->rhoZ[ic]);
                        Him     = oorhoim * (INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]) + p0_c);
                        
#ifdef STANDARD_STENCIL_PROJ1
                        tmpx = - dto2dx * hplusx[n] * cstencil * (dp2[ic] - dp2[icm]);  
#else
                        tmpx = - dto2dx * hplusx[n] * cstencil *
                        (  36.0 *  (dp2[ic] - dp2[icm])  
                         +  6.0 * (  (dp2[ic+diy] - dp2[icm+diy]) 
                                   + (dp2[ic-diy] - dp2[icm-diy])  
                                   + (dp2[ic+diz] - dp2[icm+diz])  
                                   + (dp2[ic-diz] - dp2[icm-diz])
                                   )
                         +  1.0 * (  (dp2[ic+diy+diz] - dp2[icm+diy+diz]) 
                                   + (dp2[ic-diy+diz] - dp2[icm-diy+diz])  
                                   + (dp2[ic+diy-diz] - dp2[icm+diy-diz])  
                                   + (dp2[ic-diy-diz] - dp2[icm-diy-diz])
                                   )
                         );
#endif
                        
                        frhoY = fx->rhoY[n] + tmpx;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * uim + (1.0 - upwind) * ui;
                        vs = upwind * vim + (1.0 - upwind) * vi;  
                        ws = upwind * wim + (1.0 - upwind) * wi;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xim[nsp] + (1.0 - upwind) * Xi[nsp];
                        }
                        Ys = upwind * Yim + (1.0 - upwind) * Yi;
                        Zs = upwind * Zim + (1.0 - upwind) * Zi;
                        Hs = upwind * Him + (1.0 - upwind) * Hi;
                        
                        fx->rho[n]  = Ys * tmpx;
                        fx->rhou[n] = us * tmpx + us * tmpx;
                        fx->rhov[n] = vs * tmpx;  
                        fx->rhow[n] = ws * tmpx;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fx->rhoX[nsp][n] = Xs[nsp] * tmpx; 
                        }
                        fx->rhoe[n] = 0.0 * Hs * tmpx;
                        fx->rhoY[n] = tmpx;
                        fx->rhoZ[n] = Zs * tmpx;
                    }
                } 
            }
            
            
            /* fluxes in the y-direction */
            for(i = igx; i < icx - igx; i++) {l = i * ify*icz;
                for(k = igz; k < icz - igz; k++) {m = l + k * ify;
                    for(j = igy; j < ify - igy; j++) {n = m + j;
                        jc   = k*diz + j*diy + i*dix;
                        jcp  = jc + diy;  
                        jcm  = jc - diy;
                        jcmm = jc - 2*diy; 
                        
                        p0_c = mpv->HydroState->p0[j];
                        p0_s = mpv->HydroState->p0[j-1];
                        
                        buoy->y[jc] = 0.0;
                        
                        oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
                        uj = oorhoj  * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
                        vj = oorhoj  * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
                        wj = oorhoj  * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                        }
                        Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
                        Zj = oorhoj * INTERPOL(Sol->rhoZ[jcp], Sol->rhoZ[jc], Sol->rhoZ[jcm]);
                        Hj = oorhoj * (INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]) + p0_c);
                        
                        oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                        ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                        vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                        wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                        }
                        Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
                        Zjm = oorhojm * INTERPOL(Sol->rhoZ[jcmm], Sol->rhoZ[jcm], Sol->rhoZ[jc]);
                        Hjm = oorhojm * (INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]) + p0_s);
                        
#ifdef STANDARD_STENCIL_PROJ1
                        tmpy = - dto2dy * hplusy[n] * cstencil * (dp2[jc] - dp2[jcm]);  
#else
                        tmpy = - dto2dy * hplusy[n] * cstencil *
                        (  36.0 *    (dp2[jc] - dp2[jcm])  
                         +  6.0 * (  (dp2[jc+diz] - dp2[jcm+diz]) 
                                   + (dp2[jc-diz] - dp2[jcm-diz])  
                                   + (dp2[jc+dix] - dp2[jcm+dix])  
                                   + (dp2[jc-dix] - dp2[jcm-dix])
                                   )
                         +  1.0 * (  (dp2[jc+diz+dix] - dp2[jcm+diz+dix]) 
                                   + (dp2[jc-diz+dix] - dp2[jcm-diz+dix])  
                                   + (dp2[jc+diz-dix] - dp2[jcm+diz-dix])  
                                   + (dp2[jc-diz-dix] - dp2[jcm-diz-dix])
                                   )
                         );
#endif                                                
                        frhoY       = fy->rhoY[n] + tmpy;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * ujm + (1.0 - upwind) * uj;
                        vs = upwind * vjm + (1.0 - upwind) * vj;
                        ws = upwind * wjm + (1.0 - upwind) * wj;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xjm[nsp] + (1.0 - upwind) * Xj[nsp];
                        }
                        Ys = upwind * Yjm + (1.0 - upwind) * Yj;
                        Zs = upwind * Zjm + (1.0 - upwind) * Zj;
                        Hs = upwind * Hjm + (1.0 - upwind) * Hj;
                        
                        fy->rho[n]  = Ys * tmpy;
                        fy->rhou[n] = us * tmpy;  
                        fy->rhov[n] = vs * tmpy + vs * tmpy; 
                        fy->rhow[n] = ws * tmpy;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fy->rhoX[nsp][n] = Xs[nsp] * tmpy; 
                        }
                        fy->rhoe[n] = 0.0 * Hs * tmpy;
                        fy->rhoY[n] = tmpy;
                        fy->rhoZ[n] = Zs * tmpy;
                    }   
                }
            }
            
            /* fluxes in the z-direction */
            for(j = igy; j < icy - igy; j++) {l = j * ifz*icx;
                for(i = igx; i < icx - igx; i++) {m = l + i * ifz;
                    for(k = igz; k < ifz - igz; k++) {n = m + k;
                        kc   = k*diz + j*diy + i*dix;
                        kcp  = kc + diz;   
                        kcm  = kc - diz;
                        kcmm = kc - 2*diz;  
                        
                        p0_c = mpv->HydroState->p0[j];
                        
                        oorhok = 1.0 /INTERPOL(Sol->rhoY[kcp], Sol->rhoY[kc], Sol->rhoY[kcm]);
                        uk = oorhok * INTERPOL(Sol->rhou[kcp], Sol->rhou[kc], Sol->rhou[kcm]);
                        vk = oorhok * INTERPOL(Sol->rhov[kcp], Sol->rhov[kc], Sol->rhov[kcm]);
                        wk = oorhok * INTERPOL(Sol->rhow[kcp], Sol->rhow[kc], Sol->rhow[kcm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xk[nsp]      = oorhok * INTERPOL(Sol->rhoX[nsp][kcp], Sol->rhoX[nsp][kc], Sol->rhoX[nsp][kcm]);
                        }
                        Yk = oorhok * INTERPOL(Sol->rho[kcp], Sol->rho[kc], Sol->rho[kcm]);
                        Zk = oorhok * INTERPOL(Sol->rhoZ[kcp], Sol->rhoZ[kc], Sol->rhoZ[kcm]);
                        Hk = oorhok * (INTERPOL(Sol->rhoe[kcp], Sol->rhoe[kc], Sol->rhoe[kcm]) + p0_c);
                        
                        oorhokm = 1.0 / INTERPOL(Sol->rhoY[kcmm], Sol->rhoY[kcm], Sol->rhoY[kc]);
                        ukm = oorhokm * INTERPOL(Sol->rhou[kcmm], Sol->rhou[kcm], Sol->rhou[kc]);
                        vkm = oorhokm * INTERPOL(Sol->rhov[kcmm], Sol->rhov[kcm], Sol->rhov[kc]);
                        wkm = oorhokm * INTERPOL(Sol->rhow[kcmm], Sol->rhow[kcm], Sol->rhow[kc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xkm[nsp]      = oorhokm * INTERPOL(Sol->rhoX[nsp][kcmm], Sol->rhoX[nsp][kcm], Sol->rhoX[nsp][kc]);
                        }
                        Ykm = oorhokm * INTERPOL(Sol->rho[kcmm], Sol->rho[kcm], Sol->rho[kc]);
                        Zkm = oorhokm * INTERPOL(Sol->rhoZ[kcmm], Sol->rhoZ[kcm], Sol->rhoZ[kc]);
                        Hkm = oorhokm * (INTERPOL(Sol->rhoe[kcmm], Sol->rhoe[kcm], Sol->rhoe[kc]) + p0_c);
                        
#ifdef STANDARD_STENCIL_PROJ1
                        tmpz = - dto2dz * hplusz[n] * cstencil * (dp2[kc] - dp2[kcm]);  
#else
                        tmpz = - dto2dz * hplusz[n] * cstencil *
                        (  36.0 *    (dp2[kc] - dp2[kcm])  
                         +  6.0 * (  (dp2[kc+diz] - dp2[kcm+diz]) 
                                   + (dp2[kc-diz] - dp2[kcm-diz])  
                                   + (dp2[kc+dix] - dp2[kcm+dix])  
                                   + (dp2[kc-dix] - dp2[kcm-dix])
                                   )
                         +  1.0 * (  (dp2[kc+diz+dix] - dp2[kcm+diz+dix]) 
                                   + (dp2[kc-diz+dix] - dp2[kcm-diz+dix])  
                                   + (dp2[kc+diz-dix] - dp2[kcm+diz-dix])  
                                   + (dp2[kc-diz-dix] - dp2[kcm-diz-dix])
                                   )
                         );
#endif                        
                        frhoY       = fz->rhoY[n] + tmpz;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * ukm + (1.0 - upwind) * uk;
                        vs = upwind * vkm + (1.0 - upwind) * vk;
                        ws = upwind * wkm + (1.0 - upwind) * wk;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xkm[nsp] + (1.0 - upwind) * Xk[nsp];
                        }
                        Ys = upwind * Ykm + (1.0 - upwind) * Yk;
                        Zs = upwind * Zkm + (1.0 - upwind) * Zk;
                        Hs = upwind * Hkm + (1.0 - upwind) * Hk;
                        
                        fz->rho[n]  = Ys * tmpz;
                        fz->rhou[n] = us * tmpz;  
                        fz->rhov[n] = vs * tmpz; 
                        fz->rhow[n] = ws * tmpz + ws * tmpz;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fz->rhoX[nsp][n] = Xs[nsp] * tmpz; 
                        }
                        fz->rhoe[n] = 0.0 * Hs * tmpz;
                        fz->rhoY[n] = tmpz;
                        fz->rhoZ[n] = Zs * tmpz;
                    }   
                }
            }			
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
}

#else /* SHIFTED_COEFFICIENTS_PROJ1 */
static void flux_correction_due_to_pressure_gradients(
													  ConsVars* flux[3],
													  VectorField* buoy,
													  const ElemSpaceDiscr* elem,
													  ConsVars* Sol,
													  ConsVars* Sol0,
													  const MPV* mpv, 
													  double* hplus[3],
													  double* hgrav,
													  double* dp2,
													  const double t,
													  const double dt,
													  const double implicitness) {
	
	extern User_Data ud;
	
	const int ndim = elem->ndim;
        
    int nsp;
	
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
			
			const double dto2dx = implicitness * 0.5 * dt / elem->dx;
			const double dto2dy = implicitness * 0.5 * dt / elem->dy;
            
			ConsVars* f = flux[0];
			ConsVars* g = flux[1];
			
			const double* hplusx   = hplus[0];
			const double* hplusy   = hplus[1];
			
			double oorhoi, ui, vi, wi, Yi, Zi, Hi, oorhoim, uim, vim, wim, Yim, Zim, Him;
            double Xi[NSPEC], Xim[NSPEC];
			double oorhoj, uj, vj, wj, Yj, Zj, Hj, oorhojm, ujm, vjm, wjm, Yjm, Zjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];
			double us, vs, ws, Hs, Ys, Zs;
            double Xs[NSPEC];
			double tmpx, tmpy, frhoY, grhoY;
			
			double p0_c, p0_s;
			double upwind;
            
            int i, j, m, n, ic, icm, icmm, icp, jc, jcm, jcmm, jcp;
			
			for(j = igy; j < icy - igy; j++) {
                m    = j * ifx;
				p0_c = mpv->HydroState->p0[j];
								
				for(i = igx; i < ifx - igx; i++) {
                    n    = m + i;
					ic   = j*icx + i;
					icp  = ic + 1;
					icm  = ic - 1;
					icmm = ic - 2;
                    
                    buoy->x[ic] = 0.0;
										
					/* compute interface values of u, v, Y, Z, H */
					oorhoi  = 1.0 / INTERPOL(Sol->rhoY[icp],  Sol->rhoY[ic],  Sol->rhoY[icm]); 
					ui      = oorhoi * INTERPOL(Sol->rhou[icp], Sol->rhou[ic], Sol->rhou[icm]);
					vi      = oorhoi * INTERPOL(Sol->rhov[icp], Sol->rhov[ic], Sol->rhov[icm]);
					wi      = oorhoi * INTERPOL(Sol->rhow[icp], Sol->rhow[ic], Sol->rhow[icm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xi[nsp]      = oorhoi * INTERPOL(Sol->rhoX[nsp][icp], Sol->rhoX[nsp][ic], Sol->rhoX[nsp][icm]);
                    }
					Yi      = oorhoi * INTERPOL(Sol->rho[icp], Sol->rho[ic], Sol->rho[icm]);
					Zi      = oorhoi * INTERPOL(Sol->rhoZ[icp], Sol->rhoZ[ic], Sol->rhoZ[icm]);
					Hi      = oorhoi * (INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]) + p0_c);
					
					oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
					uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
					vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
					wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                    }
					Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
					Zim     = oorhoim * INTERPOL(Sol->rhoZ[icmm], Sol->rhoZ[icm], Sol->rhoZ[ic]);
					Him     = oorhoim * (INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]) + p0_c);
						
#ifdef STANDARD_STENCIL_PROJ1
                    tmpx = - dto2dx * ( hplusx[n] * (dp2[ic] - dp2[icm]) );
#else
                    tmpx = - dto2dx * (  hplusx[n] * (dp2[ic] - dp2[icm] )  
                                       + 0.125 * ( hplusx[n] * ((dp2[ic+icx] - dp2[icm+icx] ) - (dp2[ic] - dp2[icm]))  
                                                  - hplusx[n] * ((dp2[ic] - dp2[icm] ) - (dp2[ic-icx] - dp2[icm-icx]))  
                                                  ) 
                                       );
#endif
                    
					frhoY = f->rhoY[n] + tmpx;
                    
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
					
					us = upwind * uim + (1.0 - upwind) * ui;
					vs = upwind * vim + (1.0 - upwind) * vi;  
					ws = upwind * wim + (1.0 - upwind) * wi;  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xs[nsp]      = upwind * Xim[nsp] + (1.0 - upwind) * Xi[nsp];
                    }
					Ys = upwind * Yim + (1.0 - upwind) * Yi;
					Zs = upwind * Zim + (1.0 - upwind) * Zi;
					Hs = upwind * Him + (1.0 - upwind) * Hi;
					
					f->rho[n]  = Ys * tmpx;
					f->rhou[n] = us * tmpx + us * tmpx;
					f->rhov[n] = vs * tmpx;  
					f->rhow[n] = ws * tmpx;  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        f->rhoX[nsp][n] = Xs[nsp] * tmpx; 
                    }
					f->rhoe[n] = 0.0 * Hs * tmpx;
					f->rhoY[n] = tmpx;
					f->rhoZ[n] = Zs * tmpx;
				}
			}  
			
			
			/* fluxes in the y-direction */
			for(i = igx; i < icx - igx; i++) {n = i * ify;
				for(j = igy; j < ify - igy; j++) {m = n + j;
					jc   = j * icx + i;
					jcp  = jc + icx;
					jcm  = jc - icx;
					jcmm = jc - 2*icx;
										
					p0_c = mpv->HydroState->p0[j];
					p0_s = mpv->HydroState->p0[j-1];
					
                    buoy->y[jc] = 0.0;
                    
					oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
					uj = oorhoj * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
					vj = oorhoj * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
					wj = oorhoj * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                    }
					Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
					Zj = oorhoj * INTERPOL(Sol->rhoZ[jcp], Sol->rhoZ[jc], Sol->rhoZ[jcm]);
					Hj = oorhoj * (INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]) + p0_c);
					
					oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
					ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
					vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
					wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                    }
					Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
					Zjm = oorhojm * INTERPOL(Sol->rhoZ[jcmm], Sol->rhoZ[jcm], Sol->rhoZ[jc]);
					Hjm = oorhojm * (INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]) + p0_s);
					
#ifdef STANDARD_STENCIL_PROJ1
					tmpy = - dto2dy * ( hplusy[m] * (dp2[jc] - dp2[jcm] ) );
#else
					tmpy = - dto2dy * (   hplusy[m] * (dp2[jc]   - dp2[jcm] ) 
									   + 0.125 * ( hplusy[m] * ((dp2[jc+1] - dp2[jcm+1] ) - (dp2[jc  ] - dp2[jcm  ] )) 
                                                 - hplusy[m] * ((dp2[jc  ] - dp2[jcm  ] ) - (dp2[jc-1] - dp2[jcm-1] )) 
												  ) 
									   );
#endif					
                    
                    tmpy += - 0.5 * dt * hgrav[m] * 0.5 * (dp2[jc]+dp2[jcm]); /* implicit gravity contribution */                    
					grhoY       = g->rhoY[m] + tmpy;
                    
                    
                    
#ifdef NO_UPWIND_PROJ1
					upwind = 0.5; 
#else
					upwind = 0.5 * ( 1.0 + SMOOTHSIGN(grhoY, 0.01) ); 
#endif
					
					us = upwind * ujm + (1.0 - upwind) * uj;
					vs = upwind * vjm + (1.0 - upwind) * vj;
					ws = upwind * wjm + (1.0 - upwind) * wj;
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xs[nsp]      = upwind * Xjm[nsp] + (1.0 - upwind) * Xj[nsp];
                    }
					Ys = upwind * Yjm + (1.0 - upwind) * Yj;
					Zs = upwind * Zjm + (1.0 - upwind) * Zj;
					Hs = upwind * Hjm + (1.0 - upwind) * Hj;
					
					g->rho[m]  = Ys * tmpy;
					g->rhou[m] = us * tmpy;  
					g->rhov[m] = vs * tmpy + vs * tmpy; 
					g->rhow[m] = ws * tmpy; 
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        g->rhoX[nsp][m] = Xs[nsp] * tmpy; 
                    }
					g->rhoe[m] = 0.0 * Hs * tmpy;
					g->rhoY[m] = tmpy;
					g->rhoZ[m] = Zs * tmpy;
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
			
			const double dto2dx = implicitness * 0.5 * dt / elem->dx;
			const double dto2dy = implicitness * 0.5 * dt / elem->dy;
			const double dto2dz = implicitness * 0.5 * dt / elem->dz;

            const int dix = 1;
            const int diy = icx; 
            const int diz = icx*icy;            
            
			
#ifdef STANDARD_STENCIL_PROJ1
            const double cstencil  = 1.0;
#else
            const double cstencil  = 1.0/64.0;
#endif

			ConsVars* fx = flux[0];
			ConsVars* fy = flux[1];
			ConsVars* fz = flux[2];
			
			const double* hplusx   = hplus[0];
			const double* hplusy   = hplus[1];
			const double* hplusz   = hplus[2];
			
			double oorhoi, ui, vi, wi, Yi, Zi, Hi; 
            double oorhoim, uim, vim, wim, Yim, Zim, Him;
            double Xi[NSPEC], Xim[NSPEC];

            double oorhoj, uj, vj, wj, Yj, Zj, Hj; 
            double oorhojm, ujm, vjm, wjm, Yjm, Zjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];

            double oorhok, uk, vk, wk, Yk, Zk, Hk; 
            double oorhokm, ukm, vkm, wkm, Ykm, Zkm, Hkm;
            double Xk[NSPEC], Xkm[NSPEC];

            double us, vs, ws, Hs, Ys, Zs;
            double Xs[NSPEC];
            
			double tmpx, tmpy, tmpz, frhoY;
			
			double p0_c, p0_s;
			double upwind;
			
            int i, j, k, l, m, n;
            
            int ic, icm, icmm, icp;
            int jc, jcm, jcmm, jcp;
            int kc, kcm, kcmm, kcp;

			for(k = igz; k < icz - igz; k++) {l = k * ifx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
                    p0_c = mpv->HydroState->p0[j];
                    for(i = igx; i < ifx - igx; i++) {n = m + i;
                        ic   = k*diz + j*diy + i*dix;
                        icp  = ic + dix; /* these appear unused as they only are accessed in macro INTERPOL() */
                        icm  = ic - dix;
                        icmm = ic - 2*dix; /* these appear unused as they only are accessed in macro INTERPOL() */
                                                
                        /* compute interface values of u, v, Y, Z, H */
                        
                        oorhoi  = 1.0 / INTERPOL(Sol->rhoY[icp],  Sol->rhoY[ic],  Sol->rhoY[icm]); 
                        ui      = oorhoi * INTERPOL(Sol->rhou[icp], Sol->rhou[ic], Sol->rhou[icm]);
                        vi      = oorhoi * INTERPOL(Sol->rhov[icp], Sol->rhov[ic], Sol->rhov[icm]);
                        wi      = oorhoi * INTERPOL(Sol->rhow[icp], Sol->rhow[ic], Sol->rhow[icm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xi[nsp]      = oorhoi * INTERPOL(Sol->rhoX[nsp][icp], Sol->rhoX[nsp][ic], Sol->rhoX[nsp][icm]);
                        }
                        Yi      = oorhoi * INTERPOL(Sol->rho[icp], Sol->rho[ic], Sol->rho[icm]);
                        Zi      = oorhoi * INTERPOL(Sol->rhoZ[icp], Sol->rhoZ[ic], Sol->rhoZ[icm]);
                        Hi      = oorhoi * (INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]) + p0_c);
                        
                        oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                        uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                        vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                        wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                        }
                        Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
                        Zim     = oorhoim * INTERPOL(Sol->rhoZ[icmm], Sol->rhoZ[icm], Sol->rhoZ[ic]);
                        Him     = oorhoim * (INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]) + p0_c);
                                        
#ifdef STANDARD_STENCIL_PROJ1
                        tmpx = - dto2dx * hplusx[n] * cstencil * (dp2[ic] - dp2[icm]);  
#else
                        tmpx = - dto2dx * hplusx[n] * cstencil *
                        (  36.0 *  (dp2[ic] - dp2[icm])  
                         +  6.0 * (  (dp2[ic+diy] - dp2[icm+diy]) 
                                   + (dp2[ic-diy] - dp2[icm-diy])  
                                   + (dp2[ic+diz] - dp2[icm+diz])  
                                   + (dp2[ic-diz] - dp2[icm-diz])
                                   )
                         +  1.0 * (  (dp2[ic+diy+diz] - dp2[icm+diy+diz]) 
                                   + (dp2[ic-diy+diz] - dp2[icm-diy+diz])  
                                   + (dp2[ic+diy-diz] - dp2[icm+diy-diz])  
                                   + (dp2[ic-diy-diz] - dp2[icm-diy-diz])
                                   )
                         );
#endif
                                                
                        frhoY = fx->rhoY[n] + tmpx;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * uim + (1.0 - upwind) * ui;
                        vs = upwind * vim + (1.0 - upwind) * vi;  
                        ws = upwind * wim + (1.0 - upwind) * wi;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xim[nsp] + (1.0 - upwind) * Xi[nsp];
                        }
                        Ys = upwind * Yim + (1.0 - upwind) * Yi;
                        Zs = upwind * Zim + (1.0 - upwind) * Zi;
                        Hs = upwind * Him + (1.0 - upwind) * Hi;
                        
                        fx->rho[n]  = Ys * tmpx;
                        fx->rhou[n] = us * tmpx + us * tmpx;
                        fx->rhov[n] = vs * tmpx;  
                        fx->rhow[n] = ws * tmpx;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fx->rhoX[nsp][n] = Xs[nsp] * tmpx; 
                        }
                        fx->rhoe[n] = 0.0 * Hs * tmpx;
                        fx->rhoY[n] = tmpx;
                        fx->rhoZ[n] = Zs * tmpx;
                    }
                } 
            }
			
			
			/* fluxes in the y-direction */
            for(i = igx; i < icx - igx; i++) {l = i * ify*icz;
                for(k = igz; k < icz - igz; k++) {m = l + k * ify;
                    for(j = igy; j < ify - igy; j++) {n = m + j;
                        jc   = k*diz + j*diy + i*dix;
                        jcp  = jc + diy;  
                        jcm  = jc - diy;
                        jcmm = jc - 2*diy; 
                                                
                        p0_c = mpv->HydroState->p0[j];
                        p0_s = mpv->HydroState->p0[j-1];
                        
                        buoy->y[jc] = 0.0;

                        oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
                        uj = oorhoj  * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
                        vj = oorhoj  * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
                        wj = oorhoj  * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                        }
                        Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
                        Zj = oorhoj * INTERPOL(Sol->rhoZ[jcp], Sol->rhoZ[jc], Sol->rhoZ[jcm]);
                        Hj = oorhoj * (INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]) + p0_c);
                        
                        oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                        ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                        vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                        wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                        }
                        Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
                        Zjm = oorhojm * INTERPOL(Sol->rhoZ[jcmm], Sol->rhoZ[jcm], Sol->rhoZ[jc]);
                        Hjm = oorhojm * (INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]) + p0_s);
                       
#ifdef STANDARD_STENCIL_PROJ1
                        tmpy = - dto2dy * hplusy[n] * cstencil * (dp2[jc] - dp2[jcm]);  
#else
                        tmpy = - dto2dy * hplusy[n] * cstencil *
                        (  36.0 *    (dp2[jc] - dp2[jcm])  
                         +  6.0 * (  (dp2[jc+diz] - dp2[jcm+diz]) 
                                   + (dp2[jc-diz] - dp2[jcm-diz])  
                                   + (dp2[jc+dix] - dp2[jcm+dix])  
                                   + (dp2[jc-dix] - dp2[jcm-dix])
                                   )
                         +  1.0 * (  (dp2[jc+diz+dix] - dp2[jcm+diz+dix]) 
                                   + (dp2[jc-diz+dix] - dp2[jcm-diz+dix])  
                                   + (dp2[jc+diz-dix] - dp2[jcm+diz-dix])  
                                   + (dp2[jc-diz-dix] - dp2[jcm-diz-dix])
                                   )
                         );
#endif                                                
                        frhoY       = fy->rhoY[n] + tmpy;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * ujm + (1.0 - upwind) * uj;
                        vs = upwind * vjm + (1.0 - upwind) * vj;
                        ws = upwind * wjm + (1.0 - upwind) * wj;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xjm[nsp] + (1.0 - upwind) * Xj[nsp];
                        }
                        Ys = upwind * Yjm + (1.0 - upwind) * Yj;
                        Zs = upwind * Zjm + (1.0 - upwind) * Zj;
                        Hs = upwind * Hjm + (1.0 - upwind) * Hj;
                        
                        fy->rho[n]  = Ys * tmpy;
                        fy->rhou[n] = us * tmpy;  
                        fy->rhov[n] = vs * tmpy + vs * tmpy; 
                        fy->rhow[n] = ws * tmpy;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fy->rhoX[nsp][n] = Xs[nsp] * tmpy; 
                        }
                        fy->rhoe[n] = 0.0 * Hs * tmpy;
                        fy->rhoY[n] = tmpy;
                        fy->rhoZ[n] = Zs * tmpy;
                    }   
                }
			}
             
			/* fluxes in the z-direction */
			for(j = igy; j < icy - igy; j++) {l = j * ifz*icx;
                for(i = igx; i < icx - igx; i++) {m = l + i * ifz;
                    for(k = igz; k < ifz - igz; k++) {n = m + k;
                        kc   = k*diz + j*diy + i*dix;
                        kcp  = kc + diz;   
                        kcm  = kc - diz;
                        kcmm = kc - 2*diz;  
                        
                        p0_c = mpv->HydroState->p0[j];
                       
                        oorhok = 1.0 /INTERPOL(Sol->rhoY[kcp], Sol->rhoY[kc], Sol->rhoY[kcm]);
                        uk = oorhok * INTERPOL(Sol->rhou[kcp], Sol->rhou[kc], Sol->rhou[kcm]);
                        vk = oorhok * INTERPOL(Sol->rhov[kcp], Sol->rhov[kc], Sol->rhov[kcm]);
                        wk = oorhok * INTERPOL(Sol->rhow[kcp], Sol->rhow[kc], Sol->rhow[kcm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xk[nsp]      = oorhok * INTERPOL(Sol->rhoX[nsp][kcp], Sol->rhoX[nsp][kc], Sol->rhoX[nsp][kcm]);
                        }
                        Yk = oorhok * INTERPOL(Sol->rho[kcp], Sol->rho[kc], Sol->rho[kcm]);
                        Zk = oorhok * INTERPOL(Sol->rhoZ[kcp], Sol->rhoZ[kc], Sol->rhoZ[kcm]);
                        Hk = oorhok * (INTERPOL(Sol->rhoe[kcp], Sol->rhoe[kc], Sol->rhoe[kcm]) + p0_c);
                        
                        oorhokm = 1.0 / INTERPOL(Sol->rhoY[kcmm], Sol->rhoY[kcm], Sol->rhoY[kc]);
                        ukm = oorhokm * INTERPOL(Sol->rhou[kcmm], Sol->rhou[kcm], Sol->rhou[kc]);
                        vkm = oorhokm * INTERPOL(Sol->rhov[kcmm], Sol->rhov[kcm], Sol->rhov[kc]);
                        wkm = oorhokm * INTERPOL(Sol->rhow[kcmm], Sol->rhow[kcm], Sol->rhow[kc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xkm[nsp]      = oorhokm * INTERPOL(Sol->rhoX[nsp][kcmm], Sol->rhoX[nsp][kcm], Sol->rhoX[nsp][kc]);
                        }
                        Ykm = oorhokm * INTERPOL(Sol->rho[kcmm], Sol->rho[kcm], Sol->rho[kc]);
                        Zkm = oorhokm * INTERPOL(Sol->rhoZ[kcmm], Sol->rhoZ[kcm], Sol->rhoZ[kc]);
                        Hkm = oorhokm * (INTERPOL(Sol->rhoe[kcmm], Sol->rhoe[kcm], Sol->rhoe[kc]) + p0_c);
                        
#ifdef STANDARD_STENCIL_PROJ1
                        tmpz = - dto2dz * hplusz[n] * cstencil * (dp2[kc] - dp2[kcm]);  
#else
                        tmpz = - dto2dz * hplusz[n] * cstencil *
                        (  36.0 *    (dp2[kc] - dp2[kcm])  
                         +  6.0 * (  (dp2[kc+diz] - dp2[kcm+diz]) 
                                   + (dp2[kc-diz] - dp2[kcm-diz])  
                                   + (dp2[kc+dix] - dp2[kcm+dix])  
                                   + (dp2[kc-dix] - dp2[kcm-dix])
                                   )
                         +  1.0 * (  (dp2[kc+diz+dix] - dp2[kcm+diz+dix]) 
                                   + (dp2[kc-diz+dix] - dp2[kcm-diz+dix])  
                                   + (dp2[kc+diz-dix] - dp2[kcm+diz-dix])  
                                   + (dp2[kc-diz-dix] - dp2[kcm-diz-dix])
                                   )
                         );
#endif                        
                        frhoY       = fz->rhoY[n] + tmpz;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * ukm + (1.0 - upwind) * uk;
                        vs = upwind * vkm + (1.0 - upwind) * vk;
                        ws = upwind * wkm + (1.0 - upwind) * wk;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xkm[nsp] + (1.0 - upwind) * Xk[nsp];
                        }
                        Ys = upwind * Ykm + (1.0 - upwind) * Yk;
                        Zs = upwind * Zkm + (1.0 - upwind) * Zk;
                        Hs = upwind * Hkm + (1.0 - upwind) * Hk;
                        
                        fz->rho[n]  = Ys * tmpz;
                        fz->rhou[n] = us * tmpz;  
                        fz->rhov[n] = vs * tmpz; 
                        fz->rhow[n] = ws * tmpz + ws * tmpz;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fz->rhoX[nsp][n] = Xs[nsp] * tmpz; 
                        }
                        fz->rhoe[n] = 0.0 * Hs * tmpz;
                        fz->rhoY[n] = tmpz;
                        fz->rhoZ[n] = Zs * tmpz;
                    }   
                }
			}			
            break;
		}
		default: ERROR("ndim not in {1, 2, 3}");
    }
}
#endif /* SHIFTED_COEFFICIENTS_PROJ1 */

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
	const double implicitness = 0.5*ud.implicitness;
	
	double** hplus       = mpv->Level[0]->wplus;
	
	int ndim = elem->ndim;
	
	if (elem->left == OPEN) {
		
		assert(elem->left == elem->right); /* if you get thrown out here: one-sided open domain not yet impltd */
		
		switch(ndim) {
			case 1: {
				ERROR("function not available");
				break;
			}
				
			case 2: {
				
				const double factor = 1.0 / (implicitness * elem->dx * dt);
				
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
				ERROR("ndim = 3 - option not implemented (kinetic_energy_change_explicit() )");
			}
			default: ERROR("ndim not in {1, 2, 3}");
		}
	}
}

/* ========================================================================== */

static void flux_correction_due_to_pressure_values(
												   ConsVars* flux[3],
												   VectorField* buoy,
												   const ElemSpaceDiscr* elem,
												   ConsVars* Sol,
												   double* dp2, 
												   const double dt) {
	
	extern User_Data ud;
	
	const int ndim = elem->ndim;
	    
	switch(ndim) {
		case 1: {
			ERROR("function not available");
			break;
		}
			
		case 2: {
			int i, j, m, n, ic, icm, jc, jcm;
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int ifx = elem->ifx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int ify = elem->ify;
			
			ConsVars* f = flux[0];
			ConsVars* g = flux[1];
			
            /* flux_post_correction    change from   flux ... +=  corr    to    flux ... = corr */
			
			/* Momentum flux corrections */
			for(j = igy; j < icy - igy; j++) {
				
				m = j * ifx;
				
				for(i = igx; i < ifx - igx; i++) {
					
					n   = m + i;
					ic  = n - j;
					icm = ic - 1;
					
					f->rhou[n] += ud.p_extrapol*0.5*(
                                              (ud.latw[0]*dp2[ic+icx]  + ud.latw[1]*dp2[ic]  + ud.latw[2]*dp2[ic-icx]) 
                                            + (ud.latw[0]*dp2[icm+icx] + ud.latw[1]*dp2[icm] + ud.latw[2]*dp2[icm-icx]) 
                                               );
				}
			}  
			
			
            /* flux_post_correction    change from   flux ... +=  corr    to    flux ... = corr */
			
			/* fluxes in the y-direction */
			for(i = igx; i < icx - igx; i++) {
				
				n = i * ify;
				
				for(j = igy; j < ify - igy; j++) {
					
					m   = n + j;
					jc  = j * icx + i;
					jcm = jc - icx;
					
					g->rhov[m] += ud.p_extrapol*0.5*(  (ud.latw[0]*dp2[jc+1]  + ud.latw[1]*dp2[jc]  + ud.latw[2]*dp2[jc-1])  
                                                     + (ud.latw[0]*dp2[jcm+1] + ud.latw[1]*dp2[jcm] + ud.latw[2]*dp2[jcm-1])
                                                  );
				}   
			}
			
			break;
		}
		case 3: {
			ERROR("function not available");
			break;
		}
		default: ERROR("ndim not in {1, 2,3}");
    }
#if 0
    printf("leaving flux_correction_due_to_pressure_values()\n");
#endif
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
				ERROR("ndim = 3 - option not implemented (kinetic_energy_change_explicit() )");
			}
			default: ERROR("ndim not in {1, 2, 3}");
		}
	}
}

#endif /* SOLVER_1_HYPRE */


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
