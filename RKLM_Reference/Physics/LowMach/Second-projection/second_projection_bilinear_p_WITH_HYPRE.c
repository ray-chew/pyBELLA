
/*******************************************************************************
 File:   second_projection.c
 Author: Nicola
 Date:   Fri Mar 13 14:58:07 WET 1998
 *******************************************************************************/
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <tgmath.h>
#include "Common.h"
#include "BiCGSTAB.h"
#include "set_ghostnodes_p.h"
#include "variable.h"
#include "mpv.h"
#include "error.h"
#include "warning.h"
#include "userdata.h"
#include "variable_coefficient_poisson_nodes.h"
#include "second_projection_bilinear_p.h" 
#include "io.h"
#include "userdata.h"
#include "enumerator.h"
#include "Eos.h"
#include "thermodynamic.h"
#include "boundary.h"
#include "limiter.h"
#include "math_own.h"
#include "memory.h"
#include "laplacian_nodes.h"
#include "numerical_flux.h"
#include "flux_correction.h"
#include "set_ghostcells_p.h"


#ifdef SOLVER_2_HYPRE

#ifdef SOLVER_2_HYPRE_RUPE

static double divergence_nodes(
                               double* rhs,
                               const ElemSpaceDiscr* elem,
                               const NodeSpaceDiscr* node,
                               const ConsVars* Sol,
                               double* eta,
                               const MPV* mpv,
                               const BDRY* bdry,
                               const double dt,
                               const double weight);


static enum Constraint integral_condition_nodes(
                                                double* rhs,
                                                const NodeSpaceDiscr* node,
                                                const int x_periodic,
                                                const int y_periodic,
                                                const int z_periodic);

static void	catch_periodic_directions(
                                      double* rhs,
                                      const NodeSpaceDiscr* node,
                                      const ElemSpaceDiscr* elem,
                                      const int x_periodic,
                                      const int y_periodic,
                                      const int z_periodic);

static void operator_coefficients_nodes(
                                        double* Splus[3],
                                        double* wcenter,
                                        double* wgrav,
                                        const ElemSpaceDiscr* elem,
                                        const NodeSpaceDiscr* node,
                                        ConsVars* Sol,
                                        const ConsVars* Sol0,
                                        const MPV* mpv,
                                        const double dt);

void correction_nodes(
                      ConsVars* Sol,
                      const ElemSpaceDiscr* elem,
                      const NodeSpaceDiscr* node,
                      const double* hplus[3], 
                      const double* hgrav,
                      double* p2,
                      const double t,
                      const double dt);


/* ========================================================================== */

void second_projection(
                       ConsVars* Sol,
                       MPV* mpv,
                       const ConsVars* Sol0,
                       const ElemSpaceDiscr* elem,
                       const NodeSpaceDiscr* node,
                       const double p_update,
                       const double t,
                       const double dt) {
    
    extern User_Data ud;
    extern BDRY* bdry;
    
    const int nc = node->nc;
    const double dx = node->dx;
    const double dy = node->dy;
    const double dz = node->dz;
    const double dV = dx*dy*dz;

    
    double** hplus   = mpv->Level[0]->wplus;
    double*  hcenter = mpv->Level[0]->wcenter;
    double*  hgrav   = mpv->Level[0]->wgrav;
    
    double* rhs      = mpv->Level[0]->rhs;
    double* p2       = mpv->Level[0]->p;
    
    double rhs_max;
    
    int x_periodic, y_periodic, z_periodic;
    int ii;
    
    printf("\n\n====================================================");
    printf("\nSecond Projection");
    printf("\n====================================================\n");

    /* switch for stopping loops in solver for periodic data before
     last grid line in each direction */
    x_periodic = 0;
    y_periodic = 0;
    z_periodic = 0;
    if(ud.bdrytype_min[0] == PERIODIC) x_periodic = 1;
    if(ud.bdrytype_min[1] == PERIODIC) y_periodic = 1;
    if(ud.bdrytype_min[2] == PERIODIC) z_periodic = 1;
    
    /* KEEP_OLD_POISSON_SOLUTIONS */
    for(ii=0; ii<nc; ii++){
        p2[ii] = mpv->dp2_nodes[ii];
        rhs[ii] = 0.0;
    }
    
    double rhs_weight_new = 2.0;
    double rhs_weight_old = 2.0*ud.compressibility;
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv->eta, mpv, bdry, dt, rhs_weight_new);
    if (ud.compressibility) {
        divergence_nodes(rhs, elem, node, Sol0, mpv->eta0, mpv, bdry, dt, rhs_weight_old); /*  for psinc, this can be commented out */
    }
    
    catch_periodic_directions(rhs, node, elem, x_periodic, y_periodic, z_periodic);    
    /* test */
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "%s/rhs_nodes/rhs_nodes_000.hdf", ud.file_name);
    sprintf(fieldname, "rhs-nodes");
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
#endif
    
    
    assert(integral_condition_nodes(rhs, node, x_periodic, y_periodic, z_periodic) != VIOLATED);
    operator_coefficients_nodes(hplus, hcenter, hgrav, elem, node, Sol, Sol0, mpv, dt);
    
    variable_coefficient_poisson_nodes(p2, hplus, hcenter, hgrav, rhs, x_periodic, y_periodic, z_periodic, dt);
        
    for(ii=0; ii<nc; ii++) {
        /* hunt down factor of 0.5 in the next line !! */
        mpv->p2_nodes[ii]  = SCND_PROJ_OLDP_WEIGHT*mpv->p2_nodes[ii] + p_update*SCND_PROJ_DELP_WEIGHT*0.5*(p2[ii]*dV);
#ifdef NODAL_PROJECTION_ONLY
        mpv->dp2_nodes[ii] += 0.5*(p2[ii]*dV);
#else /* NODAL_PROJECTION_ONLY */ 
        mpv->dp2_nodes[ii] = (p2[ii]*dV);
#endif /* NODAL_PROJECTION_ONLY */ 
    }

#if 0
    extern double *W0;
    int imax, jmax;
    rhs_max = 0.0;
    for (int jn=node->igy+1; jn<node->icy-node->igy-1; jn++) {
        int mn = jn*node->icx;
        for (int in=node->igx+1; in<node->icx-node->igx-1; in++) {
            int nn = mn+in;
            if (rhs_max < fabs(rhs[nn])) {
                rhs_max = fabs(rhs[nn]);
                imax = in;
                jmax = jn;
            } 
        }
    }
    printf("rhs_max before = %e\n", rhs_max);
    correction_nodes(Sol, elem, node, hplus, hgrav, mpv->dp2_nodes, t, dt);
    for(ii=0; ii<nc; ii++) rhs[ii] = 0.0;
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv->eta, mpv, bdry, dt, rhs_weight_new);
    catch_periodic_directions(rhs, node, elem, x_periodic, y_periodic, z_periodic);

#if 0
    FILE *prhs2file = NULL;
    char fn2[120], fieldname2[90];
    sprintf(fn2, "%s/rhs_nodes/rhs_nodes_post_000.hdf", ud.file_name);
    sprintf(fieldname2, "rhs-nodes-post");
    WriteHDF(prhs2file, node->icx, node->icy, node->icz, node->ndim, rhs, fn2, fieldname2);
    /*
    sprintf(fn2, "%s/rhs_nodes/lap_final.hdf", ud.file_name);
    sprintf(fieldname2, "lap-final");
    WriteHDF(prhs2file, node->icx, node->icy, node->icz, node->ndim, W0, fn2, fieldname2);
     */
#endif
    
    rhs_max = 0.0;
    for (int jn=node->igy+1; jn<node->icy-node->igy-1; jn++) {
        int mn = jn*node->icx;
        for (int in=node->igx+1; in<node->icx-node->igx-1; in++) {
            int nn = mn+in;
            if (rhs_max < fabs(rhs[nn])) {
                rhs_max = fabs(rhs[nn]);
                imax = in;
                jmax = jn;
            } 
        }
    }
    printf("rhs_max after  = %e\n", rhs_max);
#else
    correction_nodes(Sol, elem, node, hplus, hgrav, mpv->dp2_nodes, t, dt);
#endif

#ifndef NO_BDRYCONDS_PROJ2
    Bound_p_nodes(mpv, Sol, elem, node, 1);
#endif   
}

/* ========================================================================== */

static double divergence_nodes(
                               double* rhs,
                               const ElemSpaceDiscr* elem,
                               const NodeSpaceDiscr* node,
                               const ConsVars* Sol,
                               double* eta,
                               const MPV* mpv,
                               const BDRY* bdry,
                               const double dt,
                               const double weight) {
    
    extern User_Data ud;
    
    const int ndim = node->ndim;
    
    double div_max = 0.0;
    
    int i, j, k, mn;
    
    switch(ndim) {
        case 1: {
            ERROR("divergence_nodes() not implemented for 1D\n");
            break;
        }
        case 2: {
            
            /*
             const int igxn = node->igx;
             const int igyn = node->igy;
             const int icyn = node->icy;
             */
            const int icxn = node->icx;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double oodxdt = 1.0 / (dx * dt);
            const double oodydt = 1.0 / (dy * dt);
            const double oow = 1.0 / 2.0;
            const double todt = 2.0 / dt;
            const double oowdxdt = weight * oow * oodxdt;
            const double oowdydt = weight * oow * oodydt;
            
            double Y;
            
            /* predicted time level divergence via scattering */
            for(j = igye; j < icye - igye; j++) {
                const int me = j * icxe;
                const int mn = j * icxn;
                for(i = igxe; i < icxe - igxe; i++) {
                    const int n     = mn + i;
                    const int nicx  = n  + icxn;
                    const int n1    = n  + 1;
                    const int n1icx = n1 + icxn;
                    
                    const int ne = me + i;
                    double tmpfx, tmpfy, tmpfe;
                    
                    Y = Sol->rhoY[ne] / Sol->rho[ne];
                    tmpfx = oowdxdt * Y * Sol->rhou[ne];
                    tmpfy = oowdydt * Y * Sol->rhov[ne];
                    tmpfe = todt * eta[ne];
                    
                    rhs[n]           += + tmpfx + tmpfy - tmpfe;
                    rhs[n1]          += - tmpfx + tmpfy + tmpfe;
                    rhs[n1icx]       += - tmpfx - tmpfy - tmpfe;
                    rhs[nicx]        += + tmpfx - tmpfy + tmpfe;
                }
            }
            
            
            /* account for influx bottom boundary */
            j = igye;
            mn = j * icxn;
            Y  = mpv->HydroState->Y0[j];
            for(i = igxe; i < icxe - igxe; i++) {
                const int n     = mn + i;
                const int n1    = n  + 1;
                
                double rhov_wall = bdry->wall_massflux[i];
                double tmpy = oowdydt * Y * rhov_wall;
                
                rhs[n]  += - tmpy;
                rhs[n1] += - tmpy;
            }
            
            break;
        }
        case 3: {
            
            /*
             const int igxn = node->igx;
             const int igyn = node->igy;
             const int igzn = node->igz;
             const int iczn = node->icz;
             */
            const int icxn = node->icx;
            const int icyn = node->icy;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            const int igze = elem->igz;
            const int icze = elem->icz;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double dz = node->dz;
            
            const double oodxdt = 1.0 / (dx * dt);
            const double oodydt = 1.0 / (dy * dt);
            const double oodzdt = 1.0 / (dz * dt);
            
            const double oow = 1.0 / 4.0;
            const double oowdxdt = weight * oow * oodxdt;
            const double oowdydt = weight * oow * oodydt;
            const double oowdzdt = weight * oow * oodzdt;
            
            const int dixn = 1;
            const int diyn = icxn;
            const int dizn = icxn*icyn;
            
            double Y;
            
            /* predicted time level divergence via scattering */
            for(k = igze; k < icze - igze; k++) {
                const int le = k * icye * icxe;
                const int ln = k * icyn * icxn;
                for(j = igye; j < icye - igye; j++) {
                    const int me = le + j * icxe;
                    const int mn = ln + j * icxn;
                    for(i = igxe; i < icxe - igxe; i++) {
                        const int ne = me + i;
                        const int nn000  = mn + i;
                        const int nn010  = nn000 + diyn;
                        const int nn011  = nn000 + diyn + dizn;
                        const int nn001  = nn000        + dizn;
                        const int nn100  = nn000 + dixn;
                        const int nn110  = nn000 + dixn + diyn;
                        const int nn111  = nn000 + dixn + diyn + dizn;
                        const int nn101  = nn000 + dixn        + dizn;
                        
                        double tmpfx, tmpfy, tmpfz;
                        
                        Y = Sol->rhoY[ne] / Sol->rho[ne];
                        tmpfx = oowdxdt * Y * Sol->rhou[ne]; /* (rhou*Y) * 0.5 * / (dx*dt) */
                        tmpfy = oowdydt * Y * Sol->rhov[ne];
                        tmpfz = oowdzdt * Y * Sol->rhow[ne];
                        
                        rhs[nn000]       += + tmpfx + tmpfy + tmpfz;
                        rhs[nn010]       += + tmpfx - tmpfy + tmpfz;
                        rhs[nn011]       += + tmpfx - tmpfy - tmpfz;
                        rhs[nn001]       += + tmpfx + tmpfy - tmpfz;
                        rhs[nn100]       += - tmpfx + tmpfy + tmpfz;
                        rhs[nn110]       += - tmpfx - tmpfy + tmpfz;
                        rhs[nn111]       += - tmpfx - tmpfy - tmpfz;
                        rhs[nn101]       += - tmpfx + tmpfy - tmpfz;
                    }
                }
            }
            
            /* account for influx bottom boundary */
            j = igye;
            Y  = mpv->HydroState->Y0[j];
            for(k = igze; k < icze - igze; k++) {
                int ln = k * icyn*icxn;
                int mn = ln + j*icxn;
                for(i = igxe; i < icxe - igxe; i++) {
                    int nn   = mn + i;
                    int nn00 = nn;
                    int nn10 = nn + dixn;
                    int nn11 = nn + dixn + dizn;
                    int nn01 = nn +      + dizn;
                    
                    double rhov_wall = bdry->wall_massflux[i];
                    double tmpy = oowdydt * Y * rhov_wall;
                    
                    rhs[nn00] += - tmpy;
                    rhs[nn10] += - tmpy;
                    rhs[nn01] += - tmpy;
                    rhs[nn11] += - tmpy;
                }
            }
            
            /* set ghost values of rhs to zero
             j = igyn - 1;
             l = j * icxn;
             for(i = 0; i < icxn; i++)
             rhs[l + i] = 0.0;
             
             j = icyn - igyn;
             l = j * icxn;
             for(i = 0; i < icxn; i++)
             rhs[l + i] = 0.0;
             
             i = igxn - 1;
             l = i;
             for(j = 0; j < icyn; j++)
             rhs[l + j * icxn] = 0.0;
             
             i = icxn - igxn;
             l = i;
             for(j = 0; j < icyn; j++)
             rhs[l + j * icxn] = 0.0;
             */
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
    
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "rhs_nodes.hdf");
    sprintf(fieldname, "rhs-nodes");
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
#endif
    
    for (int in = 0; in<node->nc; in++) {
        div_max = MAX_own(div_max, fabs(rhs[in]));
    }
    return div_max;
    
}

/* ========================================================================== */

#ifdef P2_EXACT_PROJECTION_VK
void evaluate_eta(double* eta,
                  const ConsVars* Sol,
                  const ElemSpaceDiscr* elem)
{
    assert(elem->ndim == 2);
    
    const int icx = elem->icx;
    const int icy = elem->icz;
    
    const int igx = elem->igx;
    const int igy = elem->igy;
    
    const double dx = elem->dx;
    const double dy = elem->dy;
    
    for (int j=igy; j<icy-igy; j++) {int m = j*icx;
        for (int i=igx; i<icx-igx; i++) {int n = m+i;
            int nn = n+icx;
            int ns = n-icx;
            int ne = n+1;
            int nw = n-1;
            eta[n] =  (dy/dx) * (Sol->rhou[nn]*Sol->rhoY[nn]/Sol->rho[nn] - Sol->rhou[ns]*Sol->rhoY[ns]/Sol->rho[ns]) / dy;
            eta[n] += (dx/dy) * (Sol->rhov[ne]*Sol->rhoY[ne]/Sol->rho[ne] - Sol->rhov[nw]*Sol->rhoY[nw]/Sol->rho[nw]) / dx;
            eta[n] *= 0.125;
        }
    }
    
}

void eta_second_correction(MPV* mpv,
                           const double* hplus[3],
                           const ElemSpaceDiscr* elem,
                           const NodeSpaceDiscr* node,
                           const double dt)
{
    
    assert(elem->ndim == 2);
    
    const double* hplusx   = hplus[0];
    const double* hplusy   = hplus[1];
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    
    const int inx = node->icx;
    
    const int igx = elem->igx;
    const int igy = elem->igy;
    
    const double dx = elem->dx;
    const double dy = elem->dy;
    
    const double oodxsq = 1.0/(dx*dx);
    const double oodysq = 1.0/(dy*dy);
    
    double *eta0 = mpv->eta0;
    double *eta  = mpv->eta;
    double *p    = mpv->dp2_nodes;
    
    for (int j=igy; j<icy; j++) {
        int mc = j*icx;
        int mn = j*inx;
        for (int i=igx; i<icx; i++) {
            int nc  = mc+i;
            int nn  = mn+i;
            
            int nn1    = nn+1;
            int nninx  = nn+inx;
            int nn1inx = nn+1+inx;
            eta0[nc] = eta[nc];
            /* eta[nc]  -= 0.5*dt*0.125*(dy/dx+dx/dy)*(p[nn1inx] - p[nn1] - p[nninx] + p[nn]); */
            eta[nc]  -= 0.5*dt*0.125*(hplusx[nc]*oodxsq+hplusy[nc]*oodysq)*(p[nn1inx] - p[nn1] - p[nninx] + p[nn]);
        }
    }
}

#endif

/* ========================================================================== */

#ifdef SECOND_CORRECTION_IS_A_PROJECTION
static void controlled_variable_change_nodes(
                                             double* rhs,
                                             const ElemSpaceDiscr* elem,
                                             const NodeSpaceDiscr* node,
                                             const ConsVars* Sol,
                                             const ConsVars* Sol0,
                                             const double dt)
{
    extern User_Data ud;
    
    const int ndim = elem->ndim;
    int i, j;
    
    switch(ndim) {
        case 1: {
            ERROR("divergence_nodes() not implemented for 1D\n");
            break;
        }
        case 2: {
            
            const int icxn = node->icx;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            
            const double factor = 0.25 * ud.compressibility * 4.0 / (dt * dt);
            
            /* predicted time level divergence via scattering */
            for(j = igye; j < icye - igye; j++) {
                const int me = j * icxe;
                const int mn = j * icxn;
                for(i = igxe; i < icxe - igxe; i++) {
                    const int n     = mn + i;
                    const int nicx  = n  + icxn;
                    const int n1    = n  + 1;
                    const int n1icx = n1 + icxn;
                    const int ne = me + i;
                    
                    double d = factor*(Sol->rhoY[ne] - Sol0->rhoY[ne]);
                    
                    rhs[n]     += d;
                    rhs[n1]    += d;
                    rhs[n1icx] += d;
                    rhs[nicx]  += d;
                }
            }
            
            break;
        }
        case 3: {
            
            ERROR("controlled_variable_change_nodes() not implemented for 3D\n");
            
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
    
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "rhs_nodes.hdf");
    sprintf(fieldname, "rhs-nodes");
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
#endif
    
}
#endif

/* ========================================================================== */

static enum Constraint integral_condition_nodes(
                                                double* rhs,
                                                const NodeSpaceDiscr* node,
                                                const int x_periodic,
                                                const int y_periodic,
                                                const int z_periodic) {
    
    /* Rupert, this is valid only for periodic data so far!!! */
    
    int i, j, k, l, m, n;
    
    const int igx = node->igx;
    const int icx = node->icx;
    const int igy = node->igy;
    const int icy = node->icy;
    const int igz = node->igz;
    const int icz = node->icz;
    const int icxicy = icx * icy;
    
    const double del = DBL_EPSILON * node->nc;
    
    double tmp = 0.0, rhs_max=0.0;
    
    int periodic_shift_off = 1;
    
    for(k = igz; k < icz - igz - z_periodic*periodic_shift_off; k++) {l = k * icxicy;
        for(j = igy; j < icy - igy - y_periodic*periodic_shift_off; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx - x_periodic*periodic_shift_off; i++) {n = m + i;
                tmp += rhs[n];
                if(fabs(rhs[n])>rhs_max) rhs_max = fabs(rhs[n]);
            }
        }
    }
    
    /* if(fabs(tmp) > 100*del) { */
    if(fabs(tmp) > 10000*del*rhs_max) {
        printf("CHEATING:  tmp = %e,  rhs_max = %e, (tmp/rhs_max)/del = %e\n", tmp/rhs_max, rhs_max, (tmp/rhs_max)/del);
        /* return VIOLATED; */
        return VIOLATED;
    }
    else if(fabs(tmp) > 100*del*rhs_max) {
        printf("integral_condition_nodes = barely OK; tmp = %e,  rhs_max = %e, (tmp/rhs_max)/del = %e\n", tmp/rhs_max, rhs_max, (tmp/rhs_max)/del);
        return SATISFIED;
    }
    else {
        /* printf("integral_condition_nodes = OK;    tmp = %e,  rhs_max = %e\n", tmp, rhs_max); */
        printf("integral_condition_nodes = OK; tmp = %e\n", tmp);
        return SATISFIED;
    }
}

/* ========================================================================== */

static void	catch_periodic_directions(
                                      double* rhs,
                                      const NodeSpaceDiscr* node,
                                      const ElemSpaceDiscr* elem,
                                      const int x_periodic,
                                      const int y_periodic,
                                      const int z_periodic) {
    
    /* in any of the periodicity directions, gather boundary values of
     rhs on first nodes, set last nodes to zero
     */
    
    int i, j, k, l, m, l0, l1, m0, m1, n0, n1;
    
    const int igx = node->igx;
    const int icx = node->icx;
    const int igy = node->igy;
    const int icy = node->icy;
    const int igz = node->igz;
    const int icz = node->icz;
    const int icxicy = icx * icy;
    
    if(x_periodic == 1){
        for(k = igz; k < icz - igz; k++) {l = k * icxicy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                n0 = m + igx;
                n1 = m + icx-igx-1;
                rhs[n0] += rhs[n1];
                rhs[n1] = 0.0;
            }
        }
    }
    
    if(y_periodic == 1){
        for(k = igz; k < icz - igz; k++) {l = k * icxicy;
            m0 = l + igy * icx;
            m1 = l + (icy-igy-1) * icx;
            for(i = igx; i < icx - igx; i++) {
                n0 = m0 + i;
                n1 = m1 + i;
                rhs[n0] += rhs[n1];
                rhs[n1] = 0.0;
            }
        }
    }
    
    if(z_periodic == 1){
        l0 = igz * icxicy;
        l1 = (icz-igz-1) * icxicy;
        for(j = igy; j < icy - igy; j++) {
            m0 = l0 + j * icx;
            m1 = l1 + j * icx;
            for(i = igx; i < icx - igx; i++) {
                n0 = m0 + i;
                n1 = m1 + i;
                rhs[n0] += rhs[n1];
                rhs[n1] = 0.0;
            }
        }
    }
}

/* ========================================================================== */

static void operator_coefficients_nodes(
                                        double* hplus[3],
                                        double* hcenter,
                                        double* hgrav,
                                        const ElemSpaceDiscr* elem,
                                        const NodeSpaceDiscr* node,
                                        ConsVars* Sol,
                                        const ConsVars* Sol0,
                                        const MPV* mpv,
                                        const double dt) {
    
    extern User_Data ud;
    extern Thermodynamic th;
    
    const int ndim = node->ndim;
    
    const int impl_grav_th2 = ud.implicit_gravity_theta2;
    const int impl_grav_pr = ud.implicit_gravity_press2;
    
    /* #define CHECK_MINMAX */
#ifdef CHECK_MINMAX
    double coeffmin = 1000000.0;
    double coeffmax = -1000000.0;
#endif
    
    switch(ndim) {
        case 1: {
            ERROR("interface_enthalpy_nodes() not implemented for 1D\n");
            break;
        }
        case 2: {
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            
            const double dx = elem->dx;
            const double dy = elem->dy;
            
            double Msq = ud.Msq;
            
            double* hplusx  = hplus[0];
            double* hplusy  = hplus[1];
            double* hc      = hcenter;
            double* hg      = hgrav;
#if 1
#ifdef SECOND_CORRECTION_IS_A_PROJECTION
            assert(ud.p_flux_correction);
            const double ccenter = 0.0;
#else
            const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt);
#endif
#else
            const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt);
#endif            
            const double cexp    = 1.0-th.gamm;
            int i, j, m, n;
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = 0.0;
            
            for(j = igy; j < icy - igy; j++) {
                m = j * icx;
                for(i = igx; i < icx - igx; i++) {
                    n = m + i;
                    
                    double theta   = Sol->rhoY[n] / Sol->rho[n] ;
                    double dthetax = ((Sol->rhoY[n+1]  / Sol->rho[n+1]) - (Sol->rhoY[n-1]  / Sol->rho[n-1])) / (2.0*dx);
                    double gimpx   = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[0]/Msq)*dthetax/theta);
                    
                    double dthetay = ((Sol->rhoY[n+icx]  / Sol->rho[n+icx]) - (Sol->rhoY[n-icx]  / Sol->rho[n-icx])) / (2.0*dy);
                    double gimpy    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetay/theta);
                                        
                    hplusx[n]  = theta * gimpx;
                    hplusy[n]  = theta * gimpy;
                    
                    hg[n]      = impl_grav_pr2 * th.gamminv * pow(Sol->rhoY[n],cexp) * ud.gravity_strength[1] * gimpy * 0.5; This right?
                    
                    /* hg[n]      = impl_grav_pr2 * th.gamminv * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp) * ud.gravity_strength[1] * gimpy * 0.5; */
                    /* hg[n]      = impl_grav_pr2 * th.gamminv * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp) * (ud.gravity_strength[1]/Msq) * gimpy; */
                    /* hg[n]      = 0.0 * th.gamminv * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp) * ud.gravity_strength[1] * gimpy; */
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
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int igz = elem->igz;
            const int icz = elem->icz;
            
            const double dx = elem->dx;
            const double dy = elem->dy;
            const double dz = elem->dz;
            
            const int dix = 1;
            const int diy = icx;
            const int diz = icx*icy;
            
            double Msq = ud.Msq;
            
            double* hplusx  = hplus[0];
            double* hplusy  = hplus[1];
            double* hplusz  = hplus[2];
            double* hc      = hcenter;
            
            const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt);            
            const double cexp    = 1.0-th.gamm;
            
            int i, j, k, l, m, n;
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = hplusz[i] = 0.0;
            
            for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        {
                            double theta   = Sol->rhoY[n] / Sol->rho[n] ;
                            
                            double dthetax = ((Sol->rhoY[n+dix]  / Sol->rho[n+dix]) - (Sol->rhoY[n-dix]  / Sol->rho[n-dix])) / (2.0*dx);
                            double gimpx   = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[0]/Msq)*dthetax/theta);
                            
                            double dthetay = ((Sol->rhoY[n+diy]  / Sol->rho[n+diy]) - (Sol->rhoY[n-diy]  / Sol->rho[n-diy])) / (2.0*dy);
                            double gimpy    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetay/theta);
                            
                            double dthetaz = ((Sol->rhoY[n+diz]  / Sol->rho[n+diz]) - (Sol->rhoY[n-diz]  / Sol->rho[n-diz])) / (2.0*dz);
                            double gimpz    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetaz/theta);
                                                        
                            hplusx[n]  = theta * gimpx;
                            hplusy[n]  = theta * gimpy;
                            hplusz[n]  = theta * gimpz;
                        }
                    }
                }
            }
            
            for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;
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


#if 0
/* ========================================================================== */

static void diagonal_preconditioner(
                                    double* diaginv,
                                    const ElemSpaceDiscr* elem,
                                    const NodeSpaceDiscr* node,
                                    ConsVars* Sol) {
    
    extern User_Data ud;
    
    const int ndim = node->ndim;
    
    switch(ndim) {
        case 1: {
            ERROR("interface_enthalpy_nodes() not implemented for 1D\n");
            break;
        }
        case 2: {
            
            int n;
            
#ifdef PRECON
            int i, j, m;
            const int igxn = node->igx;
            const int icxn = node->icx;
            const int igyn = node->igy;
            const int icyn = node->icy;
            
            const int icxe = elem->icx;
            
            for(j = igyn; j < icyn - igyn; j++) {m = j * icxn;
                for(i = igxn; i < icxn - igxn; i++) {n = m + i;
                    {
                        int nne, nnw, nse, nsw;
                        
                        nne  =   i   +   j   * icxe;
                        nnw  = (i-1) +   j   * icxe;
                        nsw  = (i-1) + (j-1) * icxe;
                        nse  =   i   + (j-1) * icxe;
                        
                        /* diaginv[n] = 4.0 / (Sol->rhoY[nne] + Sol->rhoY[nnw] + Sol->rhoY[nsw] + Sol->rhoY[nse]); */
                        diaginv[n] = 4.0 / (  Sol->rhoY[nne] / Sol->rho[nne]
                                            + Sol->rhoY[nnw] / Sol->rho[nnw]
                                            + Sol->rhoY[nsw] / Sol->rho[nsw]
                                            + Sol->rhoY[nse] / Sol->rho[nse]);
                    }
                }
            }
#else
            for(n=0; n<node->nc; n++) {
                diaginv[n] = 1.0;
            }
#endif
            
            break;
        }
        case 3: {
            
#ifdef PRECON
            int i, j, k, l, m, n;
            const int igxn = node->igx;
            const int icxn = node->icx;
            const int igyn = node->igy;
            const int icyn = node->icy;
            const int igzn = node->igz;
            const int iczn = node->icz;
            
            const int icxe = elem->icx;
            const int icye = elem->icy;
            
            const int dixe = 1;
            const int diye = icxe;
            const int dize = icxe*icye;
            
            for(k = igzn; k < iczn - igzn; k++) {l = k * icxn*icyn;
                for(j = igyn; j < icyn - igyn; j++) {m = l + j * icxn;
                    for(i = igxn; i < icxn - igxn; i++) {n = m + i;
                        int nc     = (k-1)*dize + (j-1)*diye + (i-1)*dixe;
                        int nc000  = nc;
                        int nc010  = nc        + diye;
                        int nc011  = nc        + diye + dize;
                        int nc001  = nc               + dize;
                        int nc100  = nc + dixe;
                        int nc110  = nc + dixe + diye;
                        int nc111  = nc + dixe + diye + dize;
                        int nc101  = nc + dixe        + dize;
                        
                        /* diaginv[n] = 4.0 / (Sol->rhoY[nne] + Sol->rhoY[nnw] + Sol->rhoY[nsw] + Sol->rhoY[nse]); */
                        diaginv[n] = 8.0 / (  Sol->rhoY[nc000] / Sol->rho[nc000]
                                            + Sol->rhoY[nc010] / Sol->rho[nc010]
                                            + Sol->rhoY[nc011] / Sol->rho[nc011]
                                            + Sol->rhoY[nc001] / Sol->rho[nc001]
                                            + Sol->rhoY[nc100] / Sol->rho[nc100]
                                            + Sol->rhoY[nc110] / Sol->rho[nc110]
                                            + Sol->rhoY[nc111] / Sol->rho[nc111]
                                            + Sol->rhoY[nc101] / Sol->rho[nc101]);
                    }
                }
            }
#else /* PRECON */
            int n;
            
            for(n=0; n<node->nc; n++) {
                diaginv[n] = 1.0;
            }
#endif
            
            break;
        }
        default: ERROR("ndim not in {1,2,3}");
    }
}
#endif

/* ========================================================================== */

void correction_nodes(
                      ConsVars* Sol,
                      const ElemSpaceDiscr* elem,
                      const NodeSpaceDiscr* node,
                      const double* hplus[3], 
                      const double* hgrav,
                      double* p,
                      const double t,
                      const double dt) {
    
    extern User_Data ud;
    
    const int ndim = elem->ndim;
    
    
#ifndef NO_BDRYCONDS_PROJ2
    set_ghostnodes_p(p, node, 1);   
#endif 
    
    switch(ndim) {
        case 1: {
            ERROR("function not available");
            break;
        }
            
        case 2: {
            
            const int igx = node->igx;
            const int icx = node->icx;
            const int igy = node->igy;
            const int icy = node->icy;
            
            const int icxe = node->icx - 1;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double oodx = 1.0 / dx;
            const double oody = 1.0 / dy;
            const double dtowdx = 0.5*dt * oodx;
            const double dtowdy = 0.5*dt * oody;
            
            int i, j, m, me;
            
            for(j = igy; j < icy - igy - 1; j++) {
                m = j * icx; 
                me = j * icxe;
                
                for(i = igx; i < icx - igx - 1; i++) {
                    const int n = m + i;
                    const int nicx = n + icx;
                    const int n1 = n + 1;
                    const int n1icx = n1 + icx;
                    
                    const int ne = me + i; 
                    
                    const double Dpx = 0.5 * (p[n1]   - p[n] + p[n1icx] - p[nicx]);
                    
                    const double Dpy = 0.5 * (p[nicx] - p[n] + p[n1icx] - p[n1]);
                    
                    Sol->rhou[ne] += - dtowdx * Dpx;
                    Sol->rhov[ne] += - dtowdy * Dpy;
                    Sol->rhov[ne] += - 0.5*dt * hgrav[ne] * (Sol->rho[ne]/Sol->rhoY[ne]) * 0.25 * (p[n] + p[n1] + p[n1icx] + p[nicx]);
                }
            } 
            
            break;
        }
        case 3: {
            
            extern User_Data ud;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            const int igze = elem->igz;
            const int icze = elem->icz;
            
            const int dixe = 1;
            const int diye = icxe;
            const int dize = icxe*icye;
            
            const int icxn = node->icx;
            const int icyn = node->icy;
            /*
             const int iczn = node->icz;    
             */
            const int dixn = 1;
            const int diyn = icxn;
            const int dizn = icxn*icyn;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double dz = node->dz;
            const double oodx = 1.0 / dx;
            const double oody = 1.0 / dy;
            const double oodz = 1.0 / dz;
            const double dtowdx = 0.5*dt * oodx;
            const double dtowdy = 0.5*dt * oody;
            const double dtowdz = 0.5*dt * oodz;
            
            int i, j, k, ln, mn, nn, le, me, ne;
            
            for(k = igze; k < icze - igze; k++) {
                le = k * dize; 
                ln = k * dizn;
                
                for(j = igye; j < icye - igye; j++) {
                    me = le + j * diye; 
                    mn = ln + j * diyn;
                    
                    for(i = igxe; i < icxe - igxe; i++) {
                        ne = me + i * dixe;
                        nn = mn + i * dixn;
                        
                        const int nn000 = nn;
                        const int nn010 = nn        + diyn;
                        const int nn011 = nn        + diyn + dizn;
                        const int nn001 = nn               + dizn;
                        const int nn100 = nn + dixn;
                        const int nn110 = nn + dixn + diyn;
                        const int nn111 = nn + dixn + diyn + dizn;
                        const int nn101 = nn + dixn        + dizn;
                        
                        const double Dpx = 0.25 * (  p[nn100] - p[nn000]
                                                   + p[nn110] - p[nn010]
                                                   + p[nn111] - p[nn011]
                                                   + p[nn101] - p[nn001]
                                                   );
                        
                        const double Dpy = 0.25 * (  p[nn010] - p[nn000]
                                                   + p[nn110] - p[nn100]
                                                   + p[nn111] - p[nn101]
                                                   + p[nn011] - p[nn001]
                                                   );
                        
                        const double Dpz = 0.25 * (  p[nn001] - p[nn000]
                                                   + p[nn101] - p[nn100]
                                                   + p[nn111] - p[nn110]
                                                   + p[nn011] - p[nn010]
                                                   );

                        Sol->rhou[ne] += - dtowdx * Dpx;
                        Sol->rhov[ne] += - dtowdy * Dpy;
                        Sol->rhow[ne] += - dtowdz * Dpz;
                    }
                } 
            } 
            break;
        }
        default: ERROR("ndim not in {1, 2,3}");
    }
}

#else /* SOLVER_2_HYPRE_RUPE */

#include "/Users/rupert/Documents/Computation/Walls/wallPoissonTest/Version_2012.05.27/rupert_latest_patched/nodePoisson.h"
#include "/Users/rupert/Documents/Computation/Walls/wallPoissonTest/Version_2012.05.27/rupert_latest_patched/cellPoisson.h"
#include "enumerator.h"

static double* beta[3];
static double* b[3];
static double grad_bdry[3];

static enum Constraint integral_condition(
                                          const double* rhs,
                                          const NodeSpaceDiscr* node,
                                          const double* oovol);

static enum Constraint integral_condition_inn(
                                              double* rhs,
                                              const NodeSpaceDiscr* node,
                                              const int bdry_switch);

static double IntegralRHS_2(
                            double* rhs,
                            const ElemSpaceDiscr* elem,
                            const NodeSpaceDiscr* node,
                            ConsVars* Sol,
                            const MPV* mpv,
                            const double dt);

#ifdef MICHAEL_CORRECTION
static void correction_michael(
                               ConsVars *Sol,
                               const ElemSpaceDiscr *elem,
                               const NodeSpaceDiscr *node,
                               const double *int_p,
                               const double t,
                               const double dt);
#else /* MICHAEL_CORRECTION */
static void correction(
                       ConsVars* Sol,
                       const ElemSpaceDiscr* elem,
                       const NodeSpaceDiscr* node,
                       double* p2,
                       const double t,
                       const double dt);
#endif /* MICHAEL_CORRECTION */

#define CORRECTION_FROM_MG_CODE
#ifdef CORRECTION_FROM_MG_CODE
static void operator_coefficients_nodes(
                                        double* hplus[3],
                                        double* hcenter,
                                        double* hgrav,
                                        const ElemSpaceDiscr* elem,
                                        const NodeSpaceDiscr* node,
                                        ConsVars* Sol,
                                        const ConsVars* Sol0,
                                        const MPV* mpv,
                                        const double dt);
#else
static void interface_coefficients(
                                   double* h[3],
                                   const ConsVars* Sol,
                                   const ElemSpaceDiscr* elem,
                                   const NodeSpaceDiscr* node);
#endif

void periodic(
              double test[5],
              double *rhs,
              double *h0,
              double *h1,
              const NodeSpaceDiscr *node);


/* #define FIX_PERIODIC_BCS_2ND_PROJ */
#ifdef FIX_PERIODIC_BCS_2ND_PROJ
void fix_periodic_boundaries_2nd_projection_inn(
                                                double *rhs,
                                                double* h[3],
                                                const NodeSpaceDiscr* node);
#endif

/* =========================================================== */

void initSecondProjection(
                          const ElemSpaceDiscr *elem,
                          const NodeSpaceDiscr *node)
{
    extern double grad_bdry[3];
    
    static enum Boolean beta_allocated = WRONG;
    
    int i, j, n, nf[3];
    int idim;
    
    double dx[3] = {node->dx, node->dy, node->dz};
    
    int inx = node->icx - 2*node->igx;
    int iny = node->icy - 2*node->igy;
    int inz = node->icz - 2*node->igz;
    
    int ifx = node->ifx - 2*node->igx;
    int ify = node->ify - 2*node->igy;
    int ifz = node->ifz - 2*node->igz;
    
    nf[0] = ifx*iny*inz;
    nf[1] = inx*ify*inz;
    nf[2] = inx*iny*ifz;
    
    grad_bdry[0] = grad_bdry[1] = 0.0;
    
    if (beta_allocated == WRONG) {
        for (idim = 0; idim < node->ndim; idim++) {
            beta[idim] = (double*)malloc((2*nf[idim])*sizeof(double));
            memset(beta[idim], 1.0, (2*nf[idim])*sizeof(double));
            b[idim]    = (double*)malloc((2*nf[idim])*sizeof(double));
            memset(b[idim], 1.0, (2*nf[idim])*sizeof(double));
        }
        beta_allocated = CORRECT;
    }
    
    /* Note: the betas live inside the flow domain only, they don't exist on the dummy cells */
    getNodeIntegralsScalar(dx[0], dx[1], inx, iny, b[0], b[1], beta[0], beta[1]);
    
    /*
     make betas into area fractions instead of actual areas and
     flip beta[0] in the x-y-coordinates to comply with Rupert's memory scheme
     */
    {
        const double dxinv = 1.0/node->dx;
        const double dyinv = 1.0/node->dy;
        
        for (n=0; n<2*nf[0]; n++) {
            b[0][n] = beta[0][n]*dyinv;
        }
        for (n=0; n<2*nf[1]; n++) {
            beta[1][n] *= dxinv;
        }
        
        for (j=0; j<iny; j++) {
            int m  = j*ifx;
            int mf = j;
            for (i=0; i<ifx; i++) {
                int n  = m + i;
                int nf = mf + i*iny;
                beta[0][2*n]   = b[0][2*nf];
                beta[0][2*n+1] = b[0][2*nf+1];
            }
        }
    }
}

/* =========================================================== */

static enum Boolean second_projection_is_initialized = WRONG;


void second_projection(ConsVars* Sol,
                       MPV* mpv,
                       const ConsVars* Sol0,
                       const ElemSpaceDiscr* elem,
                       const NodeSpaceDiscr* node,
                       const double p_update,
                       const double t,
                       const double dt)
{
    
    extern double grad_bdry[3];
    extern double *W0, *W1;
    
    const double dx = node->dx;
    const double dy = node->dy;
    const double dz = node->dz;
    const double dV = dx*dy*dz;
    
    const int nc = node->nc;
    const int igx = node->igx;
    const int igy = node->igy;
    const int igz = node->igz;
    const int icx = node->icx;
    const int icy = node->icy;
    const int icz = node->icz;
    
    const int icx_inn = node->icx-2*node->igx;
    const int icy_inn = node->icy-2*node->igy;
    
    const int nfx = node->nfx;
    const int nfy = node->nfy;
    const int nfz = node->nfz;
    
    double test[5];
    
    int i, j, k, l, m, n;
    double rhs_max_before;
    
    double* p2    = (double*)malloc(nc * sizeof(double));
    double* rhs   = (double*)malloc(nc * sizeof(double));
    double* oovol = (double*)malloc(nc * sizeof(double));

    double** hplus    = mpv->Level[0]->wplus;
    double*  hcenter  = mpv->Level[0]->wcenter;
    double*  hgrav    = mpv->Level[0]->wgrav;

    double* h[3];
    
    printf("\n\n====================================================");
    printf("\nSecond Projection");
    printf("\n====================================================\n");
    
    if (second_projection_is_initialized == WRONG) {

        double* fakeG = (double*)malloc(nc * sizeof(double));
        
        for (int jc = 0; jc < icy; jc++) {int mc = jc*icx;
            for (int ic = 0; ic < icx; ic++) {int nc = mc + ic;
                fakeG[nc] = -1.0 - sqrt((dx*(ic-0.5*icx))*(dx*(ic-0.5*icx)) + (dy*(jc-0.5*icy))*(dy*(jc-0.5*icy)));
            }
        }
        
        initNodePoisson(node->x[node->igx], node->y[node->igy], node->dx, node->dy, icx_inn, icy_inn, fakeG);
        initSecondProjection(elem, node);
        second_projection_is_initialized = CORRECT;
        free(fakeG);
    }
    
    h[0] = (double*)malloc(2*node->nfx * sizeof(double));
    h[1] = (double*)malloc(2*node->nfy * sizeof(double));
    h[2] = (double*)malloc(2*node->nfz * sizeof(double));
    
    double *int_xn, *int_yn, *int_zn, *int_p;
    int_xn = (double*)malloc(2*node->nfx * sizeof(double));
    int_yn = (double*)malloc(2*node->nfy * sizeof(double));
    int_zn = (double*)malloc(2*node->nfz * sizeof(double));
    
    int_p = (double*)malloc(2*node->nc * sizeof(double));
    
    memset(p2,   0, nc  * sizeof(double));
    memset(rhs,  0, nc  * sizeof(double));
    memset(h[0], 0, 2*nfx * sizeof(double));
    memset(h[1], 0, 2*nfy * sizeof(double));
    memset(h[2], 0, 2*nfz * sizeof(double));
    
    for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                oovol[n] = 1.0 / (dx*dy*dz);
            }
        }
    }
    
#ifdef CORRECTION_FROM_MG_CODE
    operator_coefficients_nodes(hplus, hcenter, hgrav, elem, node, Sol, Sol0, mpv, dt);
    map2D_to_michaels_memory_nodefaces(h, hplus, elem, node);
#else
    interface_coefficients(h, Sol, elem, node);
    map2D_to_michaels_memory_nodefaces_x(h[0],node,W0,W1);
    map2D_to_michaels_memory_nodefaces_y(h[1],node,W0,W1);
#endif
    
#if 0
    extern User_Data ud;
    FILE *phfile = NULL;
    char fn[100], fieldname[90];
    
    sprintf(fn, "%s/Tests/hx_test.hdf", ud.file_name);
    sprintf(fieldname, "hx");
    
    WriteHDF(phfile,
             2*mpv->Level[0]->node->ifx,
             mpv->Level[0]->node->icy,
             mpv->Level[0]->node->icz,
             mpv->Level[0]->elem->ndim,
             h[0],
             fn,
             fieldname);
    
    sprintf(fn, "%s/Tests/hy_test.hdf", ud.file_name);
    sprintf(fieldname, "hy");
    
    WriteHDF(phfile,
             2*mpv->Level[0]->node->ify,
             mpv->Level[0]->node->icx,
             mpv->Level[0]->node->icz,
             mpv->Level[0]->node->ndim,
             h[1],
             fn,
             fieldname);
    
#endif

    rhs_max_before = IntegralRHS_2(rhs, elem, node, Sol, mpv, dt);
    
    assert(integral_condition(rhs, node, 0) != VIOLATED);
    
    for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                rhs[n] *= oovol[n];
                /* rhs[n] /= (dx*dy*dz);  */
            }
        }
    }
    
    map2D_to_michaels_memory_nodes(rhs,node,W0);
    
    periodic(test, rhs, h[0], h[1], node);
    
    grad_bdry[0] = 0.0;
    grad_bdry[1] = 0.0;
    grad_bdry[2] = 0.0;
    
    solveNodePoisson(dx, dy, icx_inn, icy_inn, h[0], h[1], rhs, p2);
    
#if 0
    extern User_Data ud;
    FILE *prhsfile = NULL;
    char fn2[100], fieldname2[90];
    
    sprintf(fn2, "%s/Tests/rhs_test.hdf", ud.file_name);
    sprintf(fieldname2, "rhs");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->node->icy-2*mpv->Level[0]->node->igy,
             mpv->Level[0]->node->icx-2*mpv->Level[0]->node->igx,
             mpv->Level[0]->node->icz,
             mpv->Level[0]->elem->ndim,
             rhs,
             fn2,
             fieldname2);
    
    sprintf(fn2, "%s/Tests/p2_test.hdf", ud.file_name);
    sprintf(fieldname2, "p2");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->node->icy-2*mpv->Level[0]->node->igy,
             mpv->Level[0]->node->icx-2*mpv->Level[0]->node->igx,
             mpv->Level[0]->node->icz,
             mpv->Level[0]->elem->ndim,
             p2,
             fn2,
             fieldname2);
#endif
    
    map2D_to_ruperts_memory_nodes(p2, node, W0);
    for(int ii=0; ii<nc; ii++) {
        /* hunt down factor of 0.5 in the next line !! */
        mpv->p2_nodes[ii]  = SCND_PROJ_OLDP_WEIGHT*mpv->p2_nodes[ii] + p_update*SCND_PROJ_DELP_WEIGHT*0.5*(p2[ii]*dV);
        mpv->dp2_nodes[ii] = (p2[ii]*dV);
    }
#ifndef NO_BDRYCONDS_PROJ2
    set_ghostnodes_p(mpv->p2_nodes,node,node->igx);
    set_ghostnodes_p(mpv->dp2_nodes,node,node->igx);
#endif
    correction(Sol, elem, node, mpv->dp2_nodes, t, dt);

    free(p2);
    free(rhs);
    free(oovol);
    free(h[0]);
    free(h[1]);
    free(h[2]);
    free(int_xn);
    free(int_yn);
    free(int_zn);
    free(int_p);
}


/* ========================================================== */

static double IntegralRHS_2(
                            double* rhs,
                            const ElemSpaceDiscr* elem,
                            const NodeSpaceDiscr* node,
                            ConsVars* Sol,
                            const MPV* mpv,
                            const double dt)
{
    
    extern User_Data ud;
    
    double rhs_max = 0.0;
    
    memset(rhs, 0, node->nc * sizeof(double));
    
    switch(node->ndim) {
        case 1: {
            ERROR("function not available");
            break;
        }
            
        case 2: {
            
            int i, j;
            
            const int icxe = elem->icx;
            const int icye = elem->icy;
            
            const int icxn = node->icx;
            
            const int igx = node->igx;
            const int igy = node->igy;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double dxodt = dx / dt;
            const double dyodt = dy / dt;
            
            for(j = igy; j < icye-igy; j++) {
                /* loop indices run over primary cells within domain */
                const int me      = j * icxe;
                const int mn      = j * icxn;
                
                for(i = igx; i < icxe-igx; i++) {
                    const int ne     = me + i;
                    const int nn     = mn + i;
                    const int nn1    = nn + 1;
                    const int nnicx  = nn + icxn;
                    const int nn1icx = nn + icxn + 1;
                    
                    const double Y    = Sol->rhoY[ne] / Sol->rho[ne];
                    const double tmpx = dyodt * Y * Sol->rhou[ne];
                    const double tmpy = dxodt * Y * Sol->rhov[ne];
                    
                    rhs[nn]     += + tmpx + tmpy;
                    rhs[nn1]    += - tmpx + tmpy;
                    rhs[nn1icx] += - tmpx - tmpy;
                    rhs[nnicx]  += + tmpx - tmpy;
                    
                }
            }
            
            for (i=0; i<node->nc; i++) {
                rhs_max = MAX_own(rhs_max,fabs(rhs[i]));
            }
            
            break;
        }
        case 3: { /* 1.6.99 */
            
            int i, j, k, l, m, le, me;
            const int igx = node->igx;
            const int icx = node->icx;
            const int igy = node->igy;
            const int icy = node->icy;
            const int igz = node->igz;
            const int icz = node->icz;
            const int icxe = icx - 1;
            const int icye = icy - 1;
            const double dx = node->dx;
            const double dy = node->dy;
            const double dz = node->dz;
            const double oow = 1.0 / 4.0;
            const double oowdydzodt = oow * dy * dz / dt;
            const double oowdzdxodt = oow * dz * dx / dt;
            const double oowdxdyodt = oow * dx * dy / dt;
            
            for(k = igz; k < icz - igz - 1; k++) {
                l = k * icy * icx;
                le = k * icye * icxe;
                for(j = igy; j < icy - igy - 1; j++) {
                    m = l + j * icx;
                    me = le + j * icxe;
                    for(i = igx; i < icx - igx - 1; i++) {
                        const int n           = m + i;
                        const int nicx        = n + icx;
                        const int n1          = n + 1;
                        const int n1icx       = n1 + icx;
                        const int nicxicy     = n + icx * icy;
                        const int nicxicxicy  = nicx + icx * icy;
                        const int n1icxicy    = n1 + icx * icy;
                        const int n1icxicxicy = n1icx + icx * icy;
                        const int ne          = me + i;
                        
                        const int nfx = 0;
                        const int nfy = 0;
                        const int nfz = 0;
                        
                        assert(node->ndim < 3); /* beta indices not correct yet */
                        
                        const double Y    = Sol->rhoY[ne] / Sol->rho[ne];
                        const double tmpx = oowdydzodt * Y * Sol->rhou[ne] * beta[0][nfx];
                        const double tmpy = oowdzdxodt * Y * Sol->rhov[ne] * beta[1][nfy];
                        const double tmpz = oowdxdyodt * Y * Sol->rhow[ne] * beta[2][nfz];
                        
                        rhs[n]           += + tmpx + tmpy + tmpz;
                        rhs[n1]          += - tmpx + tmpy + tmpz;
                        rhs[n1icx]       += - tmpx - tmpy + tmpz;
                        rhs[nicx]        += + tmpx - tmpy + tmpz;
                        rhs[nicxicy]     += + tmpx + tmpy - tmpz;
                        rhs[n1icxicy]    += - tmpx + tmpy - tmpz;
                        rhs[n1icxicxicy] += - tmpx - tmpy - tmpz;
                        rhs[nicxicxicy]  += + tmpx - tmpy - tmpz;
                        
                    }
                }
            }
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
    
    return(rhs_max);
}

/* ========================================================================= */

static enum Constraint integral_condition(
                                          const double* rhs,
                                          const NodeSpaceDiscr* node,
                                          const double* oovol)
{
    
    int i, j, k, l, m, n;
    
    const int igx = node->igx;
    const int icx = node->icx;
    const int igy = node->igy;
    const int icy = node->icy;
    const int igz = node->igz;
    const int icz = node->icz;
    const int icxicy = icx * icy;
    
    double tmp = 0.0;
    for(k = igz; k < icz - igz; k++) {l = k * icxicy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                tmp += rhs[n];
            }
        }
    }
    if(fabs(tmp) > 100*node->nc*DBL_EPSILON) {
        printf("tmp = %e\n", tmp);
        return VIOLATED;
    }
    else {
        return SATISFIED;
    }
}

/* ========================================================================= */

static enum Constraint integral_condition_inn(
                                              double* rhs,
                                              const NodeSpaceDiscr* node,
                                              const int bdry_switch)
{
    
    extern User_Data ud;
    
    int i, j, k, l, m, n;
    
    const int igx = node->igx;
    const int icx = node->icx-2*igx;
    const int igy = node->igy;
    const int icy = node->icy-2*igy;
    const int igz = node->igz;
    const int icz = node->icz-2*igz;
    
    double tmp = 0.0;
    double rhs_max = 0.0;
    
    for(k = 0; k < icz; k++) {l = k;
        for(j = 0; j < icy; j++) {m = l + j * icz;
            for(i = 0; i < icx; i++) {n = m + i*icz*icy;
                tmp += rhs[n];
                rhs_max = MAX_own(rhs_max, fabs(rhs[n]));
            }
        }
    }
    
    if(fabs(tmp) > 100000*DBL_EPSILON) {
        printf("tmp = %e\n", tmp);
        return VIOLATED;
        /* return SATISFIED; */
    }
    else {
        return SATISFIED;
    }
}

#ifdef CORRECTION_FROM_MG_CODE

/* ========================================================================== */

static void operator_coefficients_nodes(
                                        double* hplus[3],
                                        double* hcenter,
                                        double* hgrav,
                                        const ElemSpaceDiscr* elem,
                                        const NodeSpaceDiscr* node,
                                        ConsVars* Sol,
                                        const ConsVars* Sol0,
                                        const MPV* mpv,
                                        const double dt) {
    
    extern User_Data ud;
    extern Thermodynamic th;
    
    const int ndim = node->ndim;
    
    const int impl_grav_th2 = ud.implicit_gravity_theta2;
    const int impl_grav_pr2 = ud.implicit_gravity_press2;
    
    const int thermcon         = ud.thermcon;
    
    /* #define CHECK_MINMAX */
#ifdef CHECK_MINMAX
    double coeffmin = 1000000.0;
    double coeffmax = -1000000.0;
#endif
    
    switch(ndim) {
        case 1: {
            ERROR("interface_enthalpy_nodes() not implemented for 1D\n");
            break;
        }
        case 2: {
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            
            const double dx = elem->dx;
            const double dy = elem->dy;
            
            double Msq = ud.Msq;
            
            double* hplusx  = hplus[0];
            double* hplusy  = hplus[1];
            double* hc      = hcenter;
            double* hg      = hgrav;

#ifdef SECOND_CORRECTION_IS_A_PROJECTION
            assert(ud.p_flux_correction);
            const double ccenter = 0.0;
#else
            const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt);
#endif            
            const double cexp    = 1.0-th.gamm;
            int i, j, m, n;
            
            for(i=0; i<node->nc; i++) hc[i] = hplusx[i] = hplusy[i] = 0.0;
            
            for(j = igy; j < icy - igy; j++) {m = j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    
                    double theta   = Sol->rhoY[n] / Sol->rho[n] ;
                    double dthetax = ((Sol->rhoY[n+1]  / Sol->rho[n+1]) - (Sol->rhoY[n-1]  / Sol->rho[n-1])) / (2.0*dx);
                    double gimpx   = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[0]/Msq)*dthetax/theta);
                    
                    double dthetay = ((Sol->rhoY[n+icx]  / Sol->rho[n+icx]) - (Sol->rhoY[n-icx]  / Sol->rho[n-icx])) / (2.0*dy);
                    double gimpy    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetay/theta);
                    
                    
                    hplusx[n]  = theta * gimpx;
                    hplusy[n]  = theta * gimpy;
                    
                    hg[n]      = impl_grav_pr2 * thermcon * th.gamminv * pow(Sol->rhoY[n],cexp) * ud.gravity_strength[1] * gimpy * 0.5;
                    
                    /* hg[n]      = impl_grav_pr2 * thermcon * th.gamminv * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp) * ud.gravity_strength[1] * gimpy; */
                    /* hg[n]      = impl_grav_pr2 * th.gamminv * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp) * (ud.gravity_strength[1]/Msq) * gimpy; */
                    /* hg[n]      = 0.0 * thermcon * th.gamminv * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp) * ud.gravity_strength[1] * gimpy * 0.5; */
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
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int igz = elem->igz;
            const int icz = elem->icz;
            
            const double dx = elem->dx;
            const double dy = elem->dy;
            const double dz = elem->dz;
            
            const int dix = 1;
            const int diy = icx;
            const int diz = icx*icy;
            
            double Msq = ud.Msq;
            
            double* hplusx  = hplus[0];
            double* hplusy  = hplus[1];
            double* hplusz  = hplus[2];
            double* hc      = hcenter;
            
            const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt);
            const double cexp    = 1.0-th.gamm;
            
            int i, j, k, l, m, n;
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = hplusz[i] = 0.0;
            
            for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        {
                            double theta   = Sol->rhoY[n] / Sol->rho[n] ;
                            
                            double dthetax = ((Sol->rhoY[n+dix]  / Sol->rho[n+dix]) - (Sol->rhoY[n-dix]  / Sol->rho[n-dix])) / (2.0*dx);
                            double gimpx   = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[0]/Msq)*dthetax/theta);
                            
                            double dthetay = ((Sol->rhoY[n+diy]  / Sol->rho[n+diy]) - (Sol->rhoY[n-diy]  / Sol->rho[n-diy])) / (2.0*dy);
                            double gimpy    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetay/theta);
                            
                            double dthetaz = ((Sol->rhoY[n+diz]  / Sol->rho[n+diz]) - (Sol->rhoY[n-diz]  / Sol->rho[n-diz])) / (2.0*dz);
                            double gimpz    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetaz/theta);
                                                        
                            hplusx[n]  = theta * gimpx;
                            hplusy[n]  = theta * gimpy;
                            hplusz[n]  = theta * gimpz;
                        }
                    }
                }
            }
            
            for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;
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
#else
/* ========================================================================= */

static void interface_coefficients(
                                   double* h[3],
                                   const ConsVars* Sol,
                                   const ElemSpaceDiscr* elem,
                                   const NodeSpaceDiscr* node)
{
    
    
    /* let's do it right away in Michael's memory scheme */
    
    const int icxe = elem->icx;
    const int ifxn = node->ifx;
    const int ifyn = node->ify;
    
    int i, j;
    
    assert(node->ndim == 2);
    
    for (j = 0; j < elem->icy; j++ ) {
        const int me  = j * icxe;
        const int mfx = j * ifxn;
        const int mfy = j+1;
        
        for (i = 0; i < elem->icx; i++ ) {
            const int ne   = me + i;
            const int nfxm = 2*(mfx        + i+1) + 1;             /* 2*(i+1 +   j  *ifxn ) + 1 */
            const int nfxp = 2*(mfx + ifxn + i+1);                 /* 2*(i+1 + (j+1)*ifxn )     */
            const int nfym = 2*(mfy + i    *ifyn) + 1;             /* 2*(j+1 +   i  *ifyn ) + 1 */
            const int nfyp = 2*(mfy + (i+1)*ifyn);                 /* 2*(j+1 + (i+1)*ifyn )     */
            const double theta = Sol->rhoY[ne] / Sol->rho[ne];
            const double a = theta;

            h[0][nfxm] = a;
            h[0][nfxp] = a;
            h[1][nfym] = a;
            h[1][nfyp] = a;
        }
    }
}
#endif

/* ========================================================================= */

#ifdef MICHAEL_CORRECTION

/* ========================================================================= */

static void correction_michael(
                               ConsVars *Sol,
                               const ElemSpaceDiscr *elem,
                               const NodeSpaceDiscr *node,
                               const double *int_p,
                               const double t,
                               const double dt)
{
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int igx = elem->igx;
    const int igy = elem->igy;
    const int icx_inn = icx-2*igx;
    const int icy_inn = icy-2*igy;
    
    double fac   = 0.5*dt;
    
    int i, j;
    
    /* first run blindly over all cells */
    for (j=0; j<icy_inn; j++) {
        const int n_inn = j;
        const int n     = (j+igy)*icx;
        for (i=0; i<icx_inn; i++) {
            const int m_inn = n_inn + i*icy_inn;
            const int m     = n + (i+igx);

            Sol->rhou[m] -= fac*int_p[2*m_inn];
            Sol->rhov[m] -= fac*int_p[2*m_inn+1];
        }
    }
}
#else /* MICHAEL_CORRECTION */

static void correction(
                       ConsVars* Sol,
                       const ElemSpaceDiscr* elem,
                       const NodeSpaceDiscr* node,
                       double* p,
                       const double t,
                       const double dt)
{
    
    extern User_Data ud;
    
    const double M = ud.M;
    const int ndim = elem->ndim;
    
    if(M == 0.0) {
        
        switch(ndim) {
            case 1: {
                ERROR("function not available");
                break;
            }
                
            case 2: {
                
                int i, j, m, me;
                const int igx = node->igx;
                const int icx = node->icx;
                const int igy = node->igy;
                const int icy = node->icy;
                const int icxe = icx - 1;
                const double dx = node->dx;
                const double dy = node->dy;
                const double oodx = 1.0 / dx;
                const double oody = 1.0 / dy;
                const double oow = 1.0 / 2.0;
                const double dtowdx = dt * oow * oodx;
                const double dtowdy = dt * oow * oody;
                
                for(j = igy; j < icy - igy - 1; j++) {
                    m = j * icx;
                    me = j * icxe;
                    for(i = igx; i < icx - igx - 1; i++) {
                        const int n = m + i;
                        const int nicx = n + icx;
                        const int n1 = n + 1;
                        const int n1icx = n1 + icx;
                        const int ne = me + i;
                        const double Dpx = 0.5*(p[n1] - p[n] + p[n1icx] - p[nicx]);
                        const double Dpy = 0.5*(p[nicx] - p[n] + p[n1icx] - p[n1]);
                        
                        Sol->rhou[ne] += - dtowdx * Dpx;
                        Sol->rhov[ne] += - dtowdy * Dpy;
                    }
                }
                
                break;
            }
            case 3: {
                
                int i, j, k, l, m, le, me;
                const int igx = node->igx;
                const int icx = node->icx;
                const int igy = node->igy;
                const int icy = node->icy;
                const int igz = node->igz;
                const int icz = node->icz;
                const int icxicy = icx * icy;
                const int icxe = icx - 1;
                const int icye = icy - 1;
                const int icxicye = icxe * icye;
                const double dx = node->dx;
                const double dy = node->dy;
                const double dz = node->dz;
                const double oodx = 1.0 / dx;
                const double oody = 1.0 / dy;
                const double oodz = 1.0 / dz;
                const double oow = 1.0 / 4.0;
                const double dtowdx = dt * oow * oodx;
                const double dtowdy = dt * oow * oody;
                const double dtowdz = dt * oow * oodz;
                
                for(k = igz; k < icz - igz - 1; k++) {
                    l = k * icxicy;
                    le = k * icxicye;
                    for(j = igy; j < icy - igy - 1; j++) {
                        m = l + j * icx;
                        me = le + j * icxe;
                        for(i = igx; i < icx - igx - 1; i++) {
                            const int n = m + i;
                            const int nicx = n + icx;
                            const int nicxicy = n + icxicy;
                            const int nicxicxicy = nicx + icxicy;
                            const int n1 = n + 1;
                            const int n1icx = n1 + icx;
                            const int n1icxicy = n1 + icxicy;
                            const int n1icxicxicy = n1icx + icxicy;
                            const int ne = me + i;
                            
                            const double Dpx = p[n1]          - p[n]
                            + p[n1icx]       - p[nicx]
                            + p[n1icxicxicy] - p[nicxicxicy]
                            + p[n1icxicy]    - p[nicxicy];
                            
                            const double Dpy = p[nicx]        - p[n]
                            + p[nicxicxicy]  - p[nicxicy]
                            + p[n1icxicxicy] - p[n1icxicy]
                            + p[n1icx]       - p[n1];
                            
                            const double Dpz = p[nicxicy]     - p[n]
                            + p[n1icxicy]    - p[n1]
                            + p[n1icxicxicy] - p[n1icx]
                            + p[nicxicxicy]  - p[nicx];
                            
                            Sol->rhou[ne] += - dtowdx * Dpx;
                            Sol->rhov[ne] += - dtowdy * Dpy;
                            Sol->rhow[ne] += - dtowdz * Dpz;
                        }
                    }
                }
                break;
            }
            default: ERROR("ndim not in {1, 2, 3}");
        }
    }
    else {
        ERROR("function not available");
    }
}
#endif /* MICHEAL_CORRECTION */

/* ========================================================================= */
#ifdef FIX_PERIODIC_BCS_2ND_PROJ

void fix_periodic_boundaries_2nd_projection_inn(
                                                double *rhs,
                                                double *h[3],
                                                const NodeSpaceDiscr* node)
{
    /* note that we are in Michael's notation, i.e., the position
     in the rhs[] linear array reads (j + i*icy)
     */
    
    extern User_Data ud;
    
    int j;
    int igx = node->igx;
    int igy = node->igy;
    int icx = node->icx - 2*igx;
    int icy = node->icy - 2*igy;
    
    if (ud.bdrytype_min[0] != PERIODIC) {
        return;
        
    } else {
        
        for (j = 0; j < icy; j++) {
            int icl = j;
            int icr = j + (icx-1)*icy;
            rhs[icr] *= 1.0;
            rhs[icl] *= 1.0;
        }
    }
}
#endif


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void periodic(
              double test[5],
              double *rhs,
              double *h0,
              double *h1,
              const NodeSpaceDiscr *node)
{
    /*
     checks periodicity of rhs and coefficient arrays for the
     second projection Poisson problem in Michaels memory
     arrangement ( j first )
     */
    
    const int ifx_inn = node->ifx - 2*node->igx;
    const int icx_inn = node->icx - 2*node->igx;
    const int ify_inn = node->ify - 2*node->igy;
    const int icy_inn = node->icy - 2*node->igy;
    
    double test_rhs = 0.0;
    double test_h00 = 0.0;
    double test_h01 = 0.0;
    double test_h10 = 0.0;
    double test_h11 = 0.0;
    
    int j;
    
    /* periodicity check for rhs living on nodes */
    for (j=0; j<icy_inn; j++) {
        int nl = j;
        int nr = j + (icx_inn-1)*icy_inn;
        test_rhs = MAX_own(test_rhs, fabs(rhs[nl] - rhs[nr]));
    }
    
    /* periodicity check for coefficients h00, h01 living on dual cell edges facing in the x-direction */
    for (j=0; j<icy_inn; j++) {
        int nl00 = 2*j;
        int nl01 = 2*j + 1;
        int nr00 = nl00 + 2*(ifx_inn-2)*icy_inn;
        int nr01 = nl01 + 2*(ifx_inn-2)*icy_inn;
        
        int nl10 = 2*(j + icy_inn);
        int nl11 = 2*(j + icy_inn) + 1;
        int nr10 = nl10 + 2*(ifx_inn-2)*icy_inn;
        int nr11 = nl11 + 2*(ifx_inn-2)*icy_inn;
        
        test_h00 = MAX_own(test_h00, fabs(h0[nl00] - h0[nr00]));
        test_h00 = MAX_own(test_h00, fabs(h0[nl01] - h0[nr01]));
        
        test_h01 = MAX_own(test_h01, fabs(h0[nl10] - h0[nr10]));
        test_h01 = MAX_own(test_h01, fabs(h0[nl11] - h0[nr11]));
    }
    
    /* periodicity check for coefficients h10, h11 living on dual cell edges facing in the y-direction */
    for (j=0; j<ify_inn; j++) {
        int nl00 = 2*j;
        int nl01 = 2*j + 1;
        int nr00 = nl00 + 2*(icx_inn-1)*ify_inn;
        int nr01 = nl01 + 2*(icx_inn-1)*ify_inn;
        
        test_h10 = MAX_own(test_h10, fabs(h1[nl00] - h1[nr00]));
        test_h10 = MAX_own(test_h10, fabs(h1[nl01] - h1[nr01]));
    }
    
    test[0] = test_rhs;
    test[1] = test_h00;
    test[2] = test_h01;
    test[3] = test_h10;
    test[4] = test_h11;
    
}

#endif /*SOLVER_2_HYPRE_RUPE */

#else /* SOLVER_2_HYPRE */

static double divergence_nodes(
							 double* rhs,
							 const ElemSpaceDiscr* elem,
							 const NodeSpaceDiscr* node,
							 const ConsVars* Sol,
                             double* eta,
							 const MPV* mpv,
							 const BDRY* bdry,
							 const double dt,
							 const double weight);


#ifdef P2_EXACT_PROJECTION_VK
void evaluate_eta(double* eta_hyp0,
                  const ConsVars* Sol,
                  const ElemSpaceDiscr* elem);

void eta_second_correction(MPV* mpv,
                           const double* hplus[3],
                           const ElemSpaceDiscr* elem,
                           const NodeSpaceDiscr *node,
                           const double dt);
#endif

#ifdef SECOND_CORRECTION_IS_A_PROJECTION
static void controlled_variable_change_nodes(
                                             double* rhs, 
                                             const ElemSpaceDiscr* elem, 
                                             const NodeSpaceDiscr* node, 
                                             const ConsVars* Sol, 
                                             const ConsVars* Sol0, 
                                             const double dt);
#endif

static enum Constraint integral_condition_nodes(
												double* rhs,
												const NodeSpaceDiscr* node,
												const int x_periodic,
												const int y_periodic,
												const int z_periodic);

static void	catch_periodic_directions(
									  double* rhs,  
									  const NodeSpaceDiscr* node, 
									  const ElemSpaceDiscr* elem,
									  const int x_periodic,
									  const int y_periodic,
									  const int z_periodic);

static void operator_coefficients_nodes(
										double* Splus[3], 
										double* wcenter,
										double* wgrav,
										const ElemSpaceDiscr* elem,
										const NodeSpaceDiscr* node,
										ConsVars* Sol,
										const ConsVars* Sol0,
										const MPV* mpv,
                                        const double dt);

void correction_nodes(
					  ConsVars* Sol,
					  const ElemSpaceDiscr* elem,
					  const NodeSpaceDiscr* node,
                      const double* hplus[3], 
                      const double* hgrav,
					  double* p2, 
					  const double t,
					  const double dt);

/* ========================================================================== */

void second_projection(
                       ConsVars* Sol,
                       MPV* mpv,
                       const ConsVars* Sol0,
                       const ElemSpaceDiscr* elem,
                       const NodeSpaceDiscr* node,
                       const double p_update,
                       const double t,
                       const double dt) 
{
	
	extern User_Data ud;
	extern BDRY* bdry;
    
#ifdef P2_EXACT_PROJECTION_VK
    static int first_call = 1;
#endif
    
	const int nc = node->nc;
	
	double** hplus   = mpv->Level[0]->wplus;
	double*  hcenter = mpv->Level[0]->wcenter;
	double*  hgrav   = mpv->Level[0]->wgrav;
	
	double* rhs      = mpv->Level[0]->rhs;
	double* p2       = mpv->Level[0]->p;
	
    double rhs_max;
    
	int x_periodic, y_periodic, z_periodic;
	int ii;
	
    
    printf("\n\n====================================================");
    printf("\nSecond Projection");
    printf("\n====================================================\n");

	/* switch for stopping loops in solver for periodic data before
	 last grid line in each direction */
	x_periodic = 0;
	y_periodic = 0;
	z_periodic = 0;
	if(ud.bdrytype_min[0] == PERIODIC) x_periodic = 1;
	if(ud.bdrytype_min[1] == PERIODIC) y_periodic = 1;
	if(ud.bdrytype_min[2] == PERIODIC) z_periodic = 1;
	
    /* KEEP_OLD_POISSON_SOLUTIONS */
	for(ii=0; ii<nc; ii++){
		p2[ii] = mpv->dp2_nodes[ii];
		rhs[ii] = 0.0;
	}
	
#ifdef P2_EXACT_PROJECTION_VK
    /* eta-update or start from hyperbolic predictor */
    if (first_call) {
        memcpy(mpv->eta_hyp0, mpv->eta, elem->nc*sizeof(double));
        evaluate_eta(mpv->eta_hyp, Sol, elem);
        for (int ic=0; ic<elem->nc; ic++) {
            mpv->eta[ic] += mpv->eta_hyp[nc] - mpv->eta_hyp0[nc];
        }
        first_call = HYP_UPDATE_ETA;
    }
#endif
    
    double rhs_weight_new = 2.0;
    double rhs_weight_old = 2.0*ud.compressibility;
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv->eta, mpv, bdry, dt, rhs_weight_new);
    if (ud.compressibility) {
        divergence_nodes(rhs, elem, node, Sol0, mpv->eta0, mpv, bdry, dt, rhs_weight_old); /*  for psinc, this can be commented out */
    }
    
#ifdef SECOND_CORRECTION_IS_A_PROJECTION
    controlled_variable_change_nodes(rhs, elem, node, Sol, Sol0, dt);
#endif
    
	catch_periodic_directions(rhs, node, elem, x_periodic, y_periodic, z_periodic);

    /* test */
#if 1
    FILE *prhsfile = NULL;
    char fn[120], fieldname[90];
    sprintf(fn, "%s/rhs_nodes/rhs_nodes_000.hdf", ud.file_name);
    sprintf(fieldname, "rhs-nodes");
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
#endif
     
    
	assert(integral_condition_nodes(rhs, node, x_periodic, y_periodic, z_periodic) != VIOLATED); 
	operator_coefficients_nodes(hplus, hcenter, hgrav, elem, node, Sol, Sol0, mpv, dt);
	
	variable_coefficient_poisson_nodes(p2, (const double **)hplus, hcenter, hgrav, rhs, x_periodic, y_periodic, z_periodic, dt);
    
#if 1
    rhs_max = 0.0;
    double rhs_min = 100000.0;
    double hx_max  = 0.0;
    double hy_max  = 0.0;
    double p2_max  = 0.0;
    double hx_min  = 100000.0;
    double hy_min  = 100000.0;
    double p2_min  = 100000.0;
    for (int in = 0; in<node->nc; in++) {
        rhs_max = MAX_own(rhs_max, fabs(rhs[in]));
        rhs_min = MIN_own(rhs_min, fabs(rhs[in]));
    }
    for (int jn = elem->igy; jn < elem->icy-elem->igy; jn++) {
        for (int in = elem->igx; in < elem->icx-elem->igx; in++) {
            int nn = jn*elem->icx + in;
            hx_max = MAX_own(hx_max, fabs(hplus[0][nn]));
            hx_min = MIN_own(hx_min, fabs(hplus[0][nn]));
            hy_max = MAX_own(hy_max, fabs(hplus[1][nn]));
            hy_min = MIN_own(hy_min, fabs(hplus[1][nn]));
        }
    }
    for (int jn = node->igy; jn < node->icy-node->igy; jn++) {
        for (int in = node->igx; in < node->icx-node->igx; in++) {
            int nn = jn*node->icx + in;
            p2_max = MAX_own(p2_max, p2[nn]);
            p2_min = MIN_own(p2_min, p2[nn]);
        }
    }
    printf("rhs_max, rhs_min, h1_max, h1_min, h2_max, h2_max, p2_min, p2_max = %e, %e, %e, %e, %e, %e, %e, %e\n", rhs_max, rhs_min, hx_max, hx_min, hy_max, hy_min, p2_min, p2_max);

    sprintf(fn, "%s/rhs_nodes/rhs_nodes_001.hdf", ud.file_name);
    sprintf(fieldname, "rhs-nodes");
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);

#endif
        
#if 1
    extern double *W0;
    int imax, jmax;
    rhs_max = 0.0;
    for (int jn=node->igy; jn<node->icy-node->igy; jn++) {
        int mn = jn*node->icx;
        for (int in=node->igx; in<node->icx-node->igx; in++) {
            int nn = mn+in;
            if (rhs_max < fabs(rhs[nn])) {
                rhs_max = fabs(rhs[nn]);
                imax = in;
                jmax = jn;
            } 
        }
    }
    printf("rhs_max before = %e\n", rhs_max);
    correction_nodes(Sol, elem, node, (const double **)hplus, hgrav, p2, t, dt);
    EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, p2, (const double **)hplus, hcenter, hgrav, x_periodic, y_periodic, z_periodic, W0);
    precon_apply(W0,W0,node);
    for(ii=0; ii<nc; ii++) rhs[ii] = 0.0;
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv->eta, mpv, bdry, dt, rhs_weight_new);
    catch_periodic_directions(rhs, node, elem, x_periodic, y_periodic, z_periodic);

#if 1
    FILE *prhs2file = NULL;
    char fn2[120], fieldname2[90];
    sprintf(fn2, "%s/rhs_nodes/rhs_nodes_002.hdf", ud.file_name);
    sprintf(fieldname2, "rhs-nodes-post");
    WriteHDF(prhs2file, node->icx, node->icy, node->icz, node->ndim, rhs, fn2, fieldname2);
    sprintf(fn2, "%s/rhs_nodes/rhs_nodes_003.hdf", ud.file_name);
    sprintf(fieldname2, "lap-final");
    WriteHDF(prhs2file, node->icx, node->icy, node->icz, node->ndim, W0, fn2, fieldname2);
#endif

    rhs_max = 0.0;
    for (int jn=node->igy+1; jn<node->icy-node->igy-1; jn++) {
        int mn = jn*node->icx;
        for (int in=node->igx+1; in<node->icx-node->igx-1; in++) {
            int nn = mn+in;
            if (rhs_max < fabs(rhs[nn])) {
                rhs_max = fabs(rhs[nn]);
                imax = in;
                jmax = jn;
            } 
        }
    }
    printf("rhs_max after  = %e\n", rhs_max);
#else
    correction_nodes(Sol, elem, node, (const double**)hplus, hgrav, p2, t, dt);
#endif
    
    for(ii=0; ii<nc; ii++) {
        /* hunt down factor of 0.5 in the next line !! */
        mpv->p2_nodes[ii]  = SCND_PROJ_OLDP_WEIGHT*mpv->p2_nodes[ii] + p_update*SCND_PROJ_DELP_WEIGHT*0.5*p2[ii];

#ifdef NODAL_PROJECTION_ONLY
        mpv->dp2_nodes[ii] = 0.5*p2[ii];
#else /* NODAL_PROJECTION_ONLY */ 
        mpv->dp2_nodes[ii] = p2[ii];
#endif /* NODAL_PROJECTION_ONLY */ 
    }

#ifdef NODAL_PROJECTION_ONLY
    for(ii=0; ii<nc; ii++) {
        /* hunt down factor of 0.5 in the next line !! */
        mpv->p2_nodes[ii]  = SCND_PROJ_OLDP_WEIGHT*mpv->p2_nodes[ii] + p_update*SCND_PROJ_DELP_WEIGHT*0.5*p2[ii];        
        mpv->dp2_nodes[ii] = 0.5*p2[ii];
    }
    Pressure_node_to_elem(mpv, (const ConsVars*)Sol, (const ConsVars*)Sol0, elem, node);
#else /* NODAL_PROJECTION_ONLY */
    for(ii=0; ii<nc; ii++) {
        /* hunt down factor of 0.5 in the next line !! */
        mpv->p2_nodes[ii]  = SCND_PROJ_OLDP_WEIGHT*mpv->p2_nodes[ii] + p_update*SCND_PROJ_DELP_WEIGHT*0.5*p2[ii];        
        mpv->dp2_nodes[ii] = p2[ii];
    }
#endif /* NODAL_PROJECTION_ONLY */ 
   
#ifndef NO_BDRYCONDS_PROJ2
    Bound_p_nodes(mpv, Sol, elem, node, 1);
#endif
    
#ifdef P2_EXACT_PROJECTION_VK
    eta_second_correction(mpv, hplus, elem, node, dt);
#endif
            
}

/* ========================================================================== */

static double divergence_nodes(
							 double* rhs,
							 const ElemSpaceDiscr* elem,
							 const NodeSpaceDiscr* node,
							 const ConsVars* Sol,
                             double* eta,
							 const MPV* mpv,
							 const BDRY* bdry,
							 const double dt,
							 const double weight) {
	
	extern User_Data ud;
	
	const int ndim = node->ndim;

    double div_max = 0.0;
    
    int i, j, k, mn;
    
	switch(ndim) {
		case 1: {
			ERROR("divergence_nodes() not implemented for 1D\n");
			break;
		}
		case 2: {
				
            /*
			const int igxn = node->igx;
			const int igyn = node->igy;
			const int icyn = node->icy;
             */
			const int icxn = node->icx;
			
			const int igxe = elem->igx;
			const int icxe = elem->icx;
			const int igye = elem->igy;
			const int icye = elem->icy;
			
			const double dx = node->dx;
			const double dy = node->dy;
			const double oodxdt = 1.0 / (dx * dt);
			const double oodydt = 1.0 / (dy * dt);
			const double oow = 1.0 / 2.0;
            const double todt = 2.0 / dt;
			const double oowdxdt = weight * oow * oodxdt;
			const double oowdydt = weight * oow * oodydt;
			
			double Y;
			
			/* predicted time level divergence via scattering */
			for(j = igye; j < icye - igye; j++) {
				const int me = j * icxe;
				const int mn = j * icxn; 
				for(i = igxe; i < icxe - igxe; i++) {
					const int n     = mn + i;
					const int nicx  = n  + icxn;
					const int n1    = n  + 1;
					const int n1icx = n1 + icxn;
					
					const int ne = me + i; 
					double tmpfx, tmpfy, tmpfe;
					
					Y = Sol->rhoY[ne] / Sol->rho[ne];
					tmpfx = oowdxdt * Y * Sol->rhou[ne];
					tmpfy = oowdydt * Y * Sol->rhov[ne];
                    tmpfe = todt * eta[ne];
					
					rhs[n]           += + tmpfx + tmpfy - tmpfe;
					rhs[n1]          += - tmpfx + tmpfy + tmpfe;
					rhs[n1icx]       += - tmpfx - tmpfy - tmpfe;
					rhs[nicx]        += + tmpfx - tmpfy + tmpfe;
				}
			} 
			
            
			/* account for influx bottom boundary */
			j = igye;
			mn = j * icxn; 
			Y  = mpv->HydroState->Y0[j];
			for(i = igxe; i < icxe - igxe; i++) {
				const int n     = mn + i;
				const int n1    = n  + 1;
				
				double rhov_wall = bdry->wall_massflux[i]; 
				double tmpy = oowdydt * Y * rhov_wall;
				
				rhs[n]  += - tmpy;
				rhs[n1] += - tmpy;
			}
						
			break;
		}
		case 3: {
            
            /*
			const int igxn = node->igx;
			const int igyn = node->igy;
			const int igzn = node->igz;
			const int iczn = node->icz;
             */
			const int icxn = node->icx;
			const int icyn = node->icy;
			
			const int igxe = elem->igx;
			const int icxe = elem->icx;
			const int igye = elem->igy;
			const int icye = elem->icy;
			const int igze = elem->igz;
			const int icze = elem->icz;
			
			const double dx = node->dx;
			const double dy = node->dy;
			const double dz = node->dz;

			const double oodxdt = 1.0 / (dx * dt);
			const double oodydt = 1.0 / (dy * dt);
			const double oodzdt = 1.0 / (dz * dt);
            
			const double oow = 1.0 / 4.0;
			const double oowdxdt = weight * oow * oodxdt;
			const double oowdydt = weight * oow * oodydt;
			const double oowdzdt = weight * oow * oodzdt;
            
            const int dixn = 1;
            const int diyn = icxn;
            const int dizn = icxn*icyn;
			
			double Y;
			
			/* predicted time level divergence via scattering */
			for(k = igze; k < icze - igze; k++) {
				const int le = k * icye * icxe;
				const int ln = k * icyn * icxn; 
                for(j = igye; j < icye - igye; j++) {
                    const int me = le + j * icxe;
                    const int mn = ln + j * icxn; 
                    for(i = igxe; i < icxe - igxe; i++) {
                        const int ne = me + i; 
                        const int nn000  = mn + i;
                        const int nn010  = nn000 + diyn;
                        const int nn011  = nn000 + diyn + dizn;
                        const int nn001  = nn000        + dizn;
                        const int nn100  = nn000 + dixn;
                        const int nn110  = nn000 + dixn + diyn;
                        const int nn111  = nn000 + dixn + diyn + dizn;
                        const int nn101  = nn000 + dixn        + dizn;
                        
                        double tmpfx, tmpfy, tmpfz;
                        
                        Y = Sol->rhoY[ne] / Sol->rho[ne];
                        tmpfx = oowdxdt * Y * Sol->rhou[ne]; /* (rhou*Y) * 0.5 * / (dx*dt) */
                        tmpfy = oowdydt * Y * Sol->rhov[ne];
                        tmpfz = oowdzdt * Y * Sol->rhow[ne];
                        
                        rhs[nn000]       += + tmpfx + tmpfy + tmpfz;
                        rhs[nn010]       += + tmpfx - tmpfy + tmpfz;
                        rhs[nn011]       += + tmpfx - tmpfy - tmpfz;
                        rhs[nn001]       += + tmpfx + tmpfy - tmpfz;
                        rhs[nn100]       += - tmpfx + tmpfy + tmpfz;
                        rhs[nn110]       += - tmpfx - tmpfy + tmpfz;
                        rhs[nn111]       += - tmpfx - tmpfy - tmpfz;
                        rhs[nn101]       += - tmpfx + tmpfy - tmpfz;
                    }
                } 
			} 
			
			/* account for influx bottom boundary */
			j = igye;
			Y  = mpv->HydroState->Y0[j];
			for(k = igze; k < icze - igze; k++) {
                int ln = k * icyn*icxn;
                int mn = ln + j*icxn;
                for(i = igxe; i < icxe - igxe; i++) {
                    int nn   = mn + i;
                    int nn00 = nn;
                    int nn10 = nn + dixn;
                    int nn11 = nn + dixn + dizn;
                    int nn01 = nn +      + dizn;
                    
                    double rhov_wall = bdry->wall_massflux[i]; 
                    double tmpy = oowdydt * Y * rhov_wall;
                    
                    rhs[nn00] += - tmpy;
                    rhs[nn10] += - tmpy;
                    rhs[nn01] += - tmpy;
                    rhs[nn11] += - tmpy;
                }
            }
			
			/* set ghost values of rhs to zero 
			j = igyn - 1;
			l = j * icxn;
			for(i = 0; i < icxn; i++) 
				rhs[l + i] = 0.0;
			
			j = icyn - igyn;
			l = j * icxn;
			for(i = 0; i < icxn; i++)
				rhs[l + i] = 0.0;
			
			i = igxn - 1;
			l = i;
			for(j = 0; j < icyn; j++)
				rhs[l + j * icxn] = 0.0;
			
			i = icxn - igxn;
			l = i;
			for(j = 0; j < icyn; j++) 
				rhs[l + j * icxn] = 0.0;
			*/
			break;
		}
		default: ERROR("ndim not in {1, 2, 3}");
	}
    
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "rhs_nodes.hdf");
    sprintf(fieldname, "rhs-nodes");    
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
#endif
    
    for (int in = 0; in<node->nc; in++) {
        div_max = MAX_own(div_max, fabs(rhs[in]));
    }
    return div_max;
    
}

/* ========================================================================== */

#ifdef P2_EXACT_PROJECTION_VK
void evaluate_eta(double* eta,
                  const ConsVars* Sol,
                  const ElemSpaceDiscr* elem)
{
    assert(elem->ndim == 2);
    
    const int icx = elem->icx;
    const int icy = elem->icz;

    const int igx = elem->igx;
    const int igy = elem->igy;

    const double dx = elem->dx;
    const double dy = elem->dy;
    
    for (int j=igy; j<icy-igy; j++) {int m = j*icx;
        for (int i=igx; i<icx-igx; i++) {int n = m+i;
            int nn = n+icx;
            int ns = n-icx;
            int ne = n+1;
            int nw = n-1;
            eta[n] =  (dy/dx) * (Sol->rhou[nn]*Sol->rhoY[nn]/Sol->rho[nn] - Sol->rhou[ns]*Sol->rhoY[ns]/Sol->rho[ns]) / dy;
            eta[n] += (dx/dy) * (Sol->rhov[ne]*Sol->rhoY[ne]/Sol->rho[ne] - Sol->rhov[nw]*Sol->rhoY[nw]/Sol->rho[nw]) / dx;
            eta[n] *= 0.125;
        }
    }
    
}

void eta_second_correction(MPV* mpv,
                           const double* hplus[3],
                           const ElemSpaceDiscr* elem,
                           const NodeSpaceDiscr* node,
                           const double dt)
{
    
    assert(elem->ndim == 2);
    
    const double* hplusx   = hplus[0];
    const double* hplusy   = hplus[1];

    const int icx = elem->icx;
    const int icy = elem->icy;

    const int inx = node->icx;

    const int igx = elem->igx;
    const int igy = elem->igy;
    
    const double dx = elem->dx;
    const double dy = elem->dy;
    
    const double oodxsq = 1.0/(dx*dx);
    const double oodysq = 1.0/(dy*dy);
    
    double *eta0 = mpv->eta0;
    double *eta  = mpv->eta;
    double *p    = mpv->dp2_nodes;
    
    for (int j=igy; j<icy; j++) {
        int mc = j*icx;
        int mn = j*inx;
        for (int i=igx; i<icx; i++) {
            int nc  = mc+i;
            int nn  = mn+i;
            
            int nn1    = nn+1;
            int nninx  = nn+inx;
            int nn1inx = nn+1+inx;
            eta0[nc] = eta[nc];
            /* eta[nc]  -= 0.5*dt*0.125*(dy/dx+dx/dy)*(p[nn1inx] - p[nn1] - p[nninx] + p[nn]); */
            eta[nc]  -= 0.5*dt*0.125*(hplusx[nc]*oodxsq+hplusy[nc]*oodysq)*(p[nn1inx] - p[nn1] - p[nninx] + p[nn]);
        }
    }
}

#endif

/* ========================================================================== */

#ifdef SECOND_CORRECTION_IS_A_PROJECTION
static void controlled_variable_change_nodes(
                                             double* rhs, 
                                             const ElemSpaceDiscr* elem, 
                                             const NodeSpaceDiscr* node, 
                                             const ConsVars* Sol, 
                                             const ConsVars* Sol0, 
                                             const double dt)
{
	extern User_Data ud;
	
	const int ndim = elem->ndim;
	int i, j;
    
	switch(ndim) {
		case 1: {
			ERROR("divergence_nodes() not implemented for 1D\n");
			break;
		}
		case 2: {
            
			const int icxn = node->icx;
			
			const int igxe = elem->igx;
			const int icxe = elem->icx;
			const int igye = elem->igy;
			const int icye = elem->icy;
			
			const double factor = 0.25 * ud.compressibility * 4.0 / (dt * dt);
						
			/* predicted time level divergence via scattering */
			for(j = igye; j < icye - igye; j++) {
				const int me = j * icxe;
				const int mn = j * icxn; 
				for(i = igxe; i < icxe - igxe; i++) {
					const int n     = mn + i;
					const int nicx  = n  + icxn;
					const int n1    = n  + 1;
					const int n1icx = n1 + icxn;
					const int ne = me + i; 
                    					
					double d = factor*(Sol->rhoY[ne] - Sol0->rhoY[ne]);
					
					rhs[n]     += d;
					rhs[n1]    += d;
					rhs[n1icx] += d;
					rhs[nicx]  += d;
				}
			} 
			            
			break;
		}
		case 3: {
            
			ERROR("controlled_variable_change_nodes() not implemented for 3D\n");
						
			break;
		}
		default: ERROR("ndim not in {1, 2, 3}");
	}
    
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "rhs_nodes.hdf");
    sprintf(fieldname, "rhs-nodes");    
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs, fn, fieldname);
#endif
    
}
#endif

/* ========================================================================== */

static enum Constraint integral_condition_nodes(
												double* rhs,
												const NodeSpaceDiscr* node,
												const int x_periodic,
												const int y_periodic,
												const int z_periodic) {
	
	/* Rupert, this is valid only for periodic data so far!!! */
	
	int i, j, k, l, m, n;
	
	const int igx = node->igx;
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int igz = node->igz;
	const int icz = node->icz;
	const int icxicy = icx * icy;
    
    const double del = DBL_EPSILON * node->nc;
    
	double tmp = 0.0, rhs_max=0.0;
	
	for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icxicy;
		for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx; 
			for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
				tmp += rhs[n];
                if(fabs(rhs[n])>rhs_max) rhs_max = fabs(rhs[n]);
			}
		}
	}
    
    /* if(fabs(tmp) > 100*del) { */
    if(fabs(tmp) > 10000*del*rhs_max) {
        printf("CHEATING:  tmp = %e,  rhs_max = %e, (tmp/rhs_max)/del = %e\n", tmp/rhs_max, rhs_max, (tmp/rhs_max)/del);
        /* return VIOLATED; */
        return VIOLATED;
    }
    else if(fabs(tmp) > 100*del*rhs_max) {
        printf("integral_condition_nodes = barely OK; tmp = %e,  rhs_max = %e, (tmp/rhs_max)/del = %e\n", tmp/rhs_max, rhs_max, (tmp/rhs_max)/del);
        return SATISFIED;
    }
    else {
        /* printf("integral_condition_nodes = OK;    tmp = %e,  rhs_max = %e\n", tmp, rhs_max); */
        printf("integral_condition_nodes = OK; tmp = %e\n", tmp);
        return SATISFIED;
    }
}

/* ========================================================================== */

static void	catch_periodic_directions(
									  double* rhs,  
									  const NodeSpaceDiscr* node, 
									  const ElemSpaceDiscr* elem,
									  const int x_periodic,
									  const int y_periodic,
									  const int z_periodic) {
	
	/* in any of the periodicity directions, gather boundary values of
	 rhs on first nodes, set last nodes to zero
	 */
	
	int i, j, k, l, m, l0, l1, m0, m1, n0, n1;
    
	const int igx = node->igx;
	const int icx = node->icx;
	const int igy = node->igy;
	const int icy = node->icy;
	const int igz = node->igz;
	const int icz = node->icz;
	const int icxicy = icx * icy;
	
	if(x_periodic == 1){
		for(k = igz; k < icz - igz; k++) {l = k * icxicy; 
			for(j = igy; j < icy - igy; j++) {m = l + j * icx; 
				n0 = m + igx;
				n1 = m + icx-igx-1;
				rhs[n0] += rhs[n1];
				rhs[n1] = 0.0;
			}
		}
	}
	
	if(y_periodic == 1){
		for(k = igz; k < icz - igz; k++) {l = k * icxicy; 
			m0 = l + igy * icx;
			m1 = l + (icy-igy-1) * icx;
			for(i = igx; i < icx - igx; i++) {
				n0 = m0 + i;
				n1 = m1 + i;
				rhs[n0] += rhs[n1];
				rhs[n1] = 0.0;
			}
		}
	}
	
	if(z_periodic == 1){
		l0 = igz * icxicy;
		l1 = (icz-igz-1) * icxicy; 
		for(j = igy; j < icy - igy; j++) {
			m0 = l0 + j * icx; 
			m1 = l1 + j * icx; 
			for(i = igx; i < icx - igx; i++) {
				n0 = m0 + i;
				n1 = m1 + i;
				rhs[n0] += rhs[n1];
				rhs[n1] = 0.0;
			}
		}
	}
}

/* ========================================================================== */

static void operator_coefficients_nodes(
										double* hplus[3], 
										double* hcenter,
										double* hgrav,
										const ElemSpaceDiscr* elem,
										const NodeSpaceDiscr* node,
										ConsVars* Sol,
										const ConsVars* Sol0,
										const MPV* mpv,
                                        const double dt) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	
	const int ndim = node->ndim;
    
    const int impl_grav_th2 = ud.implicit_gravity_theta2;
    const int impl_grav_pr2 = ud.implicit_gravity_press2;
        
    /* #define CHECK_MINMAX */ 
#ifdef CHECK_MINMAX
    double coeffmin = 1000000.0;
    double coeffmax = -1000000.0;
#endif
    
	switch(ndim) {
		case 1: {    
			ERROR("interface_enthalpy_nodes() not implemented for 1D\n");
			break;
		}
		case 2: {			
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int igy = elem->igy;
			const int icy = elem->icy;

            const double dx = elem->dx;
            const double dy = elem->dy;

            double Msq = ud.Msq;
			
			double* hplusx  = hplus[0];
			double* hplusy  = hplus[1];
			double* hc      = hcenter;
			double* hg      = hgrav;

#ifdef SECOND_CORRECTION_IS_A_PROJECTION
            assert(ud.p_flux_correction);
			const double ccenter = 0.0;
#else
			const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt);
#endif
			const double cexp    = 1.0-th.gamm;
            int i, j, m, n;
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = 0.0;
						            
			for(j = igy; j < icy - igy; j++) {
                m = j * icx;
				
                for(i = igx; i < icx - igx; i++) {
                    n = m + i;     
                    
                    double theta   = mpv->HydroState->Y0[j];
                    
                    /*
                     double dthetax = ((Sol->rhoY[n+1]  / Sol->rho[n+1]) - (Sol->rhoY[n-1]  / Sol->rho[n-1])) / (2.0*dx);
                     */
                    double dthetax = (mpv->HydroState->Y0[j] - mpv->HydroState->Y0[j]) / (2.0*dx);
                    double gimpx   = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[0]/Msq)*dthetax/theta);
                    
                    /*
                    double dthetay = ((Sol->rhoY[n+icx]  / Sol->rho[n+icx]) - (Sol->rhoY[n-icx]  / Sol->rho[n-icx])) / (2.0*dy);
                     */
                    double dthetay = (mpv->HydroState->Y0[j+1] - mpv->HydroState->Y0[j-1]) / (2.0*dy);

                    double gimpy    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetay/theta);
                                        
                    hplusx[n]  = theta * gimpx;
                    hplusy[n]  = theta * gimpy;
                    /* hg[n]      = impl_grav_pr2 * th.gamminv * pow(Sol->rhoY[n],cexp) * ud.gravity_strength[1] * gimpy * 0.5; */
                    hg[n]      = impl_grav_pr2 * th.gamminv * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp) * ud.gravity_strength[1] * gimpy * 0.5;
                    /* hg[n]      = impl_grav_pr2 * th.gamminv * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp) * (ud.gravity_strength[1]/Msq) * gimpy; */
                    /* hg[n]      = 0.0 * th.gamminv * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp) * ud.gravity_strength[1] * gimpy; */
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
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int igz = elem->igz;
			const int icz = elem->icz;
            
            const double dx = elem->dx;
            const double dy = elem->dy;
            const double dz = elem->dz;

            const int dix = 1;
			const int diy = icx;
			const int diz = icx*icy;

            double Msq = ud.Msq;
			
			double* hplusx  = hplus[0];
			double* hplusy  = hplus[1];
			double* hplusz  = hplus[2];
			double* hc      = hcenter;
            
			const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt);
            
			const double cexp    = 1.0-th.gamm;
			int i, j, k, l, m, n;
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = hplusz[i] = 0.0;
            
			for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        {
                            double theta   = Sol->rhoY[n] / Sol->rho[n] ;
                            
                            double dthetax = ((Sol->rhoY[n+dix]  / Sol->rho[n+dix]) - (Sol->rhoY[n-dix]  / Sol->rho[n-dix])) / (2.0*dx);
                            double gimpx   = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[0]/Msq)*dthetax/theta);
                            
                            double dthetay = ((Sol->rhoY[n+diy]  / Sol->rho[n+diy]) - (Sol->rhoY[n-diy]  / Sol->rho[n-diy])) / (2.0*dy);
                            double gimpy    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetay/theta);

                            double dthetaz = ((Sol->rhoY[n+diz]  / Sol->rho[n+diz]) - (Sol->rhoY[n-diz]  / Sol->rho[n-diz])) / (2.0*dz);
                            double gimpz    = 1.0 / (1.0 + impl_grav_th2*0.25*dt*dt*(ud.gravity_strength[1]/Msq)*dthetaz/theta);

                            hplusx[n]  = theta * gimpx;
                            hplusy[n]  = theta * gimpy;
                            hplusz[n]  = theta * gimpz;
                        }
                    }
                }
			}
            
			for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;
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


#if 0
/* ========================================================================== */

static void diagonal_preconditioner(
									double* diaginv, 
									const ElemSpaceDiscr* elem,
									const NodeSpaceDiscr* node,
									ConsVars* Sol) {
	
	extern User_Data ud;
	
	const int ndim = node->ndim;
	
	switch(ndim) {
		case 1: {    
			ERROR("interface_enthalpy_nodes() not implemented for 1D\n");
			break;
		}
		case 2: {
			
			int n;
			
#ifdef PRECON
			int i, j, m;
			const int igxn = node->igx;
			const int icxn = node->icx;
			const int igyn = node->igy;
			const int icyn = node->icy;
			
			const int icxe = elem->icx;
			
			for(j = igyn; j < icyn - igyn; j++) {m = j * icxn;
				for(i = igxn; i < icxn - igxn; i++) {n = m + i;
					{
						int nne, nnw, nse, nsw;
						
						nne  =   i   +   j   * icxe;
						nnw  = (i-1) +   j   * icxe;
						nsw  = (i-1) + (j-1) * icxe;
						nse  =   i   + (j-1) * icxe;
						
						/* diaginv[n] = 4.0 / (Sol->rhoY[nne] + Sol->rhoY[nnw] + Sol->rhoY[nsw] + Sol->rhoY[nse]); */

                        diaginv[n] = 4.0 / (  Sol->rhoY[nne] / Sol->rho[nne]
											+ Sol->rhoY[nnw] / Sol->rho[nnw] 
											+ Sol->rhoY[nsw] / Sol->rho[nsw]
											+ Sol->rhoY[nse] / Sol->rho[nse]);
					}
				}
			}
#else
			for(n=0; n<node->nc; n++) {
				diaginv[n] = 1.0;
			}	
#endif
			
			break;
		}
		case 3: {
						
#ifdef PRECON
			int i, j, k, l, m, n;
			const int igxn = node->igx;
			const int icxn = node->icx;
			const int igyn = node->igy;
			const int icyn = node->icy;
			const int igzn = node->igz;
			const int iczn = node->icz;
			
			const int icxe = elem->icx;
			const int icye = elem->icy;
            
            const int dixe = 1;
            const int diye = icxe;
            const int dize = icxe*icye;
			
			for(k = igzn; k < iczn - igzn; k++) {l = k * icxn*icyn;
                for(j = igyn; j < icyn - igyn; j++) {m = l + j * icxn;
                    for(i = igxn; i < icxn - igxn; i++) {n = m + i;
						int nc     = (k-1)*dize + (j-1)*diye + (i-1)*dixe;
						int nc000  = nc;
						int nc010  = nc        + diye;
						int nc011  = nc        + diye + dize;
						int nc001  = nc               + dize;
						int nc100  = nc + dixe;
						int nc110  = nc + dixe + diye;
						int nc111  = nc + dixe + diye + dize;
						int nc101  = nc + dixe        + dize;
						
						/* diaginv[n] = 4.0 / (Sol->rhoY[nne] + Sol->rhoY[nnw] + Sol->rhoY[nsw] + Sol->rhoY[nse]); */
						diaginv[n] = 8.0 / (  Sol->rhoY[nc000] / Sol->rho[nc000]
											+ Sol->rhoY[nc010] / Sol->rho[nc010] 
											+ Sol->rhoY[nc011] / Sol->rho[nc011]
											+ Sol->rhoY[nc001] / Sol->rho[nc001]
                                            + Sol->rhoY[nc100] / Sol->rho[nc100]
											+ Sol->rhoY[nc110] / Sol->rho[nc110] 
											+ Sol->rhoY[nc111] / Sol->rho[nc111]
											+ Sol->rhoY[nc101] / Sol->rho[nc101]);
                    }
                }
			}
#else /* PRECON */
            int n;

			for(n=0; n<node->nc; n++) {
				diaginv[n] = 1.0;
			}	
#endif
			
			break;
		}
		default: ERROR("ndim not in {1,2,3}");
	}
}
#endif

/* ========================================================================== */

void correction_nodes(
					  ConsVars* Sol,
					  const ElemSpaceDiscr* elem,
					  const NodeSpaceDiscr* node,
                      const double* hplus[3], 
                      const double* hgrav,
					  double* p, 
					  const double t,
					  const double dt) {
	
	extern User_Data ud;
	
	const int ndim = elem->ndim;
	
#ifndef NO_BDRYCONDS_PROJ2
	set_ghostnodes_p(p, node, 1);   
#endif
    
	switch(ndim) {
		case 1: {
			ERROR("function not available");
			break;
		}
			
		case 2: {
			
			const int igx = node->igx;
			const int icx = node->icx;
			const int igy = node->igy;
			const int icy = node->icy;
			
			const int icxe = node->icx - 1;
			
			const double dx = node->dx;
			const double dy = node->dy;
			const double oodx = 1.0 / dx;
			const double oody = 1.0 / dy;
			const double dtowdx = 0.5*dt * oodx;
			const double dtowdy = 0.5*dt * oody;
            
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
			
			int i, j, m, me;
			
			for(j = igy; j < icy - igy - 1; j++) {
				m = j * icx; 
				me = j * icxe;
				
				for(i = igx; i < icx - igx - 1; i++) {
					const int n = m + i;
					const int nicx = n + icx;
					const int n1 = n + 1;
					const int n1icx = n1 + icx;
					
					const int ne = me + i; 
					
					const double Dpx   = 0.5 * (p[n1]   - p[n] + p[n1icx] - p[nicx]);
					const double Dpy   = 0.5 * (p[nicx] - p[n] + p[n1icx] - p[n1]);
                    const double thinv = Sol->rho[ne] / Sol->rhoY[ne] ;

					Sol->rhou[ne] += - dtowdx * thinv * hplusx[ne] * Dpx;
					Sol->rhov[ne] += - dtowdy * thinv * hplusy[ne] * Dpy;
                    Sol->rhov[ne] += - 0.5*dt * hgrav[ne] * (Sol->rho[ne]/Sol->rhoY[ne]) * 0.25 * (p[n] + p[n1] + p[n1icx] + p[nicx]);
				}
			} 
			
			break;
		}
		case 3: {
			
			extern User_Data ud;
			
			const int igxe = elem->igx;
			const int icxe = elem->icx;
			const int igye = elem->igy;
			const int icye = elem->icy;
			const int igze = elem->igz;
			const int icze = elem->icz;
            
			const int dixe = 1;
			const int diye = icxe;
			const int dize = icxe*icye;
			
			const int icxn = node->icx;
			const int icyn = node->icy;
            /*
			const int iczn = node->icz;    
             */
            const int dixn = 1;
			const int diyn = icxn;
			const int dizn = icxn*icyn;
            
			const double dx = node->dx;
			const double dy = node->dy;
			const double dz = node->dz;
			const double oodx = 1.0 / dx;
			const double oody = 1.0 / dy;
			const double oodz = 1.0 / dz;
			const double dtowdx = 0.5*dt * oodx;
			const double dtowdy = 0.5*dt * oody;
			const double dtowdz = 0.5*dt * oodz;
			
			int i, j, k, ln, mn, nn, le, me, ne;
			
			for(k = igze; k < icze - igze; k++) {
				le = k * dize; 
				ln = k * dizn;
								
				for(j = igye; j < icye - igye; j++) {
					me = le + j * diye; 
					mn = ln + j * diyn;
                    
					for(i = igxe; i < icxe - igxe; i++) {
						ne = me + i * dixe;
						nn = mn + i * dixn;
                        
						const int nn000 = nn;
						const int nn010 = nn        + diyn;
						const int nn011 = nn        + diyn + dizn;
						const int nn001 = nn               + dizn;
						const int nn100 = nn + dixn;
						const int nn110 = nn + dixn + diyn;
						const int nn111 = nn + dixn + diyn + dizn;
						const int nn101 = nn + dixn        + dizn;
						
						const double Dpx = 0.25 * (  p[nn100] - p[nn000]
                                                   + p[nn110] - p[nn010]
                                                   + p[nn111] - p[nn011]
                                                   + p[nn101] - p[nn001]
                                                   );
						
						const double Dpy = 0.25 * (  p[nn010] - p[nn000]
                                                   + p[nn110] - p[nn100]
                                                   + p[nn111] - p[nn101]
                                                   + p[nn011] - p[nn001]
                                                   );
						
						const double Dpz = 0.25 * (  p[nn001] - p[nn000]
                                                   + p[nn101] - p[nn100]
						                           + p[nn111] - p[nn110]
						                           + p[nn011] - p[nn010]
                                                   );
						Sol->rhou[ne] += - dtowdx * Dpx;
						Sol->rhov[ne] += - dtowdy * Dpy;
						Sol->rhow[ne] += - dtowdz * Dpz;
					}
				} 
			} 
			break;
		}
		default: ERROR("ndim not in {1, 2,3}");
	}
}

/* ========================================================================== */

#ifdef NODAL_PROJECTION_ONLY

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

static void Advection_Flux_Correction(ConsVars* flux[3],
                                      VectorField* adv_flux1,
                                      VectorField* adv_flux0,
                                      const ConsVars* Sol,
                                      const ElemSpaceDiscr* elem) 
{
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
                        
            ConsVars* f = flux[0];
            ConsVars* g = flux[1];
                        
            double oorhoi, ui, vi, wi, Yi, Zi, Hi, oorhoim, uim, vim, wim, Yim, Zim, Him;
            double Xi[NSPEC], Xim[NSPEC];
            double oorhoj, uj, vj, wj, Yj, Zj, Hj, oorhojm, ujm, vjm, wjm, Yjm, Zjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];
            double us, vs, ws, Hs, Ys, Zs;
            double Xs[NSPEC];
            double frhoY, grhoY;
            
            double upwind;
            
            int i, j, mc, me, mw, nc, nn, ns, ic, icm, icmm, icp, jc, jcm, jcmm, jcp;
            
            for(j = igy; j < icy - igy; j++) {
                mc    = j * ifx;
                
                for(i = igx; i < ifx - igx; i++) {
                    nn   = mc + i + ifx;
                    nc   = mc + i;
                    ns   = mc + i - ifx;
                    ic   = j*icx + i;
                    icp  = ic + 1;
                    icm  = ic - 1;
                    icmm = ic - 2;
                                        
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
                    Hi      = oorhoi * INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]); /* pressure missing from enthalpy */
                    
                    oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                    uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                    vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                    wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                    }
                    Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
                    Zim     = oorhoim * INTERPOL(Sol->rhoZ[icmm], Sol->rhoZ[icm], Sol->rhoZ[ic]);
                    Him     = oorhoim * INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]); /* pressure missing from enthalpy */
                                        
                    frhoY = 0.5*(adv_flux1->x[nc] - adv_flux0->x[nc]);
                                        
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); /* What about clean upwinding ? */
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
                    
                    f->rho[nc]  = Ys * frhoY;
                    f->rhou[nc] = us * frhoY + us * frhoY;
                    f->rhov[nc] = vs * frhoY;  
                    f->rhow[nc] = ws * frhoY;  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        f->rhoX[nsp][nc] = Xs[nsp] * frhoY; 
                    }
                    f->rhoe[nc] = 0.0 * Hs * frhoY;
                    f->rhoY[nc] = frhoY;
                    f->rhoZ[nc] = Zs * frhoY;
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
                                                            
                    oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
                    uj = oorhoj * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
                    vj = oorhoj * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
                    wj = oorhoj * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                    }
                    Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
                    Zj = oorhoj * INTERPOL(Sol->rhoZ[jcp], Sol->rhoZ[jc], Sol->rhoZ[jcm]);
                    Hj = oorhoj * INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]); /* pressure missing from enthalpy */
                    
                    oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                    ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                    vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                    wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                    }
                    Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
                    Zjm = oorhojm * INTERPOL(Sol->rhoZ[jcmm], Sol->rhoZ[jcm], Sol->rhoZ[jc]);
                    Hjm = oorhojm * INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]);  /* pressure missing from enthalpy */
                    
                    grhoY = 0.5*(adv_flux1->y[mc] - adv_flux0->y[mc]);
                    
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(grhoY, 0.01) ); /* What about clean upwinding ? */
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
                    
                    g->rho[mc]  = Ys * grhoY;
                    g->rhou[mc] = us * grhoY;  
                    g->rhov[mc] = vs * grhoY + vs * grhoY; 
                    g->rhow[mc] = ws * grhoY; 
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        g->rhoX[nsp][mc] = Xs[nsp] * grhoY; 
                    }
                    g->rhoe[mc] = 0.0 * Hs * grhoY;
                    g->rhoY[mc] = grhoY;
                    g->rhoZ[mc] = Zs * grhoY;
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
            
            const int dix = 1;
            const int diy = icx; 
            const int diz = icx*icy;            
            
            ConsVars* fx = flux[0];
            ConsVars* fy = flux[1];
            ConsVars* fz = flux[2];
            
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
            
            double frhoY;
            
            double upwind;
            
            int i, j, k, l, m, n;
            
            int ic, icm, icmm, icp;
            int jc, jcm, jcmm, jcp;
            int kc, kcm, kcmm, kcp;
            
            for(k = igz; k < icz - igz; k++) {l = k * ifx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
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
                        Hi      = oorhoi * INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]); /* pressure missing from enthalpy */
                        
                        oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                        uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                        vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                        wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                        }
                        Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
                        Zim     = oorhoim * INTERPOL(Sol->rhoZ[icmm], Sol->rhoZ[icm], Sol->rhoZ[ic]);
                        Him     = oorhoim * INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]); /* pressure missing from enthalpy */
                        
                        frhoY = 0.5*(adv_flux1->x[n] - adv_flux0->x[n]);
                        
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
                        
                        fx->rho[n]  = Ys * frhoY;
                        fx->rhou[n] = us * frhoY + us * frhoY;
                        fx->rhov[n] = vs * frhoY;  
                        fx->rhow[n] = ws * frhoY;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fx->rhoX[nsp][n] = Xs[nsp] * frhoY; 
                        }
                        fx->rhoe[n] = 0.0 * Hs * frhoY;
                        fx->rhoY[n] = frhoY;
                        fx->rhoZ[n] = Zs * frhoY;
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
                        
                        oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
                        uj = oorhoj  * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
                        vj = oorhoj  * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
                        wj = oorhoj  * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                        }
                        Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
                        Zj = oorhoj * INTERPOL(Sol->rhoZ[jcp], Sol->rhoZ[jc], Sol->rhoZ[jcm]);
                        Hj = oorhoj * INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]); /* pressure missing from enthalpy */
                        
                        oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                        ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                        vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                        wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                        }
                        Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
                        Zjm = oorhojm * INTERPOL(Sol->rhoZ[jcmm], Sol->rhoZ[jcm], Sol->rhoZ[jc]);
                        Hjm = oorhojm * INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]); /* pressure missing from enthalpy */
                        
                        frhoY = 0.5*(adv_flux1->y[n] - adv_flux0->y[n]);
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); /* What about clean upwinding ? */
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
                        
                        fy->rho[n]  = Ys * frhoY;
                        fy->rhou[n] = us * frhoY;  
                        fy->rhov[n] = vs * frhoY + vs * frhoY; 
                        fy->rhow[n] = ws * frhoY;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fy->rhoX[nsp][n] = Xs[nsp] * frhoY; 
                        }
                        fy->rhoe[n] = 0.0 * Hs * frhoY;
                        fy->rhoY[n] = frhoY;
                        fy->rhoZ[n] = Zs * frhoY;
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
                                                
                        oorhok = 1.0 /INTERPOL(Sol->rhoY[kcp], Sol->rhoY[kc], Sol->rhoY[kcm]);
                        uk = oorhok * INTERPOL(Sol->rhou[kcp], Sol->rhou[kc], Sol->rhou[kcm]);
                        vk = oorhok * INTERPOL(Sol->rhov[kcp], Sol->rhov[kc], Sol->rhov[kcm]);
                        wk = oorhok * INTERPOL(Sol->rhow[kcp], Sol->rhow[kc], Sol->rhow[kcm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xk[nsp]      = oorhok * INTERPOL(Sol->rhoX[nsp][kcp], Sol->rhoX[nsp][kc], Sol->rhoX[nsp][kcm]);
                        }
                        Yk = oorhok * INTERPOL(Sol->rho[kcp], Sol->rho[kc], Sol->rho[kcm]);
                        Zk = oorhok * INTERPOL(Sol->rhoZ[kcp], Sol->rhoZ[kc], Sol->rhoZ[kcm]);
                        Hk = oorhok * INTERPOL(Sol->rhoe[kcp], Sol->rhoe[kc], Sol->rhoe[kcm]); /* pressure missing from enthalpy */
                        
                        oorhokm = 1.0 / INTERPOL(Sol->rhoY[kcmm], Sol->rhoY[kcm], Sol->rhoY[kc]);
                        ukm = oorhokm * INTERPOL(Sol->rhou[kcmm], Sol->rhou[kcm], Sol->rhou[kc]);
                        vkm = oorhokm * INTERPOL(Sol->rhov[kcmm], Sol->rhov[kcm], Sol->rhov[kc]);
                        wkm = oorhokm * INTERPOL(Sol->rhow[kcmm], Sol->rhow[kcm], Sol->rhow[kc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xkm[nsp]      = oorhokm * INTERPOL(Sol->rhoX[nsp][kcmm], Sol->rhoX[nsp][kcm], Sol->rhoX[nsp][kc]);
                        }
                        Ykm = oorhokm * INTERPOL(Sol->rho[kcmm], Sol->rho[kcm], Sol->rho[kc]);
                        Zkm = oorhokm * INTERPOL(Sol->rhoZ[kcmm], Sol->rhoZ[kcm], Sol->rhoZ[kc]);
                        Hkm = oorhokm * INTERPOL(Sol->rhoe[kcmm], Sol->rhoe[kcm], Sol->rhoe[kc]); /* pressure missing from enthalpy */
                        
                        frhoY = 0.5*(adv_flux1->z[n] - adv_flux0->z[n]);
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); /* What about clean upwinding ? */
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
                        
                        fz->rho[n]  = Ys * frhoY;
                        fz->rhou[n] = us * frhoY;  
                        fz->rhov[n] = vs * frhoY; 
                        fz->rhow[n] = ws * frhoY + ws * frhoY;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fz->rhoX[nsp][n] = Xs[nsp] * frhoY; 
                        }
                        fz->rhoe[n] = 0.0 * Hs * frhoY;
                        fz->rhoY[n] = frhoY;
                        fz->rhoZ[n] = Zs * frhoY;
                    }   
                }
            }			
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
}

/* ========================================================================== */

void Pressure_node_to_elem(MPV* mpv, 
                           const ConsVars* Sol, 
                           const ConsVars* Sol0, 
                           const ElemSpaceDiscr* elem, 
                           const NodeSpaceDiscr* node)
{    
    extern User_Data ud;
    
    const int icxe = elem->icx;
    const int icye = elem->icy;
    const int icze = elem->icz;

    const int stride_ej = icxe;
    const int stride_ek = icye*icxe;

    const int icxn = node->icx;
    const int icyn = node->icy;
    
    const int stride_nj = icxn;
    const int stride_nk = icyn*icxn;

    const int di = 2;
    const int dj = (elem->ndim > 1 ? 2 : 1);
    const int dk = (elem->ndim > 2 ? 2 : 1  );
    
    const double Msq = ud.Msq;    
    const double *pn   = mpv->p2_nodes;
    const double *dpn  = mpv->dp2_nodes;
    const double *phye = mpv->HydroState->p0;
    const double *phyn = mpv->HydroState_n->p0;
    
    double** hplus       = mpv->Level[0]->wplus;
    double*  hcenter     = mpv->Level[0]->wcenter;
    double*  hgrav       = mpv->Level[0]->wgrav;

    double *pc  = mpv->p2_cells;
    double *dpc = mpv->dp2_cells;

    /* 
     straight-forward interpolation of nodal to cell-centered pressures, 
     modulo the hydrostatic background
     */
    for (int k=0; k<icze; k++) {
        int nek = k*stride_ek;
        int nnk = k*stride_nk;
        
        for (int j=0; j<icye; j++) {
            int nekj = nek + j*stride_ej;
            int nnkj = nnk + j*stride_nj;

            for (int i=0; i<icxe; i++) {
                int nekji = nekj + i;
                int nnkji = nnkj + i;
                int    we   = 0;
                double pcc  = 0.0;
                double dpcc = 0.0;
                
                for (int kk=0; kk<dk; kk++) {
                    for (int jj=0; jj<dj; jj++) {
                        for (int ii=0; ii<di; ii++) {
                            dpcc += dpn[nnkji + kk*stride_nk + jj*stride_nj + ii];
                            pcc  += pn[nnkji + kk*stride_nk + jj*stride_nj + ii] - phyn[j+jj]/Msq;
                            we  += 1;
                        }
                    }
                }
                dpcc /= we;
                pcc  /= we;
                pcc  += phye[j]/Msq;
                dpc[nekji] = dpcc;
                pc[nekji]  = pcc;
            }
        }
    }
    
    /* set dummy cell data for the cell centered pressure data */
    operator_coefficients(hplus, hcenter, hgrav, elem, Sol, Sol0, mpv, mpv->dt);
    set_ghostcells_p2(mpv->dp2_cells, (const double**)hplus, hgrav, elem, elem->igx);
    set_ghostcells_p2(mpv->p2_cells, (const double**)hplus, hgrav, elem, elem->igx);

}

/* ========================================================================== */

void Divergence_Control(ConsVars* Sol, 
                        ConsVars* flux[3], 
                        VectorField* adv_flux1,
                        VectorField* buoy,
                        MPV* mpv,
                        const ConsVars* Sol0, 
                        const ElemSpaceDiscr* elem, 
                        const NodeSpaceDiscr* node, 
                        const double t, 
                        const double dt,
                        const double tout,
                        const int step)
{
    extern User_Data ud;
    extern BDRY* bdry;
    extern VectorField* adv_flux0;

    extern double *W0;
    double *rhs = W0;
    
    const int nfx = elem->nfx;
    const int nfy = elem->nfy;
    const int nfz = elem->nfz;

    const double tol = ud.second_projection_precision / (0.5*dt);
    double delta;

    /* buoy needed here only so that I am able to reuse the update() function */
    VectorField_setzero(buoy, elem->nc);
    
    /* Define and initialize adv_flux0 to rhoY-components of flux */
    
    /* Iterative divergence control algorithm:
     
     1) initial setup:
     - nu = 0
     - adv_flux1_nu    holds averaged flux from previous time level   adv_fluxn
     - adv_flux0_nu    holds flux from explicit advection predictor   adv_flux*
     
     2) after the subsequent loop, 
     - adv_flux0_nu    holds    2 * adv_flux0_nu  -  adv_flux1_nu  (=  2 * adv_flux*  -  adv_fluxn   at nu=0 )
     
     3) after second projection and recomputation:
     - adv_flux1_nu+1  holds averaged flux from nu-th correction of present time level
     
     4) advective flux correction:
     - flux            holds   0.5 * (adv_flux1_nu+1 - adv_flux0_nu)
     - Sol             gets corrected with flux
     
     5) reset fluxes for next iterate:
     - adv_flux0_nu+1  holds  adv_flux1_nu+1
     - nu              holds  nu+1
     
     6) branchpoint:
     - tolerance not satisfied -> 3)
     - tolerance satisfied     -> move on
     
     */
    
    /* 
     VectorField* adv_flux0 = (VectorField*)malloc(sizeof(VectorField)); 
     adv_flux0->x = flux[0]->rhoY;
     */
    for (int i=0; i<nfx; i++) adv_flux0->x[i] = 2.0*flux[0]->rhoY[i] - 1.0*adv_flux1->x[i];
    
    if (elem->ndim > 1) {
        adv_flux0->y = flux[1]->rhoY;
        for (int i=0; i<nfy; i++) adv_flux0->y[i] = 2.0*flux[1]->rhoY[i] - 1.0*adv_flux1->y[i];        
    }
    
    if (elem->ndim > 2) {
        adv_flux0->z = flux[2]->rhoY;
        for (int i=0; i<nfz; i++) adv_flux0->z[i] = 2.0*flux[2]->rhoY[i] - 1.0*adv_flux1->z[i];        
    }
    
    /* iteration to control both the node- and cell-centered divergences */    
    delta = 10000.0;
    
    /* while (delta > 4.0*tol) { */
    for (int cnt=0; cnt<2; cnt++) {
        double rhs_weight = 2.0;  /* why that number and not another one? */
        second_projection(Sol, mpv, Sol0, elem, node, 1.0, t, dt);
        Set_Explicit_Boundary_Data(Sol, elem, mpv);
#if OUTPUT_SUBSTEPS
        putout(Sol, t, tout , step, 0, ud.file_name, "Sol", 1);
#endif
        Advective_Fluxes(adv_flux1, Sol, elem);
        Advection_Flux_Correction(flux, adv_flux1, adv_flux0, (const ConsVars*)Sol, elem);
        update(Sol, (const ConsVars**)flux, buoy, elem, dt);
        Set_Explicit_Boundary_Data(Sol, elem, mpv);
#if OUTPUT_SUBSTEPS
        putout(Sol, t, tout , step, 0, ud.file_name, "Sol", 1);
#endif
        VectorField_set(adv_flux0, (const VectorField*)adv_flux1, (const int)node->nc, (const int)node->ndim);    
        delta = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv->eta, mpv, bdry, dt, rhs_weight);
    }
    
    /* free(adv_flux0); */
}

#endif /* NODAL_PROJECTION_ONLY */


#endif /* SOLVER_2_HYPRE */

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
