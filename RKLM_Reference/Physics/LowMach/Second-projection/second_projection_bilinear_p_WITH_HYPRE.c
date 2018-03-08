
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
        
	const int nc = node->nc;
	
	double** hplus   = mpv->Level[0]->wplus;
	double*  hcenter = mpv->Level[0]->wcenter;
	
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
		p2[ii] = mpv->p2_nodes[ii];
		rhs[ii] = 0.0;
	}
    
    rhs_max = divergence_nodes(rhs, elem, node, (const ConsVars*)Sol, mpv->eta, mpv, bdry, dt, 1.0);
    printf("\nrhsmax = %e\n", rhs_max);

    catch_periodic_directions(rhs, node, elem, x_periodic, y_periodic, z_periodic);
    assert(integral_condition_nodes(rhs, node, x_periodic, y_periodic, z_periodic) != VIOLATED); 
    if (ud.is_compressible) {
        for (int nn=0; nn<node->nc; nn++) {
            rhs[nn] += hcenter[nn]*mpv->p2_nodes[nn];
        }
    }
         
	operator_coefficients_nodes(hplus, hcenter, elem, node, Sol, Sol0, mpv, dt);
	variable_coefficient_poisson_nodes(p2, (const double **)hplus, hcenter, rhs, x_periodic, y_periodic, z_periodic, dt);
    correction_nodes(Sol, elem, node, (const double**)hplus, p2, t, dt);
    
    for(ii=0; ii<nc; ii++) {
        mpv->dp2_nodes[ii] = p2[ii] - mpv->p2_nodes[ii];
        mpv->p2_nodes[ii]  = p2[ii];
    }

    Set_Explicit_Boundary_Data(Sol, elem);
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
	
    /* with weight = 1.0, this routine computes 
          rhs_out = rhs_in + (2/dt) * div(rhoY\vec{v}) 
       from the current Sol. 
     */

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
            
            const int icxn = node->icx;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double oodxdt = 1.0 / (dx * dt);
            const double oodydt = 1.0 / (dy * dt);
            const double todt = 2.0 / dt;
            const double oowdxdt = weight * oodxdt;
            const double oowdydt = weight * oodydt;
            
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
            
            for(j = igye+1; j < icye - igye; j++) {
                const int mn = j * icxn; 
                for(i = igxe+1; i < icxe - igxe; i++) {
                    int nn    = mn + i;
                    div_max = MAX_own(div_max, fabs(rhs[nn]));
                }
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
    
    return div_max;
    
}

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
    
    const int imax = MAX_own(1, icx - igx - x_periodic);
    const int jmax = MAX_own(1, icy - igy - y_periodic);
    const int kmax = MAX_own(1, icz - igz - z_periodic);
    
    const double del = DBL_EPSILON * node->nc;
    
	double tmp = 0.0, rhs_max=0.0;
	
	for(k = igz; k < kmax; k++) {l = k * icxicy;
		for(j = igy; j < jmax; j++) {m = l + j * icx; 
			for(i = igx; i < imax; i++) {n = m + i;
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
	
	if(node->ndim > 1 && y_periodic == 1){
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
	
	if(node->ndim > 2 && z_periodic == 1){
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
										const ElemSpaceDiscr* elem,
										const NodeSpaceDiscr* node,
										ConsVars* Sol,
										const ConsVars* Sol0,
										const MPV* mpv,
                                        const double dt) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	
    const double g   = ud.gravity_strength[1];
    const double Msq = ud.Msq;
    const double Gammainv = th.Gammainv;
    
    const int ndim = node->ndim;
                    
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

            const double dy = elem->dy;
			
			double* hplusx  = hplus[0];
			double* hplusy  = hplus[1];
			double* hc      = hcenter;

            /*  const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt); 
                const double cexp    = 1.0-th.gamm;
             */
			const double ccenter = - 4.0*(ud.compressibility*ud.Msq)*th.gm1inv/(mpv->dt*mpv->dt);
            const double cexp    = 2.0-th.gamm;
            
            int i, j, m, n;
            
            for(i=0; i<elem->nc; i++) hc[i] = hplusx[i] = hplusy[i] = 0.0;
						                        
			for(j = igy; j < icy - igy; j++) {
                m = j * icx;
				
                double strat = 2.0 * (mpv->HydroState_n->Y0[j+1]-mpv->HydroState_n->Y0[j]) \
                                    /(mpv->HydroState_n->Y0[j+1]+mpv->HydroState_n->Y0[j])/dy;

                for(i = igx; i < icx - igx; i++) {
                    n = m + i;     
                    
                    double Y     = 0.5 * (Sol->rhoY[n]/Sol->rho[n] + Sol0->rhoY[n]/Sol0->rho[n]); 
                    double coeff = Gammainv * Sol->rhoY[n] * Y;
                    double Nsqsc = 0.25*dt*dt * (g/Msq) * strat;                    
                    double gimpy = 1.0 / (1.0 + Nsqsc);
                                        
                    hplusx[n]    = coeff;
                    hplusy[n]    = coeff * gimpy;
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
            
            const double dy = elem->dy;

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
                    
                    double strat = 2.0 * (mpv->HydroState_n->Y0[j+1]-mpv->HydroState_n->Y0[j]) \
                                        /(mpv->HydroState_n->Y0[j+1]+mpv->HydroState_n->Y0[j])/dy;

                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        {             
                            double Y     = 0.5 * (Sol->rhoY[n]/Sol->rho[n] + Sol0->rhoY[n]/Sol0->rho[n]); 
                            double coeff = Gammainv * Sol->rhoY[n] * Y;
                            double Nsqsc = 0.25*dt*dt * (g/Msq) * strat;                    
                            double gimpy = 1.0 / (1.0 + Nsqsc);

                            hplusx[n]  = coeff;
                            hplusy[n]  = coeff * gimpy;
                            hplusz[n]  = coeff;
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
					  double* p, 
					  const double t,
					  const double dt) {
	
	extern User_Data ud;
    extern MPV* mpv;

	const int ndim = elem->ndim;
	    
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
            
            double dtowdx;
            double dtowdy;
            if (ud.time_integrator == SI_MIDPT) {
                dtowdx = 0.5*dt * oodx; /* 1.0*dt * oodx; */
                dtowdy = 0.5*dt * oody; /* 1.0*dt * oody; */
            } else {
                dtowdx = 0.5*dt * oodx;
                dtowdy = 0.5*dt * oody;                
            }
            
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
			
			int i, j, m, me;
			
			for(j = igy; j < icy - igy - 1; j++) {
				m = j * icx; 
				me = j * icxe;
				for(i = igx; i < icx - igx - 1; i++) {
					const int n     = m + i;
					const int nicx  = n + icx;
					const int n1    = n + 1;
					const int n1icx = n1 + icx;
					const int ne    = me + i; 
					
					const double Dpx   = 0.5 * (p[n1]   - p[n] + p[n1icx] - p[nicx]);
					const double Dpy   = 0.5 * (p[nicx] - p[n] + p[n1icx] - p[n1]);
                    const double thinv = Sol->rho[ne] / Sol->rhoY[ne];
					                    
                    Sol->rhou[ne] += - dtowdx * thinv * hplusx[ne] * Dpx;
					Sol->rhov[ne] += - dtowdy * thinv * hplusy[ne] * Dpy;
				}
			} 
			
			break;
		}
        case 3: {
            
            const int igx = node->igx;
            const int icx = node->icx;
            const int igy = node->igy;
            const int icy = node->icy;
            const int igz = node->igz;
            const int icz = node->icz;
            
            const int icxe = node->icx - 1;
            const int icye = node->icy - 1;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double dz = node->dz;
            const double oodx = 1.0 / dx;
            const double oody = 1.0 / dy;
            const double oodz = 1.0 / dz;
            
            double dtowdx;
            double dtowdy;
            double dtowdz;
            
            if (ud.time_integrator == SI_MIDPT) {
                dtowdx = 0.5*dt * oodx; /* 1.0*dt * oodx; */
                dtowdy = 0.5*dt * oody; /* 1.0*dt * oody; */
                dtowdz = 0.5*dt * oodz; /* 1.0*dt * oodz; */
            } else {
                dtowdx = 0.5*dt * oodx;
                dtowdy = 0.5*dt * oody;                
                dtowdz = 0.5*dt * oodz;                
            }
            
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
            const double* hplusz = hplus[2];
                        
            for(int k = igz; k < icz-igz-1; k++) {
                int l  = k*icx*icy;
                int le = k*icxe*icye;
                for(int j = igy; j < icy - igy - 1; j++) {
                    int m  = l  + j*icx; 
                    int me = le + j*icxe;
                    for(int i = igx; i < icx - igx - 1; i++) {
                        int n000 = m + i;
                        int n010 = n000 + icx;
                        int n001 = n000 + 1;
                        int n011 = n000 + 1 + icx;
                        int n100 = m + i + icx*icy;
                        int n110 = n100 + icx;
                        int n101 = n100 + 1;
                        int n111 = n100 + 1 + icx;

                        int ne   = me + i; 
                        
                        double Dpx   = 0.25 * (p[n001] - p[n000] + p[n011] - p[n010] + p[n101] - p[n100] + p[n111] - p[n110]);
                        double Dpy   = 0.25 * (p[n010] - p[n000] + p[n011] - p[n001] + p[n110] - p[n100] + p[n111] - p[n101]);
                        double Dpz   = 0.25 * (p[n100] - p[n000] + p[n110] - p[n010] + p[n101] - p[n001] + p[n111] - p[n011]);
                        double thinv = Sol->rho[ne] / Sol->rhoY[ne];
                        
                        Sol->rhou[ne] += - dtowdx * thinv * hplusx[ne] * Dpx;
                        Sol->rhov[ne] += - dtowdy * thinv * hplusy[ne] * Dpy;
                        Sol->rhow[ne] += - dtowdz * thinv * hplusz[ne] * Dpz;
                    }
                } 
            }
            
            break;
        }
		default: ERROR("ndim not in {1, 2,3}");
	}
}

/* ========================================================================== */

void euler_backward_gravity(ConsVars* Sol,
                            const MPV* mpv,
                            const ElemSpaceDiscr* elem,
                            const double dt)
{
    /* 
     evaluates Euler backward for the gravity term in a semi-implicit 
     trapezoidal or midpoint discretization of the pressure gradient and 
     gravity terms. 
     */
    extern User_Data ud;
    
    const double g   = ud.gravity_strength[1];
    const double Msq = ud.Msq;
    const double dy  = elem->dy;
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
        
    for (int k=0; k<icz; k++) {
        int l = k*icy*icx;
        for (int j=0; j<icy; j++) {
            int m = l + j*icx;
            double strat = 2.0 * (mpv->HydroState_n->Y0[j+1]-mpv->HydroState_n->Y0[j])/(mpv->HydroState_n->Y0[j+1]+mpv->HydroState_n->Y0[j])/dy;

            for (int i=0; i<icx; i++) {
                int n        = m + i;
                double Nsqsc = dt*dt * (g/Msq) * strat;
                double v     = Sol->rhov[n]/Sol->rho[n];
                double dchi  = Sol->rhoX[BUOY][n]/Sol->rho[n];
                double chi   = Sol->rho[n]/Sol->rhoY[n];
                double dbuoy = -dchi/chi;    /* -dchi/chibar; */
                
                Sol->rhov[n] = Sol->rho[n] * (v + dt * (g/Msq) * dbuoy) / (1.0 + Nsqsc);
            }
        }
    }
    Set_Explicit_Boundary_Data(Sol, elem);
}

/* ========================================================================== */

void euler_forward_non_advective(ConsVars* Sol,
                                 const MPV* mpv,
                                 const ElemSpaceDiscr* elem,
                                 const NodeSpaceDiscr* node,
                                 const double dt)
{
    /* 
     evaluates Euler forward for the pressure gradient, gravity, 
     and background stratification advection terms based on cell-
     centered data. 
     */
    extern User_Data ud;
    extern Thermodynamic th;
    
    double *p2n       = mpv->p2_nodes;
    
    const double g    = ud.gravity_strength[1];
    const double Msq  = ud.Msq;
    const double Ginv = th.Gammainv; 
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;

    const int igx = elem->igx;
    const int igy = elem->igy;
    const int igz = elem->igz;

    const int inx = node->icx;
    const int iny = node->icy;
    
    const double dx = node->dx;
    const double dy = node->dy;
    const double dz = node->dz;
        
    for (int k=igz; k<icz-igz; k++) {
        int lc = k*icy*icx;
        int ln = k*iny*inx;
        for (int j=igy; j<icy-igy; j++) {
            int mc = lc + j*icx;
            int mn = ln + j*inx;
            double S0p    = mpv->HydroState_n->S0[j+1];
            double S0m    = mpv->HydroState_n->S0[j];
            
            for (int i=igx; i<icx-igx; i++) {
                int nc        = mc + i;

                int n000 = mn   + i;
                int n010 = n000 + inx;
                int n001 = n000 + 1;
                int n011 = n000 + 1 + inx;
                int n100 = mn   + i + inx*iny;
                int n110 = n100 + inx;
                int n101 = n100 + 1;
                int n111 = n100 + 1 + inx;

                
                double dpdx   = 0.25*(p2n[n001]-p2n[n000]+p2n[n011]-p2n[n010]+p2n[n101]-p2n[n100]+p2n[n111]-p2n[n110])/dx;
                double dpdy   = 0.25*(p2n[n010]-p2n[n000]+p2n[n011]-p2n[n001]+p2n[n110]-p2n[n100]+p2n[n111]-p2n[n101])/dy;
                double dpdz   = 0.25*(p2n[n100]-p2n[n000]+p2n[n110]-p2n[n010]+p2n[n101]-p2n[n001]+p2n[n111]-p2n[n011])/dz;
                double dSdy   = (S0p-S0m) / dy;
                
                double rhoYovG = Ginv*Sol->rhoY[nc];
                double v       = Sol->rhov[nc]/Sol->rho[nc];
                double dchi    = Sol->rhoX[BUOY][nc]/Sol->rho[nc];
                double chi     = Sol->rho[nc]/Sol->rhoY[nc];
                double dbuoy   = -Sol->rho[nc]*dchi/chi;  /* -dchi/chibar; */
                
                Sol->rhou[nc]       += dt * ( - rhoYovG * dpdx);
                Sol->rhov[nc]       += dt * ( - rhoYovG * dpdy + (g/Msq) * dbuoy);
                Sol->rhow[nc]       += dt * ( - rhoYovG * dpdz);
                Sol->rhoX[BUOY][nc] += dt * ( - v * dSdy) * Sol->rho[nc];
            }
        }
    }
    
    /* last half Euler backward step equals first half Euler forward step */
    if (ud.is_compressible) {
        for (int nn=0; nn<node->nc; nn++) {
            mpv->p2_nodes[nn] += mpv->dp2_nodes[nn];
        }
    }
    Set_Explicit_Boundary_Data(Sol, elem);
}

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
