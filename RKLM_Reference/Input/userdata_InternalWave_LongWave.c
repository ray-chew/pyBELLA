/*******************************************************************************
 File:   userdata.c
 *******************************************************************************/
#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "enumerator.h"
#include "Common.h"
#include "userdata.h"
#include "time_discretization.h"
#include "error.h"
#include "variable.h"
#include "enum_bdry.h"
#include "math_own.h"
#include "space_discretization.h"
#include "thermodynamic.h"
#include "Eos.h"
#include "set_ghostcells_p.h"
#include "set_ghostnodes_p.h"
#include "boundary.h"
#include "memory.h"
#include "Hydrostatics.h"

double pressure_function(double r, double p0, double S0, double u_theta, double Msq, double Gamma);

double rho_function(double psi);

double molly(double x);

/* horizontal stretch for S&K94 IGWs: planetary -> 8*20.0;  long-wave -> 20.0;  standard -> 1.0; */
static double scalefactor = 1.0;   

void User_Data_init(User_Data* ud) {
    
    int i, max_no_of_levels; 
    
    /* ================================================================================== */
    /* =====  PROBLEM SET UP  =========================================================== */
    /* ================================================================================== */
    
    /* Earth */
    double grav     = 9.81;             /* [m/s^2]                */
    double omega    = 0.0*0.0001;     /* 1.454 */
    /* double omega  = sin(0.5*PI) * 2.0 * 0.00007272205217;   */
    
    /* thermodynamics and chemistry */
    double R_gas    = 287.4;            /* [J/kg/K]               */
    double R_vap    = 461.00;           /* [J/kg/K]               */
    double Q_vap    = 2.53e+06;         /* [J]                    */
    double gamma    = 1.4;              /* dimensionless          */

    /* references for non-dimensionalization */
    double h_ref    = 10000;                 /* [m]               */
    double t_ref    = 100;                   /* [s]               */
    double T_ref    = 300.00;                /* [K]               */
    double p_ref    = 10e+5;                 /* [Pa]              */
    double u_ref    = h_ref/t_ref;           /* [m/s]; Sr = 1     */
    /* double rho_ref = p_ref / (R_gas*T_ref); [kg/m^3]          */
    
    /* reference stratification */
    double Nsq_ref  = 1.0e-4;           /* [1/s^2]                */
    
    ud->h_ref       = h_ref;
    ud->t_ref       = t_ref;
    ud->T_ref       = T_ref;
    ud->p_ref       = p_ref;
    ud->u_ref       = u_ref;
    ud->Nsq_ref     = Nsq_ref;
    ud->g_ref       = grav;
    ud->gamm        = gamma;
    ud->Rg_over_Rv  = R_gas/R_vap;
    ud->Q           = Q_vap/(R_gas*T_ref);
    
    /* number of advected species */
    ud->nspec       = NSPEC;

    /* Low Mach */
    ud->is_compressible   = 1;   /* 0: psinc;  1: comp;  -1: psinc-comp-transition (see compressibility()) */
    ud->acoustic_timestep = 0; /* 0;  1; */
    ud->Msq =  u_ref*u_ref / (R_gas*T_ref);
    
    /* geo-stuff */
    for(i=0; i<3; i++) {
        ud->gravity_strength[i]  = 0.0;   /* corresponds to  M^2 / Fr^2  =  g href / R Tref */
        ud->coriolis_strength[i] = 0.0;
    }
    ud->gravity_strength[1]  = grav * h_ref / (R_gas * T_ref);
    ud->coriolis_strength[0] = omega * t_ref;
    ud->coriolis_strength[2] = omega * t_ref;
    
    /* integer gravity indicator */
    for (i=0; i<3; i++){
        ud->i_gravity[i]  = 0;
        ud->i_coriolis[i] = 0;
        if (ud->gravity_strength[i] > ud->eps_Machine || i==1) {
            ud->i_gravity[i] = 1;
            ud->gravity_direction = i;
        }
        if (ud->coriolis_strength[i] > ud->eps_Machine) {
            ud->i_coriolis[i] = 1;
        }
    }
    
    ud->implicit_gravity_theta  = 1;
    ud->implicit_gravity_theta2 = 1;
    
    /* flow domain */
    ud->xmin = -15.0 * scalefactor;
    ud->xmax =  15.0 * scalefactor;
    ud->ymin =   0.0;
    ud->ymax =   1.0;
    ud->zmin = - 1.0;
    ud->zmax =   1.0;
    
    /* boundary/initial conditions */
    ud->wind_speed  =  1.0 * 20.0/u_ref;
    ud->wind_shear  = -0.0;              /* velocity in [u_ref/h_ref] */
    ud->hill_height = 0.0 * 0.096447; 
    ud->hill_length_scale = 0.1535;   /* hill_length * l_ref = 1.0 km */
    
    ud->bdrytype_min[0] = PERIODIC; /* DIRICHLET; PERIODIC; WALL; */
    ud->bdrytype_min[1] = WALL; /* SLANTED_WALL; */
    ud->bdrytype_min[2] = WALL;
    ud->bdrytype_max[0] = PERIODIC; /* DIRICHLET; PERIODIC; WALL; */
    ud->bdrytype_max[1] = WALL;
    ud->bdrytype_max[2] = WALL;
    
    ud->absorber = WRONG; /* CORRECT;  WRONG; */ /*  BREAKING WAVE CHANGE */
    
    /* ================================================================================== */
    /* =====  NUMERICS  ================================================================= */
    /* ================================================================================== */
    
    /* time discretization */
    ud->time_integrator        = SI_MIDPT; /*OP_SPLIT, OP_SPLIT_MD_UPDATE, SI_MIDPT */
    ud->CFL                    = 0.96; /* 0.45; 0.9; 0.8; */
    ud->dtfixed0               = scalefactor*30.0 / ud->t_ref;
    ud->dtfixed                = scalefactor*30.0 / ud->t_ref;
    ud->no_of_steps_to_CFL     = 1;
    ud->no_of_steps_to_dtfixed = 1;
    
    set_time_integrator_parameters(ud);
    
    /* Grid and space discretization */
    ud->inx =  300+1; /* 641; 321; 161; 129; 81; */
    ud->iny =   20+1; /* 321; 161;  81;  65; 41;  */
    ud->inz = 1;
    ud->h   = MIN_own((ud->xmax-ud->xmin)/(ud->inx),MIN_own((ud->ymax-ud->ymin)/(ud->iny),(ud->zmax-ud->zmin)/(ud->inz)));
    
    /* explicit predictor step */
    /* Recovery */
    ud->recovery_order = SECOND; /* FIRST, SECOND */ 
    ud->limiter_type_scalars  = NONE; 
    ud->limiter_type_velocity = NONE; 
    /*  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
    
    /* parameters for SWEBY_MUNZ limiter family */
    ud->kp = 1.4;
    ud->kz = 1.4; /* Entro abused for velocity in split-step-aligned velocity ! */
    ud->km = 1.4;
    ud->kY = 1.4;
    ud->kZ = 1.4; /* 2.0 */
        
    ud->ncache = 27; /* 71+4; 304*44; 604*44; (ud->inx+3); (ud->inx+3)*(ud->iny+3);*/
    
    /* linear solver-stuff */
    ud->which_projection_first = 1;
    ud->Solver = BICGSTAB_PRECON;        /* options:   JACOBI, BICGSTAB, BICGSTAB_PRECON */
    ud->Solver_Node = BICGSTAB_PRECON;   /* options:   JACOBI, BICGSTAB, BICGSTAB_PRECON */
    ud->precondition = CORRECT;
    double tol = 1.e-10;
    ud->flux_correction_precision = tol;
    ud->flux_correction_local_precision = tol;   /* 1.e-05 should be enough */
    ud->second_projection_precision = tol;
    ud->second_projection_local_precision = tol;   /* 1.e-05 should be enough */
    ud->implicitness = 1.0;
    ud->flux_correction_max_MG_cycles = 100;
    ud->flux_correction_output_period = 50;
    ud->max_projection_iterations = 1;
    ud->flux_correction_max_iterations = 6000;
    ud->second_projection_max_iterations = 6000;
    
    
    /* Code could do geometric multigrid some time ago; that's no longer true. Therefore ... */
    max_no_of_levels = 1;
    ud->max_no_of_multigrid_levels = max_no_of_levels;
    ud->no_of_multigrid_levels     = max_no_of_levels-1; 
    
    /* numerics parameters */
    ud->eps_Machine = sqrt(DBL_EPSILON);
    
    /* ================================================================================== */
    /* =====  CODE FLOW CONTROL  ======================================================== */
    /* ================================================================================== */
    
    ud->tout[0] = scalefactor * 3000.0 / ud->t_ref;
    ud->tout[1] = -1.0;

    ud->stepmax = 10000;
    
    ud->write_stdout = ON;
    ud->write_stdout_period = 1;
    ud->write_file = ON;
    ud->write_file_period = 10;
    ud->file_format = HDF;
    
    {
        char *OutputBaseFolder      = "/Users/rupert/Documents/Computation/RKLM_Reference/";
        char *OutputFolderNamePsinc = "low_Mach_gravity_psinc";
        char *OutputFolderNameComp  = "low_Mach_gravity_comp";
        if (ud->is_compressible == 0) {
            sprintf(ud->file_name, "%s%s", OutputBaseFolder, OutputFolderNamePsinc);
        } else {
            sprintf(ud->file_name, "%s%s", OutputBaseFolder, OutputFolderNameComp);
        }
    }
}

/* ================================================================================== */

void Sol_initial(ConsVars* Sol, const ElemSpaceDiscr* elem, const NodeSpaceDiscr* node) {
    
    extern Thermodynamic th;
    extern User_Data ud;
    extern MPV* mpv;
    extern double *Sbg;
    
    const int compressible = (ud.is_compressible == 0 ? 0 : 1);
    
    const double u0 = ud.wind_speed;
    const double v0 = 0.0;
    const double w0 = 0.0;
    const double delth = 0.01 / ud.T_ref;  /* standard:  0.01 / ud.T_ref */
    const double xc    = -1.0*scalefactor*60.0e+03/ud.h_ref; /* -1000.0e+03/ud.h_ref;  -50.0e+03/ud.h_ref; 0.0; */
    /* const double a     = scalefactor * 5.0e+03/ud.h_ref; */
    const double a     = scalefactor*5.0e+03/ud.h_ref;   /* 1.0e+05/ud.h_ref;  1.0e+05/ud.h_ref/20; */
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int igy = elem->igy;
    
    const int icxn = node->icx;
    
    int i, j, m, n, nm, nn;
    double x, y, xn, yn, ym;
    double rho, u, v, w, p, rhoY;
    
    double g;
    
    States *HySt, *HyStn;
    double *Y, *Yn;
    
    g = ud.gravity_strength[1];
    
    HySt  = States_new(node->icy);
    Y     = (double*)malloc(node->icy*sizeof(double));
    HyStn = States_new(node->icy);
    Yn    = (double*)malloc(node->icy*sizeof(double));
    
    Hydrostatics_State(mpv, Sbg, elem, node);
    
    for(i = 0; i < icx; i++) {
        
        /* set potential temperature stratification in the column */
        for(j = 0; j < elem->icy; j++) {
            x     = elem->x[i];
            y     = elem->y[j];
            Y[j]  = stratification(y)  + delth * molly(x) * sin(PI*y)  / (1.0 + (x-xc)*(x-xc) / (a*a));
        }        
        for(j = 0; j < node->icy; j++) {
            xn    = node->x[i];
            yn    = node->y[j];
            Yn[j] = stratification(yn)  + delth * molly(xn) * sin(PI*yn)  / (1.0 + (xn-xc)*(xn-xc) / (a*a));
        }        
        Hydrostatics_Column(HySt, HyStn, Y, Yn, elem, node);
        
        /* initialization of field variables */
        for(j = igy; j < icy - igy + 1; j++) {
            
            n  = j*icx + i;
            nm = n-icx;
            
            x  = elem->x[i];
            y  = elem->y[j];
            ym = elem->y[j-1];
            
            u   = u0;
            v   = v0;
            w   = w0;
            
            p    = (compressible ? HySt->p0[j] : mpv->HydroState->p0[j]);            
            rhoY = (compressible ? HySt->rhoY0[j] : mpv->HydroState->rhoY0[j]);
            rho  = rhoY/Y[j];
            
            Sol->rho[n]    = rho;
            Sol->rhou[n]   = rho * u;
            Sol->rhov[n]   = rho * v;
            Sol->rhow[n]   = rho * w;
            Sol->rhoe[n]   = rhoe(rho, u, v, w, p, g*y);
            Sol->rhoY[n]   = rhoY;
            Sol->geopot[n] = g * y;
            
            mpv->p2_cells[n]   = HySt->p20[j];
            Sol->rhoX[BUOY][n] = Sol->rho[n] * ( Sol->rho[n]/Sol->rhoY[n] - mpv->HydroState->S0[j]);
        
            /* nodal pressure */
            nn   = j*icxn+i;
            mpv->p2_nodes[nn] = HyStn->p20[j];
        }
    }
    
    /* Pressure distribution for the hydrostatic model */
    Hydrostatic_Initial_Pressure(Sol, mpv, elem, node);
    
    /* set all dummy cells */
    /* geopotential in bottom and top dummy cells */
    for(j = 0; j < igy; j++) {m = j * icx;
        y = elem->y[j];
        for(i = 0; i < icx; i++) {n = m + i;
            Sol->geopot[n] = g * y;
        }
    }
    
    for(j = icy-igy; j < icy; j++) {m = j * icx;
        y = elem->y[j];
        for(i = 0; i < icx; i++) {n = m + i;
            Sol->geopot[n] = g * y;
        }
    }
}


/* ================================================================================== */

double pressure_function(double r, double p0, double S0, double u_theta, double Msq, double Gamma){
    return pow((pow(p0,Gamma) + Gamma*S0*Msq*u_theta*u_theta*(1.0 - pow((1.0-r),5.0)*(5.0*r+1.0))/30.0), (1.0/Gamma));
}

/* ================================================================================== */

double rho_function(double psi){
    
    double rhomax = 1.0;
    double percent = 0.0;
    double delta   = 0.1;
    
    double smooth, cosine, sgncos;
    
    cosine = cos(PI*(fabs(psi) - (0.5 - 0.5*delta)) / delta);
    sgncos = SIGNnull(cosine);
    
    smooth = 0.5*(1.0 - sgncos*sqrt(sqrt(sqrt(sgncos*cosine))) );
    
    /* return(1.0 - 0.2 * (psi*psi)); */
    /* return(rhomax * (1.0 - percent * (fabs(psi) < 0.5 ? 0.0 : 1.0 )) ); */
    return(rhomax * (1.0 - percent * (fabs(psi) < 0.5 - 0.5*delta ?
                                      0.0 :
                                      (fabs(psi) > 0.5 + 0.5*delta ?
                                       1.0 :
                                       smooth
                                       )
                                      )
                     )
           );
}

/* ================================================================================== */

double stratification(
                      double y) {
    
    extern User_Data ud;
    
    double Nsq = ud.Nsq_ref * ud.t_ref * ud.t_ref;
    double g   = ud.gravity_strength[1] / ud.Msq;
    
    return( exp(Nsq*y/g) );
}


/* ================================================================================== */

double molly(
             double x) {
    
    extern User_Data ud;
    
    const double del  = 0.25;
    const double L    = ud.xmax-ud.xmin;
    const double xi_l = MIN_own(1.0, (x-ud.xmin)/(del*L));
    const double xi_r = MIN_own(1.0, (ud.xmax-x)/(del*L));
    
    return(0.5* MIN_own(1.0-cos(PI*xi_l), 1.0-cos(PI*xi_r)));
}



/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: userdata.c,v $
 Revision 1.2  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:50:44  nicola
 Added csr.c and sod1d.c (user data for 3d code)
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
