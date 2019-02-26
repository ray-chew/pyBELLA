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
#include "boundary.h"
#include "memory.h"
#include "Hydrostatics.h"

/* mollification function used in the Skamarock-Klemp 1994 tests for theta perturbations */
double molly(double x);

/* horizontal stretch for S&K94 IGWs: planetary -> 160.0;  long-wave -> 20.0;  standard -> 1.0; */
static double scalefactor = 1.0;   

void User_Data_init(User_Data* ud) {
    
    int i; 
    
    /* ================================================================================== */
    /* =====  PROBLEM SET UP  =========================================================== */
    /* ================================================================================== */
    
    /* Earth */
    double grav     = 0.0;             /* gravitational acceleration [m/s^2]    */
    double omega    = 0.0*0.0001;       /* Coriolis parameter [1/s]              */
                                        /* sin(0.5*PI) * 2.0 * 0.00007272205217; */
    
    /* thermodynamics and chemistry */
    double R_gas    = 287.0;            /* [J/kg/K]               */
    double R_vap    = 461.00;           /* [J/kg/K]               */
    double Q_vap    = 2.53e+06;         /* [J]                    */
    double gamma    = 2.0;              /* dimensionless; 5.0/3.0       */

    double viscm  = 0.0;            /* [m^2/s]                         */
    double viscbm = 0.0;             /* [m^2/s]                         */
    double visct  = 0.0;             /* [m^2/s]                         */
    double viscbt = 0.0;             /* [m^2/s]                         */
    double cond   = 0.0;            /* [m^2/s]                         */

    /* references for non-dimensionalization */
    double h_ref    = 1.0;                   /* [m]               */
    double t_ref    = 1.0;                   /* [s]               */
    double T_ref    = 353.048780488;                /* [K]               */
    double p_ref    = 101325;                /* [Pa]              */
    double u_ref    = h_ref/t_ref;           /* [m/s]; Sr = 1     */
    double rho_ref  = p_ref / (R_gas*T_ref); /* [kg/m^3]          */
    
    /* reference stratification as (buoyancy frequency)^2 */
    double Nsq_ref  = 0.0e-4;                /* [1/s^2]           */
    
    ud->h_ref       = h_ref;
    ud->t_ref       = t_ref;
    ud->T_ref       = T_ref;
    ud->p_ref       = p_ref;
    ud->u_ref       = u_ref;
    ud->rho_ref     = rho_ref;
    ud->Nsq_ref     = Nsq_ref;
    ud->g_ref       = grav;
    ud->gamm        = gamma;
    ud->Rg_over_Rv  = R_gas/R_vap;
    ud->Q           = Q_vap/(R_gas*T_ref);
    
    /* number of advected species */
    ud->nspec       = NSPEC;

    /*FULL_MOLECULAR_TRANSPORT, STRAKA_DIFFUSION_MODEL, NO_MOLECULAR_TRANSPORT */
    ud->mol_trans   = NO_MOLECULAR_TRANSPORT; 
    ud->viscm       = viscm  * t_ref/(h_ref*h_ref);
    ud->viscbm      = viscbm * t_ref/(h_ref*h_ref);
    ud->visct       = visct  * t_ref/(h_ref*h_ref);
    ud->viscbt      = viscbt * t_ref/(h_ref*h_ref);
    ud->cond        = cond * t_ref/(h_ref*h_ref*R_gas);

    /* Low Mach */
    ud->is_nonhydrostatic = 1;    /* 0: hydrostatic;  1: nonhydrostatic;  -1: transition (see nonhydrostasy()) */
    ud->is_compressible   = 1;    /* 0: psinc;  1: comp;  -1: psinc-comp-transition (see compressibility()) */
    ud->acoustic_timestep = 0;    /* advective time step -> 0;  acoustic time step -> 1; */
    ud->Msq =  u_ref*u_ref / (R_gas*T_ref);
    
    /* geo-stuff */
    for(i=0; i<3; i++) {
        ud->gravity_strength[i]  = 0.0;  /* corresponds to  M^2 / Fr^2  =  g href / R Tref */
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
        
    /* flow domain; all lengths in units of  href  */
    ud->xmin =   0.0 * scalefactor;
    ud->xmax =   1.0 * scalefactor;
    ud->ymin =   0.0 * scalefactor;
    ud->ymax =   1.0 * scalefactor;
    ud->zmin = - 1.0;
    ud->zmax =   1.0;
    
    /* boundary/initial conditions */
    ud->wind_speed  =  1.0/u_ref;
    ud->wind_shear  = -0.0;              /* velocity in [u_ref/h_ref] */
    ud->hill_shape  = AGNESI;            /* AGNESI, SCHLUTOW */
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
    ud->time_integrator       = SI_MIDPT; /* this code version has only one option */
    ud->advec_time_integrator = STRANG; /* HEUN; EXPL_MIDPT;   best tested: STRANG; */
    ud->CFL                   = 77.0; /* 0.45; 0.9; 0.8; */
    ud->dtfixed0              = 0.0000668205;
    ud->dtfixed               = 0.0000668205;
    
    set_time_integrator_parameters(ud);
    
    /* Grid and space discretization */
    ud->inx =  256+1; /* 641; 321; 161; 129; 81; */
    ud->iny =  10+1; /* 321; 161;  81;  65; 41;  */
    ud->inz =      1;
    
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
        
    /* al explicit predictor operations are done on ncache-size data worms to save memory */ 
    ud->ncache = 175; /* 71+4; 304*44; 604*44; (ud->inx+3); (ud->inx+3)*(ud->iny+3);*/
    
    /* linear solver-stuff */
    double tol                            = 1.e-8;
    ud->flux_correction_precision         = tol;
    ud->flux_correction_local_precision   = tol;    /* 1.e-05 should be enough */
    ud->second_projection_precision       = tol;
    ud->second_projection_local_precision = tol;  /* 1.e-05 should be enough */
    ud->flux_correction_max_iterations    = 6000;
    ud->second_projection_max_iterations  = 6000;
    ud->initial_projection                = WRONG;   /* WRONG;  CORRECT; */
    ud->initial_impl_Euler                = WRONG;   /* WRONG;  CORRECT; */
    
    ud->column_preconditioner             = CORRECT; /* WRONG; CORRECT; */
    ud->synchronize_nodal_pressure        = WRONG;   /* WRONG; CORRECT; */
    ud->synchronize_weight                = 0.0;    /* relevant only when prev. option is "CORRECT"
                                                     Should ultimately be a function of dt . */  
    
    /* numerics parameters */
    ud->eps_Machine = sqrt(DBL_EPSILON);
    
    /* ================================================================================== */
    /* =====  CODE FLOW CONTROL  ======================================================== */
    /* ================================================================================== */
    
    // double t_period = sqrt(ud->Msq)*(ud->xmax-ud->xmin)/sqrt(ud->gamm);
    ud->tout[0] = 0.00267282;
    ud->tout[1] = -1.0;
    // ud->tout[2] = 3.0*t_period;
    // ud->tout[3] = 4.0*t_period;
    // ud->tout[4] = 5.0*t_period;
    // ud->tout[6] = -1.0;

    ud->stepmax = 40;
    
    ud->write_stdout = ON;
    ud->write_stdout_period = 1;
    ud->write_file = ON;
    ud->write_file_period = 40;
    ud->file_format = HDF;
    
    {
#ifdef RUPERT
        char *OutputBaseFolder      = "/Users/rupert/Documents/Computation/RKLM_Reference/";
#else
        char *OutputBaseFolder      = "/home/tommaso/work/repos/RKLM_Reference/";
#endif
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

void Sol_initial(ConsVars* Sol, 
                 ConsVars* Sol0, 
                 MPV* mpv,
                 BDRY* bdry,
                 const ElemSpaceDiscr* elem,
                 const NodeSpaceDiscr* node) {
    
    extern Thermodynamic th;
    extern User_Data ud;
        
    /* acoustic wave test */
    const double u0  = ud.wind_speed;
    const double v0  = 0.0;
    const double w0  = 0.0;
    const double del = 0.05;                          /* perturbation amplitude; standard:  0.01 */
    const double xc  = 0.0;                           /* center of perturbation */
    const double a   = 0.125 * (ud.xmax - ud.xmin);    /* characteristic width of perturbation */
    const double wn  = 2.0*PI/((ud.xmax - ud.xmin));
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
    const int igy = elem->igy;
    
    const int icxn = node->icx;
    const int icyn = node->icy;
    
    const double Ma = sqrt(ud.Msq);
    
    int i, j, n, nn;
    double x, y;
    double rho, u, v, w, p, rhoY, c;
            
    /* plain acoustic wave test */
        
    /* store background hydrostatic state in mpv auxiliary struct */
    Hydrostatics_State(mpv, elem, node);
    
    /* computations for the vertical slice at  k=0 */
    for(i = 0; i < icx; i++) {
        
        /* initialization of field variables */
        for(j = igy; j < icy - igy + 1; j++) {
            
            n  = j*icx + i;
                        
            x   = elem->x[i];
            y   = elem->y[j];

            u   = u0;
            v   = v0;
            w   = w0;
            
            /* p    = mpv->HydroState->p0[j] * (1.0 + del * molly(x) / (1.0 + (x-xc)*(x-xc) / (a*a))); */         
            p    = mpv->HydroState->p0[j] * pow(1.0 + del * sin(wn*x), 2.0*th.gamm*th.gm1inv);       
            rhoY = pow(p,th.gamminv);
            rho  = rhoY;
            c    = sqrt(th.gamm * p / rho);
            u   += (p - mpv->HydroState->p0[j])/(rho*c)/Ma;
            
            Sol->rho[n]    = rho;
            Sol->rhou[n]   = rho * u;
            Sol->rhov[n]   = rho * v;
            Sol->rhow[n]   = rho * w;
            Sol->rhoe[n]   = rhoe(rho, u, v, w, p);
            Sol->rhoY[n]   = rhoY;
            
            mpv->p2_cells[n]   = (pow(p,th.Gamma) - 1.0) / ud.Msq;
            Sol->rhoX[BUOY][n] = Sol->rho[n] * (Sol->rho[n]/Sol->rhoY[n] - mpv->HydroState->S0[j]);
        
            /* nodal pressure */
            x  = node->x[i];
            nn = j*icxn+i;
            
            p  = mpv->HydroState_n->p0[j] * pow(1.0 + del * sin(wn*x), 2.0*th.gamm*th.gm1inv);
            mpv->p2_nodes[nn] = (pow(p,th.Gamma) - 1.0) / ud.Msq;
        }
    }
    
    ud.initial_projection = WRONG;
    
    /* Copy the data just generated into all the parallel vertical slices */
    for (int j=0; j<icy; j++) {
        int mc = j*icx;
        int mn = j*icxn;
        for (int i=0; i<icx; i++) {
            int nc = mc + i;
            int nn = mn + i;
            for (int k=1; k<icz; k++) {
                int lc = nc + k*icx*icy;
                int ln = nn + k*icxn*icyn;
                
                Sol->rho[lc]    = Sol->rho[nc] ;
                Sol->rhou[lc]   = Sol->rhou[nc];
                Sol->rhov[lc]   = Sol->rhov[nc];
                Sol->rhow[lc]   = Sol->rhow[nc];
                Sol->rhoe[lc]   = Sol->rhoe[nc];
                Sol->rhoY[lc]   = Sol->rhoY[nc];
                
                Sol->rhoX[BUOY][lc] = Sol->rhoX[BUOY][nc];
                mpv->p2_cells[lc]   = mpv->p2_cells[nc];
                mpv->p2_nodes[ln]   = mpv->p2_nodes[nn];
            }
        }
    }

    set_ghostcells_p2(mpv->p2_cells, elem, elem->igx);
    set_ghostnodes_p2(mpv->p2_nodes, node, 2);   
    
    ud.nonhydrostasy   = nonhydrostasy(0);
    ud.compressibility = compressibility(0);
    
    //set_wall_massflux(bdry, Sol, elem);
    set_wall_rhoYflux(bdry, Sol, mpv, elem);
    Set_Explicit_Boundary_Data(Sol, elem);
    
    ConsVars_set(Sol0, Sol, elem->nc);
    
    /* the initial projection should ensure the velocity field is discretely
     divergence-controlled when sound-free initial data are required.
     This can mean vanishing divergence in a zero-Mach flow or the 
     pseudo-incompressible divergence constraint in an atmospheric flow setting
     */ 
    if (ud.initial_projection == CORRECT) {
        /* initial velocity div-control not implemented for the test case yet. */
        assert(0);
    }

}

/* ================================================================================== */

double stratification(
                      double y) {
        
    return( 1.0 );
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
