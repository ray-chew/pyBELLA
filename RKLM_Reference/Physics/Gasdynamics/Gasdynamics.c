/*******************************************************************************
 File:   Gasdynamics.c
 Author: Rupert 
 Date:   
 *******************************************************************************/
#include <math.h>
#include "Common.h"
#include "math_own.h"
#include "Eos.h"
#include "thermodynamic.h"
#include "variable.h"
/* #include "space_discretization.h" */#include "enumerator.h"
#include "mpv.h"
#include "error.h"
#include "Gasdynamics.h"

/* -------------------------------------------------------------------------- */

double maxspeed(const ConsVars* Sol, const int n)
{
	extern User_Data ud;
	extern ElemSpaceDiscr* elem;
	extern Thermodynamic th;
	
	const double gamm = th.gamm;
	const int igx = elem->igx;
	const int igy = elem->igy;
	const int igz = elem->igz;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	
	const double Minv = 1.0 / ud.Msq;
	
	double invrho, c, umax, amax; 
	int i, j, k, nk, njk, nijk;
		
	amax = umax = 0.0;
	
	for(k = igz; k < icz - igz; k++) {
		nk = k * icx * icy;
		for(j = igy; j < icy - igy; j++) {
			njk = nk + j * icx;
			for(i = igx; i < icx - igx; i++) {
				nijk = njk + i; 
				
				invrho = 1.0 / Sol->rho[nijk]; 
                double p = pow(Sol->rhoY[nijk],th.gamm);
				c = sqrt(gamm * p * invrho) * Minv;
				umax = ABS(Sol->rhou[nijk] * invrho);
				umax = MAX_own(umax, ABS(Sol->rhov[nijk] * invrho));
				umax = MAX_own(umax, ABS(Sol->rhow[nijk] * invrho));
				amax = MAX_own(amax, umax + c);
			}
		}
	}  
	
	return(amax);
}


/* -------------------------------------------------------------------------- */

Speeds maxspeeds(const ConsVars* Sol, const int n)
{
	extern User_Data ud;
	extern ElemSpaceDiscr* elem;
	extern Thermodynamic th;
	
	const double gamm = th.gamm;
	const int igx = elem->igx;
	const int igy = elem->igy;
	const int igz = elem->igz;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	
	const double Minv = 1.0 / sqrt(ud.Msq);
	
	double invrho, c, umax; 
	int i, j, k, nk, njk, nijk;
		
	Speeds speeds;
	
	speeds.u = 0.0;
	speeds.u_plus_c = 0.0;
	
	for(k = igz; k < icz - igz; k++) {
		nk = k * icx * icy;
		for(j = igy; j < icy - igy; j++) {
			njk = nk + j * icx;
			for(i = igx; i < icx - igx; i++) {
				nijk = njk + i; 
				
				invrho = 1.0 / Sol->rho[nijk]; 
                double p = pow(Sol->rhoY[nijk],th.gamm);
				c = sqrt(gamm * p * invrho) * Minv;
				umax = ABS(Sol->rhou[nijk] * invrho);
				umax = MAX_own(umax, ABS(Sol->rhov[nijk] * invrho));
				umax = MAX_own(umax, ABS(Sol->rhow[nijk] * invrho));
				speeds.u = MAX_own(speeds.u, umax);
				speeds.u_plus_c = MAX_own(speeds.u_plus_c, umax + c);
			}
		}
	}  
	
	return(speeds);
}

/* -------------------------------------------------------------------------- */

void dynamic_timestep(TimeStepInfo* TSI,
                      MPV* mpv,
                      const ConsVars* Sol,
                      const double time,
                      const double time_output,
                      const ElemSpaceDiscr* elem,
                      const int step) {

    extern User_Data ud;
    extern Thermodynamic th;
    
    const double gamm = th.gamm;
    const int igx = elem->igx;
    const int igy = elem->igy;
    const int igz = elem->igz;
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;

    const double Minv = 1.0 / sqrt(ud.Msq);
    const double CFL  = ud.CFL;
    
    double umax = ud.eps_Machine;
    double vmax = ud.eps_Machine;
    double wmax = ud.eps_Machine;

    double upcmax = ud.eps_Machine;
    double vpcmax = ud.eps_Machine;
    double wpcmax = ud.eps_Machine;

    double dt;
    double cfl, cfl_ac, cfl_adv;
    
    if (TSI->time_step_switch) {
        TSI->time_step_switch = 1;
        return;
    };
    
    for(int k = igz; k < icz - igz; k++) {
        int nk = k * icx * icy;
        for(int j = igy; j < icy - igy; j++) {
            int njk = nk + j * icx;
            for(int i = igx; i < icx - igx; i++) {
                int nijk = njk + i;
                
                double invrho = 1.0 / Sol->rho[nijk];
                double p = pow(Sol->rhoY[nijk],th.gamm);
                double c = sqrt(gamm * p * invrho) * Minv;
                double u = ABS(Sol->rhou[nijk] * invrho);
                double v = ABS(Sol->rhov[nijk] * invrho);
                double w = ABS(Sol->rhow[nijk] * invrho);
                
                umax = MAX_own(umax, u);
                vmax = MAX_own(vmax, v);
                wmax = MAX_own(wmax, w);

                upcmax = MAX_own(upcmax, u+c);
                vpcmax = MAX_own(vpcmax, v+c);
                wpcmax = MAX_own(wpcmax, w+c);
            }
        }
    }  

    if (ud.acoustic_timestep == 1) {
        double dtx = CFL * elem->dx / upcmax;
        double dty = CFL * elem->dy / vpcmax;
        double dtz = CFL * elem->dz / wpcmax;
        double dt_cfl;
        
        dt_cfl = dtx;
        dt_cfl = MIN_own(dt_cfl, dty);
        dt_cfl = MIN_own(dt_cfl, dtz);
        
        dt  = MIN_own(dt_cfl, ud.dtfixed0 + MIN_own(step, 1) * (ud.dtfixed - ud.dtfixed0));
        
        if (2.0*dt > time_output - time) {
            dt = 0.5 * (time_output - time) + ud.eps_Machine;
            TSI->time_step_switch = 1;
        }
        
        cfl = cfl_ac = CFL * dt / dt_cfl;
        
        cfl_adv = dt * umax / elem->dx;
        cfl_adv = MAX_own(cfl_adv, dt * vmax / elem->dy);
        cfl_adv = MAX_own(cfl_adv, dt * wmax / elem->dz);

    } else {
        double dtx = CFL * elem->dx / umax;
        double dty = CFL * elem->dy / vmax;
        double dtz = CFL * elem->dz / wmax;
        double dt_cfl;

        dt_cfl  = dtx;
        dt_cfl  = MIN_own(dt_cfl, dty);
        dt_cfl  = MIN_own(dt_cfl, dtz);
        
        dt  = MIN_own(dt_cfl, ud.dtfixed0 + MIN_own(step, 1) * (ud.dtfixed - ud.dtfixed0));
        
        dt *= MIN_own((float)(step+1), 1.0);
        
        if (2.0*dt > time_output - time) {
            dt = 0.5 * (time_output - time) + ud.eps_Machine;
            TSI->time_step_switch = 1;
        }

        cfl = cfl_adv = CFL * dt / dt_cfl;
        
        cfl_ac = dt * upcmax / elem->dx;
        cfl_ac = MAX_own(cfl_ac, dt * vpcmax / elem->dy);
        cfl_ac = MAX_own(cfl_ac, dt * wpcmax / elem->dz);
    }
    
    TSI->flow_speed.x           = umax;
    TSI->flow_speed.y           = vmax;
    TSI->flow_speed.z           = wmax;
    TSI->flow_and_sound_speed.x = upcmax;
    TSI->flow_and_sound_speed.y = vpcmax;
    TSI->flow_and_sound_speed.z = wpcmax;
    TSI->cfl                    = cfl;
    TSI->cfl_ac                 = cfl_ac;
    TSI->cfl_adv                = cfl_adv;
    TSI->cfl_gravity            = 99999999;
    TSI->time_step              = dt;
    mpv->dt                     = dt;
}


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: Gasdynamics.c,v $
 Revision 1.1  1998/03/01 18:43:32  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
