/*******************************************************************************
 File:   Eos.h
 Author: Nicola 
 Date:   Tue Feb 24 10:45:02 CET 1998
 *******************************************************************************/
#ifndef EOS_H
#define EOS_H

#include <stdio.h> 
#include "Common.h"
#include "enumerator.h"
#include "variable.h"
#include "mpv.h"


/*------------------------------------------------------------------------------
 primitive from conservative variables
 ------------------------------------------------------------------------------*/
void primitives(
				States* U, 
				const int nstart, 
				const int nende);

/*------------------------------------------------------------------------------
 conservative from primitive variables
 ------------------------------------------------------------------------------*/
void conservatives(
				   States* U, 
				   const int nstart, 
				   const int nende);

/*------------------------------------------------------------------------------
 conservative from entropy variables
 ------------------------------------------------------------------------------*/
void conservatives_from_entro_vars(
								   States* U, 
								   const int nstart, 
								   const int nende);

/*------------------------------------------------------------------------------
 conservative from primitive variables (rho, u, v, w, p)
 ------------------------------------------------------------------------------*/
void conservatives_from_primitives(
								   States* U, 
								   const int nstart, 
								   const int nende);

/*------------------------------------------------------------------------------
 conservative from (rho(Y), u, v, w, Y, Z);
 ------------------------------------------------------------------------------*/
void conservatives_from_uvwYZ(
							  States* U, 
							  const int nstart, 
							  const int nende);

/*------------------------------------------------------------------------------
 pressure from conservative variables
 ------------------------------------------------------------------------------*/
void pressure(
			  double *p, 
			  const ConsVars* U, 
			  const int nstart_p, 
			  const int nstart, 
			  const int nende);



/*------------------------------------------------------------------------------
 pressure increment computed in the first projection
 ------------------------------------------------------------------------------*/
void dp2_first_projection(
                          double *dp2, 
                          const ConsVars* U, 
                          const MPV *mpv,
                          const ElemSpaceDiscr *elem);

/*------------------------------------------------------------------------------
 deviation of Exner pressure from background state values
 ------------------------------------------------------------------------------*/
void dp_exner(
              double *dpi, 
              const ConsVars* U, 
              const MPV *mpv,
              const ElemSpaceDiscr *elem);

/*------------------------------------------------------------------------------
 deviation of thermodynamic pressure from background state (dimensional)
 ------------------------------------------------------------------------------*/
void dpress_dim(
              double *dpdim, 
              const ConsVars* U, 
              const MPV *mpv,
              const ElemSpaceDiscr *elem);

/*------------------------------------------------------------------------------
 inverse of potential temperature from conservative variables
 ------------------------------------------------------------------------------*/
void pot_temp_inv(
				  double *S, 
				  const ConsVars* U, 
				  const int nstart_p, 
				  const int nstart, 
				  const int nende);

/*------------------------------------------------------------------------------
 pressure from asymptotic expressions
 ------------------------------------------------------------------------------*/
void pressure_asymptotic(
						 double *p, 
						 const ConsVars* U, 
						 const MPV* mpv,
						 const ElemSpaceDiscr* elem);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void velox(double *u, const ConsVars* U, const int nstart, const int nende );


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void veloy(double *v, const ConsVars* U, const int nstart, const int nende );


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void veloz(double *w, const ConsVars* U, const int nstart, const int nende );


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/

double T_from_p_rho(const double p, const double rho);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
double T_from_p_theta(const double p, const double theta);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void delta_rhoY(double *drhoY, const ConsVars* U, const int nstart, const int nende);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void MassFrac_Y(
				double *Y, 
				const ConsVars* U, 
				const int nstart, 
				const int nende);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void MassFrac_Y_perturbation(
				double *Y, 
				const ConsVars* U, 
				const int nstart, 
				const int nende);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void MassFrac_Z(
				double *Z, 
				const ConsVars* U, 
				const int nstart, 
				const int nende );

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void MassFrac_q(
				double *Z, 
				const ConsVars* U, 
				const int nstart, 
				const int nende, 
                const int which_q);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void specific_enthalpy(
					   double* h, 
					   const ConsVars* U, 
					   const int nstart, 
					   const int nende);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
double rhoe(
			const double rho, 
			const double u, 
			const double v, 
			const double w,
			const double p,
			const double geopot);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void adjust_pi(ConsVars* Sol,
               MPV* mpv,
			   const ConsVars* Sol0,
			   const ElemSpaceDiscr* elem, 
               const double weight);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void reset_Y_perturbation(ConsVars* Sol,
                          const MPV* mpv,
                          const ElemSpaceDiscr* elem);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/

double qv_sat_from_p_T(const double p, const double T);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void fluctuation(double *var, const double *var0, const ElemSpaceDiscr *elem);
    
/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
double compressibility(const double t);

#endif /* EOS_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: Eos.h,v $
 Revision 1.1  1998/03/01 18:43:30  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
