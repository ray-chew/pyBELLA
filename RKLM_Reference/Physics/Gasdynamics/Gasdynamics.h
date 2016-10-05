/*******************************************************************************
 File:   Gasdynamics.h
 Author: Nicola
 Date:   Wed Feb 25 07:52:29 CET 1998
 *******************************************************************************/
#ifndef GASDYNAMICS_H
#define GASDYNAMICS_H

#include "variable.h"


/*------------------------------------------------------------------------------
 Computes the maxmum signal speed for evaluation of the CFL-condition
 ------------------------------------------------------------------------------*/
double maxspeed(const ConsVars* Sol, const int n);

Speeds maxspeeds(const ConsVars* Sol, const int n);

TimeStepInfo dynamic_timestep(
                              const ConsVars* Sol,
                              const double time,
                              const double time_output,
                              const ElemSpaceDiscr* elem,
                              const int step);

#endif /* GASDYNAMICS_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: Gasdynamics.h,v $
 Revision 1.1  1998/03/01 18:43:32  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
