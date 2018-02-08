/*******************************************************************************
 File:   second_projection.h
 Author: Nicola
 Date:   Mon Mar 30 13:09:15 MET DST 1998
 *******************************************************************************/
#ifndef SECOND_PROJECTION_H
#define SECOND_PROJECTION_H

#include "Common.h"
#include "variable.h"


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/

void second_projection(
                       ConsVars* Sol,
                       MPV* mpv,
                       const ConsVars* Sol0,
                       const ElemSpaceDiscr* elem,
                       const NodeSpaceDiscr* node,
                       const double p_update,
                       const double t,
                       const double dt);

#ifdef GRAVITY_IMPLICIT
void euler_gravity(ConsVars* Sol,
                   VectorField* buoy,
                   const MPV* mpv,
                   const ElemSpaceDiscr* elem,
                   const enum GravityTimeIntegrator GTI,
                   const double dt);
#endif


#endif /* SECOND_PROJECTION_H */


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
