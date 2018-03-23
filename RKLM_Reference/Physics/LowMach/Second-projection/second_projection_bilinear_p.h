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
                       const double t,
                       const double dt);

void euler_backward_gravity(ConsVars* Sol,
                            const MPV* mpv,
                            const ElemSpaceDiscr* elem,
                            const double dt);

void euler_forward_non_advective(ConsVars* Sol,
                                 const MPV* mpv,
                                 const ElemSpaceDiscr* elem,
                                 const NodeSpaceDiscr* node,
                                 const double dt,
                                 const enum EXPLICIT_PRESSURE with_pressure);

#endif /* SECOND_PROJECTION_H */


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
