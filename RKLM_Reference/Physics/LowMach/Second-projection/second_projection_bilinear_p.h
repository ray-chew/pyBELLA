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

void euler_backward_non_advective_impl_part(
                       ConsVars* Sol,
                       MPV* mpv,
                       const ConsVars* Sol0,
                       const ElemSpaceDiscr* elem,
                       const NodeSpaceDiscr* node,
                       const double t,
                       const double dt);

void momentum_increments(MPV* mpv, 
                         const ConsVars *Sol, 
                         const ElemSpaceDiscr *elem);

void euler_backward_non_advective_expl_part(ConsVars* Sol,
                            const MPV* mpv,
                            const ElemSpaceDiscr* elem,
                            const double dt);

void euler_forward_non_advective(ConsVars* Sol,
                                 MPV* mpv,
                                 const ConsVars* Sol0,
                                 const ElemSpaceDiscr* elem,
                                 const NodeSpaceDiscr* node,
                                 const double dt,
                                 const enum EXPLICIT_PRESSURE with_pressure);

void pressure_gradient_forces(
                              double* force[3], 
                              const ConsVars *Sol, 
                              const MPV *mpv, 
                              const ElemSpaceDiscr *elem, 
                              const NodeSpaceDiscr *node);

#endif /* SECOND_PROJECTION_H */


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
