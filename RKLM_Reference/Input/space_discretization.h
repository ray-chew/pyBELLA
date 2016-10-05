/*******************************************************************************
 File:   space_discretization.h
 Author: Nicola,						Rupert
 Date:   Fri Feb 27 19:42:36 CET 1998, Feb. 2004
 *******************************************************************************/
#ifndef SPACE_DISCRETIZATION_H
#define SPACE_DISCRETIZATION_H

#include "kgrid.h"
#include "enumerator.h"
#include "variable.h"


void ElemSpaceDiscr_ghost(
						  ConsVars* Sol, 
						  const ElemSpaceDiscr* elem, 
						  const int ig);

void ElemSpaceDiscr_ghost_rho(
							  double* rho, 
							  const ElemSpaceDiscr* elem, 
							  const int ig);

void ElemSpaceDiscr_ghost_rhou(
							   double* rhou, 
							   const ElemSpaceDiscr* elem, 
							   const int ig);

void ElemSpaceDiscr_ghost_rhov(
							   double* rhov, 
							   const ElemSpaceDiscr* elem, 
							   const int ig);

void ElemSpaceDiscr_ghost_rhow(
							   double* rhow, 
							   const ElemSpaceDiscr* elem, 
							   const int ig);

void ElemSpaceDiscr_ghost_rhoe(
							   double* rhoe, 
							   const ElemSpaceDiscr* elem, 
							   const int ig);

void ElemSpaceDiscr_ghost_rhoY(
							   double* rhoY, 
							   const ElemSpaceDiscr* elem, 
							   const int ig);

void ElemSpaceDiscr_ghost_rhoZ(
							   double* rhoZ, 
							   const ElemSpaceDiscr* elem, 
							   const int ig);

/* ================================================================== */

void NodeSpaceDiscr_ghost(
						  ConsVars* Sol, 
						  const NodeSpaceDiscr* node, 
						  const int ig);

void NodeSpaceDiscr_ghost_rho(
							  double* rho, 
							  const NodeSpaceDiscr* node, 
							  const int ig);

void NodeSpaceDiscr_ghost_rhou(
							   double* rhou, 
							   const NodeSpaceDiscr* node, 
							   const int ig);

void NodeSpaceDiscr_ghost_rhov(
							   double* rhov, 
							   const NodeSpaceDiscr* node, 
							   const int ig);

void NodeSpaceDiscr_ghost_rhow(
							   double* rhow, 
							   const NodeSpaceDiscr* node, 
							   const int ig);

void NodeSpaceDiscr_ghost_rhoe(
							   double* rhoe, 
							   const NodeSpaceDiscr* node, 
							   const int ig);

void NodeSpaceDiscr_ghost_rhoY(
							   double* rhoY, 
							   const NodeSpaceDiscr* node, 
							   const int ig);

void NodeSpaceDiscr_ghost_rhoZ(
							   double* rhoZ, 
							   const NodeSpaceDiscr* node, 
							   const int ig);

#endif /* SPACE_DISCRETIZATION_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: discretization.h,v $
 Revision 1.2  1998/03/07 09:56:45  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:33  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
