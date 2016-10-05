/*******************************************************************************
 File:   thermodynamic.h
 Author: Nicola
 Date:   Thu Feb 26 16:17:40 WET 1998 
 *******************************************************************************/
#ifndef THERMODYNAMIC_H
#define THERMODYNAMIC_H

#include "userdata.h"
#include "variable.h"
#include "space_discretization.h"

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
typedef struct {
	
	double gamm;
	double gamminv;
	double gm1;
	double gm1inv;
	double Gamma;
	double Gammainv;
    
    double Rg_over_Rv;
    double Q;
	
} Thermodynamic;


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void Thermodynamic_init(Thermodynamic* th, const User_Data* ud);

void set_Y_to_unity(ConsVars* Sol, ElemSpaceDiscr* elem);



#endif /* THERMODYNAMIC_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: thermodynamic.h,v $
 Revision 1.1  1998/03/01 18:43:36  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
