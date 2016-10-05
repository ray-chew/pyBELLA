/*******************************************************************************
 File:   variable.h
 Author: Nicola
 Date:   Thu Feb 26 11:29:23 CET 1998
 *******************************************************************************/
#ifndef VARIABLE_H
#define VARIABLE_H

#include "Common.h"
#include "kgrid.h"

struct ConsVar;


/*------------------------------------------------------------------------------
 ConsVars
 ------------------------------------------------------------------------------*/
typedef struct {
	double* rho;
	double* rhou;
	double* rhov;
	double* rhow;
	double* rhoe;
	double* rhoY;
	double* rhoZ;
    double* rhoX[NSPEC];
	double* geopot;
} ConsVars;

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
ConsVars* ConsVars_new(const int size);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_free(ConsVars* var);


/*------------------------------------------------------------------------------
 3D Vector
 ------------------------------------------------------------------------------*/
typedef struct {
    double x;
    double y;
    double z;
} Vector;

/*------------------------------------------------------------------------------
 VectorField
 ------------------------------------------------------------------------------*/
typedef struct {
	double* x;
	double* y;
	double* z;
} VectorField;


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void VectorField_setzero(VectorField* obj, const int n);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
VectorField* VectorField_new(const int size);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void VectorField_free(VectorField* var);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void VectorField_set(VectorField* obj, 
                     const VectorField* src, 
                     const int n, 
                     const int ndim);


/*------------------------------------------------------------------------------
 Sets the pointers of obj to the n component of the arrays of src
 ------------------------------------------------------------------------------*/
void ConsVars_setp(ConsVars* obj, const ConsVars* src, const int i);


/*------------------------------------------------------------------------------
 Adds n to the pointers of obj 
 ------------------------------------------------------------------------------*/
void ConsVars_addp(ConsVars* obj, const int n);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_setzero(ConsVars* obj, const int n);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_set(ConsVars* obj, const ConsVars* src, const int n);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void scalar_set(double* obj, const double* src, const int n);
void scalar_add(double* obj, const double* src, const int nstart, const int nend);

double L1_norm(double* u, const int nstart, const int nend);


/*------------------------------------------------------------------------------
 States 
 ------------------------------------------------------------------------------*/
typedef struct {
	double* rho;
	double* rhou;
	double* rhov;
	double* rhow;
	double* rhoe;
	double* rhoY;
	double* rhoZ;
    double* rhoX[NSPEC];
	double* geopot;
	double* u;
	double* v;
	double* w;
	double* q;
	double* p;
	double* c;
	double* entro;
	double* H;
	double* Y;
	double* Z;
    double* X[NSPEC];
	double* p0;
    double* pi0;
	double* p20;
	double* rho0;
	double* S0;
	double* S10;
	double* Y0;
} States;

/*------------------------------------------------------------------------------
 Speeds 
 ------------------------------------------------------------------------------*/
typedef struct {
	double u;
	double u_plus_c;
} Speeds;


/*------------------------------------------------------------------------------
 TimeStepInfo
 ------------------------------------------------------------------------------*/
typedef struct {
    Vector flow_speed;
    Vector flow_and_sound_speed;
    double cfl;
    double cfl_ac;
    double cfl_adv;
    double cfl_gravity;
    double time_step;
    int    time_step_switch;
} TimeStepInfo;

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
States* States_new(const int size);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void States_free(States* var);


/*------------------------------------------------------------------------------
 Allocates memory only for u, v, w, q, p, c, Y, Z
 ------------------------------------------------------------------------------*/
States* States_small_new(const int size);


/*------------------------------------------------------------------------------
 Releases memory only for u, v, w, q, p, c, Y, Z 
 ------------------------------------------------------------------------------*/
void States_small_free(States* var);


/*------------------------------------------------------------------------------
 Sets the pointers of obj to the n component of the arrays of src 
 ------------------------------------------------------------------------------*/
void States_setp(States* obj, const ConsVars* src, const int i);


/*------------------------------------------------------------------------------
 Adds n to the pointers of obj (only rho, rhou, rhov, rhow, rhoe, rhoY, rhoZ)
 ------------------------------------------------------------------------------*/
void States_addp(States* obj, const int n);


/*------------------------------------------------------------------------------
 Characters
 ------------------------------------------------------------------------------*/
typedef struct {
	double* plus;
	double* entro;
	double* minus;
	double* v;
	double* w;
	double* Y;
	double* Z;
    double* X[NSPEC];
} Characters;


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void States_HydroState(
					   States *Solk, 
					   const States *HydroState, 
					   const ElemSpaceDiscr *elem, 
					   const int nstart_Solk,
					   const int nmax_Solk,
					   const int nstart_full_field, 
					   const int SplitStep);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
Characters* Characters_new(const int size);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void Characters_free(Characters* var);





/*------------------------------------------------------------------------------
 Characters
 ------------------------------------------------------------------------------*/

typedef struct {
    double Y[5];
    double Yinv[5];
    double p[5];
    double p2[5];
    double p2c[5];
    double rho[5];
    double rhoY[5];
} Hydro;


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
Hydro* Hydro_new(const int size);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void Hydro_free(Hydro* var);






#endif /* VARIABLE_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: variable.h,v $
 Revision 1.2  1998/03/07 09:56:48  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:37  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
