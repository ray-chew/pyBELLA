/*******************************************************************************
 File:   variable.c
 Author: Nicola
 Date:   Thu Feb 26 11:32:41 CET 1998
 *******************************************************************************/
#include "variable.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "userdata.h"
#include "math_own.h"
#include "Common.h"
/* #include "space_discretization.h" */

#if OUTPUT_SUBSTEPS /* 5 */
#include "io.h"
#endif

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
ConsVars* ConsVars_new(const int size) {
    extern User_Data ud;
	ConsVars* var = (ConsVars*)malloc(sizeof(ConsVars));
	var->rho  = (double*)malloc(size * sizeof(double));
	var->rhou = (double*)malloc(size * sizeof(double));
	var->rhov = (double*)malloc(size * sizeof(double));
	var->rhow = (double*)malloc(size * sizeof(double));
	var->rhoe = (double*)malloc(size * sizeof(double));
	var->rhoY = (double*)malloc(size * sizeof(double));
    for (int nsp=0; nsp<ud.nspec; nsp++) {
        var->rhoX[nsp] = (double*)malloc(size * sizeof(double));
    }
	return var;
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_free(ConsVars* var) {
    extern User_Data ud;
    int nsp;
	free(var->rho);
	free(var->rhou);
	free(var->rhov);
	free(var->rhow);
	free(var->rhoe);
	free(var->rhoY);
    for (nsp=0; nsp<ud.nspec; nsp++) {
        free(var->rhoX[nsp]);
    }
	free(var); 
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_setp(ConsVars* obj, const ConsVars* src, const int i) {
    extern User_Data ud;
    int nsp;
	obj->rho  = &src->rho[i];
	obj->rhou = &src->rhou[i];
	obj->rhov = &src->rhov[i];
	obj->rhow = &src->rhow[i];
	obj->rhoe = &src->rhoe[i];
	obj->rhoY = &src->rhoY[i];
    for (nsp=0; nsp<ud.nspec; nsp++) {
        obj->rhoX[nsp] = &src->rhoX[nsp][i];
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_addp(ConsVars* obj, const int n) {
    extern User_Data ud;
	obj->rho  += n;
	obj->rhou += n;
	obj->rhov += n;
	obj->rhow += n;
	obj->rhoe += n;
	obj->rhoY += n;
    for (int nsp=0; nsp<ud.nspec; nsp++) {
        obj->rhoX[nsp] += n;
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_setzero(ConsVars* obj, const int n) {
    extern User_Data ud;
	int i, nsp;
	for(i = 0; i < n; i++) {
		obj->rho[i]  = 0.0;
		obj->rhou[i] = 0.0;
		obj->rhov[i] = 0.0;
		obj->rhow[i] = 0.0;
		obj->rhoe[i] = 0.0;
		obj->rhoY[i] = 0.0;
        for (nsp=0; nsp<ud.nspec; nsp++) {
            obj->rhoX[nsp][i] = 0.0;
        }
	}
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_add(ConsVars* obj,
                  const ConsVars* src,
                  const int n)
{
    extern User_Data ud;
    int i, nsp;
    for(i = 0; i < n; i++) {
        obj->rho[i]  += src->rho[i] ;
        obj->rhou[i] += src->rhou[i];
        obj->rhov[i] += src->rhov[i];
        obj->rhow[i] += src->rhow[i];
        obj->rhoe[i] += src->rhoe[i];
        obj->rhoY[i] += src->rhoY[i];
        for (nsp=0; nsp<ud.nspec; nsp++) {
            obj->rhoX[nsp][i] += src->rhoX[nsp][i];
        }
    }    
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_average(ConsVars* obj, 
                      const ConsVars* src, 
                      const int n)
{
    extern User_Data ud;
    int i, nsp;
    for(i = 0; i < n; i++) {
        obj->rho[i]  = 0.5 * (obj->rho[i]  + src->rho[i] );
        obj->rhou[i] = 0.5 * (obj->rhou[i] + src->rhou[i]);
        obj->rhov[i] = 0.5 * (obj->rhov[i] + src->rhov[i]);
        obj->rhow[i] = 0.5 * (obj->rhow[i] + src->rhow[i]);
        obj->rhoe[i] = 0.5 * (obj->rhoe[i] + src->rhoe[i]);
        obj->rhoY[i] = 0.5 * (obj->rhoY[i] + src->rhoY[i]);
        for (nsp=0; nsp<ud.nspec; nsp++) {
            obj->rhoX[nsp][i] = 0.5 * (obj->rhoX[nsp][i] + src->rhoX[nsp][i]);
        }
    }    
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_multiply_by_scalar(ConsVars* obj,
                                 const double factor,
                                 const int n)
{
    extern User_Data ud;
    int i, nsp;
    for(i = 0; i < n; i++) {
        obj->rho[i]  *= factor;
        obj->rhou[i] *= factor;
        obj->rhov[i] *= factor;
        obj->rhow[i] *= factor;
        obj->rhoe[i] *= factor;
        obj->rhoY[i] *= factor;
        for (nsp=0; nsp<ud.nspec; nsp++) {
            obj->rhoX[nsp][i] *= factor;
        }
    }    
}

/* void memcpy(double *a,  double *b,  int n_length); */

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ConsVars_set(ConsVars* obj, const ConsVars* src, const int n) {
    extern User_Data ud;
	int nsp;
	memcpy(obj->rho,  src->rho,  n * sizeof(double));
	memcpy(obj->rhou, src->rhou, n * sizeof(double));
	memcpy(obj->rhov, src->rhov, n * sizeof(double));
	memcpy(obj->rhow, src->rhow, n * sizeof(double));
	memcpy(obj->rhoe, src->rhoe, n * sizeof(double));
	memcpy(obj->rhoY, src->rhoY, n * sizeof(double));
    for (nsp=0; nsp<ud.nspec; nsp++) {
        memcpy(obj->rhoX[nsp], src->rhoX[nsp], n * sizeof(double));
    }
    
#if OUTPUT_SUBSTEPS /* 5 */
    extern User_Data ud;
    extern ElemSpaceDiscr* elem;
    extern NodeSpaceDiscr* node;
    putout(obj, ud.file_name, "Sol", elem, node, 1);
#endif

}


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void VectorField_setzero(VectorField* obj, const int n) {

    for(int i = 0; i < n; i++) {
        obj->x[i] = 0.0;
        obj->y[i] = 0.0;
        obj->z[i] = 0.0;
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
VectorField* VectorField_new(const int size) {
	VectorField* var = (VectorField*)malloc(sizeof(VectorField));
	var->x = (double*)malloc(size * sizeof(double));
	var->y = (double*)malloc(size * sizeof(double));
	var->z = (double*)malloc(size * sizeof(double));
	return var;
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void VectorField_free(VectorField* var) {      
	free(var->x);
	free(var->y);
	free(var->z);
	free(var); 
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void VectorField_set(VectorField* obj, const VectorField* src, const int n, const int ndim) {
    memcpy(obj->x, src->x, n * sizeof(double));
    if (ndim > 1) memcpy(obj->y, src->y, n * sizeof(double));
    if (ndim > 2) memcpy(obj->z, src->z, n * sizeof(double));
}




/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void scalar_set(double* obj, const double* src, const int n){
	memcpy(obj,  src,  n * sizeof(double));
}

void scalar_add(double* obj, const double* src, const int nstart, const int nend){
	int i;
	for(i=nstart; i<nend; i++) obj[i] += src[i];
}




/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
double L1_norm(double* u, const int nstart, const int nend){
	double norm = 0.0;
	int i;
	for(i=nstart; i<nend; i++) norm += ABS(u[i]);
	return(norm);
}


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
States* States_new(const int size) {
    extern User_Data ud;
	States* var = (States*)malloc(sizeof(States));
	var->rho  = (double*)malloc(size * sizeof(double));
	var->rhou = (double*)malloc(size * sizeof(double));
	var->rhov = (double*)malloc(size * sizeof(double));
	var->rhow = (double*)malloc(size * sizeof(double));
	var->rhoe = (double*)malloc(size * sizeof(double));
	var->rhoY = (double*)malloc(size * sizeof(double));
    for (int nsp=0; nsp<ud.nspec; nsp++) {
        var->rhoX[nsp] = (double*)malloc(size * sizeof(double));
    }
	var->u    = (double*)malloc(size * sizeof(double));
	var->v    = (double*)malloc(size * sizeof(double));
	var->w    = (double*)malloc(size * sizeof(double));
	var->q    = (double*)malloc(size * sizeof(double));
	var->p    = (double*)malloc(size * sizeof(double));
	var->c    = (double*)malloc(size * sizeof(double));
	var->entro= (double*)malloc(size * sizeof(double));
	var->H    = (double*)malloc(size * sizeof(double));
	var->Y    = (double*)malloc(size * sizeof(double));
    for (int nsp=0; nsp<ud.nspec; nsp++) {
        var->X[nsp] = (double*)malloc(size * sizeof(double));
    }
	var->p0   = (double*)malloc(size * sizeof(double));
    var->pi0  = (double*)malloc(size * sizeof(double));
	var->p20  = (double*)malloc(size * sizeof(double));
	var->rho0 = (double*)malloc(size * sizeof(double));
    var->rhoY0 = (double*)malloc(size * sizeof(double));
	var->S0   = (double*)malloc(size * sizeof(double));
	var->S10  = (double*)malloc(size * sizeof(double));
	var->Y0  = (double*)malloc(size * sizeof(double));
	return var;
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void States_free(States* var) {     
    extern User_Data ud;
	free(var->rho);
	free(var->rhou);
	free(var->rhov);
	free(var->rhow);
	free(var->rhoe);
	free(var->rhoY);
    for (int nsp=0; nsp<ud.nspec; nsp++) {
        free(var->rhoX[nsp]);
    }
	free(var->u);
	free(var->v);
	free(var->w);
	free(var->q);
	free(var->p);
	free(var->c);
	free(var->entro);
	free(var->H);
	free(var->Y);
    for (int nsp=0; nsp<ud.nspec; nsp++) {
        free(var->X[nsp]);
    }
	free(var->p0);
    free(var->pi0);
	free(var->p20);
	free(var->rho0);
	free(var->S0);
	free(var->S10);
	free(var->Y0);
	free(var);
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
States* States_small_new(const int size) {
    extern User_Data ud;
	int nsp;
	States* var = (States*)malloc(sizeof(States));
	var->u     = (double*)malloc(size * sizeof(double));
	var->v     = (double*)malloc(size * sizeof(double));
	var->w     = (double*)malloc(size * sizeof(double));
	var->q     = (double*)malloc(size * sizeof(double));
	var->p     = (double*)malloc(size * sizeof(double));
	var->entro = (double*)malloc(size * sizeof(double));
	var->H     = (double*)malloc(size * sizeof(double));
	var->c     = (double*)malloc(size * sizeof(double));
	var->Y     = (double*)malloc(size * sizeof(double));
    for (nsp=0; nsp<ud.nspec; nsp++) {
        var->X[nsp] = (double*)malloc(size * sizeof(double));
    }
	var->p0    = (double*)malloc(size * sizeof(double));
	var->p20   = (double*)malloc(size * sizeof(double));
	var->rho0  = (double*)malloc(size * sizeof(double));
	var->S0    = (double*)malloc(size * sizeof(double));
	var->S10   = (double*)malloc(size * sizeof(double));
	var->Y0   = (double*)malloc(size * sizeof(double));
	return var;
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void States_small_free(States* var) {      
    extern User_Data ud;
	int nsp;
	free(var->u);
	free(var->v);
	free(var->w);
	free(var->q);
	free(var->p);
	free(var->entro);
	free(var->H);
	free(var->c);
	free(var->Y);
    for (nsp=0; nsp<ud.nspec; nsp++) {
        free(var->X[nsp]);
    }
	free(var->p0);
	free(var->p20);
	free(var->rho0);
	free(var->S0);
	free(var->S10);
	free(var->Y0);
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void States_setp(States* obj, const ConsVars* src, const int i) {
    extern User_Data ud;
	obj->rho  = &src->rho[i];
	obj->rhou = &src->rhou[i];
	obj->rhov = &src->rhov[i];
	obj->rhow = &src->rhow[i];
	obj->rhoe = &src->rhoe[i];
	obj->rhoY = &src->rhoY[i];
    for (int nsp=0; nsp<ud.nspec; nsp++) {
        obj->rhoX[nsp] = &src->rhoX[nsp][i];
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void States_addp(States* obj, const int n) {
    extern User_Data ud;
	obj->rho  += n;
	obj->rhou += n;
	obj->rhov += n;
	obj->rhow += n;
	obj->rhoe += n;
	obj->rhoY += n;
    for (int nsp=0; nsp<ud.nspec; nsp++) {
        obj->rhoX[nsp] += n;
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
Characters* Characters_new(const int size) {
    extern User_Data ud;
	int nsp;
	Characters* var = (Characters*)malloc(sizeof(Characters));
	var->plus  = (double*)malloc(size * sizeof(double));
	var->entro = (double*)malloc(size * sizeof(double));
	var->minus = (double*)malloc(size * sizeof(double));
	var->v     = (double*)malloc(size * sizeof(double));
	var->w     = (double*)malloc(size * sizeof(double));
	var->Y     = (double*)malloc(size * sizeof(double));
    for (nsp=0; nsp<ud.nspec; nsp++) {
        var->X[nsp] = (double*)malloc(size * sizeof(double));
    }
	return var;
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void Characters_free(Characters* var) {   
    extern User_Data ud;
	int nsp;
	free(var->plus);
	free(var->entro);
	free(var->minus);
	free(var->v);
	free(var->w);
	free(var->Y);
    for (nsp=0; nsp<ud.nspec; nsp++) {
        free(var->X[nsp]);
    }
}


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
Hydro* Hydro_new(const int size) {
	Hydro* var = (Hydro*)malloc(size * sizeof(Hydro));
	return var;
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void Hydro_free(Hydro* var) {   
	free(var);
}


/* LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL */
