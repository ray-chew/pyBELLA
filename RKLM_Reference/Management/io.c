/*******************************************************************************
 File:   io.c
 Author: Nicola 
 Date:   Wed Feb 18 08:53:35 CET 1998
 
 Notice: the functions WriteASCII, WriteSILO and WriteGNU should be re-designed 
 according to the following guidelines (priority according to order):
 
 - avoid the use of global structures
 
 - use the state equations made available in Eos (instead of implementing 
 their own version of the state equation!)
 
 - grid and solution in distinct files!
 
 - a different file for each variable (as in Rupert's WriteHDF)? 
 
 - merge WriteASCII and WriteGNU (why two different routines?)
 *******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdarg.h>
#include "Common.h"
#include "math_own.h"
#include "error.h"
#include "warning.h"
#include "kgrid.h"
#include "Eos.h"
#include "io.h"
#include "userdata.h"
#include "enumerator.h"
#include "thermodynamic.h"
#include "variable.h"
/* #include "space_discretization.h" */
#include "set_ghostcells_p.h"
#include "set_ghostnodes_p.h"
#include "boundary.h"
#include "mpv.h"
#include "memory.h"

#ifdef HDFFORMAT
#ifdef __MWERKS__
#define MAC
#undef __WINDOWS__
#endif
#include  "/opt/local/include/dfsd.h" 
#endif
#ifdef SILOFORMAT
#include "silo.h"
#endif

static void (*rotate[])(ConsVars* Sol, double* rhs, double *Yinvbg, const enum Direction dir) = {NULL, rotate2D, rotate3D};

#ifdef SILOFORMAT
static void putoutSILO(char* file_name);
#endif

/* pointers to files */
static FILE *prhofile      = NULL;   
static FILE *prhoefile     = NULL;   
static FILE *prhoYfile     = NULL;   
static FILE *pdrhoYfile     = NULL;   
static FILE *pufile        = NULL;     
static FILE *pvfile        = NULL;     
static FILE *pwfile        = NULL;
static FILE *ppfile        = NULL;     
static FILE *pSfile        = NULL;     
static FILE *pYfile        = NULL;     
static FILE *pdYfile       = NULL;     
static FILE *pZfile        = NULL;     
static FILE *pqfile        = NULL;
static FILE *pp2file       = NULL;     
static FILE *dpdimfile     = NULL;     
static FILE *pgeopotfile   = NULL;     

void putout(
			ConsVars* Sol, 
			const double t, 
			const double tout, 
			const int step,
			const int SplitStep, 
			char* dir_name, 
			char* field_name,
            const int writeout) {
	
	/* User data */
	extern User_Data ud;
	
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	
	/* Arrays */
	extern MPV* mpv;
	extern double *W0, *Yinvbg;

#ifdef MATLAB_OUTPUT
	static int output_counter = 0;
#endif

	const int ndim = elem->ndim;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int nc = elem->nc; 
		
	double *var;
	char fn[200], fieldname[90], step_string[30];
	int nsp;
	
	switch(ud.file_format) {
		case ASCII: {
            assert(0);
			sprintf(fn, "%s.ascii.%d", ud.file_name, step);
			/*
            if(ud.write_stdout == ON) printf("writing %s ...\n", fn);   
			putoutASCII(Sol, mpv, elem, fn);
			if(ndim == 3) {
                const int igz = elem->igz;
                int k;
				for(k = igz; k < icz - igz; k++) {
					sprintf(fn, "%s.%d.ascii.%d", ud.file_name, k, step);
                    WriteSliceASCII(Sol, mpv, elem, fn, k); 
				}
			}
             */
			break;
		}
		case HDF: {
			int i;
			
			/* rotate forward and set boundary data */          
			for(i = 0; i < ndim; i++) { 
				const double lambda = 1.0;
				Bound(Sol, mpv->HydroState, lambda, nc, SplitStep+i); 
				if(i < ndim - 1) (*rotate[ndim - 1])(Sol, mpv->Level[0]->rhs, Yinvbg, FORWARD);
			}         
			/* rotate back */          
			for(i = ndim-1; i > 0; i--) {
				(*rotate[ndim - 1])(Sol, mpv->Level[0]->rhs, Yinvbg, BACKWARD);
			}
            
            if (writeout == 0) {
                return;
            }
			
#ifdef MATLAB_OUTPUT			
			if(output_counter<10) {
				sprintf(step_string, "00%d", output_counter);
			}
			else if(output_counter<100) {
				sprintf(step_string, "0%d", output_counter);
			}
			else {
				sprintf(step_string, "%d", output_counter);
			}
            output_counter++;
#else
			if(step<10) {
				sprintf(step_string, "000%d", step);
			}
			else if(step<100) {
				sprintf(step_string, "00%d", step);
			}
			else if(step<1000) {
				sprintf(step_string, "0%d", step);
			}
			else {
				sprintf(step_string, "%d", step);
			}
#endif
			var = W0;
			
			if(ud.write_stdout == ON) printf("\n");
			
			/* geopotential */
			sprintf(fn, "%s/geopot/geopot_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON) printf("writing %s ...\n", fn);
			sprintf(fieldname, "geopot_%s_%s", field_name, step_string);
			WriteHDF(pgeopotfile, icx, icy, icz, ndim, Sol->geopot, fn, fieldname);
			
			/* density */
			sprintf(fn, "%s/rho/rho_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON) printf("writing %s ...\n", fn);
			sprintf(fieldname, "rho_%s_%s", field_name, step_string);
			WriteHDF(prhofile, icx, icy, icz, ndim, Sol->rho, fn, fieldname);
			
			/* energy density */
			sprintf(fn, "%s/rhoe/rhoe_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON) printf("writing %s ...\n", fn);
			sprintf(fieldname, "rhoe_%s_%s", field_name, step_string);
			WriteHDF(prhoefile, icx, icy, icz, ndim, Sol->rhoe, fn, fieldname);
			
			/* u velocity component */
			velox(var, Sol, 0, nc);
			sprintf(fn, "%s/u/u_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "u_%s_%s", field_name, step_string);
			WriteHDF(pufile, icx, icy, icz, ndim, var, fn, fieldname);
			
			/* v velocity component */
			veloy(var, Sol, 0, nc);
			sprintf(fn, "%s/v/v_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "v_%s_%s", field_name, step_string);
			WriteHDF(pvfile, icx, icy, icz, ndim, var, fn, fieldname);
			
			/* w velocity component */
			 veloz(var, Sol, 0, nc);
			 sprintf(fn, "%s/w/w_%s.hdf", dir_name, step_string);
			 if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			 sprintf(fieldname, "w_%s_%s", field_name, step_string);
			 WriteHDF(pwfile, icx, icy, icz, ndim, var, fn, fieldname);
			
			/* pressure */
			pressure(var, Sol, 0, 0, nc);
			sprintf(fn, "%s/p/p_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "p_%s_%s", field_name, step_string);
			WriteHDF(ppfile, icx, icy, icz, ndim, var, fn, fieldname);

            /* pressure difference */
			dpress_dim(var, Sol, mpv, elem);
			sprintf(fn, "%s/dpdim/dpdim_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "dpdim_%s_%s", field_name, step_string);
			WriteHDF(dpdimfile, icx, icy, icz, ndim, var, fn, fieldname);
            
			/* inverse of potential temperature */
			/* pot_temp_inv_change(var, Sol, 0, 0, nc); */
			pot_temp_inv(var, Sol, 0, 0, nc);
			sprintf(fn, "%s/S/S_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "S_%s_%s", field_name, step_string);
			WriteHDF(pSfile, icx, icy, icz, ndim, var, fn, fieldname);
			
			/* rho Theta */
			sprintf(fn, "%s/rhoY/rhoY_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON) printf("writing %s ...\n", fn);
			sprintf(fieldname, "rhoY_%s_%s", field_name, step_string);
			WriteHDF(prhoYfile, icx, icy, icz, ndim, Sol->rhoY, fn, fieldname);
			
            /* rho Theta - (rho Theta)_bg */
            delta_rhoY(var, Sol, 0, nc);
            sprintf(fn, "%s/drhoY/drhoY_%s.hdf", dir_name, step_string);
            if(ud.write_stdout == ON) printf("writing %s ...\n", fn);
            sprintf(fieldname, "drhoY_%s_%s", field_name, step_string);
            WriteHDF(pdrhoYfile, icx, icy, icz, ndim, var, fn, fieldname);
            
            /* Species mass fraction */
			MassFrac_Y(var, Sol, 0, nc);
			sprintf(fn, "%s/Y/Y_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "Y_%s_%s", field_name, step_string);
			WriteHDF(pYfile, icx, icy, icz, ndim, var, fn, fieldname);
			
			/* Species mass fraction */
			MassFrac_Y_perturbation(var, Sol, 0, nc);
			sprintf(fn, "%s/dY/dY_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "dY_%s_%s", field_name, step_string);
			WriteHDF(pdYfile, icx, icy, icz, ndim, var, fn, fieldname);
			
            /* rho Z */
            sprintf(fn, "%s/rhoZ/rhoZ_%s.hdf", dir_name, step_string);
            if(ud.write_stdout == ON) printf("writing %s ...\n", fn);
            sprintf(fieldname, "rhoZ_%s_%s", field_name, step_string);
            WriteHDF(prhoYfile, icx, icy, icz, ndim, Sol->rhoZ, fn, fieldname);

            /* auxiliary Scalar Z */
			MassFrac_Z(var, Sol, 0, nc);
			sprintf(fn, "%s/Z/Z_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "Z_%s_%s", field_name, step_string);
			WriteHDF(pZfile, icx, icy, icz, ndim, var, fn, fieldname);
            
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                
                switch (nsp) {
                    case QV:
                        sprintf(fn, "%s/qv/qv_%s.hdf", dir_name, step_string);
                        sprintf(fieldname, "qv_%s_%s", field_name, step_string);
                        break;
                    case QC:
                        sprintf(fn, "%s/qc/qc_%s.hdf", dir_name, step_string);
                        sprintf(fieldname, "qc_%s_%s", field_name, step_string);
                        break;
                    case QR:
                        sprintf(fn, "%s/qr/qr_%s.hdf", dir_name, step_string);
                        sprintf(fieldname, "qr_%s_%s", field_name, step_string);
                        break;
                    default:
                        break;
                }
                
                MassFrac_q(var, Sol, 0, nc, nsp);
                if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
                WriteHDF(pqfile, icx, icy, icz, ndim, var, fn, fieldname);
            }

			sprintf(fn, "%s/dp2_nodes/dp2_n_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "dp2_nodes_%s_%s", field_name, step_string);
			WriteHDF(pp2file, 
					 mpv->Level[0]->node->icx, 
					 mpv->Level[0]->node->icy, 
					 mpv->Level[0]->node->icz, 
					 mpv->Level[0]->node->ndim, 
					 mpv->dp2_nodes,   
					 fn, 
					 fieldname);

            sprintf(fn, "%s/p2_nodes/p2_n_%s.hdf", dir_name, step_string);
            if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
            sprintf(fieldname, "p2_nodes_%s_%s", field_name, step_string);
            WriteHDF(pp2file,
                     mpv->Level[0]->node->icx,
                     mpv->Level[0]->node->icy,
                     mpv->Level[0]->node->icz,
                     mpv->Level[0]->node->ndim,
                     mpv->p2_nodes,   
                     fn, 
                     fieldname);

            /* fluctuation(var, mpv->p2_cells, elem); */
            memcpy(var, mpv->p2_cells, elem->nc*sizeof(double));
            sprintf(fn, "%s/p2_c/p2_c_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "p2_c_%s", step_string);
			WriteHDF(pp2file, 
					 mpv->Level[0]->elem->icx, 
					 mpv->Level[0]->elem->icy, 
					 mpv->Level[0]->elem->icz, 
					 mpv->Level[0]->elem->ndim, 
					 var,
					 fn, 
					 fieldname);

            /* dp_exner(var, Sol, mpv, elem); */
            dp2_first_projection(var, Sol, mpv, elem);
			sprintf(fn, "%s/dp2_c/dp2_c_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
            sprintf(fieldname, "dp2_c_%s", step_string);
			WriteHDF(pp2file, 
					 mpv->Level[0]->elem->icx, 
					 mpv->Level[0]->elem->icy, 
					 mpv->Level[0]->elem->icz, 
					 mpv->Level[0]->elem->ndim, 
					 var, 
					 fn, 
					 fieldname);
            
			break;
		}
		case SILO: {
            assert(0);
			sprintf(fn, "%s.silo.%d", ud.file_name, step);
			/* putoutSILO(fn); */
			break;
		}
		default: 
			ERROR("file format not available");
	}
	
}

/* ----------------------------------------------------------------------- */

void WriteTimeHistories(const ConsVars* Sol,
                               const ElemSpaceDiscr* elem,
                               const double time,
                               const int step,
                               const int first_running_last){
    
    extern User_Data ud;
    extern Thermodynamic th;
    
    static FILE* EnergyHistoryFile = NULL;
    
    static double e_total_0;
    static double e_internal_0;
    static double e_kinetic_0;
    static double e_potential_0;
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;

    const int igx = elem->igx;
    const int igy = elem->igy;
    const int igz = elem->igz;
    
    const double dV = elem->dxyz[0] * (elem->ndim > 1 ? elem->dxyz[1] : 1.0) * (elem->ndim > 2 ? elem->dxyz[2] : 1.0);
    
    double e_total_sum, e_internal_sum, e_kinetic_sum, e_potential_sum;
    
    int i, j, k, l, m, n;
    
    assert(ud.i_gravity[1]);
    
    if (first_running_last == -1) {
        fclose(EnergyHistoryFile);
        return;
    }
    
    /* collect energy components */
    e_internal_sum  = 0.0;
    e_kinetic_sum   = 0.0;
    e_potential_sum = 0.0;
    
    for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                e_internal_sum  += th.gm1inv * pow(Sol->rhoY[n], th.gamm);
                e_kinetic_sum   += 0.5*(  Sol->rhou[n]*Sol->rhou[n] 
                                       + Sol->rhov[n]*Sol->rhov[n] 
                                       + Sol->rhow[n]*Sol->rhow[n]
                                       ) / Sol->rho[n];
                e_potential_sum += ud.gravity_strength[1] * Sol->rho[n] * elem->y[j];
            }
        }
    }
    e_internal_sum  *= dV;
    e_kinetic_sum   *= dV*ud.Msq;
    e_potential_sum *= dV;
    e_total_sum      = e_internal_sum + e_kinetic_sum + e_potential_sum;

    if (first_running_last == 0) {
        fprintf(EnergyHistoryFile, "%d\t %f\t %e\t %e\t %e\t %e\n", step, time, e_internal_sum - e_internal_0, e_kinetic_sum - e_kinetic_0, e_potential_sum - e_potential_0, e_total_sum - e_total_0);
    } else if (first_running_last == 1) {
        EnergyHistoryFile = fopen("EnergyHistory.txt", "w+");
        fprintf(EnergyHistoryFile, "step\t time\t de_int\t de_kin\t de_pot\t de_tot\n");
        fprintf(EnergyHistoryFile, "%d\t %f\t %e\t %e\t %e\t %e\n", step, time, 0.0, 0.0, 0.0, 0.0);
        e_total_0     = e_total_sum;
        e_internal_0  = e_internal_sum;
        e_kinetic_0   = e_kinetic_sum;
        e_potential_0 = e_potential_sum;
        return;
    }
}

/* ----------------------------------------------------------------------- */

void WriteHDF(
			  FILE* pfile, 
			  int rows, 
			  int cols, 
			  int layers, 
			  int ndim,
			  double* Data, 
			  char* file_name, 
			  char* var_name) {
	
#ifdef HDFFORMAT
	
	/* needed libs
	 
	 libjpeg.62.dylib
	 libdf.a
	 libmfhdf.a
	 libz.dylib
	 */
	
	float  *image, *pimage;
	float  pmax, pmin;
	double *pData;
	int row, col, layer;
	
	int	dims[ 3 ];
	
	dims[0] = rows;
	dims[1] = cols;
	dims[2] = layers;
	
	image = (float *)malloc( (unsigned)((rows*cols*layers)*sizeof(float)) );
	
	pimage  = image;
	pData   = Data;
	pmax    = -100000.0;
	pmin    =  100000.0;
	
	
	/* throw data on float array and find min/max */
	for ( row = 0; row < rows; row++ )
	{
		for ( col = 0; col < cols; col++ )
		{
			for ( layer = 0; layer < layers; layer++ )
			{
				*pimage = (float)(Data[layer*rows*cols+col*rows+row]);
				pmin    = MIN_own(pmin, *pimage);
				pmax    = MAX_own(pmax, *pimage);
				pimage++;
				pData++;
			}
		}
	}
	pimage = image; /* reset data pointer */
	
	/* save data as HDF file	*/
	DFSDsetdims( ndim, dims );
	DFSDsetrange( &pmax, &pmin ); 
	DFSDsetdimstrs( 1, "x", "-", "F7.4" );
	DFSDsetdimstrs( 2, "y", "-", "F7.4" );
	DFSDsetdimstrs( 3, "z", "-", "F7.4" );
	DFSDsetdatastrs( var_name,"-", "F7.4", "cartesian" );
	DFSDputdata( file_name, ndim, dims, pimage );
	free( image );
	
#else
	
	WARNING("function not available");
	
#endif
}


#if 0

void ElemSpaceDiscrWriteASCII(
							  double* var, 
							  const ElemSpaceDiscr* elem,
							  const char* filename,
							  const char* varname) {
	
	int ndim = elem->ndim;
	
	FILE* file;
	
	file = fopen(filename, "w+");
	if(!file) ERROR("cannot open file");
	
	switch(ndim) {
		case 1: {
			ERROR("function not available");
			break;
		}  
		case 2: {
			const int igx = elem->igx;
			const int igy = elem->igy;
			const int icx = elem->icx;
			const int icy = elem->icy;
			const char form[] = {"%e   %e   %e\n"};
			const char head[] = {"# x y %s\n"}; 
			int i, j, m, n;
			fprintf(file, head);
			for(j = igy; j < icy - igy; j++) {m = j * icx;  
				for(i = igx; i < icx - igx; i++) {n = m + i;	       
					fprintf(file, form, elem->x[i], elem->y[j], var[n]);
				}
				fprintf(file, "\n"); /* block separation */
			}
			fclose(file);
			break;
		}
		case 3: {
			ERROR("function not available");
			break;
		}
		default: {
			fclose(file); 
			ERROR("bad dimension");
		}
	}
}


void NodeSpaceDiscrWriteASCII(
							  double* var, 
							  const NodeSpaceDiscr* node,
							  const char* filename,
							  const char* varname) {
	
	int ndim = node->ndim;
	
	FILE* file;
	
	file = fopen(filename, "w+");
	if(!file) ERROR("cannot open file");
	
	switch(ndim) {
		case 1: {
			ERROR("function not available");
			break;
		}  
		case 2: {
			const int igx = node->igx;
			const int igy = node->igy;
			const int icx = node->icx;
			const int icy = node->icy;
			const char form[] = {"%e   %e   %e\n"};
			const char head[] = {"# x y %s\n"}; 
			int i, j, m, n;
			fprintf(file, head);
			for(j = igy; j < icy - igy; j++) {m = j * icx;  
				for(i = igx; i < icx - igx; i++) {n = m + i;	       
					fprintf(file, form, node->x[i], node->y[j], var[n]);
				}
				fprintf(file, "\n"); /* block separation */
			}
			fclose(file);
			break;
		}
		case 3: {
			ERROR("function not available");
			break;
		}
		default: {
			fclose(file); 
			ERROR("bad dimension");
		}
	}
}


static void putoutASCII(
						const ConsVars* Sol,
						const MPV* mpv,
						const ElemSpaceDiscr* elem, 
						char* filename) {
	
	extern double* W0;
	
	const int ndim = elem->ndim;
	const int nc = elem->nc;
	FILE* file;
	
	pressure(W0, Sol, 0, 0, nc);
	
	file = fopen(filename, "w+");
	if(!file) ERROR("cannot open file");
	
	switch(ndim) {
		case 1: {
			const int igx = elem->igx;
			const int icx = elem->icx;
			const char head[] = {"#  x   rho   u    p   Y   Z\n"};
			const char form[] = {"%e   %e   %e   %e   %e   %e\n"};
			int i;
			double x, rho, u, p, Y, Z;
			fprintf(file, head); 
			for(i = igx; i < icx - igx; i++) {
				x   = elem->x[i];
				rho = Sol->rho[i];
				u   = Sol->rhou[i] / rho;
				p   = W0[i];
				Y   = Sol->rhoY[i] / rho;
				Z   = Sol->rhoZ[i] / rho;
				fprintf(file, form, x, rho, u, p, Y, Z);
			}
			fclose(file);
			break;
		}  
		case 2: {
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const char head[] = {"#  x   y   rho   u   v   p   Y   Z\n"};
			const char form[] = {"%e   %e   %e   %e   %e   %e   %e   %e\n"};
			int i, j, m, n;
			double x, y, rho, u, v, p, Y, Z;
			fprintf(file, head);
			for(j = igy; j < icy - igy; j++) {m = j * icx;
				for(i = igx; i < icx - igx; i++) {n = m + i;	       
					x    = elem->x[i];
					y    = elem->y[j];
					rho  = Sol->rho[n];
					u    = Sol->rhou[n] / rho;
					v    = Sol->rhov[n] / rho;
					p    = mpv->p2_cells[n];
					/*p    = W0[n];*/ 
					Y    = Sol->rhoY[n] / rho;
					Z    = Sol->rhoZ[n] / rho;	
					fprintf(file, form, x, y, rho, u, v, p, Y, Z);
				}
				fprintf(file, "\n");
			}
			fclose(file);
			break;
		}
		case 3: {
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int igz = elem->igz;
			const int icz = elem->icz;
			const char head[] = {"#  x   y   z   rho   u    v   w   p   Y   Z\n"};
			const char form[] = {"%e   %e   %e   %e   %e   %e   %e   %e   %e   %e\n"};
			int i, j, k, l, m, n;
			double x, y, z, rho, oorho, u, v, w, p, Y, Z;
			fprintf(file, head);
			for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
				for(j = igy; j < icy - igy; j++) {m = l + j * icx;
					for(i = igx; i < icx - igx; i++) {n = m + i;
						x    = elem->x[i];
						y    = elem->y[j];
						z    = elem->z[k];
						rho  = Sol->rho[n];
						oorho = 1.0 / rho;
						u    = Sol->rhou[n] * oorho;
						v    = Sol->rhov[n] * oorho;
						w    = Sol->rhow[n] * oorho;
						p    = W0[n];
						Y    = Sol->rhoY[n] * oorho;
						Z    = Sol->rhoZ[n] * oorho;
						fprintf(file, form, x, y, z, rho, u, v, w, p, Y, Z);
					}
					fprintf(file, "\n"); 
				}
				fprintf(file, "\n");
			}
			fclose(file);
			break;
		}
		default: {
			fclose(file); 
			ERROR("bad dimension");
		}
	}
}


static void WriteSliceASCII(
							const ConsVars* Sol,
							const MPV* mpv,
							const ElemSpaceDiscr* elem, 
							char* filename,
							const int k) {
	
	extern double* W0;
	
	const int ndim = elem->ndim;
	const int nc = elem->nc;
	FILE* file;
	
	pressure(W0, Sol, 0, 0, nc);
	
	file = fopen(filename, "w+");
	if(!file) ERROR("cannot open file");
	
	switch(ndim) {
		case 1: {
			ERROR("bad dimension");
			break;
		}  
		case 2: {
			ERROR("bad dimension");
			break;
		}
		case 3: {
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const char head[] = {"#  x   y   rho   u    v   w   p   Y   Z\n"};
			const char form[] = {"%e   %e   %e   %e   %e   %e   %e   %e   %e\n"};
			int i, j, l, m, n;
			double x, y, rho, oorho, u, v, w, p, Y, Z;
			fprintf(file, head);
			l = k * icx * icy;
			for(j = igy; j < icy - igy; j++) {m = l + j * icx;
				for(i = igx; i < icx - igx; i++) {n = m + i;
					x    = elem->x[i];
					y    = elem->y[j];
					rho  = Sol->rho[n];
					oorho = 1.0 / rho;
					u    = Sol->rhou[n] * oorho;
					v    = Sol->rhov[n] * oorho;
					w    = Sol->rhow[n] * oorho;
					p    = mpv->p2_cells[n];
					Y    = Sol->rhoY[n] * oorho;
					Z    = Sol->rhoZ[n] * oorho;
					fprintf(file, form, x, y, rho, u, v, w, p, Y, Z);
				}
				fprintf(file, "\n"); 
			}
			fclose(file);
			break;
		}
		default: {
			fclose(file); 
			ERROR("bad dimension");
		}
	}
}


#ifdef SILOFORMAT
static void putoutSILO(char* file_name) {
	
#ifdef SILOFORMAT
	
	/* User data */
	extern User_Data ud;
	
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	
	/* Global variables */
	extern ConsVars* Sol;
	
	const int nws = (ix > 1 ? ix - 4 : 1) * (iy > 1 ? iy - 4 : 1) * (iz > 1 ? iz - 4 : 1);
	
	const int ix = elem->icx;
	const int iy = elem->icy;
	const int iz = elem->icz;
	const int ndim = elem->ndim;
	const double dx = elem->dx;
	
	FILE* file = NULL; 
	
	switch(ndim) {
		case 1: {
			int i;
			
			
			float* x = (float*)malloc((unsigned)((ix - 4)*sizeof(float)));
			char** pn = (char**)malloc((unsigned)(ndim*sizeof(char)));
			float** pc = (float**)malloc((unsigned)(ndim*sizeof(float)));
			int* pd = (int*)malloc((unsigned)(ndim*sizeof(int)));
			float* var = (float*)malloc((unsigned)(nws*sizeof(float)));
			
			for(i = 2; i < ix - 2; i++) x[i - 2] = xmin + (0.5 + (i - 2)) * dx;
			for(i = 2; i < ix - 2; i++) var[i - 2] = Sol->rho[i];
			pn[0] = strdup("AX");
			pc[0] = x;
			pd[0] = ix - 4;
			
			file = DBCreate(file_name, DB_CLOBBER, DB_SGI, NULL, DB_PDB);
			
			i = DBPutQuadmesh(file, "grid", pn, pc, pd, ndim, DB_FLOAT, DB_COLLINEAR, NULL);
			
			if(!DBPutQuadvar1(
							  file, 
							  "density", 
							  "grid", 
							  var, 
							  pd, 
							  ndim, 
							  NULL, 
							  0, 
							  DB_FLOAT, 
							  DB_NODECENT, 
							  NULL))
				ERROR("troubles now");
			
			DBClose(file);
			free(x);
			free(pn);
			free(pc);
			free(pd);
			free(var);
			break;
		}  
		case 2: {
			ERROR("function not available");
			break;
		}
		case 3: {
			ERROR("function not available");
			break;
		}
		default: 
			ERROR("bad dimension");
	}
	
#else
	
	ERROR("function not available");
	
#endif
}
#endif

#endif /* 0 */

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: io.c,v $
 Revision 1.2  1998/03/07 09:56:46  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:34  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/

