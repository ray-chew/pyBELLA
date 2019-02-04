/*******************************************************************************
 File:   io.c
 WriteSILO() has not been tested!
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
#include "boundary.h"
#include "mpv.h"
#include "memory.h"

#ifdef HDFFORMAT
#ifdef __MWERKS__
#define MAC
#undef __WINDOWS__
#endif
#include  "../hdf/dfsd.h" 
#endif

#ifdef SILOFORMAT
#include "silo.h"
#endif

#ifdef SILOFORMAT
static void putoutSILO(char* file_name);
#endif

/* pointers to files */
static FILE *prhofile      = NULL;   
static FILE *prhoefile     = NULL;   
static FILE *prhoYfile     = NULL;   
static FILE *pdrhoYfile    = NULL;   
static FILE *pufile        = NULL;     
static FILE *pvfile        = NULL;     
static FILE *pwfile        = NULL;
static FILE *pvortfile     = NULL;
static FILE *ppfile        = NULL;     
static FILE *pSfile        = NULL;     
static FILE *pYfile        = NULL;     
static FILE *pdYfile       = NULL;     
static FILE *pqfile        = NULL;
static FILE *pp2file       = NULL;     
static FILE *dpdimfile     = NULL;     

/* ============================================================================= */

void putout(ConsVars* Sol, 
			char* dir_name, 
			char* field_name,
            const ElemSpaceDiscr *elem,
            const NodeSpaceDiscr *node,
            const int writeout) {
	
    /* TODO:  Each call to putout() seems to increase the
     occupied memory by a little bit. I have not found the
     allocation/release mismatch yet.
     */
    
	/* User data */
	extern User_Data ud;
	    	
	/* Arrays */
	extern MPV* mpv;
	extern double *W0;
    extern enum Boolean W0_in_use;
    
	static int output_counter = 0;

	const int ndim = elem->ndim;
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int nc = elem->nc; 

    const int icxn = node->icx;
    const int icyn = node->icy;
    const int iczn = node->icz;
    const int nn = node->nc; 

	double *var;
	char fn[200], fieldname[90], step_string[30];
	int nsp;
	
	switch(ud.file_format) {
		case HDF: {
			            
            if (writeout == 0) {
                return;
            }
			
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
            
#if OUTPUT_FLUXES
            extern ConsVars* flux[3];
            extern User_Data ud;
            FILE *prhs2file = NULL;
            char fn2[120], fieldname2[90];
            sprintf(fn2, "%s/fluxes/flux_rhou_%s.hdf", dir_name, step_string);
            sprintf(fieldname2, "flux_rhou");
            WriteHDF(prhs2file, elem->ifx, elem->icy, elem->icz, elem->ndim, flux[0]->rhou, fn2, fieldname2);
            sprintf(fn2, "%s/fluxes/flux_rhov_%s.hdf", dir_name, step_string);
            sprintf(fieldname2, "flux_rhov");
            WriteHDF(prhs2file, elem->ify, elem->icx, elem->icz, elem->ndim, flux[1]->rhou, fn2, fieldname2);
#endif

            assert(W0_in_use == WRONG);
            W0_in_use = CORRECT;
            var = W0;
			
			if(ud.write_stdout == ON) printf("\n");
						
			/* density */
			sprintf(fn, "%s/rho/rho_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON) printf("writing %s ...\n", fn);
			sprintf(fieldname, "rho_%s_%s", field_name, step_string);
			WriteHDF(prhofile, icx, icy, icz, ndim, Sol->rho, fn, fieldname);
			
            /* momentum x */
            sprintf(fn, "%s/rhou/rhou_%s.hdf", dir_name, step_string);
            if(ud.write_stdout == ON) printf("writing %s ...\n", fn);
            sprintf(fieldname, "rhou_%s_%s", field_name, step_string);
            WriteHDF(prhoefile, icx, icy, icz, ndim, Sol->rhou, fn, fieldname);
            
            /* momentum y */
            sprintf(fn, "%s/rhov/rhov_%s.hdf", dir_name, step_string);
            if(ud.write_stdout == ON) printf("writing %s ...\n", fn);
            sprintf(fieldname, "rhov_%s_%s", field_name, step_string);
            WriteHDF(prhoefile, icx, icy, icz, ndim, Sol->rhov, fn, fieldname);
            
            /* momentum z */
            sprintf(fn, "%s/rhow/rhow_%s.hdf", dir_name, step_string);
            if(ud.write_stdout == ON) printf("writing %s ...\n", fn);
            sprintf(fieldname, "rhow_%s_%s", field_name, step_string);
            WriteHDF(prhoefile, icx, icy, icz, ndim, Sol->rhow, fn, fieldname);

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

            /* u-v vorticity */
            vortz(var, Sol, elem, node, 0, nc);
            sprintf(fn, "%s/vortz/vortz_%s.hdf", dir_name, step_string);
            if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
            sprintf(fieldname, "vortz_%s_%s", field_name, step_string);
            WriteHDF(pvortfile, icxn, icyn, iczn, ndim, var, fn, fieldname);

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
            
            /* temperature */
            temperature(var, Sol, 0, 0, nc);
            sprintf(fn, "%s/T/T_%s.hdf", dir_name, step_string);
            if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
            sprintf(fieldname, "T_%s_%s", field_name, step_string);
            WriteHDF(ppfile, icx, icy, icz, ndim, var, fn, fieldname);
            
            /* temperature difference */
            dtemperature(var, Sol, mpv, elem);
            sprintf(fn, "%s/dT/dT_%s.hdf", dir_name, step_string);
            if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
            sprintf(fieldname, "dT_%s_%s", field_name, step_string);
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
			
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                
                switch (nsp) {
                    case BUOY:
                        sprintf(fn, "%s/buoy/buoy_%s.hdf", dir_name, step_string);
                        sprintf(fieldname, "buoy_%s_%s", field_name, step_string);
                        break;
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
			WriteHDF(pp2file, icxn, icyn, iczn, ndim, mpv->dp2_nodes, fn, fieldname);

            sprintf(fn, "%s/p2_nodes/p2_n_%s.hdf", dir_name, step_string);
            if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
            sprintf(fieldname, "p2_nodes_%s_%s", field_name, step_string);
            WriteHDF(pp2file,icxn, icyn, iczn, ndim, mpv->p2_nodes, fn, fieldname);

            /* fluctuation(var, mpv->p2_cells, elem); */
            memcpy(var, mpv->p2_cells, elem->nc*sizeof(double));
            sprintf(fn, "%s/p2_c/p2_c_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
			sprintf(fieldname, "p2_c_%s", step_string);
			WriteHDF(pp2file, icx, icy, icz, ndim, var, fn, fieldname);

            /* dp_exner(var, Sol, mpv, elem); */
            dp2_first_projection(var, Sol, mpv, elem);
			sprintf(fn, "%s/dp2_c/dp2_c_%s.hdf", dir_name, step_string);
			if(ud.write_stdout == ON ) printf("writing %s ...\n", fn);
            sprintf(fieldname, "dp2_c_%s", step_string);
			WriteHDF(pp2file, icx, icy, icz, ndim, var, fn, fieldname);
            
			break;
		}
		default: 
			ERROR("file format not available");
	}
	
    W0_in_use = WRONG;
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

#ifdef SILOFORMAT
static void putoutSILO(char* file_name) {
		
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
}
#endif



/* ============================================================================= */
static ConsVars *time_series;

void initialize_time_series(void)
{
    extern User_Data ud;
    
    time_series = ConsVars_new(ud.n_time_series);
}

/* ============================================================================= */

void store_time_series_entry(const ConsVars *Sol,
                             const ElemSpaceDiscr *elem,
                             const int step)
{
    extern User_Data ud;
    
    int i_ts  = 2*elem->icx/3;
    int j_ts  = 2*elem->icy/3; 
    int n_ts = j_ts*elem->icx + i_ts;
    int its  = MIN_own(ud.n_time_series-1, step);
    time_series->rho[its]  = Sol->rho[n_ts];
    time_series->rhou[its] = Sol->rhou[n_ts];
    time_series->rhov[its] = Sol->rhov[n_ts];
    time_series->rhow[its] = Sol->rhow[n_ts];
    time_series->rhoe[its] = Sol->rhoe[n_ts];
    time_series->rhoY[its] = Sol->rhoY[n_ts];
    for (int nsp = 0; nsp < ud.nspec; nsp++) {
        time_series->rhoX[nsp][its] = Sol->rhoX[nsp][n_ts];
    }
}

/* ============================================================================= */

void close_time_series()
{
    extern User_Data ud;
    
    char fn[200];
    sprintf(fn, "%s/time_series.txt", ud.file_name);
    FILE *TSfile = fopen(fn, "w+");
    fprintf(TSfile, "rho_ts rhou_ts rhov_ts rhow_ts rhoe_ts rhoY_ts \n");
    for(int i=0; i<ud.n_time_series; i++) {
        fprintf(TSfile, "%e %e %e %e %e %e \n",  
                time_series->rho[i], 
                time_series->rhou[i], 
                time_series->rhov[i], 
                time_series->rhow[i], 
                time_series->rhoe[i], 
                time_series->rhoY[i]-time_series->rhoY[0]);
    }
    fclose(TSfile);
    ConsVars_free(time_series);
}

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: io.c,v $
 Revision 1.2  1998/03/07 09:56:46  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:34  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/

