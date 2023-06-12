/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: manage_print.c
* BRANCH NAME: main
* 
* CONTRIBUTORS: Shuaitao WANG, Lauriane VILMIN, Aurélien BORDET, Masihullah HASANYAR, Thomas ROMARY, Nicolas FLIPO
*
* PROJECT MANAGER: Nicolas FLIPO
* 
* SOFTWARE BRIEF DESCRIPTION: The ProSe-PA is a software for simulating the hydro-biogeochemical functioning of rivers, particularly heavily urbanised rivers, and streams. The sofware can
* operate in two modes: direct calculation or data assimilation.
*
* In direct calculation mode, based on a semi-implicit Eulerian numerical scheme, the software simulates the functioning of the water column in contact with a benthic compartment made up of unconsolidated
* sediments and periphyton (librive library). It can be used to simulate the anthropisation of environments, through the explicit representation of developments such as navigation dams, sluice
* gates and river navigation, as well as discharges into the environment, such as those from wastewater treatment plants or combined sewer overflows.
* The software explicitly simulates the growth of micro-organisms in the water column and in the benthic compartment, enabling the carbon, oxygen and nutrient (nitrogen, phosphrus, silica) cycles
* associated with these biological processes to be quantified (librive library). Water temperature is also simulated by the software (libseb library), as are particulate and dissolved exchanges
* between the water column and the benthic compartment. The software can simulate 1D, pseudo-2D hydraulics of river and streams (discharge, water height) using the libhyd library. The advection-dispersion 
* process is simulated using libttc library.
* 
* In data assimilation mode, ProSe-PA includes two filters for assimilating high frequency dissolved oxygen data. These two filters are a particle filter and the ensemble Kalman filter.
*
* ANSI C software developed at the Geosciences and geoengineering Department, joint research center of Mines Paris-PSL and ARMINES, Fontainebleau, France. The code is based on the coupling 
* of 12 libraries developed also at the Geosciences and geoengineering Department, mostly in ANSI C: libprint, libts, libpc, libchronos, libio, libhyd, libtube, libttc, librive, libseb, libmb, scripts.
*
* CITATION: 
* Wang, S., Flipo, N., Romary, T.. (2019). Oxygen data assimilation for estimating micro-organism communities' parameters in river systems. Water Research, 165, 115021. doi:10.1016/j.watres.2019.115021
* Flipo, N., Even, S., Poulin, M., Tusseau-Vuillemin, M-H., Ameziane, T., Dauta, A. (2004). Biogeochemical modelling at the river scale: plankton and periphyton dynamics (Grand Morin case study, France).
*    Ecol. Model. 176(3-4), 333-347. doi:10.1016/j.ecolmodel.2004.01.012. 
* Even S, Poulin M, Garnier J, Billen G, Servais P, Chesterikoff A (1998). River ecosystem modelling: application of the PROSE model to the Seine River (France). Hydrobiologia, 373, pp. 27-45.
*    doi: 10.1023/A:1017045522336
* Vilmin, L, Aissa, N., Garnier, J., Billen, G., Mouchel, J-M., Poulin, M., Flipo, N. (2015). Impact of hydro-sedimentary processes on the dynamics of soluble reactive phosphorus in the Seine River.
*    Biogeochemistry, 122, 229-251. doi:10.1007/s10533-014-0038-3
* Wang, S., Flipo, N., Romary, T. (2023). Which filter for data assimilation in water quality models? Focus on oxygen reaeration and heterotrophic bacteria activity. Journal of Hydrology, 620, 129423. 
*    doi:10.1016/j.jhydrol.2023.129423
* 
* COPYRIGHT: (c) 2023 Contributors to the ProSe-PA software. 
* CONTACT: Nicolas FLIPO <nicolas.flipo@minesparis.psl.eu>
*          
* 
* All rights reserved. This software and the accompanying materials
* are made available under the terms of the Eclipse Public License v2.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v20.html
* 
*------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "spmatrix.h"
#include <string.h>
#ifdef OMP
#include <omp.h>
#endif
//#include <libprint.h>
#include <libprint.h>
#include "time_series.h"
#include "libpc.h"
#include "IO.h"
#include "GC.h"
#include "CHR.h"
//#ifdef COUPLED
//#include "spa.h"
//#include "reservoir.h"
//#include "FP.h"
//#endif
#include "HYD.h"
#include "TTC.h"
#include "RIVE.h"
#include "SEB.h"
#include "TUB.h"
#include "MB.h"
#include "LA.h"
#include "PROSE.h"
//#include "global_PROSE.h"
#include "ext_PROSE.h"

void PROSE_print_mats(s_species_ttc *pspecies,FILE *fpout)
{

int i;
 s_gc *pgc;
 pgc=pspecies->pgc;
  LP_printf(fpout,"\nmat4\n");
  for (i = 0; i < pgc->ndl + 1; i++) LP_printf(fpout,"%d ",pgc->mat4[i]);
  LP_printf(fpout,"\nmat5\n");
  for (i = 0; i < pgc->lmat5; i++) LP_printf(fpout,"%d ",pgc->mat5[i]);
  LP_printf(fpout,"\nmat6\n");
  for (i = 0; i < pgc->lmat5; i++) LP_printf(fpout,"%f ",pgc->mat6[i]);
  LP_printf(fpout,"\nb\n");
  for (i = 0; i < pgc->ndl; i++) LP_printf(fpout,"%f ",pgc->b[i]);
  LP_printf(fpout,"\n");
  print_LHS_RHS(pgc);
}

void print_LHS_RHS(s_gc *pgc)
{
	int nrow = 0, ncol; // row and column number
	int i,j;
	int ind_mat6 = 0;
	int ndl;
	double **mat_tab;
	ndl = pgc->ndl; 
	
	mat_tab = (double **)malloc(sizeof(double *)*ndl);
	
	for(i = 0; i < ndl; i++)
	{
		mat_tab[i] = (double *)malloc(sizeof(double)*ndl);
		memset(mat_tab[i],0,sizeof(double)*ndl);
	}
	for(i = 1; i < pgc->ndl + 1; i++) // mat4 dimension
	{

		for(j = pgc->mat4[i-1]; j < pgc->mat4[i]; j++)
		{
			ncol = pgc->mat5[j] - 1;
			mat_tab[nrow][ncol] = pgc->mat6[ind_mat6++];
		}
		nrow++;
	}
	LP_printf(Simul->poutputs,"ind_mat6 = %d\n",ind_mat6);
	for(nrow = 0; nrow < ndl; nrow++)
	{
		for(ncol = 0; ncol < ndl; ncol++)
			LP_printf(Simul->poutputs,"%10.6f\t",mat_tab[nrow][ncol]); // print LHS
		LP_printf(Simul->poutputs,"%10.6f\n",pgc->b[nrow]); //print RHS
	}
		
	for(i = 0; i < ndl; i++) // free memery
	{
		free(*(mat_tab + i));
	}
	mat_tab = NULL;
}

/*void PROSE_fill_sparse(s_gc *pgc,double *RHS_b,	void *mat)
{
	int nrow = 0, ncol;
	int i, j;
	int ind_mat6 = 0;
    spREAL *pelement;
	for(i = 1; i < pgc->ndl + 1; i++) // mat4 dimension
	{

		for(j = pgc->mat4[i-1]; j < pgc->mat4[i]; j++)
		{
			ncol = pgc->mat5[j] - 1;
			pelement = spGetElement(mat,nrow+1,ncol+1);
			*pelement = pgc->mat6[ind_mat6];
			ind_mat6++;
		}
		RHS_b[nrow+1] = pgc->b[nrow];
		nrow++;		
	} 
    //RHS_b[0] = 0;	
}*/

void PROSE_fill_sparse(s_gc *pgc,double *RHS_b,	void *mat)
{
	int i, j;
	//int ind_mat6 = 0;
	//int num_threads;
	//num_threads = Simul->psmp->nthreads;
	//omp_set_num_threads(num_threads);
	//#pragma omp parallel for shared(pgc,mat,RHS_b)     
	for(i = 1; i < pgc->ndl + 1; i++) // mat4 dimension
	{
        spREAL *pelement;
		int ncol;
		for(j = pgc->mat4[i-1]; j < pgc->mat4[i]; j++)
		{
			ncol = pgc->mat5[j] - 1;
			pelement = spGetElement(mat,i,ncol+1);
			*pelement = pgc->mat6[j];
			//ind_mat6++;
		}
		RHS_b[i] = pgc->b[i-1];		
	}
}

void PROSE_allocate_spmatrix(int nspecies, int nele, int nparticules, FILE *fp)
{
	int ns,error,np;
	//int num_threads;
	//double *rhs_b;
	//char *mat;
	Simul->RHS_b = (spREAL ***)calloc(nparticules,sizeof(double **));
	
	Simul->mat_adv = (char ***)calloc(nparticules,sizeof(char **));
	
    //num_threads = Simul->num_threads_par;
	//omp_set_num_threads(num_threads);
	//#pragma omp parallel for shared(nele,nspecies,nparticules) private(np,ns)
	for(np = 0; np < nparticules; np++)
	{	
      Simul->RHS_b[np] = (spREAL **)calloc(nspecies,sizeof(double *));
	  Simul->mat_adv[np] = (char **)calloc(nspecies,sizeof(char *));
	for(ns = 0; ns < nspecies; ns++)
	{
		Simul->RHS_b[np][ns] = (spREAL *)calloc(nele+1,sizeof(double));
		Simul->mat_adv[np][ns] = spCreate(nele,0,&error);
        bzero((char *) Simul->RHS_b[np][ns],(nele+1)*sizeof(double));
		
	}
	}
}

void PROSE_allocate_spmatrix_phy(int nele, int nparticules, FILE *fp)
{
	int ns,error,np;
	int nspecies;
	//int num_threads;
	nspecies = 3*Simul->counter_bio->nsubspecies[PHY];
	Simul->RHS_b_phy = (spREAL ***)calloc(nparticules,sizeof(double **));
	Simul->mat_adv_phy = (char ***)calloc(nparticules,sizeof(char **));

    //num_threads = Simul->psmp->nthreads;
	//omp_set_num_threads(num_threads);
	//#pragma omp parallel for shared(nele,nspecies) private(np,ns)
	for(np = 0; np < nparticules; np++)
	{
		Simul->RHS_b_phy[np] = (spREAL **)calloc(nspecies,sizeof(double *));
		Simul->mat_adv_phy[np] = (char **)calloc(nspecies,sizeof(char *));
	for(ns = 0; ns < nspecies; ns++)
	{
		Simul->RHS_b_phy[np][ns] = (spREAL *)calloc(nele+1,sizeof(double));
		Simul->mat_adv_phy[np][ns] = spCreate(nele,0,&error);
        bzero((char *) Simul->RHS_b_phy[np][ns],(nele+1)*sizeof(double));
		
	}
	}
}


