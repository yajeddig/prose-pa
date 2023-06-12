/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: solve_prose.c
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
#include <strings.h>
//#include <libprint.h>
#include <libprint.h>
#ifdef OMP
#include <omp.h>
#endif
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
//boucle temporelle

void PROSE_ttc_one_species(s_carac_ttc *pcarac_ttc,s_param_calc_ttc *pparam_calc_ttc,s_species_ttc *pspecies,double dt,double t,FILE *fpout)
{
  //nt nombre d'iteration temporelle le temps final de simulation = dt*nt
  //dt discretisation temporelle (s)
  s_gc *pgc;
  int nele;
  int i,regime,type;
  int ts_count = 0; // SW 15/06/2021 add ts_count for pgc configuration

  nele=pcarac_ttc->count[NELE_TTC];
  regime=pcarac_ttc->regime;
  type=pspecies->type;
  pgc=pspecies->pgc;

  //initialisation //TTC_init_tab_gc
  printf("%lf %lf\n",t-dt,Simul->chronos->t[BEGINNING]);
  //if(t - dt <= Simul->chronos->t[BEGINNING]) //SW 02/02/2021
  if(t - dt - Simul->chronos->t[BEGINNING] < EPS_TS) //SW 02/02/2021
    TTC_init_mat5_mat4(pspecies,nele,fpout);
  else
    ts_count = 2; // SW 15/06/2021 add ts_count for pgc configuration
  
  //initialisation //TTC_init_tab_gc
  TTC_init_RHS(pgc);
  TTC_init_sol_vec(pgc);
  TTC_free_mat6(pgc);
  //
  //TTC_init_mat5_mat4(pspecies,nele,fpout);
  pgc->mat6=GC_create_mat_double(pspecies->pgc->lmat5);
  //On remplit mat6 qui sera utile pour le solveur libgc
  TTC_fill_mat6_and_b_rive(pcarac_ttc,pparam_calc_ttc,pspecies,dt,fpout);
  
  //PROSE_print_mats(pspecies,fpout);
  //print_LHS_RHS(pgc); //SW
  //GC_solve_gc(pgc);

  GC_configure_gc(pgc, &ts_count, fpout); // SW 15/06/2021 set a specific pgc configuration
  GC_solve_gc(pgc, TTC_GC, fpout); // SW 10/06/2021 update libgc (preconditionnement)

  //On met la solution de gc dans le tableau var a linstant t qui servira d initialisation pour literation temporelle suivante
  for(i=0;i<nele;i++)
    {
      pspecies->plink->pvar_ttc->var[i]=pgc->x[i];
      //LP_printf(fpno3,"%f\t",pspecies->plink->pvar_ttc->var[i]); // SW
      //LP_printf(fpout,"t = %f,id = %d, var = %f\n",t,i+1,pspecies->plink->pvar_ttc->var[i]); // SW
    }
   /* SW 02/02/2021 free mat4 mat5*/
   /* free(pgc->mat4);
   pgc->mat4 = NULL;
   free(pgc->mat5);
   pgc->mat5 = NULL; */

  //LP_printf(fpno3,"\n");
}

void PROSE_ttc_one_species_sparse(s_carac_ttc *pcarac_ttc,s_param_calc_ttc *pparam_calc_ttc,s_species_ttc *pspecies,double dt,double t,double *RHS_b, void *mat,int num_threads, FILE *fpout)
{
  //nt nombre d'iteration temporelle le temps final de simulation = dt*nt
  //dt discretisation temporelle (s)
  s_gc *pgc;
  int nele;
  int i,regime,type;
  double *sol_x;
  int error;
  
  nele=pcarac_ttc->count[NELE_TTC];
  regime=pcarac_ttc->regime;
  type=pspecies->type;
  pgc=pspecies->pgc;
  sol_x = (spREAL *)calloc(nele+1,sizeof(double));
  bzero((char *) sol_x,(nele+1)*sizeof(double));
  //initialisation //TTC_init_tab_gc
  //if(t - dt <= Simul->chronos->t[BEGINNING])  
  if(t - dt - Simul->chronos->t[BEGINNING] < EPS_TS) //SW 02/02/2021
     TTC_init_mat5_mat4(pspecies,nele,fpout);

  
  //initialisation //TTC_init_tab_gc
  TTC_init_RHS(pgc);
  //TTC_init_sol_vec(pgc);
  TTC_free_mat6(pgc);
  //
  //TTC_init_mat5_mat4(pspecies,nele,fpout);
  pgc->mat6=GC_create_mat_double(pspecies->pgc->lmat5);
  //On remplit mat6 qui sera utile pour le solveur libgc
  TTC_fill_mat6_and_b_rive(pcarac_ttc,pparam_calc_ttc,pspecies,dt,fpout);
  
  //PROSE_print_mats(pspecies,fpout);
  //print_LHS_RHS(pgc); //SW
  //GC_solve_gc(pgc);
  spClear(mat);
  PROSE_fill_sparse(pgc,RHS_b,mat);
  TTC_free_mat6(pgc);
  //spPrint(mat,0,1,1);
  
  /* SW 02/02/2021, This routine is the companion routine to spOrderAndFactor(). Unlike spOrderAndFactor(), spFactor() cannot change the ordering. 
  It is also faster than spOrderAndFactor(). The standard way of using these two routines is to first use spOrderAndFactor() for the initial factorization. 
  For subsequent factorizations, spFactor() is used if there is some assurance that little growth will occur (say for example, that the matrix is diagonally dominant). 
  If spFactor() is called for the initial factorization of the matrix, then spOrderAndFactor() is automatically called with the default threshold. This routine uses "row at a time" LU factorization. 
  Pivots are associated with the lower triangular matrix and the diagonals of the upper triangular matrix are ones. */

  /* SW 02/02/2021, bring some problem for heat transport, very strange, so comment spFactor for heat transport*/
  //if(t - dt == Simul->chronos->t[BEGINNING])
  switch(type) {
      case SOLUTE_TTC :
          if(fabs(t - dt - Simul->chronos->t[BEGINNING]) < EPS_TS) 
             error = spOrderAndFactor(mat,RHS_b,0.01,0,0);
          else
             error = spFactor(mat); 
          break;
      case HEAT_TTC: 
          error = spOrderAndFactor(mat,RHS_b,0.01,0,0);
          break;
      default: LP_printf(fpout,"Unknown species type, SOLUTE_TTC or HEAT_TTC. Species name is %s\n", pspecies->name);
  }
   spSolve(mat,RHS_b,sol_x);
  //On met la solution de gc dans le tableau var a linstant t qui servira d initialisation pour literation temporelle suivante
  //#ifndef CDA
  //#ifdef OMP
  //omp_set_num_threads(num_threads);
  //#pragma omp parallel for schedule(dynamic,1) shared(pspecies,sol_x,nele) 
  //#endif
  //#endif
  for(i=0;i<nele;i++)
  {
  //pspecies->plink->pvar_ttc->var[i]=pgc->x[i];
  pspecies->plink->pvar_ttc->var_old[i] = pspecies->plink->pvar_ttc->var[i];
  pspecies->plink->pvar_ttc->var[i]=sol_x[i+1];
  //LP_printf(fpno3,"%f\t",pspecies->plink->pvar_ttc->var[i]); // SW
  //LP_printf(fpout,"t = %f,id = %d, var = %f\n",t,i+1,pspecies->plink->pvar_ttc->var[i]); // SW
  }
  free(sol_x);
  sol_x = NULL; 

   /* SW 02/02/2021 free mat4 mat5*/
   /*free(pgc->mat4);
   pgc->mat4 = NULL;
   free(pgc->mat5);
   pgc->mat5 = NULL;*/
 
  //LP_printf(fpno3,"\n");
}

//void PROSE_steady_solve(s_carac_ttc *pcarac_ttc,s_param_calc_ttc *pparam_calc_ttc,s_species_ttc *pspecies, FILE *fpout)
void PROSE_steady_solve(s_carac_ttc *pcarac_ttc,s_param_calc_ttc *pparam_calc_ttc,s_species_ttc *pspecies, int ns,FILE *fpout)
{
  s_gc *pgc;
  int nele;
  int i,regime,type;
  double dt=0;
  //FILE* fpvar= NULL;
  //FILE* fpmb= NULL;
  nele=pcarac_ttc->count[NELE_TTC];
  regime=pcarac_ttc->regime;
  type=pspecies->type;
  pgc=pspecies->pgc;

  //TTC_header_output_steady(fpvar,type);
  //TTC_header_output_MB_steady(fpmb,type);
  //TTC_init_tab_gc
  TTC_init_mat5_mat4(pspecies,nele,fpout);
  pgc->mat6=GC_create_mat_double(pspecies->pgc->lmat5);
  TTC_init_RHS(pgc);
  pgc->x=GC_create_mat_double(pspecies->pgc->ndl);
  //On remplit mat6 qui sera utile pour le solveur libgc
    TTC_fill_mat6_and_b_rive(pcarac_ttc,pparam_calc_ttc,pspecies,dt,fpout);
	  //print_LHS_RHS(pspecies->pgc); // debug
    //GC_solve_gc(pgc);
	PROSE_solve_gc(pgc, ns);
	
  //On met la solution de gc dans le tableau var a linstant t qui servira dinitialisation pour literation temporelle suivante

  for(i=0;i<nele;i++)
	{
	  pspecies->plink->pvar_ttc->var[i]=pspecies->pgc->x[i];
	  //pspecies->plink->pvar_ttc->var[i]=sol_x[i];
	  //printf("%s valeur var a ele %i \t %lf\n",pspecies->name,i,pspecies->plink->pvar_ttc->var[i]);
	  //fprintf(fpout,"%s valeur var a ele %i \t %lf\n",pspecies->name,i,pspecies->plink->pvar_ttc->var[i]);
	}
	

  //for(i=0;i<nele;i++)
	//{
   //TTC_cal_condflux(pcarac_ttc,pparam_calc_ttc,pspecies,i);
   //TTC_cal_advflux(pcarac_ttc,pparam_calc_ttc,pspecies,i);
   //TTC_cal_totflux(pcarac_ttc,pparam_calc_ttc,pspecies,i);
	//}
   //TTC_print_output_var_steady(pcarac_ttc,pspecies->plink->pvar_ttc,fpvar);
   //TTC_print_output_MB_steady(pcarac_ttc,pparam_calc_ttc,pspecies,fpmb);
}

void PROSE_steady_solve_sparse(s_carac_ttc *pcarac_ttc,s_param_calc_ttc *pparam_calc_ttc,s_species_ttc *pspecies, double *RHS_b, void *mat,FILE *fpout)
{
  s_gc *pgc;
  int nele;
  int i,regime,type;
  double dt=0;
  //FILE* fpvar= NULL;
  //FILE* fpmb= NULL;
  double *sol_x;
  int error;
  nele=pcarac_ttc->count[NELE_TTC];
  regime=pcarac_ttc->regime;
  type=pspecies->type;
  pgc=pspecies->pgc;
  sol_x = (spREAL *)calloc(nele+1,sizeof(double));
  bzero((char *) sol_x,(nele+1)*sizeof(double));
  //TTC_header_output_steady(fpvar,type);
  //TTC_header_output_MB_steady(fpmb,type);
  //TTC_init_tab_gc
  TTC_init_mat5_mat4(pspecies,nele,fpout);
  pgc->mat6=GC_create_mat_double(pspecies->pgc->lmat5);
  TTC_init_RHS(pgc);
  //pgc->x=GC_create_mat_double(pspecies->pgc->ndl);
  //On remplit mat6 qui sera utile pour le solveur libgc
    TTC_fill_mat6_and_b_rive(pcarac_ttc,pparam_calc_ttc,pspecies,dt,fpout);
	  //print_LHS_RHS(pspecies->pgc); // debug
    //GC_solve_gc(pgc);
	//PROSE_solve_gc(pgc, ns);
  //spClear(mat);
  PROSE_fill_sparse(pgc,RHS_b,mat);
  TTC_free_mat6(pgc);
  //spPrint(mat,0,1,1);
  error = spOrderAndFactor(mat,RHS_b,0.01,0,0);
  spSolve(mat,RHS_b,sol_x);
  
  //On met la solution de gc dans le tableau var a linstant t qui servira dinitialisation pour literation temporelle suivante
	  for(i=0;i<nele;i++)
	{
	  //pspecies->plink->pvar_ttc->var[i]=pspecies->pgc->x[i];
	  pspecies->plink->pvar_ttc->var[i]=sol_x[i+1];
	  //printf("%s valeur var a ele %i \t %lf\n",pspecies->name,i,pspecies->plink->pvar_ttc->var[i]);
	  //fprintf(fpout,"%s valeur var a ele %i \t %lf\n",pspecies->name,i,pspecies->plink->pvar_ttc->var[i]);
	}	
	free(sol_x);
	sol_x = NULL;
  //for(i=0;i<nele;i++)
	//{
   //TTC_cal_condflux(pcarac_ttc,pparam_calc_ttc,pspecies,i);
   //TTC_cal_advflux(pcarac_ttc,pparam_calc_ttc,pspecies,i);
   //TTC_cal_totflux(pcarac_ttc,pparam_calc_ttc,pspecies,i);
	//}
   //TTC_print_output_var_steady(pcarac_ttc,pspecies->plink->pvar_ttc,fpvar);
   //TTC_print_output_MB_steady(pcarac_ttc,pparam_calc_ttc,pspecies,fpmb);
}


void PROSE_ttc_all_species(s_carac_ttc *pcarac_ttc, s_chyd *pchyd, double dt,double t,int np,double tempe, double Osat, FILE *fpout)
{
	int ns;
	int nspecies, nele;
	int phy, nsub;
	int num_threads,taille;
	//s_species_ttc *pspecies;
	//s_param_calc_ttc *pparam_calc_ttc;	
	
	
    nspecies = pcarac_ttc->count[NSPECIES_TTC];
    nele = pcarac_ttc->count[NELE_TTC];	
	num_threads = Simul->psmp->nthreads;
	//int tid1 = omp_get_thread_num();
	//printf("tid1 = %d np = %d\n",tid1,np); 
	#ifndef CDA
	#ifdef OMP
	Simul->psmp->chunk = PC_set_chunk_size_silent(fpout,nspecies-1,num_threads);
	taille = Simul->psmp->chunk;
	omp_set_num_threads(num_threads);
	#pragma omp parallel for schedule(dynamic,taille) shared(pcarac_ttc,pchyd,dt,t,fpout,nele,np,tempe,Osat) private(ns)
	#endif
	#endif
	for(ns = 0; ns < nspecies; ns++)
	{
		//int tid = omp_get_thread_num();
		s_species_ttc *pspecies1;
		s_param_calc_ttc *pparam_calc_ttc1;			
		pspecies1 = pcarac_ttc->p_species[ns];
		//if(ns == 0 || ns == 150) //printf("ns = %d\n",ns);
		//printf("tid = %d name = %s np = %d\n",tid,pspecies1->name,np);
		pparam_calc_ttc1 = pspecies1->plink->pparam_calc_ttc;		
		PROSE_fill_var_one_species(pspecies1, nele, np, fpout); // SW
		PROSE_fill_u_one_species(pspecies1, pchyd,tempe, Osat,fpout); // SW 24/10/2018 add temperature for reoxygenation
		PROSE_find_apport_t_one_species(pspecies1, pchyd, t,np, fpout);
	    if(Simul->calc_mode[MB_BIO] == YES_TS && ns == 0) // SW 17/10/2018 include mass balance for transport
	       PROSE_calc_mb_minit_all_sections(Simul->psimul_bio[np], t, dt, nele, np,Simul->poutputs);

		if(pcarac_ttc->regime == STEADY)
		{	
             if(Simul->solver == SP_PROSE)	
				 PROSE_steady_solve_sparse(pcarac_ttc,pparam_calc_ttc1,pspecies1, Simul->RHS_b[np][ns], Simul->mat_adv[np][ns],fpout);
			 else
			     PROSE_steady_solve(pcarac_ttc, pparam_calc_ttc1, pspecies1,ns,fpout);
			
			//LP_printf(fpout,"ok ns = %d iappli = %d\n",ns,pspecies1->iappli_gc); // debug
			//print_LHS_RHS(pspecies1->pgc); // debug
		}
		else if(pcarac_ttc->regime == TRANSIENT)
		{	
	        //PROSE_fill_surf_ttc_one_species(pspecies1, pchyd, fpout);	
			if(Simul->solver == SP_PROSE)
          	    PROSE_ttc_one_species_sparse(pcarac_ttc,pparam_calc_ttc1,pspecies1,dt, t,Simul->RHS_b[np][ns], Simul->mat_adv[np][ns],num_threads,fpout);
			else
			    PROSE_ttc_one_species(pcarac_ttc, pparam_calc_ttc1, pspecies1, dt, t, fpout);		
		}
		//SW 17/10/2018 add mass balance, only mb in total domaine validate
		//if(Simul->calc_mode[MB_BIO] == YES_TS)
		//{
		   //Prose_cal_tot_advflux(pcarac_ttc,pparam_calc_ttc1,pspecies1);
           //PROSE_calc_mb_one_species(pspecies1, Simul->total_mb[np][0],nele, np, dt, fpout);
		//}		   
	}
	

	for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
	{
	   #ifndef CDA
	   #ifdef OMP
	   omp_set_num_threads(3);
	   #pragma omp parallel for schedule(dynamic,1) shared(pcarac_ttc,pchyd,dt,t,fpout,nele,np,tempe) 
	   #endif
	   #endif
		for(phy = 0; phy < 3; phy++)
		{
		   //int tid = omp_get_num_procs();
			int num_s;
		    s_species_ttc *pspecies1;
		    s_param_calc_ttc *pparam_calc_ttc1;
            num_s = nsub*3 + phy;			
			//printf("tid = %d\n",tid);
	        pspecies1 = Simul->p_phy_species[np][nsub][phy];
			pparam_calc_ttc1 = pspecies1->plink->pparam_calc_ttc;	
			PROSE_fill_var_one_annex_species(pspecies1, nele, nsub , phy , np,fpout);
			PROSE_fill_u_one_species(pspecies1, pchyd, tempe,Osat,fpout);
			PROSE_find_apport_t_one_annex_species(pspecies1, pchyd, t, nsub, phy, fpout);
			if(pcarac_ttc->regime == STEADY)
		   {	
               if(Simul->solver == SP_PROSE)	
                  PROSE_steady_solve_sparse(pcarac_ttc,pparam_calc_ttc1,pspecies1, Simul->RHS_b_phy[np][num_s], Simul->mat_adv_phy[np][num_s],fpout);
			   else   
			      PROSE_steady_solve(pcarac_ttc, pparam_calc_ttc1, pspecies1,num_s,fpout);
		      //print_LHS_RHS(pspecies->pgc); // debug
			}
		   else if(pcarac_ttc->regime == TRANSIENT)
		   {	
	        //PROSE_fill_surf_ttc_one_species(pspecies1, pchyd, fpout);
			if(Simul->solver == SP_PROSE)
	           PROSE_ttc_one_species_sparse(pcarac_ttc,pparam_calc_ttc1,pspecies1,dt, t,Simul->RHS_b_phy[np][num_s], Simul->mat_adv_phy[np][num_s],num_threads,fpout);
			else
			   PROSE_ttc_one_species(pcarac_ttc, pparam_calc_ttc1, pspecies1, dt, t, fpout);		
		   }  	   
		//SW 17/10/2018 add mass balance, only mb in total domaine validate
		//if(Simul->calc_mode[MB_BIO] == YES_TS)
           //PROSE_calc_mb_one_species(pspecies1, Simul->total_mb[np][0],nele, np, dt, fpout);		
	}	
	}	
	
}

//void PROSE_one_section(s_simul *psimulbio, double t, double dt, double tempe,int calc_bilan, FILE *fp )
void PROSE_one_section(s_simul *psimulbio, double t, double dt, int calc_bilan, FILE *fp ) // SW 04/12/2019 
{
	   	
   double i0, tempe;;
   //clock_t start,end; 
   //PROSE_calc_volume_water(t, psimulbio);
                //if( t > 9401400 && isnan(psimulbio->section->compartments[1][0]->pspecies[BACT][0]->C))
                  //LP_printf(fp," before bio t = %f id = %d Cup_mod1 = %f, cdown_mod1 = %f cbact1 = %f vol_down = %f \n",t,psimulbio->section->id_section,psimulbio->section->compartments[WATER][0]->pspecies[7][0]->C,psimulbio->section->compartments[VASE][0]->pspecies[7][0]->C,psimulbio->section->compartments[VASE][0]->pspecies[BACT][0]->C,psimulbio->section->compartments[VASE][0]->state->volume);

	/***SW 04/12/2019 add parameter calculation using temperature calulated by libseb***/
   if((Simul->calc_mode[H_T] == YES_TS) && (Simul->calc_mode[SEB] == YES_TS))
   {	   
	   tempe = psimulbio->section->meteo->tempe_value; // must have a value
	   i0 = calc_radiation(t,psimulbio);
	   calc_O2_sat(tempe,psimulbio);
	   calc_N2O_sat(tempe,psimulbio);
	   calc_em_o2_t(tempe,psimulbio);
	   calc_param_bio_t(tempe,psimulbio);
	   if(Simul->calc_mode[DA] == NO_TS)
		   calc_param_bio_t_variable(t,psimulbio); // SW 03/02/2020
	   calc_param_bio_t(tempe,psimulbio);
   }
   else
	  i0 = calc_param_t(t,psimulbio, Simul->calc_mode[DA]); 
   //calc_hyd(t,psimulbio);
   PROSE_calc_hyd(t,psimulbio);
   //start = clock();   
   //if(isnan(psimulbio->section->compartments[0][0]->pspecies[7][0]->C) || isnan(psimulbio->section->compartments[1][0]->pspecies[7][0]->C) || isnan(psimulbio->section->compartments[1][0]->pspecies[BACT][0]->C))
	// LP_error(fp,"yes is nan in solve prose line 381 id_section = %d\n",psimulbio->section->id_section);   
   processus_bio(t,dt,i0,psimulbio,calc_bilan);       
   //end = clock(); 
   //if(isnan(psimulbio->section->compartments[0][0]->pspecies[7][0]->C) || isnan(psimulbio->section->compartments[1][0]->pspecies[7][0]->C ) || isnan(psimulbio->section->compartments[1][0]->pspecies[BACT][0]->C))
	 //LP_error(fp,"yes is nan in solve prose line 385 id_section = %d cup = %f cdown = %f cbct_down = %f vol = %f mass = %f\n",psimulbio->section->id_section,psimulbio->section->compartments[0][0]->pspecies[7][0]->C, psimulbio->section->compartments[1][0]->pspecies[7][0]->C,psimulbio->section->compartments[1][0]->pspecies[BACT][0]->C,psimulbio->section->compartments[1][0]->state->volume,psimulbio->section->compartments[1][0]->state->mass);   
   corrige_mass_sed(psimulbio);
   //if(isnan(psimulbio->section->compartments[0][0]->pspecies[7][0]->C) || isnan(psimulbio->section->compartments[1][0]->pspecies[7][0]->C) || isnan(psimulbio->section->compartments[1][0]->pspecies[BACT][0]->C))
	 //LP_error(fp,"yes is nan in solve prose line 388 id_section = %d\n",psimulbio->section->id_section);   
   //LP_printf(fp,"id == %d t = %f vol = %f\n",psimulbio->section->id_section,(double)(end-start)/CLOCKS_PER_SEC,psimulbio->section->compartments[1][0]->state->volume);
   flows_at_interfaces(t, dt, psimulbio);
   //if(isnan(psimulbio->section->compartments[0][0]->pspecies[7][0]->C) || isnan(psimulbio->section->compartments[1][0]->pspecies[7][0]->C) || isnan(psimulbio->section->compartments[1][0]->pspecies[BACT][0]->C))
	// LP_error(fp,"yes is nan in solve prose line 392 id_section = %d\n",psimulbio->section->id_section);    
}

void PROSE_all_sections(s_simul **psimul_bio, double t, double dt, int nsections, int np, FILE *fp)
{
	int ns;
	int num_threads,taille,calc_bilan = 0;
	s_simul *psimulbio;
    //clock_t start,end;
    //double tempe = 0;	
	//PROSE_calc_volume_water_all_sections(psimul_bio, nsections, t, fp);
	/*calulation of initial mass for all sections and the sum*/
   
	//nmission = (int) nsections/num_threads;
	//if(num_threads == 1)
		//num_threads = 10;
	//CHR_begin_timer();
	//start = clock();
	if(Simul->calc_mode[MB_BIO] == YES_TS)
		calc_bilan = 1;
    #ifndef CDA
	#ifdef OMP
	num_threads = Simul->psmp->nthreads;
	Simul->psmp->chunk = PC_set_chunk_size_silent(fp,nsections-1,num_threads);
	taille = Simul->psmp->chunk;	
	omp_set_num_threads(num_threads);
	#pragma omp parallel for schedule(dynamic,taille) shared(psimul_bio,dt,nsections,t,fp) private(psimulbio)
	#endif
    #endif	
	for(ns = 0; ns < nsections; ns++)
	{

		psimulbio = psimul_bio[ns];
        //if(ns == 3366 && t >= 1200)	
          //LP_printf(fp,"id == 3366\n");			
		//PROSE_one_section(psimulbio, t, dt, tempe, calc_bilan,fp);
		PROSE_one_section(psimulbio, t, dt, calc_bilan,fp); //SW 04/12/2019
		/*calculation of mass balance variation, and sum processus by processus for the period dt*nsteps*/
	    //PROSE_calc_mb_one_section(t, psimulbio, fp);
	}
	
	//Simul->clock->time_spent[SOLVE_RIVE] += CHR_end_timer();
	//#ifdef OMP
	//omp_set_num_threads(num_threads);
	//#pragma omp parallel for  shared(psimul_bio,t,fp,Simul) private(psimulbio)
	//#endif	
	//end = clock(); 
   //LP_printf(fp,"tot = %f\n",(double)(end-start)/CLOCKS_PER_SEC);
	
}

void PROSE_solve_gc(s_gc *pgc, int ns)
{
  int *mat4,*mat5,*sw_int,*ndl,*lmat5,appl_nb;
  double *epsgc,*sol,*b,*mat6;
  mat4=pgc->mat4;
  mat5=pgc->mat5;
  ndl=&pgc->ndl;
  lmat5=&pgc->lmat5;
  sw_int=pgc->sw_int;
  sol=pgc->x;
  epsgc=&pgc->epsgc;
  mat6=pgc->mat6;
  b=pgc->b;
  appl_nb=pgc->appl_nb; 
  /* SW 10/06/2021 initialization proposed by Niolcas G. should be added */
  if(ns == 0)
    gc_init_sys_(mat4,mat5,ndl);

 gc_solve_(mat6,lmat5,b,ndl,sol,sol,sw_int,epsgc);
}
