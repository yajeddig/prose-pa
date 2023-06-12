/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: manage_assimilation.c
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
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
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

void Prose_create_obs_points(s_carac_assim *passim, s_carac_obs_assim *ppobs)
{
  int i;
  int nobs;
  nobs= passim->N_obs;
  ppobs = Prose_browse_ts_obs(ppobs,BEGINNING);

  passim->pobs = (s_carac_obs_assim **)calloc(nobs,sizeof(s_carac_obs_assim *));
  for (i = 0; i < nobs; i++) {
    passim->pobs[i] = ppobs;
    ppobs = ppobs->next;
  }
}

s_carac_obs_assim *Prose_chain_ts_obs(s_carac_obs_assim *pobs1,s_carac_obs_assim *pobs2)
{
  pobs1->next = pobs2;
  pobs2->prev = pobs1;

  return pobs1;
}


s_carac_obs_assim *Prose_browse_ts_obs(s_carac_obs_assim *ppk_tmp,int beg_end)
{

  s_carac_obs_assim *ppk;
  ppk=ppk_tmp;
  if (ppk != NULL) {

    if (beg_end == BEGINNING) {
      while (ppk->prev != NULL)
	ppk = ppk->prev;
    }
    
    if (beg_end == END) {
      while (ppk->next != NULL)
	ppk = ppk->next;
    }
  }

  return ppk;
}


void Prose_find_obs_id_ele(s_carac_obs_assim *ppk,s_chyd *pchyd,FILE *fp)
{
  int r = 0;
  int ne;
  s_reach_hyd *preach;
  s_element_hyd *pele;
  double pkinit;
  int found = NO_TS;
  

  

  while ((found == NO_TS) && (r < pchyd->counter->nreaches))
    {
      
      preach = pchyd->p_reach[r];
      //LP_printf(fp,"river in preach : %s, river in ppk : %s\n",preach->river,ppk->river); // BL to debug
      // LP_printf(fp,"in river : %s pk at uptsream : %f,pk at downstream : %f \n",preach->river,preach->limits[UPSTREAM]->pk,preach->limits[DOWNSTREAM]->pk); // BL to debug
      // LP_printf(fp,"pk at output %f \n",ppk->pk); // BL to debug
      if ((strcmp(preach->river,ppk->river) == 0) && (preach->branch_nb == ppk->branch_nb))
	{
	  if(ppk->pk_type==LENGTH_HYD)
	    {
	      
	      ppk->pk=preach->limits[UPSTREAM]->pk + ppk->pk;
	      //LP_printf(fp,"pk_new at output %f \n",ppk->pk); //BL to debug
	    }
	  if((preach->limits[UPSTREAM]->pk <= ppk->pk) && (preach->limits[DOWNSTREAM]->pk >= ppk->pk)) 
	    {
	      
	      found = YES_TS;
	      ppk->reach_nb = preach->id[ABS_HYD];
	      pkinit = preach->limits[UPSTREAM]->pk;
	      
	      ne = 0;
	      pele = preach->p_ele[ne];
	      
	      while ((ne < preach->nele - 1) && (pkinit + 0.5 * pele->length < ppk->pk))
		{
		  pkinit = pele->face[X_HYD][TWO]->element[TWO]->center->pk;
		  pele = pele->face[X_HYD][TWO]->element[TWO];
		  ne++;
		}
	      ppk->id_ele_obs = pele->id[ABS_HYD];
	      LP_printf(fp,"output pk %s is linked with ele %d in reach\n",ppk->river,ppk->id_ele_obs);
	   
	    }
	}	
      r++;
    }

  if (found == NO_TS)
    LP_error(fp,"Cannot find pk %f in branch %d of %s in given river network\n",ppk->pk / 1000.,ppk->branch_nb,ppk->river);
}


int Prose_calc_difference_obs_simul(s_carac_assim *passim, s_simul ***psimul_bio, int nparticules, double t, FILE *fp)
{
   int np, id_section,nobs;
   int var, num;
   int num_obs = 0;
   int answer = NO_TS;
   double var_obs, sd_sigma, epsilon;
   double varY, sigmaY;
   s_carac_obs_assim *pobs;
   s_ft *obs_t;
   for(nobs = 0; nobs < passim->N_obs; nobs++)
   {
       pobs = passim->pobs[nobs];
       id_section = pobs->id_ele_obs;
       //var_obs = TS_function_value_t(t,pobs->obs, fp); // il faut pas faire interpolation
       pobs->obs = TS_function_t_pointed(t,pobs->obs, fp); // SW 04/02/2020
       obs_t = TS_function_t_pointed(t,pobs->obs, fp); // find the pointer the closest to t
	   
       if(fabs(obs_t->t - t) < 1.) // SW 30/04/2019 replace EPS_TS by 1. second, EPS_TS is too small
       {
           pobs->Obs = obs_t->ft;
           var_obs = obs_t->ft;
	   answer = YES_TS;
	   pobs->answer_obs = YES_TS;
	   var = pobs->var;
	   num = pobs->num;	
           num_obs++; // count how many stations which have a observation at time t;		   
	   
           for(np = 0; np < nparticules; np++)
           {
               // SW 22/11/2021 add obs pertubation for EnKF
               if(passim->method == ENKF_PROSE)
	       {
	           sd_sigma = passim->error_obs_sigma*var_obs;
	           epsilon = Prose_generate_random_norm_param(sd_sigma);
	           var_obs += epsilon; // obs value after the pertubation
                   //sigmaY = psimul_bio[np][id_section]->section->compartments[WATER][0]->pspecies[var][num]->C * 0.1;
                   //varY = psimul_bio[np][id_section]->section->compartments[WATER][0]->pspecies[var][num]->C + Prose_generate_random_norm_param(sigmaY);
                   varY = psimul_bio[np][id_section]->section->compartments[WATER][0]->pspecies[var][num]->C;
                   pobs->difference[np] = var_obs - varY; //(var_obs - varY); // (Y* - HY)
                   pobs->epsilon[np] = epsilon; // this epsilon is used to calculte error covariance matrix of observation
	       }
               else // SW 19/01/2022 add predict pertubation for EnKF
                   pobs->difference[np] = (var_obs - psimul_bio[np][id_section]->section->compartments[WATER][0]->pspecies[var][num]->C);
	       //LP_printf(fp,"nobs = %d np = %d diff = %f\n",nobs,np,pobs->difference[np]);  
	   }
       }
   }
   passim->num_t_obs = num_obs;
   return answer;   
}

void Prose_init_assimilation(s_carac_assim *passim, int nparticules, FILE *fp)
{
         int nobs,np, npd;
	 passim->omega_prev = (double *) malloc(nparticules * sizeof(double));
	 passim->omega = (double *) malloc(nparticules * sizeof(double));
	 passim->omega_normlized = (double *) malloc(nparticules * sizeof(double));
	 passim->cdf = (double *) malloc(nparticules * sizeof(double));
	 
	 passim->eliminated = (int *) malloc(nparticules * sizeof(int));
	 passim->ndupli = (int *) malloc(nparticules * sizeof(int));
	 
	 //passim->param = (double **) malloc(nparticules * sizeof(double*));
	 for(np = 0; np < nparticules; np++)
	 {
		 passim->omega_prev[np] = 1./nparticules;
		 passim->omega[np] = 1./nparticules;
		 passim->omega_normlized[np] = 1./nparticules;
		 passim->eliminated[np] = NO_TS;
		 passim->ndupli[np] = 0;
		 //passim->param[np] = (double *) malloc(NPARAMDA * sizeof(double));
	 }
	 
         passim->Neff = nparticules; // SW 03/01/2023

	 for(np = 0; np < NPARAMDA; np++)
         {
                 // SW 24/01/2022 add check if parameter np is assimilated
                 if(passim->param_yesOrno[np] == YES_TS)
		   {
		     passim->param[np] = (double *) malloc(nparticules * sizeof(double));

		      // SW 20/05/2022 allocation table for perturbation density and initialized to 1.0
                     passim->param_perturbation_density[np] = (double *) malloc(nparticules * sizeof(double));
                     for(npd = 0; npd < nparticules; npd++)
                         passim->param_perturbation_density[np][npd] = 1.;

                 }
                 // SW 14/01/2022
                 //if(passim->method == ENKF_PROSE)
                 //    passim->paramGaussian[np] = (double *) malloc(nparticules * sizeof(double)); 

         }
	 // SW 20/05/2022 allocation of perturbation densities production
         passim->param_perturbation_density_product = (double *) malloc(nparticules * sizeof(double));

	 // MH 24/01/2022 : initialization of param range to CODE_TS in order to facilitate the if clause condition in manage_link.c
	 /*for(np = 0; np < NPARAMDA; np++) {
	   passim->param_range[np] = (double *) malloc(PARAMR * sizeof(double));
	   passim->param_range[np][PARAM_DOWN] = NULL;
	   passim->param_range[np][PARAM_UP] = NULL;
	   }*/

	 
	 for(nobs = 0; nobs < passim->N_obs; nobs++)
	 {
		 //Simul->passim->pobs[nobs]->y = (s_ft **) malloc(nparticules * sizeof(s_ft *));
		 // find _id_ele_obs
		 Prose_find_obs_id_ele(Simul->passim->pobs[nobs],Simul->pchyd,Simul->poutputs);
         for(np = 0; np < nparticules; np++)
		 {
             //Simul->passim->pobs[nobs]->y[np] = new_function();
			 Simul->passim->pobs[nobs]->difference = (double *) malloc(nparticules * sizeof(double));
                         Simul->passim->pobs[nobs]->epsilon = (double *) malloc(nparticules * sizeof(double)); // SW 19/01/2022 to store epsilon for EnKF
		 }			 
	}
	passim->state = BLOOM;
}

void PROSE_init_output_simulation(int nparticules, FILE *fp)
{
  int e, layer;
  s_simul *psimulbio;
  psimulbio = Simul->psimul_bio[0][0];	
  for (layer = 0; layer < NLAYERS; layer++) {
    if (psimulbio->section->compartments[layer] != NULL) { 
	   Simul->concentrations[layer]  = (FILE ****)calloc(nparticules,sizeof(FILE ***));
	   for (e = 0; e < NSPECIES; e++)
	       Simul->mass_balances[layer][e] = (FILE ****)calloc(nparticules,sizeof(FILE ***));
	   for (e = 0; e < NANNEX_VAR; e++)
		   Simul->mass_balances[layer][NSPECIES+e] = (FILE ****)calloc(nparticules,sizeof(FILE ***));
	   /*
	   for (e = 0; e < NDISS; e++)
		   Simul->ads_mass_balances[layer][e] = (FILE ****)calloc(nparticules,sizeof(FILE ***));
	   */
	}
  }
}

// This function used for resetting answer_obs = NO_TS which is default. When there is a
// observation, it is assigned to be YES_TS
void PROSE_reset_no_answer_obs(s_carac_assim *passim, FILE *fp)
{
    int nobs;
	int num_threads,taille;
    s_carac_obs_assim *pobs;
    
    //#ifdef CDA     
    //#ifdef OMP
    //num_threads = Simul->num_threads_par;
    //taille = PC_set_chunk_size_silent(fp,passim->N_obs,num_threads);
    //omp_set_num_threads(num_threads);
    //#pragma omp parallel for schedule(dynamic,taille) private(np)
    //#endif
	//#endif	
    for(nobs = 0; nobs < passim->N_obs; nobs++)
    {
	    pobs = passim->pobs[nobs];
	    pobs->answer_obs = NO_TS;
    }
}
	
