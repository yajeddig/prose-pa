/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: manage_bio_prose.c
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
#include <strings.h>
#ifdef OMP
#include <omp.h>
#endif
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

void PROSE_init_param_reaction_state_comp(s_simul **psimul_bio, int nsections, FILE *fp)
{
	int ns;
	s_simul *psimulbio;
	for(ns = 0; ns < nsections; ns++)
	{ 
        psimulbio = psimul_bio[ns];

	determine_parameters(psimulbio);
	init_reaction_species(psimulbio);		
	init_state_compartments(psimulbio);	
        if(Simul->calc_mode[DA] == NO_TS)	// SW 03/02/2020
		calc_param_bio_t_variable(Simul->chronos->t[BEGINNING], psimulbio);
        init_annex_var(psimulbio->section, psimulbio);			
	}
	
}

void PROSE_update_phy_ctot(s_simul **psimul_bio, int nsections, int np,FILE *fp)
{
	int ns, nsub,phy;
	s_simul *psimulbio;
	s_compartment *pcomp;
	double ctot;
	for(ns = 0; ns < nsections; ns++)
	{		
        psimulbio = psimul_bio[ns];
		for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
		{
			pcomp = psimulbio->section->compartments[WATER][0];
			ctot = 0.0;
			for(phy = 0; phy < 3; phy++)
			{
				ctot += Simul->p_phy_species[np][nsub][phy]->plink->pvar_ttc->var[ns];
			}
			//LP_printf(fp,"Cphy_tot = %f, Cphy = %f\n",ctot,Simul->pcarac_ttc->p_species[0]->plink->pvar_ttc->var[ns]);
			// SW 02/10/2018 add newC
			pcomp->pspecies[PHY][nsub]->C = ctot;
			pcomp->pspecies[PHY][nsub]->newC = ctot;
			for(phy = 0; phy < 3; phy++)
			{
				psimulbio->settings->phy2[phy] = Simul->p_phy_species[np][nsub][phy]->plink->pvar_ttc->var[ns]/ctot;
				pcomp->pannex_var[phy][nsub]->C = Simul->p_phy_species[np][nsub][phy]->plink->pvar_ttc->var[ns];
				// SW 02/10/2018 add newC
				pcomp->pannex_var[phy][nsub]->newC = Simul->p_phy_species[np][nsub][phy]->plink->pvar_ttc->var[ns];
				//LP_printf(fp,"phy2 %d == %f\n",phy,psimulbio->settings->phy2[phy]);
			}
		}
	}
}

void PROSE_init_O2_conc_sat(s_simul **psimul_bio, double t, int nsections, FILE *fp)
{
	int ns;
	double tem;
	s_simul *psimulbio;
	for(ns = 0; ns < nsections; ns++)
	{
		psimulbio = psimul_bio[ns];		
        if((Simul->calc_mode[H_T] == NO_TS) && (Simul->calc_mode[SEB] == NO_TS)) /* SW 04/12/2019*/
			tem = calc_temperature(t,psimulbio);
		else if(Simul->calc_mode[H_T] == YES_TS)
			tem = psimulbio->section->meteo->tempe_value; // SW 04/12/2019 il faudrait initialiser par libseb dans main_Prose.c			

		calc_O2_sat(tem,psimulbio);
		init_O2_conc(psimulbio);
		calc_em_o2_t(tem,psimulbio);
		calc_N2O_sat(tem,psimulbio);		
	}
}

void PROSE_link_hyd_rive_hydro(s_simul **psimul_bio, double t, s_chyd *pchyd, int nthreads, FILE *fp)
{
	int ne, r, id_abs_ele;
    s_reach_hyd *preach;
    s_element_hyd *pele;
	s_simul *psimulbio;
	//int nmission;
	//double tempe = 0;
	//double hydrad;

  for(r = 0; r < pchyd->counter->nreaches; r++)
  {
	  preach = pchyd->p_reach[r];
	#ifndef CDA
	#ifdef OMP
	omp_set_num_threads(nthreads);
	#pragma omp parallel for shared(psimul_bio,preach,t,fp) private(psimulbio,pele,id_abs_ele)
	#endif
	#endif
	  for(ne = 0; ne < preach->nele; ne++)
	  {
		 pele = preach->p_ele[ne];
		 id_abs_ele = pele->id[ABS_HYD];          		  
	     psimulbio =  psimul_bio[id_abs_ele];
		 //psimulbio->section->hydro->h = TS_create_function(t,5);//TS_create_function(t,pele->center->hydro->H[T_HYD]);
		 //if(psimulbio->section->hydro->h != NULL){
			// psimulbio->section->hydro->h = TS_free_ts(psimulbio->section->hydro->h,fp);
			 //psimulbio->section->hydro->h = NULL;
		 //}
		 psimulbio->section->hydro->h->t = t;
		 //psimulbio->section->hydro->h->ft = pele->center->hydro->Surf/pele->center->hydro->Width;//pele->center->hydro->H[T_HYD];
                 /* SW 09/06/2021 use HL as wet surface, rectangular section */
                 psimulbio->section->hydro->h->ft = pele->center->hydro->H[T_HYD];

            //psimulbio->section->hydro->q = TS_create_function(t,8);//TS_create_function(t,pele->center->hydro->Q[T_HYD]);
		 //if(psimulbio->section->hydro->q != NULL){
			// psimulbio->section->hydro->q = TS_free_ts(psimulbio->section->hydro->q,fp);
			// psimulbio->section->hydro->q = NULL;
		 //}
		 psimulbio->section->hydro->q->t = t;
		 psimulbio->section->hydro->q->ft = pele->center->hydro->Q[T_HYD];
		 //psimulbio->section->hydro->vt = TS_create_function(t,0.16); //TS_create_function(t,pele->center->hydro->Vel);
		 //if(psimulbio->section->hydro->vt != NULL){
			//psimulbio->section->hydro->vt = TS_free_ts(psimulbio->section->hydro->vt,fp);
            //psimulbio->section->hydro->vt = NULL;		
		 //}
		 
		 psimulbio->section->hydro->surf->t = t;
		 //psimulbio->section->hydro->surf->ft = pele->center->hydro->Surf;

                 /* SW 09/06/2021 use HL as wet surface, rectangular section */
		 psimulbio->section->hydro->surf->ft = pele->center->hydro->Width * pele->center->hydro->H[T_HYD];

		 psimulbio->section->hydro->vt->t = t;
		 psimulbio->section->hydro->vt->ft = pele->center->hydro->Vel;
		 //hydrad = pele->center->hydro->Surf/pele->center->hydro->Peri;
		 //hydrad = 2.5;
		 //if(psimulbio->section->hydro->hydrad != NULL){
			 //psimulbio->section->hydro->hydrad = TS_free_ts(psimulbio->section->hydro->hydrad,fp);
			// psimulbio->section->hydro->hydrad = NULL;
		 //}
		 psimulbio->section->hydro->hydrad->t = t;
		 psimulbio->section->hydro->hydrad->ft = pele->center->hydro->Surf/pele->center->hydro->Peri;
		 if(psimulbio->section->hydro->hydrad->ft < EPS_TS)
			 LP_error(psimulbio->poutputs,"t = %f surf = %f peri = %f id_section = %f",t,pele->center->hydro->Surf,pele->center->hydro->Peri,pele->id[ABS_HYD]);
	  }
  }
}

/* Calculation of the volume of water in each water compartment at time t */
//void calc_volume_water(double t)
void PROSE_calc_volume_water(double t, s_simul *Simul) // SW 26/04/2018
{
  /* Loop index */
  int nl;
  int e,j,phy;
  /* Total water height in the section */
  double hwater;
  double v_old;
  double rapp_vol;
  /* Water compartment */
  s_compartment *pcomp;
  
  //hwater = 5; //Simul->section->hydro->h->ft;
  hwater = Simul->section->hydro->h->ft;
  for (nl = 0; nl < Simul->section->nsublayers[WATER]; nl++) {

    pcomp = Simul->section->compartments[WATER][nl];
    v_old = pcomp->state->volume;
    //pcomp->state->volume_old = v_old;
    //pcomp->state->volume = (*Simul->section->description->dx) * 
    //  (*Simul->section->description->width) * hwater;
	pcomp->state->volume = (*Simul->section->description->dx) * Simul->section->hydro->surf->ft;
    pcomp->state->mass = pcomp->state->volume * pcomp->state->rho;
     rapp_vol = v_old/pcomp->state->volume;
    //pcomp->state->volume_old = pcomp->state->volume;
	 //printf("t = %f rapp_vol = %f\n",t,rapp_vol);
    /*if(v_old > 0.){
	for(e = 0; e < NSPECIES; e++) {
		for (j = 0; j < Simul->counter->nsubspecies[e]; j++) {
	  pcomp->pspecies[e][j]->C *= rapp_vol;
	  if(e == PHY){
		  for(phy = 0; phy < 3; phy++)
			 pcomp->pannex_var[phy][j]->C *= rapp_vol;
	  }
	}
      }
    }*/
  }	
}


void PROSE_calc_volume_water_all_sections(s_simul **psimul_bio, int nsections, double t, int nthreads,FILE *fp)
{
	int ns;
	s_simul *psimulbio;
	#ifdef OMP
	omp_set_num_threads(nthreads);
	#pragma omp parallel for shared(psimul_bio,nsections,t,fp) private(psimulbio)
	#endif	
	for(ns = 0; ns < nsections; ns++)
	{
		psimulbio = psimul_bio[ns];
		PROSE_calc_volume_water(t, psimulbio);
	}
}

void PROSE_link_hyd_rive_geom(s_simul **psimul_bio, s_chyd *pchyd, int nthreads,FILE *fp)
{
	int ne, r, id_abs_ele;
    s_reach_hyd *preach;
    s_element_hyd *pele;
	s_simul *psimulbio;

  for(r = 0; r < pchyd->counter->nreaches; r++)
  {
	  preach = pchyd->p_reach[r];
	  //#ifndef CDA // SW 21/11/201/
	  //#ifdef OMP
	  //omp_set_num_threads(nthreads);
	  //#pragma omp parallel for shared(psimul_bio,preach,fp) private(psimulbio,pele,id_abs_ele)
	  //#endif
	  //#endif
	  for(ne = 0; ne < preach->nele; ne++)
	  {
		 pele = preach->p_ele[ne];
		 id_abs_ele = pele->id[ABS_HYD];          		  
	     psimulbio =  psimul_bio[id_abs_ele];
		 //free(psimulbio->section->description->dx);
         //free(psimulbio->section->description->width); // SW 30/04/2012 c'est alloue dans init_description		 
		 *psimulbio->section->description->dx = pele->length;
		 *psimulbio->section->description->width = pele->center->hydro->Width;//pele->center->hydro->Surf/pele->center->hydro->H[T_HYD];
		 //psimulbio->section->description->width = &pele->center->hydro->Width;		 
	     psimulbio->section->hydro->ks = preach->strickler;
	  }
  }
}

/*calculation of water-sediemnt interface surface = wet perimeter * length*/
 
void PROSE_calc_interface_water_sed(s_simul **psimul_bio,s_chyd *pchyd, int nthreads, FILE *fp)
{
	int ne, r, id_abs_ele;
    s_reach_hyd *preach;
    s_element_hyd *pele;
	s_simul *psimulbio;
    s_interface *pint;
	
    for(r = 0; r < pchyd->counter->nreaches; r++)
    {
	  preach = pchyd->p_reach[r];
	  //#ifndef CDA // SW 21/11/2018
	//#ifdef OMP
	//omp_set_num_threads(nthreads);
	//#pragma omp parallel for shared(psimul_bio,preach,fp) private(psimulbio,pele,id_abs_ele,pint)
	//#endif
     //#endif	
	  for(ne = 0; ne < preach->nele; ne++)
	  {
		 pele = preach->p_ele[ne];
		 id_abs_ele = pele->id[ABS_HYD];          		  
	     psimulbio =  psimul_bio[id_abs_ele];
		 if(psimulbio->section->compartments[WATER][0]->down != NULL){		 
		 pint = psimulbio->section->compartments[WATER][0]->down[0];
		 pint->surface = pele->center->hydro->Width * pele->length;
		 //pint->surface = pele->center->hydro->Peri * pele->length;
         psimulbio->section->compartments[WATER][0]->state->surface = pint->surface; // SW 30/04/2018 pour l'instant qu'on considere que c'est la meme 
		 psimulbio->section->compartments[VASE][0]->state->surface = pint->surface;
		 if(pint->surface == 0 || pint->surface < EPS_TS)
			 LP_warning(fp,"yes intface null\n");
	  }
		 }
	}         		 
}

void PROSE_calc_hyd(double t,s_simul *Simul) // SW 26/04/2018
{
  double Q,vt,rh,ks,J,h;

  Q = Simul->section->hydro->q->ft;
  vt = Simul->section->hydro->vt->ft;
  rh = Simul->section->hydro->hydrad->ft;
  ks = TS_function_value_t(Q,Simul->section->hydro->ks,Simul->poutputs);
  //h = TS_function_value_t(Q,Simul->section->hydro->h,Simul->poutputs);//LV 09082012
  h = Simul->section->hydro->h->ft;
  
  //J = fabs(vt / (ks * sqrt(pow(rh,4./3.))));
 J = fabs(vt / (ks * pow(rh,2./3.))); //LV 09082012 : correspond au "sqrt(J)" de la formule de Manning-Strickler // SW 24/04/2020
  //Simul->section->hydro->Ds = J * sqrt(GR / rh);
  Simul->section->hydro->Ds = J * sqrt(GR * h);//LV 09082012
  //Simul->section->hydro->pe = Simul->section->param_ero[ETA_HYD] * J * fabs(vt);
  //Simul->section->hydro->pe = Simul->section->exchange_settings->param_ero[ETA_HYD] * J * fabs(vt);//LV 07/09/2012
  
  /*correspond au "sqrt(J)" de la formule de Manning-Strickler*/ // SW 04/12/2017
  Simul->section->hydro->pe = Simul->section->exchange_settings->param_ero[ETA_HYD] * J * J * fabs(vt);  // SW 24/04/2020 no RHO_WATER, it's included in ETA_HYD
  Simul->section->hydro->pe += Simul->section->exchange_settings->param_ero[PNAVIG];
  if(isnan(Simul->section->hydro->pe))
	  LP_warning(Simul->poutputs,"pe  = %f vt = %f, J = %f, id_section = %d eta_hyd = %f, pnavig = %f rh = %f,h = %f ks = %f Q = %f\n",Simul->section->hydro->pe,vt,J,Simul->section->id_section,Simul->section->exchange_settings->param_ero[ETA_HYD],Simul->section->exchange_settings->param_ero[PNAVIG],rh,h,ks,Q);
  //fprintf(stderr,"Q = %f, vt = %f, rh = %f, ks = %f, h = %f pe = %f\n",Q,vt,rh,ks,h,Simul->section->hydro->pe); //SW
}
