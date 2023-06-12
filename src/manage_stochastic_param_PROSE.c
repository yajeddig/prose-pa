/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: manage_stochastic_param_PROSE.c
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



/* Function to convert TOC into the six fractions  */ //MH 13/09/2021

void PROSE_toc_to_mod_mop_fract(int itype, s_simul_PROSE *Simul, s_inflow_hyd *pinflow, FILE *fp) {

    double t_param, b1_param, b2_param, s1_param, s2_param; //
    double sharemod1, sharemod2, sharemod3, sharemop1, sharemop2, sharemop3;
    int i, j, k ;
    LP_printf(Simul->poutputs,"1 \n");

    // MH 29/10/2021 [0][0] added to conform with Nparticlse and Nsections
    t_param = Simul->p_macrospecies[itype]->threshold->val;
    b1_param = Simul->p_macrospecies[itype]->degradOrgMat[MACMOD][B]->val;
    b2_param =  Simul->p_macrospecies[itype]->degradOrgMat[MACMOP][B]->val;
    s1_param =  Simul->p_macrospecies[itype]->degradOrgMat[MACMOD][S]->val;
    s2_param =  Simul->p_macrospecies[itype]->degradOrgMat[MACMOP][S]->val;
    LP_printf(Simul->poutputs,"2 \n");

    sharemod1 = t_param * b1_param * s1_param ;
    sharemod2 = t_param * b1_param * (1-s1_param);
    sharemod3 = t_param * (1-b1_param);
    double sharemod[] = {sharemod1,sharemod2, sharemod3};
    sharemop1 = (1-t_param) * b2_param * s2_param;
    sharemop2 = (1-t_param) * b2_param * (1-s2_param);
    sharemop3 = (1-t_param) * (1-b2_param);
    double sharemop[] = {sharemop1,sharemop2, sharemop3};

    // loop to distribute TOC among the six fractions for all times
    // put the pointer pft at the beginning
    s_ft *pft = pinflow->pt_inflow->flow_in_macrospecies[itype];
    // now calculate the first row of MOD1..MOP3 to be used for setting the pointer
    // they will be stored inside the flow_in[species] struct
    for (j= 0; j<Simul->counter_bio->nsubspecies[MOD]; j++)
      {
	// the time 
	 pinflow->pt_inflow->app_bio[MOD][j]->t = pinflow->pt_inflow->flow_in_macrospecies[itype]->t;
	 pinflow->pt_inflow->app_bio[MOP][j]->t = pinflow->pt_inflow->flow_in_macrospecies[itype]->t;
	
	//filling the MOD123 and MOP123 with sharemo * TOC
	 pinflow->pt_inflow->app_bio[MOD][j]->ft = pinflow->pt_inflow->flow_in_macrospecies[itype]->ft * sharemod[j];
	 pinflow->pt_inflow->app_bio[MOP][j]->ft = pinflow->pt_inflow->flow_in_macrospecies[itype]->ft * sharemop[j];
	}
    // now printing the first six fractions
    LP_printf(Simul->poutputs,"TOC:  %3.2f %3.2f \n",pft->t/86400,pft->ft/83.3);
    for (k= 0; k< Simul->counter_bio->nsubspecies[MOD]; k++)
      {
	LP_printf(fp,"MOD%d: %3.2f %3.2f \n",k+1, pinflow->pt_inflow->app_bio[MOD][k]->t/86400, pinflow->pt_inflow->app_bio[MOD][k]->ft/83.3);
	LP_printf(fp,"MOP%d: %3.2f %3.2f \n",k+1, pinflow->pt_inflow->app_bio[MOP][k]->t/86400, pinflow->pt_inflow->app_bio[MOP][k]->ft/83.3);
      } 

    // assigning the first row to the pft_mo in order to attach the rest like a linked list inside WHILE loop
    s_ft **pft_mod =  pinflow->pt_inflow->app_bio[MOD];
    s_ft **pft_mop =  pinflow->pt_inflow->app_bio[MOP];
    
    //loop for distributing the TOC into MOD1...MOP3 of the remaining time steps
    while(pft->next != NULL)
      {
	 pft = pft->next;
	 
         //creating time and conc of MOD:
	 for(j=0; j<Simul->counter_bio->nsubspecies[MOD]; j++) // separate MOD and MOP
     {
	   
	  s_ft *temp_mod = TS_create_function(pft->t,pft->ft * sharemod[j]);
	  s_ft *temp_mop = TS_create_function(pft->t,pft->ft * sharemop[j]);
	   
	  pft_mod[j]->next = TS_secured_chain_fwd_ts(pft_mod[j],temp_mod);
	  pft_mop[j]->next = TS_secured_chain_fwd_ts(pft_mop[j],temp_mop);
		 
     }  

		 // reassigning pfts to the next 	
     for(i=0; i<Simul->counter_bio->nsubspecies[MOP]; i++) 
     { 
       pft_mod[i] = pft_mod[i]->next;
	   pft_mop[i] = pft_mop[i]->next;	   
	   }

	
     // LP_printf(fp,"TOC:  %3.2f %3.2f \n",pft->t/86400,pft->ft/83.3);
     // LP_printf(fp,"MOD%d: %3.2f %3.2f \n",1,pft_mod[0]->t/86400,pft_mod[0]->ft/83.3);
     // LP_printf(fp,"MOP%d: %3.2f %3.2f \n",1,pft_mop[0]->t/86400,pft_mop[0]->ft/83.3);

     // values of mod1 and mop1 through comp2
	 // LP_printf(fp,"TOC*:  %3.2f %3.2f \n", pinflow->pt_inflow->flow_in_macrospecies[itype]->t/86400, pinflow->pt_inflow->flow_in_macrospecies[itype]->ft/83.3);
	 // LP_printf(fp,"MOD%d*: %3.2f %3.2f \n",1, pinflow->pt_inflow->app_bio[MOD][0]->t/86400, pinflow->pt_inflow->app_bio[MOD][0]->ft/83.3);
	 // LP_printf(fp,"MOP%d*: %3.2f %3.2f \n",1, pinflow->pt_inflow->app_bio[MOP][0]->t/86400, pinflow->pt_inflow->app_bio[MOP][0]->ft/83.3);
      }
     // now printing the last six fractions
    LP_printf(Simul->poutputs,"TOC:  %3.2f %3.2f \n",pft->t/86400,pft->ft/83.3);
    for (k= 0; k< Simul->counter_bio->nsubspecies[MOD]; k++)
      {
	LP_printf(fp,"MOD%d: %3.2f %3.2f \n",k+1,pft_mod[k]->t/86400,pft_mod[k]->ft/83.3);
	LP_printf(fp,"MOP%d: %3.2f %3.2f \n",k+1,pft_mop[k]->t/86400,pft_mop[k]->ft/83.3);
      }
    // moving the pointer back to the head of the linked list
     for(i=0; i<Simul->counter_bio->nsubspecies[MOP]; i++) 
     {
	     pinflow->pt_inflow->app_bio[MOD][i] = TS_browse_ft( pinflow->pt_inflow->app_bio[MOD][i],BEGINNING_TS);
	     pinflow->pt_inflow->app_bio[MOP][i] = TS_browse_ft( pinflow->pt_inflow->app_bio[MOP][i],BEGINNING_TS);
	 }

     // now printing the length of time series
     LP_printf(fp,"TOC time series length =  %d \n",TS_length_ts(pinflow->pt_inflow->flow_in_macrospecies[itype]));
     LP_printf(fp,"mod1 time series length = %d \n", TS_length_ts( pinflow->pt_inflow->app_bio[MOD][0]));
 
     //printing contents of pft_mod1     
     // for (s_ft *tmp = pinflow->pt_inflow->app_bio[MOD][0];tmp!=NULL;tmp=tmp->next) //tmp->next moves to the next n in the list
       // {
       // printf("list value: %3.2f %3.2f \n",tmp->t/86400,tmp->ft/83.3);
      // }
    
     LP_printf(fp," %s is distributed for all time steps \n",RIV_name_macrospecies(itype,fp));

}



/* Function to give the share of the six fractions for DA of b1  */ //MH 04/12/2021

double PROSE_DA_mod_mop_fract(int itype, int e, int nsub, s_simul_PROSE *Simul, int np, FILE *fp) { 

    double t_param, b1_param, b2_param, s1_param, s2_param; //
    double b1_param_temp;
    double sharemod1, sharemod2, sharemod3, sharemop1, sharemop2, sharemop3, output;

    /*
    // to include other params in the function: add int inflow_type and char* inflow_name in the arguments//

    if (inflow_type == UPSTREAM_INFLOW || inflow_type == INFLUENT){
	  b1_param = Simul->passim->param[B1_RIVER_DA][np];       
    }
    else if (inflow_type == EFFLUENT && (strcmp(inflow_name,"seine_amont")==0 || strcmp(inflow_name,"seine_amont_bp" )==0)) {
	  b1_param = Simul->passim->param[B1_WWTP_DA][np];       
    }
    else if (inflow_type == EFFLUENT && (strcmp(inflow_name,"marne_aval")==0 || strcmp(inflow_name,"marne_aval_bp" )==0)) {
	  b1_param = Simul->passim->param[B1_WWTP_DA][np];       
    }
    else if (inflow_type == EFFLUENT && (strcmp(inflow_name,"seine_centre")==0 || strcmp(inflow_name,"seine_centre_bp" )==0)) {
	  b1_param = Simul->passim->param[B1_WWTP_DA][np];       
    }
    else if (inflow_type == EFFLUENT && (strcmp(inflow_name,"seine_aval")==0 || strcmp(inflow_name,"seine_aval_bp" )==0)) {
	  b1_param = Simul->passim->param[B1_WWTP_DA][np];       
    }
    else if (inflow_type == EFFLUENT && (strcmp(inflow_name,"seine_gresillons")==0 || strcmp(inflow_name,"seine_gresillons_bp" )==0)) {
	  b1_param = Simul->passim->param[B1_WWTP_DA][np];       
    }
    else {
	  b1_param = Simul->passim->param[B1_CSO_DA][np];       
    } */
     

    t_param = Simul->p_macrospecies[itype]->threshold->val;
    b1_param_temp = Simul->passim->param[B1_RIVER_DA][np]; // before epsilon,so shall not be used
    b1_param = Simul->passim->p_macro_da[itype][np]->degradOrgMat[MACMOD][B]->val;
    //LP_printf(fp,"TOC share func) np = %i ,  b1 (stoch) = %3.4f, b1 (passim_param) = %3.4f \n", np, b1_param, b1_param_temp);

    b2_param =  Simul->p_macrospecies[itype]->degradOrgMat[MACMOP][B]->val;
    s1_param =  Simul->p_macrospecies[itype]->degradOrgMat[MACMOD][S]->val;
    s2_param =  Simul->p_macrospecies[itype]->degradOrgMat[MACMOP][S]->val;

    sharemod1 = t_param * b1_param * s1_param ;
    sharemod2 = t_param * b1_param * (1-s1_param);
    sharemod3 = t_param * (1-b1_param);
    double sharemod[] = {sharemod1,sharemod2, sharemod3};
    sharemop1 = (1-t_param) * b2_param * s2_param;
    sharemop2 = (1-t_param) * b2_param * (1-s2_param);
    sharemop3 = (1-t_param) * (1-b2_param);
    double sharemop[] = {sharemop1,sharemop2, sharemop3};

    if (e == MOD) {
        output = sharemod[nsub];
    } else if (e == MOP) {
        output = sharemop[nsub];
    } else {
    	LP_printf(fp,"This species doesn't exist \n");

    }
    //LP_printf(fp,"np = %i,e=%s,nsub=%i, b1 = %3.4f, shardmod1= %3.4f,sharemod2=%3.4f,sharedmod3=%3.4f,output=%3.4f \n",np,name_species(e),nsub,b1_param,sharemod1,sharemod2,sharemod3,output);
    return output;
}




/* Function to convert TOC into the six fractions using varying b1 */ //MH 23/05/2022

void PROSE_toc_to_mod_mop_fract_var_b1(int itype, s_simul_PROSE *Simul, s_inflow_hyd *pinflow, FILE *fp) {

    double t_param, b1_param, b2_param, s1_param, s2_param; //
    double sharemod1, sharemod2, sharemod3, sharemop1, sharemop2, sharemop3;
    int i, j, k ;
    LP_printf(Simul->poutputs,"1 \n");

    // MH 29/10/2021 [0][0] added to conform with Nparticlse and Nsections
    t_param = Simul->p_macrospecies[itype]->threshold->val;
    s_ft *p_b1 = Simul->p_macrospecies[itype]->degradOrgMat[MACMOD][B]->val_variable;
    b1_param = p_b1->ft;
    b2_param =  Simul->p_macrospecies[itype]->degradOrgMat[MACMOP][B]->val;
    s1_param =  Simul->p_macrospecies[itype]->degradOrgMat[MACMOD][S]->val;
    s2_param =  Simul->p_macrospecies[itype]->degradOrgMat[MACMOP][S]->val;
    LP_printf(Simul->poutputs,"2 \n");
    LP_printf(Simul->poutputs,"t = %f,   b1_param=%f \n", p_b1->t/NSEC_DAY_TS, b1_param);


    sharemod1 = t_param * b1_param * s1_param ;
    sharemod2 = t_param * b1_param * (1-s1_param);
    sharemod3 = t_param * (1-b1_param);
    double sharemod[] = {sharemod1,sharemod2, sharemod3};
    sharemop1 = (1-t_param) * b2_param * s2_param;
    sharemop2 = (1-t_param) * b2_param * (1-s2_param);
    sharemop3 = (1-t_param) * (1-b2_param);
    double sharemop[] = {sharemop1,sharemop2, sharemop3};

    // loop to distribute TOC among the six fractions for all times
    // put the pointer pft at the beginning
    s_ft *pft = pinflow->pt_inflow->flow_in_macrospecies[itype];
    // now calculate the first row of MOD1..MOP3 to be used for setting the pointer
    // they will be stored inside the flow_in[species] struct
    for (j= 0; j<Simul->counter_bio->nsubspecies[MOD]; j++)
      {
	// the time 
	 pinflow->pt_inflow->app_bio[MOD][j]->t = pinflow->pt_inflow->flow_in_macrospecies[itype]->t;
	 pinflow->pt_inflow->app_bio[MOP][j]->t = pinflow->pt_inflow->flow_in_macrospecies[itype]->t;

	// interpolating b1 value at time t and recalculating sharemod1,2&3
	 double time = pinflow->pt_inflow->flow_in_macrospecies[itype]->t; //NSEC_DAY_TS;
    b1_param = TS_function_value_t(time, p_b1,fp);
    LP_printf(Simul->poutputs,"b1_param_1st_interpol at t(%f)=%f \n", time/NSEC_DAY_TS, b1_param);
    sharemod1 = t_param * b1_param * s1_param ;
    sharemod2 = t_param * b1_param * (1-s1_param);
    sharemod3 = t_param * (1-b1_param);
    double sharemod[] = {sharemod1,sharemod2, sharemod3};
    
	//filling the MOD123 and MOP123 with sharemo * TOC
	 pinflow->pt_inflow->app_bio[MOD][j]->ft = pinflow->pt_inflow->flow_in_macrospecies[itype]->ft * sharemod[j];
	 pinflow->pt_inflow->app_bio[MOP][j]->ft = pinflow->pt_inflow->flow_in_macrospecies[itype]->ft * sharemop[j];
	}
    // now printing the first six fractions
    LP_printf(Simul->poutputs,"TOC:  %3.2f %3.2f \n",pft->t/86400,pft->ft/83.3);
    for (k= 0; k< Simul->counter_bio->nsubspecies[MOD]; k++)
      {
	LP_printf(fp,"MOD%d: %3.2f %3.2f \n",k+1, pinflow->pt_inflow->app_bio[MOD][k]->t/86400, pinflow->pt_inflow->app_bio[MOD][k]->ft/83.3);
	LP_printf(fp,"MOP%d: %3.2f %3.2f \n",k+1, pinflow->pt_inflow->app_bio[MOP][k]->t/86400, pinflow->pt_inflow->app_bio[MOP][k]->ft/83.3);
      } 

    // assigning the first row to the pft_mo in order to attach the rest like a linked list inside WHILE loop
    s_ft **pft_mod =  pinflow->pt_inflow->app_bio[MOD];
    s_ft **pft_mop =  pinflow->pt_inflow->app_bio[MOP];
    
    //loop for distributing the TOC into MOD1...MOP3 of the remaining time steps
    while(pft->next != NULL)
      {
	 pft = pft->next;

    //reinterpolation of b1 and sharemod1,2&3
    double time_i = pft->t; //NSEC_DAY_TS;
    b1_param = TS_function_value_t(time_i, p_b1,fp);
    LP_printf(Simul->poutputs,"b1_param_new_interpol at t(%f)=%f \n",time_i/NSEC_DAY_TS, b1_param);
    sharemod1 = t_param * b1_param * s1_param ;
    sharemod2 = t_param * b1_param * (1-s1_param);
    sharemod3 = t_param * (1-b1_param);
    double sharemod[] = {sharemod1,sharemod2, sharemod3};
	 
         //creating time and conc of MOD pairs:
	 for(j=0; j<Simul->counter_bio->nsubspecies[MOD]; j++) // separate MOD and MOP
     {
	   
	  s_ft *temp_mod = TS_create_function(pft->t,pft->ft * sharemod[j]);
	  s_ft *temp_mop = TS_create_function(pft->t,pft->ft * sharemop[j]);
	   
	  pft_mod[j]->next = TS_secured_chain_fwd_ts(pft_mod[j],temp_mod);
	  pft_mop[j]->next = TS_secured_chain_fwd_ts(pft_mop[j],temp_mop);
		 
     }  

		 // reassigning pfts to the next 	
     for(i=0; i<Simul->counter_bio->nsubspecies[MOP]; i++) 
     { 
       pft_mod[i] = pft_mod[i]->next;
	   pft_mop[i] = pft_mop[i]->next;	   
	   }

	
     // LP_printf(fp,"TOC:  %3.2f %3.2f \n",pft->t/86400,pft->ft/83.3);
     // LP_printf(fp,"MOD%d: %3.2f %3.2f \n",1,pft_mod[0]->t/86400,pft_mod[0]->ft/83.3);
     // LP_printf(fp,"MOP%d: %3.2f %3.2f \n",1,pft_mop[0]->t/86400,pft_mop[0]->ft/83.3);

     // values of mod1 and mop1 through comp2
	 // LP_printf(fp,"TOC*:  %3.2f %3.2f \n", pinflow->pt_inflow->flow_in_macrospecies[itype]->t/86400, pinflow->pt_inflow->flow_in_macrospecies[itype]->ft/83.3);
	 // LP_printf(fp,"MOD%d*: %3.2f %3.2f \n",1, pinflow->pt_inflow->app_bio[MOD][0]->t/86400, pinflow->pt_inflow->app_bio[MOD][0]->ft/83.3);
	 // LP_printf(fp,"MOP%d*: %3.2f %3.2f \n",1, pinflow->pt_inflow->app_bio[MOP][0]->t/86400, pinflow->pt_inflow->app_bio[MOP][0]->ft/83.3);
      }
     // now printing the last six fractions
    LP_printf(Simul->poutputs,"TOC:  %3.2f %3.2f \n",pft->t/86400,pft->ft/83.3);
    for (k= 0; k< Simul->counter_bio->nsubspecies[MOD]; k++)
      {
	LP_printf(fp,"MOD%d: %3.2f %3.2f \n",k+1,pft_mod[k]->t/86400,pft_mod[k]->ft/83.3);
	LP_printf(fp,"MOP%d: %3.2f %3.2f \n",k+1,pft_mop[k]->t/86400,pft_mop[k]->ft/83.3);
      }
    // moving the pointer back to the head of the linked list
     for(i=0; i<Simul->counter_bio->nsubspecies[MOP]; i++) 
     {
	     pinflow->pt_inflow->app_bio[MOD][i] = TS_browse_ft( pinflow->pt_inflow->app_bio[MOD][i],BEGINNING_TS);
	     pinflow->pt_inflow->app_bio[MOP][i] = TS_browse_ft( pinflow->pt_inflow->app_bio[MOP][i],BEGINNING_TS);
	 }

     // now printing the length of time series
     LP_printf(fp,"TOC time series length =  %d \n",TS_length_ts(pinflow->pt_inflow->flow_in_macrospecies[itype]));
     LP_printf(fp,"mod1 time series length = %d \n", TS_length_ts( pinflow->pt_inflow->app_bio[MOD][0]));
 
     //printing contents of pft_mod1     
     // for (s_ft *tmp = pinflow->pt_inflow->app_bio[MOD][0];tmp!=NULL;tmp=tmp->next) //tmp->next moves to the next n in the list
       // {
       // printf("list value: %3.2f %3.2f \n",tmp->t/86400,tmp->ft/83.3);
      // }
    
     LP_printf(fp," %s is distributed for all time steps \n",RIV_name_macrospecies(itype,fp));

}
