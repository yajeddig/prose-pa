/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: calc_outputs.c
* BRANCH NAME: main
* 
* CONTRIBUTORS: Shuaitao WANG, Lauriane VILMIN, Aur�lien BORDET, Masihullah HASANYAR, Thomas ROMARY, Nicolas FLIPO
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
#include <string.h>
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
/* Creates the files in which the mass balances will be printed */
//void create_files_mb()
void PROSE_create_files_mb_bio(int nparticules, int np, FILE *fp) // SW 26/04/2018
{
  /* Loop indexes */
  int layer,e,j,i,p;
  s_simul *psimulbio;
  psimulbio = Simul->psimul_bio[np][0]; 
  
  for (layer = 0; layer < NLAYERS; layer++) {
    if (psimulbio->section->compartments[layer] != NULL) { 
      
      for (e = 0; e < NSPECIES; e++) {
	Simul->mass_balances[layer][e][np] = (FILE ***)calloc(Simul->counter_bio->nsubspecies[e],sizeof(FILE **));
	for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++)
	  Simul->mass_balances[layer][e][np][j] = (FILE **)calloc(Simul->counter_bio->nmb,sizeof(FILE *));
      }
      
      for (e = 0; e < NANNEX_VAR; e++) {
	Simul->mass_balances[layer][NSPECIES+e][np] = (FILE ***)calloc(Simul->counter_bio->nsubannex_var[e],sizeof(FILE **));
	for (j = 0; j < Simul->counter_bio->nsubannex_var[e]; j++)
	  Simul->mass_balances[layer][NSPECIES+e][np][j] = (FILE **)calloc(Simul->counter_bio->nmb,sizeof(FILE *));
      }      
      /*
      for (e = 0; e < NDISS; e++) {
	Simul->ads_mass_balances[layer][e][np] = (FILE ***)calloc(Simul->counter_bio->nsubspecies[e+NPART],sizeof(FILE **));
	for (j = 0; j < Simul->counter_bio->nsubspecies[e+NPART]; j++)
	  Simul->ads_mass_balances[layer][e][np][j] = (FILE **)calloc(Simul->counter_bio->nmb,sizeof(FILE *));
      }*/

      for (i = 0; i < Simul->counter_bio->nmb; i++) {
	char s[4];
	sprintf(s,"%d",i+1);
	
	for (e = 0; e < NSPECIES; e++) {
	  for(j = 0; j < Simul->counter_bio->nsubspecies[e]; j++) {
	    if (Simul->total_mb[np][i]->calc_mb_species[e][j] == YES_RIVE) { // SW 14/10/2018 here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done
	      char filename[MAXCHAR_PROSE];
	      char num_spe[4];
		  char *name_spe;
		  char *name_lay;
		  char np_string[5];
		  char nmb_string[5];
		  name_spe = name_species(e);
		  name_lay = name_layer(layer);
	      sprintf(num_spe,"%d",j+1);
	      sprintf(filename,"%s/mass_balance_bio/",getenv("RESULT"));
	      strcat(filename,"mass_balance_");
	      strcat(filename,s);
	      strcat(filename,"_");
	      strcat(filename,name_lay);
	      strcat(filename,"_");
	      strcat(filename,name_spe);
	      strcat(filename,"_");
	      strcat(filename,num_spe);
		  if(nparticules > 1)
		  {
		     strcat(filename,"_");
			 sprintf(np_string, "%d", np); 
			 strcat(filename,np_string); // SW 17/09/2018
		  }	
	      strcat(filename,"\0");
	      printf("%s\n",filename);
	      free(name_spe);
		  free(name_lay);
		  name_spe = name_lay = NULL;
	      if ((Simul->mass_balances[layer][e][np][j][i] = fopen(filename,"w")) == NULL) 
		printf("problem when opening the file %s\n",filename),exit(3);
	      
              fprintf(Simul->mass_balances[layer][e][np][j][i],"Date\t"); /* SW 06/01/2021 print date */
              fprintf(Simul->mass_balances[layer][e][np][j][i],"day\t");

	      for (p = HYDROLYSIS; p < NVAR_MB; p++) {
		//LV 15/11/2013 : pour avoir des bilans de N plus détaillés
		if ((e == NH4) && (p == RESP))
		fprintf(Simul->mass_balances[layer][e][np][j][i],"NITROS\t");
		else if ((e == NO2) && (p == GROWTH))
		fprintf(Simul->mass_balances[layer][e][np][j][i],"NITROS\t");
		else if ((e == NO2) && (p == RESP))
		fprintf(Simul->mass_balances[layer][e][np][j][i],"NITRAT\t");
		else if ((e == NO3) && (p == GROWTH))
		fprintf(Simul->mass_balances[layer][e][np][j][i],"NITRAT\t");
		else if ((e == NO3) && (p == RESP))
		fprintf(Simul->mass_balances[layer][e][np][j][i],"DENIT\t");
		else
		fprintf(Simul->mass_balances[layer][e][np][j][i],"%s\t",name_process(p));
	      }
	      fprintf(Simul->mass_balances[layer][e][np][j][i],"\n");
	    }
	  }
	}
	
	for (e = 0; e < NANNEX_VAR; e++) {
	  for(j = 0; j < Simul->counter_bio->nsubannex_var[e]; j++) {
	    if (Simul->total_mb[np][i]->calc_mb_annex_var[e][j] == YES_RIVE) { // SW 14/10/2018 here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done
	      char filename[MAXCHAR_PROSE];
	      char num_ann[4];
		  char np_string[5];
		  char *name_lay;
		  char *name_annex;
		  name_lay = name_layer(layer);
		  name_annex = name_annex_var(e);
	      sprintf(num_ann,"%d",j+1);
	      sprintf(filename,"%s/mass_balance_bio/",getenv("RESULT"));
	      strcat(filename,"mass_balance_");
	      strcat(filename,s);
	      strcat(filename,"_");
	      strcat(filename,name_lay);
	      strcat(filename,"_");
	      strcat(filename,name_annex);
	      strcat(filename,"_");
	      strcat(filename,num_ann);
		  if(nparticules > 1)
		  {
		     strcat(filename,"_");
			 sprintf(np_string, "%d", np); 
			 strcat(filename,np_string); // SW 17/09/2018
		  }			  
	      strcat(filename,"\0");
	      printf("%s\n",filename);
	      free(name_lay);
		  free(name_annex);
		  name_lay = name_annex = NULL;
	      if ((Simul->mass_balances[layer][NSPECIES+e][np][j][i] = fopen(filename,"w")) == NULL) 
		printf("problem when opening the file %s\n",filename),exit(3);

              fprintf(Simul->mass_balances[layer][NSPECIES+e][np][j][i],"Date\t"); /* SW 06/01/2021 print date */
	      fprintf(Simul->mass_balances[layer][NSPECIES+e][np][j][i],"day\t");
	      for (p = HYDROLYSIS; p < NVAR_MB; p++) {
		fprintf(Simul->mass_balances[layer][NSPECIES+e][np][j][i],"%s\t",name_process(p));
	      }	    
	      fprintf(Simul->mass_balances[layer][NSPECIES+e][np][j][i],"\n");
	    }
	  }
	}
    /*
	for (e = 0; e < NDISS; e++) {
	  for(j = 0; j < Simul->counter_bio->nsubspecies[e+NPART]; j++) {
	    if (Simul->total_mb[np][i]->calc_mb_adsorbed_species[e][j] == YES_TS) {
	      char filename[MAXCHAR_PROSE];
	      char num_spe[4];
	      sprintf(num_spe,"%d",j+1);
	      sprintf(filename,"%s/mass_balance_bio/",getenv("RESULT"));
	      strcat(filename,"adsorbed_mass_balance_");
	      strcat(filename,s);
	      strcat(filename,"_");
	      strcat(filename,name_layer(layer));
	      strcat(filename,"_");
	      strcat(filename,name_species(e+NPART));
	      strcat(filename,"_");
	      strcat(filename,num_spe);
		  if(nparticules > 1)
		  {
		     strcat(filename,"_");
			 //sprintf(np_string, "%d", np); 
			 strcat(filename,np_string); // SW 17/09/2018
		  }
	      strcat(filename,"\0");
	      printf("%s\n",filename);
	      
	      if ((Simul->ads_mass_balances[layer][e][np][j][i] = fopen(filename,"w")) == NULL) 
		printf("problem when opening the file %s\n",filename),exit(3);
	      fprintf(Simul->ads_mass_balances[layer][e][np][j][i],"t\tMINI\tMEND");
	      fprintf(Simul->ads_mass_balances[layer][e][np][j][i],"\n");
	    }
	  }
	}*/
      }
    } 
  }
}

void PROSE_create_files_weights(int nparticules, FILE *fp) // SW 26/04/2018
{
	int np;
	char filename[MAXCHAR_PROSE];
	char filename2[MAXCHAR_PROSE];
	char cmd[MAXCHAR_PROSE];
	
	//sprintf(cmd,"mkdir %s/weights",getenv("RESULT"));
	//system(cmd);
        sprintf(cmd,"%s/weights",getenv("RESULT"));
        IO_mkdir(cmd); // SW 19/10/2022 mkdir recurssive to avoid issues for dummies, system may has problems sometimes

	sprintf(filename,"%s/weights/",getenv("RESULT"));
	strcat(filename,"weights_all_particules");
	strcat(filename,"\0");
	sprintf(filename2,"%s/weights/",getenv("RESULT"));
	strcat(filename2,"resampling_size");
	strcat(filename2,"\0");
	if ((Simul->passim->pout_weight = fopen(filename,"w")) == NULL)
	    printf("problem when opening the file %s\n",filename),exit(3);
	if((Simul->passim->pout_sampling_size = fopen(filename2,"w")) == NULL)
		LP_error(fp,"can not open resampling size file\n");
    //fprintf(Simul->passim->pout_weight,"Date\t"); /* SW 06/01/2021 print date */
    fprintf(Simul->passim->pout_weight,"day\t");
    for(np = 0; np < nparticules; np++)
       fprintf(Simul->passim->pout_weight, "%d\t",np+1);
    fprintf(Simul->passim->pout_weight, "\n");   
}

void PROSE_print_weights(double t, int nparticules, FILE *fp)
{
	int np;
	if(t <= Simul->chronos->t[END_CHR]) // SW 03/02/2020
	{
	fprintf(Simul->passim->pout_weight,"%f\t",t);
	for(np = 0; np < nparticules; np++)
	{
		fprintf(Simul->passim->pout_weight,"%f\t",Simul->passim->omega_normlized[np]);
		//LP_printf(fp,"print weight func) np = %i: omega_prev = %3.4f, omega = %3.4f, omega_norm = %3.4f \n ",np, Simul->passim->omega_prev[np], Simul->passim->omega[np], Simul->passim->omega_normlized[np]);
	}
	fprintf(Simul->passim->pout_weight,"\n");
	fflush(Simul->passim->pout_weight);
	}
}

void PROSE_create_files_parameters(int nparticules, FILE *fp) // SW 26/04/2018
{
	int np, nparam;
	char filename[MAXCHAR_PROSE];
	char cmd[MAXCHAR_PROSE];
	char *np_string;
	//sprintf(cmd,"mkdir %s/parameters",getenv("RESULT"));
	//system(cmd);
        sprintf(cmd,"%s/parameters",getenv("RESULT"));
	IO_mkdir(cmd); // SW 19/10/2022 mkdir recurssive to avoid issues for dummies, system may has problems sometimes

	
	Simul->passim->pout_param = (FILE **) malloc(NPARAMDA * sizeof(FILE *));
	for(nparam = 0; nparam < NPARAMDA; nparam++)
	{
            // SW 25/01/2022 check if parameter is assimilated
            if(Simul->passim->param_yesOrno[nparam] == YES_TS)
            {
	        //sprintf(np_string, "%d",np);
	        sprintf(filename,"%s/parameters/",getenv("RESULT"));
	        strcat(filename,"parameters_");		
		np_string = PROSE_name_param(nparam);
		strcat(filename,np_string);
		strcat(filename,"\0");
	        if((Simul->passim->pout_param[nparam] = fopen(filename,"w")) == NULL)
	            printf("problem when opening the file %s\n",filename),exit(3);		
        
                //fprintf(Simul->passim->pout_param[nparam],"Date\t");
                fprintf(Simul->passim->pout_param[nparam],"day\t");
	        for(np = 0; np < nparticules; np++)
		    fprintf(Simul->passim->pout_param[nparam], "%d\t",np+1);
	        fprintf(Simul->passim->pout_param[nparam], "\n");
	        bzero(filename, MAXCHAR_PROSE);
	    }
        }             
}

/* Creates the files in which the concentration outputs will be printed */
//void create_files_conc()
void PROSE_transv_profile_format_bio(s_chyd *pchyd,s_output_hyd ***p_outputs, int np, FILE *fp) // SW 26/04/2018
{
  /* Loop indexes */
  int layer,j,i;
  int npk, nts;
  s_ts_pk_hyd *ppk;
  s_simul *psimulbio;
  
  psimulbio = Simul->psimul_bio[np][0]; 

  for (layer = 0; layer < NLAYERS; layer++) { 

    if (psimulbio->section->compartments[layer] != NULL) { 
	    nts = pchyd->counter->nts;
	    Simul->concentrations[layer][np] = (FILE ***)calloc(nts,sizeof(FILE **));      	
        for (j = 0; j < nts; j++) {
            npk = p_outputs[TRANSV_PROFILE][j]->npk;  
			Simul->concentrations[layer][np][j] = (FILE **)calloc(npk,sizeof(FILE *));			
            for (i = 0; i < npk; i++) {
			   ppk = p_outputs[TRANSV_PROFILE][j]->ts_pk[i];
	           char filename[MAXCHAR_PROSE];
	           sprintf(filename,"%s/time_series/ts_pk%4.2f_%s_branch%d_%s_moyenne_%d.txt",getenv("RESULT"),ppk->pk / 1000.,ppk->river,ppk->branch_nb,name_layer(layer),np+1);	      
	      if ((Simul->concentrations[layer][np][j][i] = fopen(filename,"w")) == NULL) 
		      LP_error(fp,"problem when opening the file %s\n",filename);
	      fprintf(Simul->concentrations[layer][np][j][i],"#File generated by ProSe%4.2f\n",NVERSION_PROSE);
          
          fprintf(Simul->concentrations[layer][np][j][i], "#Date "); /* SW 06/01/2021 print Date */
          fprintf(Simul->concentrations[layer][np][j][i], "day ");
          PROSE_print_variables(Simul->concentrations[layer][np][j][i],p_outputs[TRANSV_PROFILE][j]);
		  //fclose(Simul->concentrations[layer][j][i]);
          /* SW 23/01/2024 to print sediment layer variables such as volume, height.*/
          if(layer == VASE)
              PROSE_print_sediment_variables(Simul->concentrations[layer][np][j][i],p_outputs[TRANSV_PROFILE][j]);
   
          fprintf(Simul->concentrations[layer][np][j][i],"\n"); // SW 24/01/2024 add here this function
	    }
	  }
    }
  }	
}

void PROSE_print_variables(FILE *fout,s_output_hyd *pout_hyd)
{
  int e, nsub;
  
  for (e = 0; e < NSPECIES; e++){
    if (pout_hyd->pout->biovar[e] == YES_TS){
		for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++){
            fprintf(fout,"%s %d ",name_species(e), nsub+1);			
		}
	}
  }
// SW 04/06/2019 print annex_variables  
  for(e = 0; e < NANNEX_VAR; e++)
  {
	  if(Simul->calc_bio_annexvar[e] == YES_TS)
	  {
		  for(nsub = 0; nsub < Simul->counter_bio->nsubannex_var[e]; nsub++)
		  {  
              fprintf(fout,"%s %d ",name_annex_var(e), nsub+1); 
		  }
	  }
  }
  //fprintf(fout,"\n"); // SW 24/01/2024 move after call of this function
}

/* SW 23/01/2024 to print sediment layer variables such as volume, height.*/
void PROSE_print_sediment_variables(FILE *fout,s_output_hyd *pout_hyd)
{
  int e;
  
  for (e = 0; e < NSEDVAR; e++){
    if (pout_hyd->pout->biosedvar[e] == YES_TS)
        fprintf(fout,"%s ",name_sedvar(e));
  }
  //fprintf(fout,"\n"); // SW 24/01/2024 move after call of this function
}

void PROSE_calc_mb_minit_all_sections(s_simul **psimul_bio, double t, double dt, int nsections, int np, FILE *fp)
{
	int ns, layer, e, j, i,nl,ne;
    /* Porosity */
    double v_poral = 1.;
    /* Current compartment */  
    s_compartment *pcomp;
    /* Total mass balance */
    s_total_mb *pmb;
	s_simul *psimulbio;
	s_lp_pk_hyd *ppk;
	s_element_hyd *pele;
	s_reach_hyd *preach;
	s_chyd *pchyd;
	int nr;
	//int tid = omp_get_thread_num();
	//LP_printf(fp,"tid = %d np = %d\n",tid,np);
	for (i = 0; i < Simul->counter_bio->nmb; i++) {
		pmb = Simul->total_mb[np][i];
		if((t >= pmb->t[BEGINNING]) && (t <= pmb->t[END])){ // SW 30/01/2019
    //if ((t <= pmb->t0) && (t + dt >= pmb->t0)) {
		/*** SW 06/12/2019 add mass_balance output for the defined domaine by user***/
		if(Simul->lp_pk[i] != NULL) // output pk defined
		{
			ppk = Simul->lp_pk[i];
			nr = 0;
			pchyd = Simul->pchyd;
			if(Simul->lp_pk[i]->reach_nb[0] == -1)
				HYD_calculate_output_domain(ppk,MASS_BALANCE,pchyd,fp);
            preach = NULL;
            
			while((ppk->reach_nb[nr] >= 0) &&(nr < pchyd->counter->nreaches)) 
			{
				preach = pchyd->p_reach[ppk->reach_nb[nr]];				
				for (ne = 0; ne < preach->nele; ne++) 
				{
					pele = preach->p_ele[ne];					
				    if ((pele->center->pk >= ppk->pk_up) && (pele->center->pk <= ppk->pk_down)) 
					{
						ns = pele->id[ABS_HYD];
		                psimulbio = psimul_bio[ns];	
						for (layer = 0; layer < NLAYERS; layer++) 
						{
							if (psimulbio->section->compartments[layer] != NULL) 
							{
							/*reinitialization of the mass balances of species */
								for(e = 0; e < NSPECIES; e++) 
								{
									for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++) 
									{
										if ((pmb->calc_mb_species[e][j] == YES_RIVE)) 
										{ // SW 14/10/2018 here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done// on calcul MINI seulement apres printing ou reinitialisation
											for (nl = 0; nl < psimulbio->section->nsublayers[layer]; nl++) 
											{ 
												pcomp = psimulbio->section->compartments[layer][nl];
												if (layer > WATER)
													v_poral = pcomp->state->phi;
												else
													v_poral = 1.;
												if (e >= NPART)
												{
													pmb->mbspecies[layer][MINI][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume * v_poral; // SW 21/10/2018 volume_old for MINIT
													if(pmb->mbspecies[layer][MINI][e][j] < 0.)
														LP_error(fp,"error after sedimentation_erosion: ns = %d e = %d j = %d layer = %d vol = %f v_poral = %f\n",ns,e,j,layer,pcomp->state->volume,v_poral);
												}
												else
												{
													pmb->mbspecies[layer][MINI][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume; // SW 21/10/2018 volume_old for MINIT
													if(pmb->mbspecies[layer][MINI][e][j] < 0.)
													LP_error(fp,"error after sedimentation_erosion: ns = %d e = %d j = %d layer = %d c = %f vol = %f v_poral = %f\n",ns,e,j,layer,pcomp->pspecies[e][j]->C,pcomp->state->volume,v_poral);
												}
											}
										}
									} // for j = 0
								} // for e = 0

							}
						}
					}
				} // for ne = 0
				nr++; // next reach_nb
			}					
		}
		/*** SW 06/12/2019 add mass_balance output for the defined domaine by user end***/
    else {	// SW 06/12/2019 if no domaine pk defined, calculate the entire domaine				
	for(ns = 0; ns < nsections; ns++)
	{
		psimulbio = psimul_bio[ns];
		//pmb = psimulbio->total_mb[i];		
	    for (layer = 0; layer < NLAYERS; layer++) {
	       if (psimulbio->section->compartments[layer] != NULL) {
		   /*reinitialization of the mass balances of species */
		     for(e = 0; e < NSPECIES; e++) {
				 for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++) {
					if ((pmb->calc_mb_species[e][j] == YES_RIVE)) { // SW 14/10/2018 here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done// on calcul MINI seulement apres printing ou reinitialisation
		               for (nl = 0; nl < psimulbio->section->nsublayers[layer]; nl++) { 
		                   pcomp = psimulbio->section->compartments[layer][nl];

		               if (layer > WATER)
		                   v_poral = pcomp->state->phi;
		               else
		                   v_poral = 1.;
		    
		               if (e >= NPART){
		                   pmb->mbspecies[layer][MINI][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume * v_poral; // SW 21/10/2018 volume_old for MINIT
						   //pmb->mbspecies[layer][MINI][e][j] += pcomp->pspecies[e][j]->newC * pcomp->state->volume_old * v_poral;
	                       if(pmb->mbspecies[layer][MINI][e][j] < 0.)
							   LP_error(fp,"error after sedimentation_erosion: ns = %d e = %d j = %d layer = %d vol = %f v_poral = %f\n",ns,e,j,layer,pcomp->state->volume,v_poral);
					       //Simul->total_mb[i]->mbspecies[layer][MINI][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume * v_poral; // sum of all section
					   }
		               else{
		                   pmb->mbspecies[layer][MINI][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume; // SW 21/10/2018 volume_old for MINIT
						   //pmb->mbspecies[layer][MINI][e][j] += pcomp->pspecies[e][j]->newC * pcomp->state->volume_old;
	                       if(pmb->mbspecies[layer][MINI][e][j] < 0.)
							   LP_error(fp,"error after sedimentation_erosion: ns = %d e = %d j = %d layer = %d c = %f vol = %f v_poral = %f\n",ns,e,j,layer,pcomp->pspecies[e][j]->C,pcomp->state->volume,v_poral);
					       //if(e == 0 && j == 0 && layer == 0)
							   //LP_printf(fp,"ok\n");
						   //Simul->total_mb[i]->mbspecies[layer][MINI][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume; // sum of all section
					   }
					   }
					}
				 }
			 }
 	          /*for (e = 0; e < NANNEX_VAR; e++) {
	               for (j = 0; j < Simul->counter_bio->nsubannex_var[e]; j++) {
					  if (pmb->calc_mb_annex_var[e][j] == YES && pmb->mbannex[layer][MINI][e][j] < EPS) { // on calcul MINI seulement apres printing ou reinitialisation
		                 for (nl = 0; nl < psimulbio->section->nsublayers[layer]; nl++) {
		                  pcomp = psimulbio->section->compartments[layer][nl];
                          pmb->mbannex[layer][MINI][e][j] += pcomp->pannex_var[e][j]->C * pcomp->state->volume;
						 }
					  }
				   }
			  }*/
		   }
		}
		  }	
	}
	}
	}
}


void PROSE_calc_mb_sections(double t, s_simul **psimul_bio, int nele, int np, FILE *fp)
{
  /* Loop indexes */
  int i,layer,e,p,j,nl,ns,nr,ne;
  /* Porosity */
  double v_poral = 1.;
  /* Current compartment */  
  s_compartment *pcomp;
  /* Total mass balance */
  s_total_mb *pmb;
	
	s_simul *psimulbio;
    s_lp_pk_hyd *ppk;
	s_element_hyd *pele;
	s_reach_hyd *preach;
	s_chyd *pchyd;
	
  for (i = 0; i < Simul->counter_bio->nmb; i++) {
    //pmb = psimulbio->total_mb[i];
	pmb = Simul->total_mb[np][i];
    if((t >= pmb->t[BEGINNING]) && (t <= pmb->t[END])){	
  
		/*** SW 06/12/2019 add mass_balance output for the defined domaine by user***/
		if(Simul->lp_pk[i] != NULL) // output pk defined
		{
			ppk = Simul->lp_pk[i];
			nr = 0;
			pchyd = Simul->pchyd;
			if(Simul->lp_pk[i]->reach_nb[0] == -1)
				HYD_calculate_output_domain(ppk,MASS_BALANCE,pchyd,fp);
            preach = NULL;
            
			while((ppk->reach_nb[nr] >= 0) &&(nr < pchyd->counter->nreaches)) 
			{
				preach = pchyd->p_reach[ppk->reach_nb[nr]];				
				for (ne = 0; ne < preach->nele; ne++) 
				{
					pele = preach->p_ele[ne];
				    if ((pele->center->pk >= ppk->pk_up) && (pele->center->pk <= ppk->pk_down)) 
					{
						ns = pele->id[ABS_HYD];
		                psimulbio = psimul_bio[ns];
						//if(ns == 421)
							//LP_printf(fp,"debug in calc_ouputs.c id = 421\n");
						for (layer = 0; layer < NLAYERS; layer++) 
						{
							for (nl = 0; nl < psimulbio->section->nsublayers[layer]; nl++) 
							{
								pcomp = psimulbio->section->compartments[layer][nl];	  
								if (layer > WATER)
									v_poral = pcomp->state->phi;
								for(p = 0; p < MINI; p++)
								{
									for (e = 0; e < NSPECIES; e++) 
									{
										for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++)
										{
											if (pmb->calc_mb_species[e][j] == YES_RIVE) 
											{ // SW 14/10/2018 here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done
											pmb->mbspecies[layer][p][e][j] += pcomp->pspecies[e][j]->mb->deltamb[p];
											}
										}
									}
								}
								/* SW 11/05/2020 I add here EXCEPTION term*/
								for(e = 0; e < NSPECIES; e++)
								{
									for(j = 0; j < Simul->counter_bio->nsubspecies[e]; j++)
									{
										if(pmb->calc_mb_species[e][j] == NO_TS)
											pmb->mbspecies[layer][EXCEPTION][e][j] += pcomp->pspecies[e][j]->mb->deltamb[EXCEPTION];
									}
								}
										
							}
	             /* for (e = 0; e < NANNEX_VAR; e++) {
	              for (j = 0; j < Simul->counter_bio->nsubannex_var[e]; j++) {
		              if (pmb->calc_mb_annex_var[e][j] == YES_TS) {
		                 pmb->mbannex[layer][p][e][j] += pcomp->pannex_var[e][j]->mb->deltamb[p];
		                 //Simul->total_mb[i]->mbannex[layer][p][e][j] += pcomp->pannex_var[e][j]->mb->deltamb[p];
				   }
				  }
				  }*/
						}
					}
				}
				nr++;
			}
		}
	  else{			
      /* Calculation of the mass balances for entire domaine */	  
	  for(ns = 0; ns < nele; ns++){
	  psimulbio = psimul_bio[ns];
      for (layer = 0; layer < NLAYERS; layer++) {
	      for (nl = 0; nl < psimulbio->section->nsublayers[layer]; nl++) {
	          pcomp = psimulbio->section->compartments[layer][nl];	  
	          if (layer > WATER)
	              v_poral = pcomp->state->phi;
			  for(p = 0; p < MINI; p++){
				  for (e = 0; e < NSPECIES; e++) {
					  for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++){
						 if (pmb->calc_mb_species[e][j] == YES_RIVE) { // SW 14/10/2018 here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done
                             pmb->mbspecies[layer][p][e][j] += pcomp->pspecies[e][j]->mb->deltamb[p];
						     //Simul->total_mb[i]->mbspecies[layer][p][e][j] += pcomp->pspecies[e][j]->mb->deltamb[p]; // total_mb a faire
						 }
				  }
				  }
	             /* for (e = 0; e < NANNEX_VAR; e++) {
	              for (j = 0; j < Simul->counter_bio->nsubannex_var[e]; j++) {
		              if (pmb->calc_mb_annex_var[e][j] == YES_TS) {
		                 pmb->mbannex[layer][p][e][j] += pcomp->pannex_var[e][j]->mb->deltamb[p];
		                 //Simul->total_mb[i]->mbannex[layer][p][e][j] += pcomp->pannex_var[e][j]->mb->deltamb[p];
				   }
				  }
				  }*/
			  }
			/* SW 11/05/2020 I add here EXCEPTION term*/
			for(e = 0; e < NSPECIES; e++)
			{
				for(j = 0; j < Simul->counter_bio->nsubspecies[e]; j++)
				{
					if(pmb->calc_mb_species[e][j] == NO_TS)
						pmb->mbspecies[layer][EXCEPTION][e][j] += pcomp->pspecies[e][j]->mb->deltamb[EXCEPTION];
				}
			}			  
		  }
	  }
	  }
	  }
  }	
  }  
}

void PROSE_calc_mb_mend_all_sections(s_simul **psimul_bio, double t, double dt, int nsections, int np, FILE *fp)
{
	int ns, layer, e, j, i,nl,p,nr,ne;
    /* Porosity */
    double v_poral = 1.;
	//double t0;
    /* Current compartment */  
    s_compartment *pcomp;
    /* Total mass balance */
    s_total_mb *pmb;
	s_simul *psimulbio;
    s_lp_pk_hyd *ppk;
	s_element_hyd *pele;
	s_reach_hyd *preach;
	s_chyd *pchyd;
	
	for (i = 0; i < Simul->counter_bio->nmb; i++) {
	  pmb = Simul->total_mb[np][i];
	  if((t >= pmb->t[BEGINNING]) && (t <= pmb->t[END])){ // SW 30/01/2019
	// if ((t <= pmb->t0) && (t + dt >= pmb->t0)) {
	
		/*** SW 06/12/2019 add mass_balance output for the defined domaine by user***/
		if(Simul->lp_pk[i] != NULL) // output pk defined
		{
			ppk = Simul->lp_pk[i];
			nr = 0;
			pchyd = Simul->pchyd;
			if(Simul->lp_pk[i]->reach_nb[0] == -1)
				HYD_calculate_output_domain(ppk,MASS_BALANCE,pchyd,fp);
            preach = NULL;
            
			while((ppk->reach_nb[nr] >= 0) &&(nr < pchyd->counter->nreaches)) 
			{
				preach = pchyd->p_reach[ppk->reach_nb[nr]];				
				for (ne = 0; ne < preach->nele; ne++) 
				{
					pele = preach->p_ele[ne];					
				    if ((pele->center->pk >= ppk->pk_up) && (pele->center->pk <= ppk->pk_down)) 
					{
						ns = pele->id[ABS_HYD];
		                psimulbio = psimul_bio[ns];	
						for (layer = 0; layer < NLAYERS; layer++) {
							if (psimulbio->section->compartments[layer] != NULL) {
								for (e = 0; e < NSPECIES; e++) {  
									for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++) {
										if (pmb->calc_mb_species[e][j] == YES_RIVE) {	// SW 14/10/2018 here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done	 			 
											for (nl = 0; nl < psimulbio->section->nsublayers[layer]; nl++) 
											{ 
												pcomp = psimulbio->section->compartments[layer][nl];
												if (layer > WATER)
													v_poral = pcomp->state->phi;
												else
													v_poral = 1.;		    
												if (e >= NPART){
													pmb->mbspecies[layer][MEND][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume * v_poral;
												}
												else{
													pmb->mbspecies[layer][MEND][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume;
												}
											}
										}
									}
								}
	         /* Printing and reinitialization of the mass balances of annex_variables */
	          /*for (e = 0; e < NANNEX_VAR; e++) {
	              for (j = 0; j < Simul->counter_bio->nsubannex_var[e]; j++) {
                       if (pmb->calc_mb_annex_var[e][j] == YES_TS) {
		                   for (nl = 0; nl < psimulbio->section->nsublayers[layer]; nl++) {
		                        pcomp = psimulbio->section->compartments[layer][nl];
						        pmb->mbannex[layer][MEND][e][j] += pcomp->pannex_var[e][j]->C * pcomp->state->volume;
								//Simul->total_mb[i]->mbannex[layer][MEND][e][j] += pcomp->pannex_var[e][j]->C * pcomp->state->volume;
						   }
						pmb->mbannex[layer][ERROR][e][j] = pmb->mbannex[layer][MEND][e][j] - pmb->mbannex[layer][MINI][e][j];   
			            for (p = HYDROLYSIS; p < MINI; p++)
							pmb->mbannex[layer][ERROR][e][j] -= pmb->mbannex[layer][p][e][j];
						//Simul->total_mb[i]->mbannex[layer][ERROR][e][j] += pmb->mbannex[layer][ERROR][e][j];
		                if(fabs(pmb->mbannex[layer][ERROR][e][j]) > 1. && t > pmb->t[BEGINNING] ) // SW 08/03/2017 ajouter erreur liée aux concentrations negatives. Il faut aojuter processus exception voir dans prose
		                {
			               pmb->mbannex[layer][EXCEPTION][e][j] = pmb->mbannex[layer][ERROR][e][j];
			               pmb->mbannex[layer][ERROR][e][j] -= pmb->mbannex[layer][EXCEPTION][e][j];
		                   LP_warning(Simul->poutputs,"id_section = %d layer = %d bilan faut!!!!!!!annex_variables especes %d error = %e\n",ns+1,layer,e,pmb->mbannex[layer][EXCEPTION][e][j]);
		                }						
						/*reinitialization of pmb for all sections*/
						//for (p = HYDROLYSIS; p < NVAR_MB; p++) 
							//pmb->mbannex[layer][p][e][j] = 0.;					   
					  // }
				  //}
			 // }*/
							}  
						} // for layer
					}
				} // for ne
				nr++;
			}
		}
		else{				
	for(ns = 0; ns < nsections; ns++)
	{
		psimulbio = psimul_bio[ns];
		//pmb = psimulbio->total_mb[i];
	    for (layer = 0; layer < NLAYERS; layer++) {
	       if (psimulbio->section->compartments[layer] != NULL) {
		     for (e = 0; e < NSPECIES; e++) {  
                 for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++) {
                    if (pmb->calc_mb_species[e][j] == YES_RIVE) {	// SW 14/10/2018 here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done	 			 
		                for (nl = 0; nl < psimulbio->section->nsublayers[layer]; nl++) { 
		                    pcomp = psimulbio->section->compartments[layer][nl];
		                    if (layer > WATER)
		                        v_poral = pcomp->state->phi;
		                    else
		                        v_poral = 1.;		    
		                   if (e >= NPART){
		                         pmb->mbspecies[layer][MEND][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume * v_poral;
		                         //if(e == O2 && j == 0 && layer == WATER)
									 //LP_printf(fp,"id_section = %d t = %f mb end  = %f c= %f vol = %f\n",ns+1,t,pmb->mbspecies[layer][MEND][e][j],pcomp->pspecies[e][j]->C,pcomp->state->volume);							   
                               //if(np == 1 && ns == 0 && e == O2)
								   //LP_printf(fp,"np == 1\n");
							   //if(np == 0 && ns == 0 && e == O2)
								   //LP_printf(fp,"np == 0\n");
								 //Simul->total_mb[i]->mbspecies[layer][MEND][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume * v_poral;
						   }
						   else {
		                         pmb->mbspecies[layer][MEND][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume;

						         //Simul->total_mb[i]->mbspecies[layer][MEND][e][j] += pcomp->pspecies[e][j]->C * pcomp->state->volume;
						   }
						}
						  					

						
					}
				 }
			 }
	         /* Printing and reinitialization of the mass balances of annex_variables */
	          /*for (e = 0; e < NANNEX_VAR; e++) {
	              for (j = 0; j < Simul->counter_bio->nsubannex_var[e]; j++) {
                       if (pmb->calc_mb_annex_var[e][j] == YES_TS) {
		                   for (nl = 0; nl < psimulbio->section->nsublayers[layer]; nl++) {
		                        pcomp = psimulbio->section->compartments[layer][nl];
						        pmb->mbannex[layer][MEND][e][j] += pcomp->pannex_var[e][j]->C * pcomp->state->volume;
								//Simul->total_mb[i]->mbannex[layer][MEND][e][j] += pcomp->pannex_var[e][j]->C * pcomp->state->volume;
						   }
						pmb->mbannex[layer][ERROR][e][j] = pmb->mbannex[layer][MEND][e][j] - pmb->mbannex[layer][MINI][e][j];   
			            for (p = HYDROLYSIS; p < MINI; p++)
							pmb->mbannex[layer][ERROR][e][j] -= pmb->mbannex[layer][p][e][j];
						//Simul->total_mb[i]->mbannex[layer][ERROR][e][j] += pmb->mbannex[layer][ERROR][e][j];
		                if(fabs(pmb->mbannex[layer][ERROR][e][j]) > 1. && t > pmb->t[BEGINNING] ) // SW 08/03/2017 ajouter erreur liée aux concentrations negatives. Il faut aojuter processus exception voir dans prose
		                {
			               pmb->mbannex[layer][EXCEPTION][e][j] = pmb->mbannex[layer][ERROR][e][j];
			               pmb->mbannex[layer][ERROR][e][j] -= pmb->mbannex[layer][EXCEPTION][e][j];
		                   LP_warning(Simul->poutputs,"id_section = %d layer = %d bilan faut!!!!!!!annex_variables especes %d error = %e\n",ns+1,layer,e,pmb->mbannex[layer][EXCEPTION][e][j]);
		                }						
						/*reinitialization of pmb for all sections*/
						//for (p = HYDROLYSIS; p < NVAR_MB; p++) 
							//pmb->mbannex[layer][p][e][j] = 0.;					   
					  // }
				  //}
			 // }*/
			}
		}
	  }
		}	  
	  
	for (layer = 0; layer < NLAYERS; layer++) {
	  for(e = 0; e < NSPECIES; e++)
	{
	  for(j = 0; j < Simul->counter_bio->nsubspecies[e]; j++)
	 {
      if (pmb->calc_mb_species[e][j] == YES_RIVE) {	// SW 30/04/2019 add this check here YES_RIVE = 1, which is equivalent to YES_RIVE in librive where the intialization is done	 
	  pmb->mbspecies[layer][ERROR][e][j] = pmb->mbspecies[layer][MEND][e][j] - pmb->mbspecies[layer][MINI][e][j];
      for (p = HYDROLYSIS; p < MINI; p++) 
		  pmb->mbspecies[layer][ERROR][e][j] -= pmb->mbspecies[layer][p][e][j];
	  //Simul->total_mb[i]->mbspecies[layer][ERROR][e][j] += pmb->mbspecies[layer][ERROR][e][j];
	  if((fabs(pmb->mbspecies[layer][ERROR][e][j]) / pmb->unit_mb_species[e][j]) > 1. &&t > pmb->t[BEGINNING]) // SW 08/03/2017 ajouter erreur liée aux concentrations negatives. Il faut aojuter processus exception voir dans prose
      {
		  LP_warning(Simul->poutputs,"bilan faut!!!!!!! check the location of two PKs. mb number = %d, layer = %d especes %s sub = %d error = %e exception = %f np = %d\n",i+1,layer,name_species(e),j,pmb->mbspecies[layer][ERROR][e][j],pmb->mbspecies[layer][EXCEPTION][e][j],np+1);
		  //pmb->mbspecies[layer][EXCEPTION][e][j] = pmb->mbspecies[layer][ERROR][e][j];
		  //pmb->mbspecies[layer][ERROR][e][j] -= pmb->mbspecies[layer][EXCEPTION][e][j];
      }	
	 }
	 }
	}
	}
	//}
	
	//if((t <= pmb->t0) && (t + dt >= pmb->t0)) {
	//}	
}
	}
}

void PROSE_print_mb_bio_domaine(s_simul **psimul_bio, double t, double dt, int nsections, int np, FILE *fp)
{
	int ns = 0, layer, e, j, i, p;
    /* Porosity */
    //double v_poral = 1.;
	//double t0;
    /* Current compartment */  
    //s_compartment *pcomp;
    /* Total mass balance */
    s_total_mb *pmb;
	s_simul *psimulbio;
	
	for (i = 0; i < Simul->counter_bio->nmb; i++) {
	  pmb = Simul->total_mb[np][i];
	  if((t >= pmb->t[BEGINNING]) && (t <= pmb->t0) && 
	  (t + dt > pmb->t0)&& (t <= pmb->t[END])){ // SW 30/01/2019
	// if ((t <= pmb->t0) && (t + dt >= pmb->t0)) {
	//for(ns = 0; ns < nsections; ns++)
	//{
		psimulbio = psimul_bio[ns];
		//pmb = psimulbio->total_mb[i];
	
	for (layer = 0; layer < NLAYERS; layer++) {
	  if (psimulbio->section->compartments[layer] != NULL) { // n'import lequelle psimulbio va marcher, ici le dernier        
		 for(e = 0; e < NSPECIES; e++)
		 {
			 for(j = 0; j < Simul->counter_bio->nsubspecies[e]; j++)
			 {
				 if(pmb->calc_mb_species[e][j] == YES_RIVE) // SW 14/10/2018 here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done
				 {
				     /* SW 06/01/2021 print date */
                                     TS_print_date(Simul->chronos->pd[CUR_CHR],Simul->date_format,Simul->mass_balances[layer][e][np][j][i], Simul->poutputs);
                                     fprintf(Simul->mass_balances[layer][e][np][j][i], " ");
                                     TS_print_time(Simul->chronos->pd[CUR_CHR],Simul->mass_balances[layer][e][np][j][i], Simul->poutputs); 
                                     fprintf(Simul->mass_balances[layer][e][np][j][i], "\t");
                                     fprintf(Simul->mass_balances[layer][e][np][j][i], "%f\t",t / pmb->time_unit);
				     for (p = HYDROLYSIS; p < NVAR_MB; p++) {
						 fprintf(Simul->mass_balances[layer][e][np][j][i], "%e\t",pmb->mbspecies[layer][p][e][j]/pmb->unit_mb_species[e][j]);
					 pmb->mbspecies[layer][p][e][j] = 0.0;	//reinitialisation of the sum terms 
				     }
					 fprintf(Simul->mass_balances[layer][e][np][j][i], "\n");
				 }
			 }
		 }
	    /* Printing and reinitialization of the mass balances of annex_variables */
	    /*for (e = 0; e < NANNEX_VAR; e++) {
	      for (j = 0; j < Simul->counter_bio->nsubannex_var[e]; j++) {
		if (pmb->calc_mb_annex_var[e][j] == YES_TS) {

                        //SW 06/01/2021 print date
                        TS_print_date(Simul->chronos->pd[CUR_CHR],Simul->date_format,Simul->mass_balances[layer][NSPECIES+e][np][j][i], Simul->poutputs);
                        fprintf(Simul->mass_balances[layer][NSPECIES+e][np][j][i], " ");
                        TS_print_time(Simul->chronos->pd[CUR_CHR],Simul->mass_balances[layer][NSPECIES+e][np][j][i], Simul->poutputs); 
                        fprintf(Simul->mass_balances[layer][NSPECIES+e][np][j][i], "\t");

			fprintf(Simul->mass_balances[layer][NSPECIES+e][np][j][i], "%f\t",t / pmb->time_unit);
			for (p = HYDROLYSIS; p < NVAR_MB; p++)
				fprintf(Simul->mass_balances[layer][NSPECIES+e][np][j][i], "%e\t",pmb->mbannex[layer][p][e][j]/pmb->unit_mb_annex_var[e][j]);
			pmb->mbannex[layer][p][e][j] = 0; //reinitialisation of the sum terms
			fprintf(Simul->mass_balances[layer][NSPECIES+e][j][i], "\n");
		}
		  }
		}*/
		/*reinitialization of pmb for all sections*/
		//for (p = HYDROLYSIS; p < NVAR_MB; p++) 
		   //pmb->mbspecies[layer][p][e][j] = 0.;
	  }
	}
	pmb->t0 += pmb->ndt;
	}
	//}
	}	
}


void PROSE_print_outputs_formats_bio(s_output_hyd ***p_outputs,s_inout_set_io *pinout,s_chyd *pchyd,int np, FILE *fp)
{ 
  if (pinout->calc[TRANSV_PROFILE] == YES_TS) {
    PROSE_transv_profile_format_bio(pchyd,p_outputs,np,fp);
  }
}

/* SW 30/03/2021 pas seulement bio, mais aussi HYD, TUBE */
void PROSE_print_outputs_bio(double t,s_simul **psimul_bio,s_output_hyd ***p_outputs,s_inout_set_io *pinout,s_chyd *pchyd,s_chronos_CHR *chronos, int np, FILE *fp)
{
  double t0;
  s_output_hyd *result;
  s_output_tube_type *result_tube;
  int i, layer;
  s_element_hyd *pele;
  s_ts_pk_hyd *ppk;
  int ipk;

  if (pinout->calc[TRANSV_PROFILE] == YES_TS) {

    for (i = 0; i < pchyd->counter->nts; i++) {
      result = p_outputs[TRANSV_PROFILE][i];
	  //result->pout->t_out[CUR_IO] = t; // SW 30/01/2018 CUR_IO, n'est pas fait
      t0 = result->pout->t_out[CUR_IO];
	  
      if ((t >= result->pout->t_out[BEGINNING]) && (t <= t0) && 
	  (t + chronos->dt > t0) && (t <= result->pout->t_out[END])) {  
      if(np == 0) // SW 20/05/2019 only one time
      {
          /* SW 06/01/2021 print date */
          for (ipk = 0; ipk < result->npk; ipk++) {

              ppk = result->ts_pk[ipk];
              pele = ppk->pele;
    
              if (pele != NULL) {
      
              if ((ppk->pfic = fopen(ppk->file_name,"a")) == 0) 
	          LP_error(fp,"file %s\n",ppk->file_name);
            
              TS_print_date(Simul->chronos->pd[CUR_CHR],Simul->date_format,ppk->pfic, Simul->poutputs);
              fprintf(ppk->pfic, " ");
              TS_print_time(Simul->chronos->pd[CUR_CHR],ppk->pfic, Simul->poutputs); 
              fprintf(ppk->pfic, "\t");
    
              fclose(ppk->pfic);
              }
           }
           /* SW 06/01/2021 end print date */

	   HYD_transv_profile(t,result,fp);// SW 20/05/2019 print hyd variables

        }
	  if(Simul->calc_mode[RIVE] == YES_TS)
	  {
	    for (layer = 0; layer < NLAYERS ; layer++) { 
	       PROSE_transv_profile_bio(t,result,psimul_bio,layer, i, np,fp);
		}
	  }
	  //if(np == (Simul->passim->N_particules - 1) && (Simul->calc_mode[H_T] == NO_TS)) // SW 20/05/2019 only for the last simulation / SW 26/01/2021 update at end t, see PROSE_update_outputs_hyd_bio_heat_cur_io
	  	//result->pout->t_out[CUR_IO] += result->pout->deltat;	  	
	  }
      }
    }

	// SW 18/05/2019 add longitudinal profile
  if ((pinout->calc[LONG_PROFILE] == YES_TS)) {
	//if(Simul->calc_mode[DA] == NO_TS) { // SW 18/05/2019 no LONG_PROFILE for data assimilation by default. If not, too files opened. 

    for (i = 0; i < pchyd->counter->nlp; i++) {
		//LP_printf(fp,"nlp = %d\n",pchyd->counter->nlp);
      result = p_outputs[LONG_PROFILE][i];
	  //result->pout->t_out[CUR_IO] = t;
      t0 = result->pout->t_out[CUR_IO];
      //LP_printf(fp, "before entering long output t0 = %f t = %f\n", t0, t);

      if ((t >= result->pout->t_out[BEGINNING]) &&(t <= t0) && 
	  (t + chronos->dt > t0) && (t <= result->pout->t_out[END])) {
	  if(np == 0) // SW 20/05/2019 only one simulation
	     HYD_long_profile(t,result,pinout,pchyd,fp); 	// SW 20/05/2019 print hyd variables		  
	  if((Simul->calc_mode[RIVE] == YES_TS) && (Simul->calc_mode[DA] == NO_TS))
	  {
      for (layer = 0; layer < NLAYERS; layer++) {	
	  //LP_printf(fp,"layer = %d",layer);
	  PROSE_long_profile(psimul_bio,t,result,layer,pinout,pchyd,fp);
	  }
	  }
	  else if((Simul->calc_mode[DA] == YES_TS) && (t == result->pout->t_out[BEGINNING]))
	  {
		  LP_warning(fp,"LONG_PROFILE is asked by user, but Data assimilation is calculated. \nThis can open too many files.\nYou can modify the code in calc_outputs.c for LONG_PROFILE files.\n");
	  }	  
	  //if(np == (Simul->passim->N_particules - 1) && (Simul->calc_mode[H_T] == NO_TS)) // SW 20/05/2019 only for the last simulation  / SW 26/01/2021 update at end t, see PROSE_update_outputs_hyd_bio_heat_cur_io
	  	//result->pout->t_out[CUR_IO] += result->pout->deltat; 
      }
    }
	//}
  }
  
  /* SW 30/03/2021 add here tube MESH output and tube HYD output */
  if(Simul->pout_tube->calc_output == YES_TS)
  {
     for(i = 0; i < Simul->pout_tube->nout; i++) 
     {
     
         result_tube = Simul->pout_tube->poutput_tube_type[i];
         t0 = result_tube->t_out[CUR_IO];
         //LP_printf(fp, "before entering tube output t0 = %f t = %f\n", t0, t);
         if((t >= result_tube->t_out[BEGINNING]) && (t <= t0) && (t + chronos->dt > t0) && (t <= result_tube->t_out[END])) 
         {
            //LP_printf(fp, "enter tube output\n");
            if(np == 0)
                PROSE_print_tube_mesh_hyd(t, result_tube, pchyd, Simul->p3_tube, fp);  
 
            //LP_printf(fp, "exit tube output\n");
 
         }
     }
  }

  if (pinout->calc[MASS_BALANCE] == YES_TS) {

    for (i = 0; i < pchyd->counter->nmb; i++) {
      result = p_outputs[MASS_BALANCE][i];
	  //result->pout->t_out[CUR_IO] = t; // SW 30/01/2018 CUR_IO, n'est pas fait
      t0 = result->pout->t_out[CUR_IO];
      
	  if(np == 0) // SW 20/05/2019 only one simulation
	  {
      if ((t >= result->pout->t_out[BEGINNING]) && (t <= result->pout->t_out[END])) {		
	      HYD_mass_balance(t,result,pchyd,chronos,fp);
	      //result->pout->t_out[CUR_IO] += result->pout->deltat; // SW 26/01/2021 update at end t, see PROSE_update_outputs_hyd_bio_heat_cur_io
      }
	  }
    }
  }  
}

// SW 20/05/2019 update t_out[CUR_IO] at end of time loop

void PROSE_update_outputs_hyd_bio_heat_cur_io(double t,s_output_hyd ***p_outputs,s_inout_set_io *pinout,s_chyd *pchyd,s_chronos_CHR *chronos, int np, FILE *fp)
{
  double t0;
  s_output_hyd *result;
  s_output_tube_type *result_tube;
  int i, layer;
  s_element_hyd *pele;
  s_ts_pk_hyd *ppk;
  int ipk;

  if (pinout->calc[TRANSV_PROFILE] == YES_TS) {

    for (i = 0; i < pchyd->counter->nts; i++) {
      result = p_outputs[TRANSV_PROFILE][i];
	  //result->pout->t_out[CUR_IO] = t; // SW 30/01/2018 CUR_IO, n'est pas fait
      t0 = result->pout->t_out[CUR_IO];
	  
      if ((t >= result->pout->t_out[BEGINNING]) && (t <= t0) && 
	  (t + chronos->dt > t0) && (t <= result->pout->t_out[END])) {  
	  if(np == 0) 
	  	result->pout->t_out[CUR_IO] += result->pout->deltat;	  	
	  }
      }
    }

	// SW 18/05/2019 add longitudinal profile
  if ((pinout->calc[LONG_PROFILE] == YES_TS)) {
	//if(Simul->calc_mode[DA] == NO_TS) { // SW 18/05/2019 no LONG_PROFILE for data assimilation by default. If not, too files opened. 

    for (i = 0; i < pchyd->counter->nlp; i++) {
		//LP_printf(fp,"nlp = %d\n",pchyd->counter->nlp);
      result = p_outputs[LONG_PROFILE][i];
	  //result->pout->t_out[CUR_IO] = t;
      t0 = result->pout->t_out[CUR_IO];
      if ((t >= result->pout->t_out[BEGINNING]) &&(t <= t0) && 
	  (t + chronos->dt > t0) && (t <= result->pout->t_out[END])) { 
	  if(np == 0) // SW 20/05/2019 only for the last simulation
	  	result->pout->t_out[CUR_IO] += result->pout->deltat;
      }
    }
	//}
  }

  /* SW 30/03/2021 add here tube MESH output and tube HYD output */
  if(Simul->pout_tube->calc_output == YES_TS)
  {
     for(i = 0; i < Simul->pout_tube->nout; i++) 
     {
     
       result_tube = Simul->pout_tube->poutput_tube_type[i];
       t0 = result_tube->t_out[CUR_IO];
       if((t >= result_tube->t_out[BEGINNING]) && (t <= t0) && (t + chronos->dt > t0) && (t <= result_tube->t_out[END])) 
       {
         if(np == 0)
          result_tube->t_out[CUR_IO] += result_tube->dt;    
       }
     }
  }


  if (pinout->calc[MASS_BALANCE] == YES_TS) {

    for (i = 0; i < pchyd->counter->nmb; i++) {
      result = p_outputs[MASS_BALANCE][i];
	  //result->pout->t_out[CUR_IO] = t; // SW 30/01/2018 CUR_IO, n'est pas fait
      t0 = result->pout->t_out[CUR_IO];
      
	  if(np == 0) // SW 20/05/2019 only one simulation
	  {
      if ((t >= result->pout->t_out[BEGINNING]) && (t <= result->pout->t_out[END])) {		
	      result->pout->t_out[CUR_IO] += result->pout->deltat;
      }
	  }
    }
  }  
}


void PROSE_transv_profile_bio(double t,s_output_hyd *result, s_simul **psimul_bio, int layer, int i, int np, FILE *fp)
{
  s_element_hyd *pele;
  s_ts_pk_hyd *ppk;
  int j;
  //int id_abs;
  s_simul *psimulbio;
  
  psimulbio = psimul_bio[0];

  for (j = 0; j < result->npk; j++) {

    ppk = result->ts_pk[j];
    pele = ppk->pele;
    //id_abs = ppk->pele->id[ABS_HYD];
    if (pele != NULL && psimulbio->section->compartments[layer] != NULL) {

      /* SW 06/01/2021 print date */
      TS_print_date(Simul->chronos->pd[CUR_CHR],Simul->date_format,Simul->concentrations[layer][np][i][j], Simul->poutputs);
      fprintf(Simul->concentrations[layer][np][i][j], " ");
      TS_print_time(Simul->chronos->pd[CUR_CHR],Simul->concentrations[layer][np][i][j], Simul->poutputs); 
      fprintf(Simul->concentrations[layer][np][i][j], ";");
            
      fprintf(Simul->concentrations[layer][np][i][j],"%f;",t * result->pout->time_unit);
      //PROSE_print_average_result_bio(pele,psimul_bio, layer, pele->center->pk,result,Simul->concentrations[layer][i][j]);
	  PROSE_print_average_result_bio(pele,psimul_bio, Simul->pchyd, layer, ppk->pk,ppk->river,result,Simul->concentrations[layer][np][i][j]);
      
      //fclose(ppk->pfic);
    }
  }
}


void  PROSE_print_average_result_bio(s_element_hyd *pele, s_simul **psimul_bio, s_chyd *pchyd, int layer, double pk, char *river, s_output_hyd *result,FILE *pfic)
{
  int e, nsub, yes_no;
  int id_abs, id_abs_n;
  int nl = 0,n, ndown_limits; // one sublayer
  int phyfsr;
  double dx_bio;
  double val_e, val_n;
  s_element_hyd *pelen;
  double val;

  //s_face_hyd *pface;
  
  //pface = pele->face[X_HYD][TWO];
  //pelen = pele->face[X_HYD][TWO]->element[TWO];
  id_abs = pele->id[ABS_HYD];
  
  /* SW 16/12/2020 check upstream or downstream elements */
  n = HYD_find_river(river,TRANSV_PROFILE,pchyd,Simul->poutputs);
  // SW 11/02/2022 bug for Marne river
  //ndown_limits = pchyd->counter->ndownstr_limits; 
  if((fabs(pk - pchyd->upstr_reach[n]->limits[UPSTREAM]->pk) > EPS_TS) && (fabs(pk - pchyd->downstr_reach[0]->limits[DOWNSTREAM]->pk) > EPS_TS)) // pk is not the upstream and downstream pk
  {

      if(pk > pele->center->pk)
      {
	  pelen = pele->face[X_HYD][TWO]->element[TWO]; //element aval
	  if(pelen != NULL)
	      dx_bio = (pk - pele->center->pk)/(pelen->center->pk - pele->center->pk);
	  else
	  {
		  pelen = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO];
          dx_bio = (pk - pele->center->pk)/(pelen->center->pk - pele->center->pk);		  
	  }
	  id_abs_n = pelen->id[ABS_HYD];
      }
      else
      {
	  pelen = pele->face[X_HYD][ONE]->element[ONE]; //element amont
	  if(pelen != NULL) // SW 18/05/2019 face libre
	      dx_bio = (pele->center->pk - pk)/(pele->center->pk - pelen->center->pk);
	  else // SW 18/05/2019 face limits
	  {
		  pelen = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]; // element amont
		  dx_bio = (pele->center->pk - pk)/(pele->center->pk - pelen->center->pk);
	  }   
	  id_abs_n = pelen->id[ABS_HYD];
      }

      //id_abs_n = pelen->id[ABS_HYD];
      //if (pface->pk != pele->center->pk){
      // dx_bio = (pk - pele->center->pk) / (pelen->center->pk - pele->center->pk);
      //}
      //else {
	//dx_bio = 0.;
      //}

      dx_bio = dx_bio > 1.0 ? 1.0 : dx_bio;
      dx_bio = dx_bio < 0.0 ? 0.0 : dx_bio;
  }
  else
  {
     id_abs_n = id_abs; // pk is the downstream or upstream pk
     dx_bio = 0.;
  }

  /* SW 16/12/2020 check upstream or downstream elements */
  n = HYD_find_river(river,TRANSV_PROFILE,pchyd,Simul->poutputs);
  if(fabs(pk - pchyd->upstr_reach[n]->limits[UPSTREAM]->pk) < EPS_TS) // pk is the upstream pk
  {
     id_abs_n = id_abs;
     dx_bio = 0.;
  }
  if(fabs(pk - pchyd->downstr_reach[0]->limits[DOWNSTREAM]->pk) < EPS_TS) // pk is the downstream pk
  {
     id_abs_n = id_abs;
     dx_bio = 0.;
  }

  // SW 28/05/2018 pour l'instant one sublayer
if((psimul_bio[id_abs]->section->compartments[layer] != NULL) && (psimul_bio[id_abs_n]->section->compartments[layer] != NULL)) 
{	
  for(e = 0; e < NSPECIES; e++){
	  if(result->pout->biovar[e] == YES_TS){
	      for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++){
		     val_e = psimul_bio[id_abs]->section->compartments[layer][nl]->pspecies[e][nsub]->C;
		     val_n = psimul_bio[id_abs_n]->section->compartments[layer][nl]->pspecies[e][nsub]->C;
		     val = (1.0 - dx_bio) * val_e + dx_bio * val_n;
			 if(isnan(val))
				 fprintf(pfic,"%e;",val_e / result->pout->biovar_unit[e]); // SW 11/10/2021 just for consacre project add ;
				 //LP_error(Simul->poutputs,"Val is nan in PROSE_print_average_result_bio, pk_pele_center = %f pk_pelen_center = %f\n",pele->center->pk,pelen->center->pk);
			 else
		         fprintf(pfic,"%e;",val / result->pout->biovar_unit[e]); // SW 11/10/2021 just for consacre project add ;	
            //fprintf(pfic,"%g ",val_e / result->pout->biovar_unit[e]);			 
             yes_no = YES_TS;
			}
	  }
  }
  // SW 04/06/2019 print annex_variables
  for(e = 0; e < NANNEX_VAR; e++)
  {
	  if(Simul->calc_bio_annexvar[e] == YES_TS)
	  {
		  for(nsub = 0; nsub < Simul->counter_bio->nsubannex_var[e]; nsub++)
		  {
			  val_e = psimul_bio[id_abs]->section->compartments[layer][nl]->pannex_var[e][nsub]->C;
			  val_n = psimul_bio[id_abs_n]->section->compartments[layer][nl]->pannex_var[e][nsub]->C;
			  val = (1.0 - dx_bio) * val_e + dx_bio * val_n;
			  if(isnan(val))
				  fprintf(pfic,"%e;",val_e / Simul->unit_bio_annexvar[e]);
				  //LP_error(Simul->poutputs,"Val is nan in PROSE_print_average_result_bio, pk_pele_center = %f pk_pelen_center = %f\n",pele->center->pk,pelen->center->pk);
			  else
		          fprintf(pfic,"%e;",val / Simul->unit_bio_annexvar[e]);				  
		  }
	  }
  }

  //SW 23/01/2024 print sediment variables such as volume, height
  if(layer == VASE)
  {
      for(e = 0; e < NSEDVAR_IO; e++)
      {
          if(result->pout->biosedvar[e] == YES_TS)
          {
              switch(e) 
             {
                  case SED_VOL_IO : {
                     val_e = psimul_bio[id_abs]->section->compartments[VASE][0]->state->volume;
                     break;
                  }
                  case SED_H_IO : {
                     if(psimul_bio[id_abs]->section->compartments[VASE][0]->state->surface > EPS_TS)
                         val_e = psimul_bio[id_abs]->section->compartments[VASE][0]->state->volume / psimul_bio[id_abs]->section->compartments[VASE][0]->state->surface;
                     else
                         val_e = 0.;
                     break;
                  }
                  default : {
                  val_e = 0.; // unknown variables
                  break;
                  }
              }
          fprintf(pfic,"%e;",val_e / result->pout->biosedvar_unit[e]);
        }
    } // end for
  } // end if

  fprintf(pfic,"\n");
  fflush(pfic);
}
}
			  
void PROSE_print_conc_init(int nsections, int np, s_inout_set_io *pinout, FILE *fpout)
{
	FILE *fpc;
	int ns, nsub, e;
	s_species *pspecies;
	char *name_spe;
	char conc_init_path[MAXCHAR_PROSE];
	sprintf(conc_init_path,"%s/conc_init",pinout->name_outputs);
	if((fpc = fopen(conc_init_path,"w")) == 0)
		LP_error(fpout,"Impossible to open conc init state file\n");
   for (e = 0; e < NSPECIES; e++){
		for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++){
			pspecies = Simul->psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[e][nsub];
            name_spe = name_species(e);
            if(pspecies != NULL)
			 fprintf(fpc,"%s %d ",name_species(e), nsub+1);	
            free(name_spe);	 
		}
	}
	
	fprintf(fpc,"\n");
	
	for(ns = 0; ns < nsections; ns++)
	{
		fprintf(fpc,"%d ", ns+1);
		for (e = 0; e < NSPECIES; e++){
			for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++){
				pspecies = Simul->psimul_bio[np][ns]->section->compartments[WATER][0]->pspecies[e][nsub];
				if(pspecies != NULL)
					fprintf(fpc,"%f ",pspecies->C);	
			}
		}
		fprintf(fpc,"\n");
	}
	fclose(fpc);		
}


void PROSE_print_conc_volume_final(int nsections, s_inout_set_io *pinout, FILE *fpout)
{
	FILE *fpc_water,*fpc_vase, *fpc_volume_vase;
	char conc_final_path[2][MAXCHAR_PROSE];//0 pour Z et 1 pour U
	char vol_final_path[MAXCHAR_PROSE];
	int ns, nsub, e;
	s_species *pspecies_water;
	s_species *pspecies_vase;
	char *name_spe;
	sprintf(conc_final_path[0],"%s/conc_final_water",pinout->name_outputs);

        //LP_printf(fpout, "nsublayers = %d \n",Simul->psimul_bio[0][0]->section->nsublayers[VASE]);
        // SW 14/01/2022
        if(Simul->psimul_bio[0][0]->section->nsublayers[VASE] > 0)
        {
	    sprintf(conc_final_path[1],"%s/conc_final_vase",pinout->name_outputs);
	    sprintf(vol_final_path,"%s/volume_final_vase",pinout->name_outputs);
        }
	if((fpc_water = fopen(conc_final_path[0],"w")) == 0)
	    LP_error(fpout,"Impossible to open concentration water final state file\n");
        
        // SW 14/01/2022
        if(Simul->psimul_bio[0][0]->section->nsublayers[VASE] > 0)
        {
	    if((fpc_vase = fopen(conc_final_path[1],"w")) == 0)
	        LP_error(fpout,"Impossible to open concentration vase final state file\n");
	    if((fpc_volume_vase = fopen(vol_final_path,"w")) == 0)
		LP_error(fpout,"Impossible to open volume vase final state file\n");
        }
   for (e = 0; e < NSPECIES; e++){
		for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++){
			pspecies_water = Simul->psimul_bio[0][0]->section->compartments[WATER][0]->pspecies[e][nsub];
            //name_spe = name_species(e);
            if(pspecies_water != NULL)
			{
			 fprintf(fpc_water,"%s %d ",name_species(e), nsub+1);
                         // SW 14/01/2022
                         if(Simul->psimul_bio[0][0]->section->nsublayers[VASE] > 0)
	                     fprintf(fpc_vase,"%s %d ",name_species(e), nsub+1);
			}
            //free(name_spe);	 
		}
	}
	// SW 04/06/2019  // SW 30/07/2019 comment print annex variables
	//for(e = 0; e < ACTBACT; e++){
		//for(nsub = 0; nsub < Simul->counter_bio->nsubannex_var[e]; nsub++){
			//pspecies_water = Simul->psimul_bio[0][0]->section->compartments[WATER][0]->pannex_var[e][nsub];
		    //name_spe = name_annex_var(e);
            //if(pspecies_water != NULL)
			//{
			 //fprintf(fpc_water,"%s %d ",name_annex_var(e), nsub+1);	
			 //fprintf(fpc_vase,"%s %d ",name_annex_var(e), nsub+1);
			//}
	//}
	//}
	
	fprintf(fpc_water,"\n");

        // SW 14/01/2022
        if(Simul->psimul_bio[0][0]->section->nsublayers[VASE] > 0)
	    fprintf(fpc_vase,"\n");

	for(ns = 0; ns < nsections; ns++)
	{
		fprintf(fpc_water,"%d ", ns+1);

                // SW 14/01/2022
                if(Simul->psimul_bio[0][0]->section->nsublayers[VASE] > 0)
                {
		    fprintf(fpc_vase,"%d ", ns+1);
		    fprintf(fpc_volume_vase,"%d ", ns+1);
                
	    //if(Simul->psimul_bio[0][ns]->section->compartments[VASE][0] != NULL)
		    fprintf(fpc_volume_vase,"%f ",Simul->psimul_bio[0][ns]->section->compartments[VASE][0]->state->volume); // unit m3
                }
		for (e = 0; e < NSPECIES; e++){
			for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++){
				pspecies_water = Simul->psimul_bio[0][ns]->section->compartments[WATER][0]->pspecies[e][nsub];
				if(pspecies_water != NULL)
					fprintf(fpc_water,"%f ",pspecies_water->C);
				//if(Simul->psimul_bio[0][ns]->section->compartments[VASE][0] != NULL)
                                // SW 14/01/2022
                                if(Simul->psimul_bio[0][0]->section->nsublayers[VASE] > 0)
				{
					pspecies_vase = Simul->psimul_bio[0][ns]->section->compartments[VASE][0]->pspecies[e][nsub];	
					if(pspecies_vase != NULL)
					    fprintf(fpc_vase,"%f ",pspecies_vase->C);
				}
			}
		}
		// SW 04/06/2019 print PHYF PHYS PHYR  // SW 30/07/2019 comment print annex variables

		//for(e = 0; e < ACTBACT; e++)
		//{
			//for(nsub = 0; nsub < Simul->counter_bio->nsubannex_var[e]; nsub++)
			//{
				//pspecies_water = Simul->psimul_bio[0][ns]->section->compartments[WATER][0]->pannex_var[e][nsub];
				//if(pspecies_water != NULL)
					//fprintf(fpc_water,"%f ",pspecies_water->C);
				//if(Simul->psimul_bio[0][ns]->section->compartments[VASE][0] != NULL)
				//{
					//pspecies_vase = Simul->psimul_bio[0][ns]->section->compartments[VASE][0]->pannex_var[e][nsub];	
					//if(pspecies_vase != NULL)
					    //fprintf(fpc_vase,"%f ",pspecies_vase->C);
				//}				
			//}
				
		//}
		
		fprintf(fpc_water,"\n");
                // SW 14/01/2022
               if(Simul->psimul_bio[0][0]->section->nsublayers[VASE] > 0)
               {

		    fprintf(fpc_vase,"\n");	
		    fprintf(fpc_volume_vase,"\n");
               }			
	}
	fclose(fpc_water);
        // SW 14/01/2022
        if(Simul->psimul_bio[0][0]->section->nsublayers[VASE] > 0)
        {

	    fclose(fpc_vase);
	    fclose(fpc_volume_vase);
        }	
}

/*calculation of different processus for check time
  argument length of check type needed
  FILE *fp : log file
*/
#ifdef CHR_CHECK
void PROSE_calculate_check_time( int ngrid, FILE *fp)
{
	int i,e;
	double sum_time;
	double time_tot = 0.;
	for(i = 0; i < CHECK_TYPE; i++)
	{
		sum_time = 0.;
		for(e = 0; e < ngrid; e++)
		{
			sum_time += Simul->psimul_bio[e]->check_time->check_time[i];
		}
		time_tot += sum_time;
		LP_printf(fp,"check_type = %d sum_time = %f\n",i,sum_time);
	}
	fprintf(fp,"time_tot_check = %f\n",time_tot);
}
#endif


void PROSE_print_parameters(s_simul ***psimul_bio, int nparticules,int param, double t, FILE *fp)
{
	double val,unit;
	int np;
switch(param) {
   case MAINT_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   fprintf(Simul->passim->pout_param[MAINT_DA],"%f\t",t);
	   unit = Simul->passim->units_param[MAINT_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->respiration->resp[MAINT20];
	   //Simul->passim->param[MAINT_DA][np] = val; // before resampling
	   fprintf(Simul->passim->pout_param[MAINT_DA],"%f\t",val/unit);		   
	   }
        fprintf(Simul->passim->pout_param[MAINT_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(MAINT_DA));     
	   break;}
   case ALPHA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   unit = Simul->passim->units_param[ALPHA_DA];
	   fprintf(Simul->passim->pout_param[ALPHA_DA],"%f\t",t);
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[A_RIVE];
	   //Simul->passim->param[ALPHA_DA][np] = val;
	   fprintf(Simul->passim->pout_param[ALPHA_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[ALPHA_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ALPHA_DA));	   
	   break;}
   case PMAX_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   unit = Simul->passim->units_param[PMAX_DA];
	   fprintf(Simul->passim->pout_param[PMAX_DA],"%f\t",t);
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[PMAX20];
	   //Simul->passim->param[PMAX_DA][np] = val;
	   fprintf(Simul->passim->pout_param[PMAX_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[PMAX_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(PMAX_DA));	   
	   break;}
   case ETA_CHLA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   unit = Simul->passim->units_param[ETA_CHLA_DA];
	   fprintf(Simul->passim->pout_param[ETA_CHLA_DA],"%f\t",t);
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA_CHLA];
	   //Simul->passim->param[ETA_CHLA_DA][np] = val;
	   fprintf(Simul->passim->pout_param[ETA_CHLA_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[ETA_CHLA_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ETA_CHLA_DA));	   
	   break;}
   case C_CHLA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   unit = Simul->passim->units_param[C_CHLA_DA];
	   fprintf(Simul->passim->pout_param[C_CHLA_DA],"%f\t",t);
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->nut_C[NCOMP_RIVE]; //C/chla
	   //Simul->passim->param[C_CHLA_DA][np] = val;
	   fprintf(Simul->passim->pout_param[C_CHLA_DA],"%f\t",1/(val*unit)); // C/chla		   
	   }
       fprintf(Simul->passim->pout_param[C_CHLA_DA],"\n");	
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ETA_CHLA_DA));   
	   break;}
   case ETA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   fprintf(Simul->passim->pout_param[ETA_DA],"%f\t",t);
	   unit = Simul->passim->units_param[ETA_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA]; //eta_water
	   //Simul->passim->param[ETA_DA][np] = val;
	   fprintf(Simul->passim->pout_param[ETA_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[ETA_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ETA_DA)); 	   
	   break;}
   case TOPT_PHY_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   fprintf(Simul->passim->pout_param[TOPT_PHY_DA],"%f\t",t);
	   unit = Simul->passim->units_param[TOPT_PHY_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->paraml[TOPT]; // Topt
	   //Simul->passim->param[TOPT_PHY_DA][np] = val;
	   fprintf(Simul->passim->pout_param[TOPT_PHY_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[TOPT_PHY_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(TOPT_PHY_DA));	   
	   break;}
   case MU_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   fprintf(Simul->passim->pout_param[MU_BACT_DA],"%f\t",t);
	   unit = Simul->passim->units_param[MU_BACT_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[MUMAX20]; // mumax
	   //Simul->passim->param[MU_BACT_DA][np] = val;
	   fprintf(Simul->passim->pout_param[MU_BACT_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[MU_BACT_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(MU_BACT_DA));	   
	   break;}
   case Y_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   fprintf(Simul->passim->pout_param[Y_BACT_DA],"%f\t",t);
	   unit = Simul->passim->units_param[Y_BACT_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[YIELD];
	   //Simul->passim->param[Y_BACT_DA][np] = val;
	   fprintf(Simul->passim->pout_param[Y_BACT_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[Y_BACT_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(Y_BACT_DA));	   
	   break;}
   case MORT_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   fprintf(Simul->passim->pout_param[MORT_BACT_DA],"%f\t",t);
	   unit = Simul->passim->units_param[MORT_BACT_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->mortality->mort[MORT20]; // mort
	   //Simul->passim->param[MORT_BACT_DA][np] = val;
	   fprintf(Simul->passim->pout_param[MORT_BACT_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[MORT_BACT_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(MORT_BACT_DA));	   
	   break;}
   case TOPT_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   fprintf(Simul->passim->pout_param[TOPT_BACT_DA],"%f\t",t);
	   unit = Simul->passim->units_param[TOPT_BACT_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->paraml[TOPT]; // Topt
	   //Simul->passim->param[TOPT_BACT_DA][np] = val;
	   fprintf(Simul->passim->pout_param[TOPT_BACT_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[TOPT_BACT_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(TOPT_BACT_DA));		   
	   break;}
   case KNAVIG_DA : {
        if(Simul->counter_bio->nsubspecies[O2] > 0)
        {
	   fprintf(Simul->passim->pout_param[KNAVIG_DA],"%f\t",t); 
	   unit = Simul->passim->units_param[KNAVIG_DA];	   
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[O2][0]->dissolved->gas->reaeration->rea[REA_NAVIG]; // Topt
	   //Simul->passim->param[TOPT_BACT_DA][np] = val;
	   fprintf(Simul->passim->pout_param[KNAVIG_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[KNAVIG_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(KNAVIG_DA));
	   break;}
   case B1_RIVER_DA : {  // b1_river
      if ( (Simul->counter_bio->nmacrospecies > EPS_TS) && (Simul->passim->param_range[B1_RIVER_DA][PARAM_UP] > EPS_TS) ) {
	   fprintf(Simul->passim->pout_param[B1_RIVER_DA],"%f\t",t); 
	   unit = Simul->passim->units_param[B1_RIVER_DA];	   
	   for(np = 0; np < nparticules; np++)
	   {
	      val = Simul->passim->p_macro_da[TOC][np]->degradOrgMat[MACMOD][B]->val; // b1_river
	      // val =  psimul_bio[np][0]->section->compartments[WATER][0]->pmacro[TOC]->degradOrgMat[MACMOD][B]->val; // b1_river 
	   //Simul->passim->param[B1_RIVER_DA][np] = val;
	   fprintf(Simul->passim->pout_param[B1_RIVER_DA],"%f\t",val/unit);		   
	   }
        fprintf(Simul->passim->pout_param[B1_RIVER_DA],"\n");
      }
	   break;}
   default : {
	   LP_error(Simul->poutputs,"Unknown parameter for assimilation. func:PROSE_print_parameter \n");
	   break;}
   }
}

void PROSE_print_extract_parameters(s_simul ***psimul_bio, int nparticules,int param, double t, FILE *fp)
{
	double val,unit;
	int np;
if(t <= Simul->chronos->t[END_CHR])
{
switch(param) {
   case MAINT_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   fprintf(Simul->passim->pout_param[MAINT_DA],"%f\t",t);
	   unit = Simul->passim->units_param[MAINT_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->respiration->resp[MAINT20];
	   Simul->passim->param[MAINT_DA][np] = val; // before resampling
	   fprintf(Simul->passim->pout_param[MAINT_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[MAINT_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(MAINT_DA));
	   break;}
   case ALPHA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   unit = Simul->passim->units_param[ALPHA_DA];
	   fprintf(Simul->passim->pout_param[ALPHA_DA],"%f\t",t);
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[A_RIVE];
	   Simul->passim->param[ALPHA_DA][np] = val;
	   fprintf(Simul->passim->pout_param[ALPHA_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[ALPHA_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ALPHA_DA));	   
	   break;}
   case PMAX_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   unit = Simul->passim->units_param[PMAX_DA];
	   fprintf(Simul->passim->pout_param[PMAX_DA],"%f\t",t);
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[PMAX20];
	   Simul->passim->param[PMAX_DA][np] = val;
	   fprintf(Simul->passim->pout_param[PMAX_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[PMAX_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(PMAX_DA));	   
	   break;}
   case ETA_CHLA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   unit = Simul->passim->units_param[ETA_CHLA_DA];
	   fprintf(Simul->passim->pout_param[ETA_CHLA_DA],"%f\t",t);
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA_CHLA];
	   Simul->passim->param[ETA_CHLA_DA][np] = val;
	   fprintf(Simul->passim->pout_param[ETA_CHLA_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[ETA_CHLA_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ETA_CHLA_DA));	   
	   break;}
   case C_CHLA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   unit = Simul->passim->units_param[C_CHLA_DA];
	   fprintf(Simul->passim->pout_param[C_CHLA_DA],"%f\t",t);
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->nut_C[NCOMP_RIVE]; //C/chla
	   Simul->passim->param[C_CHLA_DA][np] = val;
	   fprintf(Simul->passim->pout_param[C_CHLA_DA],"%f\t",1/(val*unit)); // C/chla		   
	   }
       fprintf(Simul->passim->pout_param[C_CHLA_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(C_CHLA_DA));	   
	   break;}
   case ETA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   fprintf(Simul->passim->pout_param[ETA_DA],"%f\t",t);
	   unit = Simul->passim->units_param[ETA_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA]; //eta_water
	   Simul->passim->param[ETA_DA][np] = val;
	   fprintf(Simul->passim->pout_param[ETA_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[ETA_DA],"\n");	
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ETA_DA));   
	   break;}
   case TOPT_PHY_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   fprintf(Simul->passim->pout_param[TOPT_PHY_DA],"%f\t",t);
	   unit = Simul->passim->units_param[TOPT_PHY_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->paraml[TOPT]; // Topt
	   Simul->passim->param[TOPT_PHY_DA][np] = val;
	   fprintf(Simul->passim->pout_param[TOPT_PHY_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[TOPT_PHY_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(TOPT_PHY_DA));	   
	   break;}
   case MU_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   fprintf(Simul->passim->pout_param[MU_BACT_DA],"%f\t",t);
	   unit = Simul->passim->units_param[MU_BACT_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[MUMAX20]; // mumax
	   Simul->passim->param[MU_BACT_DA][np] = val;
	   fprintf(Simul->passim->pout_param[MU_BACT_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[MU_BACT_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(MU_BACT_DA));	   
	   break;}
   case Y_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   fprintf(Simul->passim->pout_param[Y_BACT_DA],"%f\t",t);
	   unit = Simul->passim->units_param[Y_BACT_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[YIELD];
	   Simul->passim->param[Y_BACT_DA][np] = val;
	   fprintf(Simul->passim->pout_param[Y_BACT_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[Y_BACT_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(Y_BACT_DA));	   
	   break;}
   case MORT_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   fprintf(Simul->passim->pout_param[MORT_BACT_DA],"%f\t",t);
	   unit = Simul->passim->units_param[MORT_BACT_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->mortality->mort[MORT20]; // mort
	   Simul->passim->param[MORT_BACT_DA][np] = val;
	   fprintf(Simul->passim->pout_param[MORT_BACT_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[MORT_BACT_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(MORT_BACT_DA));	   
	   break;}
   case TOPT_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   fprintf(Simul->passim->pout_param[TOPT_BACT_DA],"%f\t",t);
	   unit = Simul->passim->units_param[TOPT_BACT_DA];
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->paraml[TOPT]; // Topt
	   Simul->passim->param[TOPT_BACT_DA][np] = val;
	   fprintf(Simul->passim->pout_param[TOPT_BACT_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[TOPT_BACT_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(TOPT_BACT_DA));		   
	   break;}
   case KNAVIG_DA : {
        if(Simul->counter_bio->nsubspecies[O2] > 0)
        {
	   fprintf(Simul->passim->pout_param[KNAVIG_DA],"%f\t",t); 
	   unit = Simul->passim->units_param[KNAVIG_DA];	   
	   for(np = 0; np < nparticules; np++)
	   {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[O2][0]->dissolved->gas->reaeration->rea[REA_NAVIG]; // Topt
	   Simul->passim->param[KNAVIG_DA][np] = val;
	   fprintf(Simul->passim->pout_param[KNAVIG_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[KNAVIG_DA],"\n");
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species O2 is not defined\n", PROSE_name_param(KNAVIG_DA));
	   break;}
   case B1_RIVER_DA : {  // b1_river
     if ( (Simul->counter_bio->nmacrospecies > EPS_TS) && (Simul->passim->param_range[B1_RIVER_DA][PARAM_UP] > EPS_TS) ) {
	   fprintf(Simul->passim->pout_param[B1_RIVER_DA],"%f\t",t); 
	   unit = Simul->passim->units_param[B1_RIVER_DA];	   
	   for(np = 0; np < nparticules; np++)
	   {
	     val = Simul->passim->p_macro_da[TOC][np]->degradOrgMat[MACMOD][B]->val; // b1_river
	      //val = psimul_bio[np][0]->section->compartments[WATER][0]->pmacro[TOC]->degradOrgMat[MACMOD][B]->val;
	     //LP_printf(Simul->poutputs,"b1 val_2 = %3.4f  \n", val);
	   
	   Simul->passim->param[B1_RIVER_DA][np] = val;
	   fprintf(Simul->passim->pout_param[B1_RIVER_DA],"%f\t",val/unit);		   
	   }
       fprintf(Simul->passim->pout_param[B1_RIVER_DA],"\n");
      }
	   break;}
   default : {
	   LP_error(Simul->poutputs,"Unknown parameter for assimilation. Func: PROSE_print_extract_parameters\n");
	   break;}
   }
}
}

s_total_mb * Prose_copy_total_mb(s_total_mb *mb2)
{
	int i,j,e; // one mb
	s_total_mb *mb1;
	//for(np = 0; np < nparticules; np++)
	//{
		mb1 = init_mass_balance(Simul->psimul_bio[0][0],Simul->chronos->t[BEGINNING], Simul->chronos->t[END], Simul->chronos->dt);
		mb1->t0 = mb2->t0;
		//mb1->t0 = mb2->t[BEGINNING];
		for(i = 0; i < NEXTREMA; i++)
			mb1->t[i] = mb2->t[i];
		mb1->time_unit = mb2->time_unit;
		mb1->ndt = mb2->ndt;
		for(e = 0; e < NSPECIES; e++)
		{
			for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++) 
			{
				mb1->calc_mb_species[e][j] = mb2->calc_mb_species[e][j];
				mb1->unit_mb_species[e][j] = mb2->unit_mb_species[e][j];
			}	
		}
		for(e = 0; e < NANNEX_VAR; e++)
		{
			for (j = 0; j < Simul->counter_bio->nsubannex_var[e]; j++) 
			{
				mb1->calc_mb_annex_var[e][j] = mb2->calc_mb_annex_var[e][j];
				mb1->unit_mb_annex_var[e][j] = mb2->unit_mb_annex_var[e][j];
			}
				
		}
		for(e = 0; e < NDISS; e++)
		{
			for (j = 0; j < Simul->counter_bio->nsubspecies[e+NPART]; j++) 
			{
				mb1->calc_mb_adsorbed_species[e][j] = mb2->calc_mb_adsorbed_species[e][j];
				mb1->unit_mb_adsorbed_species[e][j] = mb2->unit_mb_adsorbed_species[e][j];
			}
		}
		return mb1;
	//}
}


void PROSE_print_resampling_size(double t, int size, FILE *fpsize)
{
	fprintf(fpsize,"%f\t%d\n", t/86400.,size);
	fflush(fpsize);
}


// void PROSE_print_outputs_long(s_simul **psimul_bio,double t,s_output_hyd ***p_outputs,s_inout_set_io *pinout,s_chyd *pchyd,s_chronos_CHR *chronos,FILE *fp)
// {
  // double t0;
  // s_output_hyd *result;
  // int i;

  // if (pinout->calc[LONG_PROFILE] == YES_TS) {

    // for (i = 0; i < pchyd->counter->nlp; i++) {
		//LP_printf(fp,"nlp = %d\n",pchyd->counter->nlp);
      // result = p_outputs[LONG_PROFILE][i];
	  // result->pout->t_out[CUR_IO] = t; // SW 30/01/2018 CUR_IO, n'est pas fait
      // t0 = result->pout->t_out[CUR_IO];
      
      // if ((t >= result->pout->t_out[BEGINNING]) &&(t <= t0) && 
	  // (t + chronos->dt > t0) && (t <= result->pout->t_out[END])) {
	
	// PROSE_long_profile(psimul_bio,t,result,pinout,pchyd,fp);
	// result->pout->t_out[CUR_IO] += result->pout->deltat;
      // }
    // }
  // }
// }


void PROSE_long_profile(s_simul **psimul_bio,double t,s_output_hyd *result,int layer, s_inout_set_io *pinout,s_chyd *pchyd,FILE *fp)
{
  s_lp_pk_hyd *ppk;
  char name[MAXCHAR_PROSE];
  //int pk0_entier,pkf_entier,ligne;
  int k,pkf_atteint,nb;
  //int premier_passage;
  int lx,lx0;
  //double t_entier;
  double pk,pk0,pkf,pklim;
  s_element_hyd *pele;
  s_reach_hyd *preach;
  s_file_io *fic_water;
  //s_file_io *fic_vase;
  int i;
 
if(psimul_bio[0]->section->compartments[layer] != NULL) // SW 20/05/2019
{	
  for (i = 0; i < result->npk; i++) {
  nb = 0; // SW 25/01/2018 il faut intialiser, sinon la ligne 573 nb++
    //premier_passage = NO_TS;
    ppk = result->lp_pk[i];
    //LP_printf(fp,"reach_nb = %d\n",ppk->reach_nb[0]);
    if (ppk->reach_nb[0] == -1) {
      HYD_calculate_output_domain(ppk,LONG_PROFILE,pchyd,fp);
    }

    /* Initialisations */
    pk0 = ppk->pk_up;
    pkf = ppk->pk_down;
    preach = NULL;
    
    if (ppk->reach_nb[0] >= 0) {
      preach = pchyd->p_reach[ppk->reach_nb[0]];
      pklim = preach->limits[DOWNSTREAM]->pk;
    }
    
    pele = ppk->pele;
    //    ligne = nb = 0;
    pkf_atteint = NO_TS;
    k = 1;
    pk = pele->center->pk;
    
    /* Ouverture des fichiers :
     * Plusieurs cas :
     * - Soit un fichier unique qui contient tout 
     * (resultat->fic->adresse != NULL) ;
     * - Soit triés en plusieurs fichiers ...
     * - Fichiers au formats GNUPLOT et fichiers au format OPENDX 
     */
    //t_entier = (int)(t*result->time_unit);
    //pk0_entier = (int)(pk0/1000.0);
    //pkf_entier = (int)(pkf/1000.0);
    
    /* Fichiers au format GNUPLOT, triés par PK et par date */
    /* Ouverture à chaque nouvelle date.
     * Contiendra, à chaque date et pour chaque PK, les informations
     * de positons curviligne (PK) et les valeurs moyennes des variables
     * demandées. */
    //ppk->fic = (s_file_io *)calloc(1,sizeof(s_file_io));
	fic_water = (s_file_io *)calloc(1,sizeof(s_file_io));
	//fic_vase = (s_file_io *)calloc(1,sizeof(s_file_io));
    
	
	sprintf(name,"%s/longitudinal_profiles/lp_%s_t%4.2f_pk%4.2f_pk%4.2f_%s_%d.txt",pinout->name_outputs,name_layer(layer),t*result->pout->time_unit,pk0/1000.0,pkf/1000.,ppk->river,nb);
    if ((fic_water->address = fopen(name,"w")) == 0) 
      LP_error(fp,"file %s\n",name);	
    fic_water->name = strdup(name);	
	
    fprintf(fic_water->address,"#File generated by %s%4.2f\n",CODE_NAME,NVERSION_PROSE);//NF 11/10/2020
    fprintf(fic_water->address,"# Longitudinal profile (mean transversal values) between PK %4.2f and %4.2f at time %f in %s\n",pk0/1000.0,pkf/1000.0,t*result->pout->time_unit,name_layer(layer));
    fprintf(fic_water->address,"#PK POINT(XC YC) "); // SW 04/06/2021 add X Y coordinates
	
    //HYD_print_variables(ppk->fic->address,result);
    PROSE_print_variables(fic_water->address,result);
    
    /* SW 23/01/2024 to print sediment layer variables such as volume, height.*/
    if(layer == VASE)
        PROSE_print_sediment_variables(fic_water->address,result);

    fprintf(fic_water->address,"\n"); // SW 24/01/2024 add here this function

	
    lx = lx0 = 0;
    
    while ((preach != NULL)&&(pk < pkf)) {
      
      if ((preach->limits[DOWNSTREAM]->pk < pklim)||(pkf_atteint == YES_TS)) {
	if (nb == 0)
	  lx0 = lx;
	nb++;
	/* Les données de chaque branche sont réparties dans différents 
	 * fichiers. On ouvre donc un nouveau fichier.
	 */
	fclose(fic_water->address);
	sprintf(name,"%s/longitudinal_profiles/lp_%s_t%4.2f_pk%4.2f_pk%4.2f_%s_%d.txt",pinout->name_outputs,name_layer(layer),t*result->pout->time_unit,pk0/1000.,pkf/1000.,ppk->river,nb);
	if ((fic_water->address = fopen(name,"w")) == 0) 
	  LP_error(fp,"file %s\n",name);	
	fic_water->name = strdup(name);
      }	  
      
      while (pele != NULL) {
	/* Pour éviter les retours en arrière */
	if (pk <= preach->limits[DOWNSTREAM]->pk) {
	  /* Format GNUPLOT */
	  //fprintf(fic_water->address,"%f ",pele->center->pk / 1000.);
          fprintf(fic_water->address,"%f;",pele->center->pk / 1000.); // SW 11/10/2021  just for porject consacre add ;
          // SW 04/06/2021 add print of X Y coordinates
          //fprintf(fic_water->address,"%f %f ",pele->center->description->xc, pele->center->description->yc);
          fprintf(fic_water->address,"POINT(%f %f);",pele->center->description->xc, pele->center->description->yc); // SW 11/10/2021  just for porject consacre add ;

	  //HYD_print_average_result(pele,pele->center->pk,result,ppk->fic->address);
	  PROSE_print_average_result_bio(pele,psimul_bio, pchyd, layer, pele->center->pk, ppk->river,result,fic_water->address);
	  lx++;
	}
	pk += pele->length;
	if (pk > pkf) {
	  pkf_atteint = YES_TS;
	  break;
	}
	pele = pele->face[X_HYD][TWO]->element[TWO];
      }
      
      if (preach->limits[DOWNSTREAM]->pk >= pklim)
	pklim = preach->limits[DOWNSTREAM]->pk;
      if ((ppk->reach_nb[k] >= 0) &&(k < pchyd->counter->nreaches)){
	preach = pchyd->p_reach[ppk->reach_nb[k++]];
	pele = preach->p_ele[0];
	pk = preach->limits[UPSTREAM]->pk;
      }
      else
	preach = NULL;
    }    
    if (nb == 0)
      lx0 = lx;
    
    
    fclose(fic_water->address);
    //free(fic_water); // SW 28/01/2021 free memory
    //fic_water = NULL;
    //if (result->graphics == GNUPLOT)
      //HYD_print_GNUPLOT_lp(result,pinout,ppk,nb,t,fp);
  }
}
}


// void PROSE_print_GNUPLOT_lp(s_output_hyd *result,s_inout_set_io *pinout,s_lp_pk_hyd *ppk,int nb,double t,FILE *fp)
// {
  // int k,j,nvar = 1;  
  // s_file_io *gnu_fic;
  
  // for (k = 0; k < NSPECIES_IO; k++) { // NSPECIES_IO == NSPECIES
    
    // if (result->pout->biovar[k] == YES_TS) {
      // nvar++;
      
      // if ((gnu_fic->address = fopen(ppk->gnu_fic[k]->name,"a")) == 0)
	// LP_error(fp,"Impossible to open time series output file %s\n",ppk->gnu_fic[k]->name);
      
      // fprintf(ppk->gnu_fic[k]->address,"set title \"Longitudinal profile PK %4.2f to %4.2f, t = %4.2f\"\n",
	      // ppk->pk_up / 1000.,ppk->pk_down / 1000.,t * result->pout->time_unit);
      // fprintf(ppk->gnu_fic[k]->address,"plot ");
      
      // for (j = 0; j <= nb; j++) {
	// fprintf(ppk->gnu_fic[k]->address,"\"%s/longitudinal_profiles/lp_t%4.2f_pk%4.2f_pk%4.2f_%s_%d.txt\" ",
		// pinout->name_outputs,t * result->pout->time_unit,ppk->pk_up / 1000.,ppk->pk_down / 1000.,ppk->river,j);
	// fprintf(ppk->gnu_fic[k]->address,"using 1:%d title \"%s\" w l",nvar,HYD_name_hyd_var(k));
	
	// if (j < nb) {
	  // fprintf(ppk->gnu_fic[k]->address,",\\");
	  // fprintf(ppk->gnu_fic[k]->address,"\n");
	// }
	// else
	  // fprintf(ppk->gnu_fic[k]->address,"\n");
      // }
      
      // fprintf(ppk->gnu_fic[k]->address,"pause -1\n\n");
      
      // fclose(ppk->gnu_fic[k]->address);
    // }
  // }
// }


void PROSE_create_output_extents_bio(int npk,s_lp_pk_hyd *ppk)
{
  int i;

  //ppk = rewind_lp_pk(ppk);
  ppk = HYD_browse_lp_pk(ppk,BEGINNING);

  Simul->lp_pk = (s_lp_pk_hyd **)calloc(npk,sizeof(s_lp_pk_hyd *));
  for (i = 0; i < npk; i++) {
    Simul->lp_pk[i] = ppk;
    ppk = ppk->next;
  }
}

void PROSE_print_tube_mesh_hyd(double t, s_output_tube_type  *result, s_chyd *pchyd, s_def_tub ****tube, FILE *fp)
{
  char dir_path[MAXCHAR_PROSE], file_path[MAXCHAR_PROSE], file_path0[MAXCHAR_PROSE];
  char cmd[MAXCHAR_PROSE];
  FILE *fp_tube_mesh, *fp_tube_hyd;
  
  int r, e, ntb, tbis, i;
  
  /* SW 02/04/2021 moved into input.y */
  //if(fabs(t - Simul->chronos->t[BEGINNING_CHR]) < EPS_TS)
  //{
   //  sprintf(cmd,"mkdir %s/tube", getenv("RESULT"));
     //system(cmd);
     //LP_printf(fp, "%s created\n");
  //}

  switch(result->output_type) {
      case TUBE_MESH:
           sprintf(file_path0,"%s/tube/coord_geotube_t%4.2f.txt", getenv("RESULT"), t*result->time_unit);
           if ((fp_tube_mesh = fopen(file_path0,"w")) == 0)
               LP_error(fp,"Impossible to open %s file\n", file_path0);
           //fprintf(fp_tube_mesh,"#File generated by %s%4.2f, tube coordinates at time %f\n",CODE_NAME,NVERSION_PROSE, t*result->time_unit);       
           fprintf(fp_tube_mesh,"ID_TUBE;CELL_CENTRE;CELL_VERTICES\n");
           break;
      case TUBE_HYD:
           sprintf(file_path,"%s/tube/hyd_geotube_t%4.2f.txt", getenv("RESULT"), t*result->time_unit);
           if((fp_tube_hyd = fopen(file_path,"w")) == 0)
               LP_error(fp,"Impossible to open %s file\n", file_path);
           //fprintf(fp_tube_hyd,"#File generated by %s%4.2f, tube hydraulics at time %f\n",CODE_NAME,NVERSION_PROSE, t*result->time_unit);
           fprintf(fp_tube_hyd,"TIME;ID_TUBE;H_WATER;U_WATER;ALTI_WATER\n");
           break;
      default: LP_error(fp,"Unknown tube output type \n");
  }

  /* print tubes in entire domain */
  for (r = 0; r < pchyd->counter->nreaches; r++)
  {
    for (e = 0; e < pchyd->p_reach[r]->nele; e++)
    {
      for (ntb = 0; ntb < pchyd->p_reach[r]->n_tube; ntb++)
	{
	  if (tube[r][e][ntb]->geocoords != NULL)
	    {
	      if (error_geocoord_TUB(tube[r][e][ntb]->geocoords) == 1)
		{
		  LP_warning(fp, "Error with tube[%d][%d][%d] of abs. id: %d, some coords are similar.\n",r,e,t,tube[r][e][ntb]->id_tube_abs);
		  /*printf("Destruction of the coresponding tubes' series.\n");
		  for (tbis=0; tbis<pchyd->p_reach[r]->n_tube; tbis++)
		  {
		  for (tbis=0; tbis<pchyd->p_reach[r]->n_tube; tbis++)
		  {
		  for (i=0; i<N_VERTEX; i++)
		  free(tube[r][e][tbis]->geocoords->vertices[i]);
		  free(tube[r][e][tbis]->geocoords->centre);
		  free(tube[r][e][tbis]->geocoords);
		  }
		  break;
		  }*/
		}
	      else
		{
                  if(result->output_type == TUBE_MESH)
		      fprintf(fp_tube_mesh,"%d;POINT(%lf %lf);POLYGON((%lf %lf, %lf %lf, %lf %lf, %lf %lf, %lf %lf))\n", tube[r][e][ntb]->id_tube_abs, tube[r][e][ntb]->geocoords->centre[0], tube[r][e][ntb]->geocoords->centre[1], tube[r][e][ntb]->geocoords->vertices[LEFT_UP][0], tube[r][e][ntb]->geocoords->vertices[LEFT_UP][1], 
                                           tube[r][e][ntb]->geocoords->vertices[RIGHT_UP][0], tube[r][e][ntb]->geocoords->vertices[RIGHT_UP][1], tube[r][e][ntb]->geocoords->vertices[RIGHT_DOWN][0], tube[r][e][ntb]->geocoords->vertices[RIGHT_DOWN][1], tube[r][e][ntb]->geocoords->vertices[LEFT_DOWN][0], tube[r][e][ntb]->geocoords->vertices[LEFT_DOWN][1], tube[r][e][ntb]->geocoords->vertices[LEFT_UP][0], tube[r][e][ntb]->geocoords->vertices[LEFT_UP][1]);
		  else if(result->output_type == TUBE_HYD)
                      fprintf(fp_tube_hyd,"%f;%d;%lf;%lf;%lf\n", t*result->time_unit, tube[r][e][ntb]->id_tube_abs, (tube[r][e][ntb]->param_amont->H + tube[r][e][ntb]->param_aval->H)/2., (tube[r][e][ntb]->param_amont->U + tube[r][e][ntb]->param_aval->U)/2., pchyd->p_reach[r]->p_ele[e]->center->hydro->Z[T_HYD]);
		}
	     }
	  }
      }
    }

  if(result->output_type == TUBE_MESH)
      fclose(fp_tube_mesh);
  if(result->output_type == TUBE_HYD)
      fclose(fp_tube_hyd); 
}


