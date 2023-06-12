/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: calc_base_assimilation.c
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

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <time.h>
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
#include "HYD.h"
#include "TTC.h"
#include "RIVE.h"
#include "SEB.h"
#include "TUB.h"
#include "MB.h"
#include "LA.h"
#include "PROSE.h"
////#include "global_PROSE.h"
#include "ext_PROSE.h"

void Prose_init_rand_param(s_simul ***psimul_bio, s_carac_assim *passim, int nele, FILE *fp)
{
	int nparam, np;
	double var, varP, varG;
	double up, down;
	for(nparam = 0; nparam < NPARAMDA; nparam++)
	{
            // SW 24/01/2022 add check if parameter is assimilated
            if(passim->param_yesOrno[nparam] == YES_TS)
            {
	        up = passim->param_range[nparam][PARAM_UP];
	        down = passim->param_range[nparam][PARAM_DOWN];
                for(np = 0; np < passim->N_particules; np++)
	        {	
	            var = Prose_generate_random_param(down, up);
	    
	            passim->param[nparam][np] = var;

                    /* SW 13/01/2022 just for EnKF */
                    //if(passim->method == ENKF_PROSE)
                    //{
                    // /* SW 13/01/2022 transform parameter values into values in range (0,1), varP, for EnKF. X -> X1 */
                    // /* This value, varP, is treated as  a cumulative propability. */
   
                    //    varP = (var - down)/(up - down);
                    //  //LP_printf(fp,"nparam = %d np = %d down = %f up = %f varP = %.10f\n",nparam,np,down,up,varP);

                    // /* SW 13/01/2022 transform varP into gaussian values, varG, which follow the standard normal distribution N(0, 1). X1 -> X2 */
                    // /* These gaussain values varG are used in EnKF */

                    //   varG = gsl_cdf_ugaussian_Pinv(varP);
                    //   //LP_printf(fp,"nparam = %d np = %d down = %f up = %f varG = %.10f\n",nparam,np,down,up,varG);
                    //   passim->paramGaussian[nparam][np] = varG; 
                    //}
                    PROSE_assign_parameters_val(psimul_bio, nele,  np, nparam, var, fp);
                    //LP_printf(fp,"nparam = %d np = %d down = %f up = %f var = %.10f\n",nparam,np,down,up,var);
                }
            }			
	}
	
}

/*calculation of normalized weights and the effective sample size used in main thread
return is the effective sample size Neff*/
int Prose_calc_normalized_weights_and_sample_size(s_carac_assim *passim, FILE *fp)
{
   int i,nobs, nele;
   double sum = 0, sumcumul = 0,Neff;
   double sum_carre = 0;
   s_carac_obs_assim *pobs;
   
   for(i = 0; i < passim->N_particules; i++)
      sum += passim->omega[i];

   LP_printf(Simul->poutputs, "sum = %.15f \n ", sum);
   
   if(fabs(sum) < EPS_TS)
  {
          
    /*
          //MH 23/12/2022 regeneration of random parameter values with 1/Np weights
    
	  nele = Simul->pcarac_ttc[0]->count[NELE_TTC];
	  Prose_init_rand_param(Simul->psimul_bio, passim, nele, Simul->poutputs); 
	  Prose_reinitialization_weights(passim, Simul->poutputs);
	  sum = 1.0;
    	  LP_warning(Simul->poutputs,"all weights are zero !!! Reset omega and param values\n");
    */

	    // old method of caring for zero weights

	  //Prose_calc_weight_all_particules(passim, Simul->poutputs);

	  LP_printf(Simul->poutputs,"np\t omega\t omega_prev\n");
	  for(i = 0; i < passim->N_particules; i++)
	  {
		  LP_printf(Simul->poutputs,"%d\t%f\t%f\n",i,passim->omega[i],passim->omega_prev[i]);
	  }
	  LP_printf(Simul->poutputs,"Differences between obs and simul\n");
	  for(nobs = 0; nobs < passim->N_obs; nobs++)
	  {
		  pobs = passim->pobs[nobs];
		  LP_printf(Simul->poutputs,"nobs = %d answer_obs = %d var at t = %f\n",nobs,pobs->answer_obs,pobs->Obs);
		  for(i = 0;i < passim->N_particules; i++)
			 LP_printf(Simul->poutputs,"%f\t",pobs->difference[i]);
          LP_printf(Simul->poutputs,"\n");		 
	  }
	  LP_warning(Simul->poutputs,"all weights are zero !!! Reset omega as omega_prev\n");
	  //SW 04/05/2019 if all weights are zero, it means that a observation is not good.

	  for(i = 0; i < passim->N_particules; i++)
      {
		  passim->omega[i] = passim->omega_prev[i];
		  sum += passim->omega[i];
	  } 
	 
  }
   for(i = 0; i < passim->N_particules; i++)
   {
      passim->omega_normlized[i] = passim->omega[i]/sum;
      //LP_printf(fp,"Calc_norm_Neff func) np = %i omega_prev = %3.4f, omega = %3.4f, omega_norm = %3.4f \n ",i, passim->omega_prev[i], passim->omega[i], passim->omega_normlized[i]);

      sumcumul += passim->omega_normlized[i];
      passim->cdf[i] = sumcumul;
      sum_carre += passim->omega_normlized[i]*passim->omega_normlized[i];
   }
   Neff = ceil(1.0/sum_carre);
   passim->Neff = (int) Neff;
   return (int) Neff;
}

/*fill the array eliminated and ndupli of a particule used after construction CDF*/
void Prose_fill_re_sample_elim_dupli(s_carac_assim *passim, FILE *fp)
{
    double u, u_init; //random variable from U[0,1/N]
    int i = 0,j; //loop index
    int num = 0;
    u_init = Prose_generate_random_param(0, 1.0/passim->N_particules);
    for(j = 1; j < passim->N_particules+1; j++)
    {
        u = u_init + (j - 1.0)/passim->N_particules;

        while(passim->cdf[i] < u) // i is the same with rank processeur or particule identity
        {
		if(passim->ndupli[i] > 0)
			i++;
        else{
			passim->eliminated[i] = YES_TS;
            //passim->ndupli[i] = 0; // initialized to 0 in reinitialization of weights
            i++;			
		}		
        }
        passim->ndupli[i]++; 
		num++;
		if(num == passim->N_particules)
		{
			while((i + 1) < passim->N_particules)
			{
				passim->eliminated[i+1] = YES_TS;
				i++;
			}
			break;
		}
    }    
}

/*reinitialization of weights after resampling, used in each noeud*/
void Prose_reinitialization_weights(s_carac_assim *passim, FILE *fp)
{
    int i,num_threads,taille;
    //#ifdef OMP
    //num_threads = Simul->num_threads_par;
    //taille = PC_set_chunk_size_silent(Simul->poutputs,passim->N_particules,num_threads);
    //omp_set_num_threads(num_threads);
    //#pragma omp parallel for schedule(dynamic,taille) shared(passim) private(i)
   // #endif   
    for(i = 0; i < passim->N_particules; i++)
    {  
        passim->omega[i] = 1.0/passim->N_particules;
		passim->omega_normlized[i] = 1.0/passim->N_particules;
		passim->ndupli[i] = 0;
		passim->omega_prev[i] = passim->omega_normlized[i];
    }
    return;        
}

/*duplicate particule used in all noeuds*/
void Prose_resample_particules(s_carac_assim *passim, FILE *fp) 
{
    int i,j;
    int dest = 0,source = 0;

    for(i = 0; i < passim->N_particules; i++)
    {
        while(passim->ndupli[i] > 1)
        {
            source = i;
            for(j = dest; j < passim->N_particules; j++)
            {
                if(passim->eliminated[j] == YES_TS)
                {
                    dest = j;
                    //function to duplicate all concentrations, all parameters (11), all masse_balances values from i to j MPI
                    //check it is in the same noeud
                    //Prose_send_receive_all_messages(psimul_rank, rank, dest, source, fp);

                    Prose_duplications_particules(source, dest, fp);
                    passim->ndupli[i]--; //one duplication
					passim->eliminated[dest] = NO_TS;
                    //passim->eliminated[dest] = NO_TS; // SW 26/01/2019 redefined in perturb_param
					dest++;
                    break;                   
                }
            } 
        }
    }           
}

void Prose_duplications_particules(int source, int dest, FILE *fp)
{
    int ns,nsections,layer,sublayer,e,esub;
    s_species *pspecies_source, *pspecies_dest;
    s_annex_var *pspecies_source_var, *pspecies_dest_var;
    s_compartment *pcomp_source, *pcomp_dest;
    nsections = Simul->pcarac_ttc[source]->count[NELE_TTC];


    // MH 10/12/2021 duplication of MACROSPECIES parameters
     if ( (Simul->counter_bio->nmacrospecies > EPS_TS) && (Simul->passim->param_range[B1_RIVER_DA][PARAM_UP] > EPS_TS) ) {
       for (e = 0; e< NMACROSPECIES; e++)
        {
            switch(e) {
                case TOC: {
		  Simul->passim->p_macro_da[TOC][dest]->degradOrgMat[MACMOD][B]->val =  Simul->passim->p_macro_da[TOC][source]->degradOrgMat[MACMOD][B]->val; //b1_river
                }
		  // other params like b1_WWTP and b1_CSO can be added here as well

		  break;  
                }
		}
     }

    for(ns = 0; ns < nsections; ns++)
    {
        for(layer = 0; layer < NLAYERS; layer++)
        {
            for(sublayer = 0; sublayer < Simul->psimul_bio[source][ns]->section->nsublayers[layer]; sublayer++)
            {
                pcomp_source = Simul->psimul_bio[source][ns]->section->compartments[layer][sublayer];
                pcomp_dest = Simul->psimul_bio[dest][ns]->section->compartments[layer][sublayer];
                
                for(e = 0; e < NSPECIES; e++)
                {
                    for(esub = 0; esub < Simul->counter_bio->nsubspecies[e]; esub++)
                    {
                        //copy concentrations of all species
                        pspecies_source = pcomp_source->pspecies[e][esub];
                        pspecies_dest = pcomp_dest->pspecies[e][esub];
                        Prose_duplication_conc_one_specie(pspecies_source, pspecies_dest, fp);

                        //copy mass balance
						//for instance only mbs with the same species
						if(Simul->calc_mode[MB_BIO] == YES_TS)
						{							
						if(Simul->total_mb[source][0]->calc_mb_species[e][esub] == YES_RIVE) // SW 30/04/2019 copy only the species asked by user here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done
						{
                          Prose_duplication_mass_balance_one_specie(pspecies_source, pspecies_dest, fp);
                        
                        //copy total mass blance, sum of all girds
                        if(ns == 0) // one time calculation
                        {
                            Prose_duplication_total_mass_balance(layer, e, esub, source, dest, fp); 
                        }
						}
						}
                        //copy all parameters to assimilate
						// SW 22/11/2018 if BLOOM, we don't copy bacterial parameters in order to maintaining the particle diversity
                        //if(e == PHY || e == BACT)
                        //{
                           // Prose_duplication_param_one_specie(pspecies_source, pspecies_dest, fp);                                            
                        //}
						if((e == BACT) || (e == O2))
							Prose_duplication_param_one_specie(pspecies_source, pspecies_dest, fp);
						//else if((e == PHY) && (Simul->passim->state == BLOOM))
						else if(e == PHY) // SW 11/06/2019 no difference bloom and horsbloom
							Prose_duplication_param_one_specie(pspecies_source, pspecies_dest, fp);
						
                     }
                }
                // all annexe variables
                for(e = 0; e < NANNEX_VAR; e++)
                {
                    for(esub = 0; esub < Simul->counter_bio->nsubannex_var[e]; esub++) 
                    {
                        //copy concentrations of all species
                        pspecies_source_var = pcomp_source->pannex_var[e][esub];
                        pspecies_dest_var = pcomp_dest->pannex_var[e][esub];
                        Prose_duplication_conc_one_annnexe_specie(pspecies_source_var, pspecies_dest_var, fp);
                      
                        //copy all mass balances
						if(Simul->calc_mode[MB_BIO] == YES_TS)
						{
						if(Simul->total_mb[source][0]->calc_mb_annex_var[e][esub] == YES_RIVE) // SW 30/04/2019 copy only the species asked by user here YES_RIVE == 1, which is equivalent to YES_RIVE in librive where the intialization is done
						{
                        Prose_duplication_mass_balance_one_annexe_specie(pspecies_source_var, pspecies_dest_var, fp);

                        //copy total mass blance, sum of all girds
                        if(ns == 0) // one time calculation
                        {
                            Prose_duplication_total_mass_balance_annexe(layer, e, esub, source, dest, fp);
                        }
						}
						}
                    }
                }
 
            }
        }
    }       

}

/*function uesd to send all concentrations, parameters, mass blance from source noeud to dest noeud*/
/*
void Prose_send_receive_all_messages(s_simul_PROSE *psimul_rank,int rank,int dest,int source,FILE *fp)
{
    int noeud_dest,noeud_source;
    int ind_dest,ind_source;
    int conc_tag = 0,param_tag = 1,spmb_tag = 2,spmb_total_tag = 3;
    int ns,nsections,layer,sublayer,e,sube,ind_nmb;
    double conc[NC_MPI]; // define C,newC NC_MPI = 2
    double param[NPARAMDA]; // NPARAMDA = 11, 11 parameters for data assimilation
    double spmb[NVAR_MB]; // all mass balances
    double spmb_total[NVAR_MB]; // total mass balance
    s_species *pspecies;


    noeud_dest = dest/(psimul_rank->passim->nparticules - 1);
    noeud_source = source/(psimul_rank->passim->nparticules - 1);
    ind_dest = dest % psimul_rank->passim->nparticules;
    ind_source = source % psimul_rank->passim->nparticules; 
 
    if(noeud_dest == noeud_source) // two simulations located in the same noeud
    {
        if(rank == noeud_dest)
        {
        //copies de concentrations, paramÃ¨tres, bilan de masse sur un noeud
        nsections = psimul_rank->pcarac_ttc->count[NELE_TTC];
        for(ns = 0; ns < nsections; ns++) //all grids or all sections
        {
            //conc[ID_SECTION_MPI] = ns;
            for(layer = 0; layer < NLAYERS; layer++)
            {
                //conc[LAYER_MPI] = layer;
                for(sublayer = 0; sublayer < psimul_rank->psimul_bio[ind_source][ns]->section->nsublayers[layer]; sublayer++)
                {
                    //conc[SUB_LAYER_MPI] = sublayer;
                    for(e = 0; e < NSPECIES; e++)
                    {
                        //conc[SPECIES_MPI] = e;
                            for(esub = 0; esub < psimul_rank->counter_bio->nsubspecies[e]; esub++)
                            {
                               // ajout des copies conc, parameters, mass balance
                    }
                }
             }
          }
        }
    }
    else
    {
        nsections = psimul_rank->pcarac_ttc->count[NELE_TTC];
        for(ns = 0; ns < nsections; ns++) //all grids or all sections
        {
            //conc[ID_SECTION_MPI] = ns;
            for(layer = 0; layer < NLAYERS; layer++)
            {
                //conc[LAYER_MPI] = layer;
                for(sublayer = 0; sublayer < psimul_rank->psimul_bio[ind_source][ns]->section->nsublayers[layer]; sublayer++)
                {
                    //conc[SUB_LAYER_MPI] = sublayer;
                    for(e = 0; e < NSPECIES; e++)
                    {
                        //conc[SPECIES_MPI] = e;
                            for(esub = 0; esub < psimul_rank->counter_bio->nsubspecies[e]; esub++)
                            {
                                 //conc[SUB_SPECIES_MPI] = esub;
                                 if(rank == noeud_source) // find the noeud source
                                 {   
                                     pspecies = psimul_rank->psimul_bio[ind_source][ns]->section->compartments[layer][sublayer]->pspecies[e][esub];
                                     conc[C_MPI] = pspecies->C;
                                     conc[NEWC_MPI] = pspecies->newC;
                                     MPI_Send(&conc, NC_MPI, MPI_DOUBLE, noeud_dest, conc_tag, MPI_COMM_WORLD); // send concentrations from noeud_source to noeud_dest
                                     if(e == PHY || e == BACT)
                                     {
                                         Prose_extraction_all_param_da(pspecies, e, esub, param, fp);
                                         MPI_Send(&param, NPARAMDA, MPI_DOUBLE, noeud_dest, param_tag, MPI_COMM_WORLD); // send parameters from noeud_source to noeud_dest
                                     }
                                     Prose_extraction_all_mass_balance(pspecies, spmb, fp);
                                     MPI_Send(&spmb, NVAR_MB, MPI_DOUBLE, noeud_dest, spmb_tag, MPI_COMM_WORLD); // send all deltamb
 
                                 }
                                 if(rank == noeud_dest) // find the noeud dest
                                 {
                                    pspecies = psimul_rank->psimul_bio[ind_dest][ns]->section->compartments[layer][sublayer]->pspecies[e][esub];
                                    MPI_Recv(&conc,NC_MPI,MPI_DOUBLE,noeud_source,conc_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);// receive concentrations 
                                    pspecies->C = conc[C_MPI];
                                    pspecies->newC = conc[NEWC_MPI];
                                    if(e == PHY || e == BACT)
                                    {
                                        MPI_Recv(&param,NPARAMDA,MPI_DOUBLE,noeud_source,param_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);// receive parameters                      
                                        Prose_assign_all_param_da(pspecies, e, esub, param, fp);   
                                    }
                                    MPI_Recv(&spmb,NVAR_MB,MPI_DOUBLE,noeud_source,spmb_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // receive all deltamb                   
                                    Prose_assign_all_mass_balance(pspecies, spmb, fp); // assign all deltamb
                                 }  
                                     // for total mass balanc, we can modify this by printing the element mass balance in output file 
                                     // only one total mb in prose-p, ensemble domaine
                                     if(ns == 0) // one time calculation
                                     {
                                     for(ind_nmb = 0; ind_nmb < psimul_rank->psimul_bio[ind_source][ns]->counter_bio->nmb; ind_nmb++)
                                     {
                                         if(rank == noeud_source) 
                                         {

                                             Prose_extraction_all_total_mass_balance(psimul_rank->psimul_bio[ind_source], ind_nmb, layer, e, esub, spmb_total, fp);
                                             MPI_Send(&spmb_total, NVAR_MB, MPI_DOUBLE, noeud_dest, spmb_total_tag, MPI_COMM_WORLD); 
                                         }
                                         if(rank == noeud_dest) 
                                         {
                                             MPI_Recv(&spmb_total,NVAR_MB,MPI_DOUBLE,noeud_source,spmb_total_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);                                             
                                             Prose_assign_all_total_mass_balance(psimul_rank->psimul_bio[ind_dest], ind_nmb,layer, e, esub, spmb_total, fp);
                                         }
                                     }
                                     } // fin if
                           } // fin esub
                   } // fin e
               } // fin sublayer
           }// fin layer
            
        } // fin ns
    } // fin else
    
}


void Prose_extraction_all_param_da(s_species *pspecies, int species, int esub, double *param, FILE *fp)
{
    switch(species) {
       case PHY: { // phytoplanktonic parameters 
            param[MAINT_DA] = pspecies->particulate->living->respiration->resp[MAINT20];
            param[ALPHA_DA] = pspecies->particulate->living->photosynthesis->phot[A_RIVE];
            param[PMAX_DA] = pspecies->particulate->living->photosynthesis->phot[PMAX20];
            param[ETA_CHLA_DA] = pspecies->particulate->living->photosynthesis->phot[ETA_CHLA];
            param[ETA_DA] = pspecies->particulate->living->photosynthesis->phot[ETA];
            param[C_CHLA_DA] = pspecies->nut_C[NCOMP_RIVE]; // C/chla
            param[TOPT_PHY_DA] = pspecies->particulate->living->paraml[TOPT];
            break;
       }
       case BACT: { // bacterial parameters
            param[MU_BACT_DA] = pspecies->particulate->living->growth->growth[MUMAX20];
            param[Y_BACT_DA] = pspecies->mb->yields;
            param[MORT_BACT_DA] = pspecies->particulate->living->mortality->mort[MORT20];
            param[TOPT_BACT_DA] = pspecies->particulate->living->paraml[TOPT];
            break;
       }
    }   
}

void Prose_assign_all_param_da(s_species *pspecies, int species, int esub, double *param, FILE *fp)
{
    switch(species) {
       case PHY: { // phytoplanktonic parameters 
            pspecies->particulate->living->respiration->resp[MAINT20] = param[MAINT_DA];
            pspecies->particulate->living->photosynthesis->phot[A_RIVE] = param[ALPHA_DA];
            pspecies->particulate->living->photosynthesis->phot[PMAX20] = param[PMAX_DA];
            pspecies->particulate->living->photosynthesis->phot[ETA_CHLA] = param[ETA_CHLA_DA];
            pspecies->particulate->living->photosynthesis->phot[ETA] = param[ETA_DA];
            pspecies->nut_C[NCOMP_RIVE] = param[C_CHLA_DA]; // C/chla
            pspecies->particulate->living->paraml[TOPT] = param[TOPT_PHY_DA];
            break;
       }
       case BACT: { // bacterial parameters
            pspecies->particulate->living->growth->growth[MUMAX20] = param[MU_BACT_DA];
            pspecies->mb->yield = param[Y_BACT_DA];
            pspecies->particulate->living->mortality->mort[MORT20] = param[MORT_BACT_DA];
            pspecies->particulate->living->paraml[TOPT] = param[TOPT_BACT_DA];
            break;
       }
    }   
}*/

void Prose_duplication_conc_one_specie(s_species *pspecies_source, s_species *pspecies_dest, FILE *fp)
{
    pspecies_dest->C = pspecies_source->C;
    pspecies_dest->newC = pspecies_source->newC;
}

void Prose_duplication_conc_one_annnexe_specie(s_annex_var *pspecies_source, s_annex_var *pspecies_dest, FILE *fp)
{
    pspecies_dest->C = pspecies_source->C;
    pspecies_dest->newC = pspecies_source->newC;
}

void Prose_duplication_param_one_specie(s_species *pspecies_source, s_species *pspecies_dest, FILE *fp)
{
    int var;
    var = pspecies_source->var;

    switch(var) {
       case PHY: {    //phytoplanktonic parameters 
            pspecies_dest->particulate->living->respiration->resp[MAINT20] = pspecies_source->particulate->living->respiration->resp[MAINT20]; //maint
            pspecies_dest->particulate->living->photosynthesis->phot[A_RIVE] = pspecies_source->particulate->living->photosynthesis->phot[A_RIVE]; //alpha
            pspecies_dest->particulate->living->photosynthesis->phot[PMAX20] = pspecies_source->particulate->living->photosynthesis->phot[PMAX20]; //pmax
            pspecies_dest->particulate->living->photosynthesis->phot[ETA_CHLA] = pspecies_source->particulate->living->photosynthesis->phot[ETA_CHLA]; //eta_chla
            pspecies_dest->particulate->living->photosynthesis->phot[ETA] = pspecies_source->particulate->living->photosynthesis->phot[ETA]; //eta_water
            pspecies_dest->nut_C[NCOMP_RIVE] = pspecies_source->nut_C[NCOMP_RIVE]; //C/chla
            pspecies_dest->particulate->living->paraml[TOPT] = pspecies_source->particulate->living->paraml[TOPT]; // Topt
            break;
       }
       case BACT: { // bacterial parameters
            pspecies_dest->particulate->living->growth->growth[MUMAX20] = pspecies_source->particulate->living->growth->growth[MUMAX20]; // mumax
            //pspecies_dest->mb->yield = pspecies_source->mb->yield; // Y SW 15/11/2018 error
            pspecies_dest->particulate->living->growth->growth[YIELD] = pspecies_source->particulate->living->growth->growth[YIELD]; //SW 15/11/2018			
            pspecies_dest->particulate->living->mortality->mort[MORT20] = pspecies_source->particulate->living->mortality->mort[MORT20]; // mort
            pspecies_dest->particulate->living->paraml[TOPT] = pspecies_source->particulate->living->paraml[TOPT]; // Topt
            break;
			}
       case O2: {  // Knavig
		    pspecies_dest->dissolved->gas->reaeration->rea[REA_NAVIG] = pspecies_source->dissolved->gas->reaeration->rea[REA_NAVIG];
			
	   }
            break;
       }
    }

void Prose_duplication_mass_balance_one_specie(s_species *pspecies_source, s_species *pspecies_dest, FILE *fp)
{
    int p;

    for(p = HYDROLYSIS; p < NVAR_MB; p++)  
    {
        pspecies_dest->mb->deltamb[p] = pspecies_source->mb->deltamb[p];
        pspecies_dest->mb->deltamb[p] = pspecies_source->mb->deltamb[p];
    }
}

void Prose_duplication_mass_balance_one_annexe_specie(s_annex_var *pspecies_source, s_annex_var *pspecies_dest, FILE *fp)
{
    int p;

    for(p = HYDROLYSIS; p < NVAR_MB; p++)  
    {
        pspecies_dest->mb->deltamb[p] = pspecies_source->mb->deltamb[p];
        pspecies_dest->mb->deltamb[p] = pspecies_source->mb->deltamb[p];
    }
}

//function to copy total mass_balance, Simul->total[N][NMB], N particules
void Prose_duplication_total_mass_balance(int layer, int e,int esub, int source, int dest, FILE *fp) 
{
    int p, ind_nmb;

    for(ind_nmb = 0; ind_nmb < Simul->counter_bio->nmb; ind_nmb++)
    {
        for(p = HYDROLYSIS; p < NVAR_MB; p++)  
        {
            Simul->total_mb[dest][ind_nmb]->mbspecies[layer][p][e][esub] = Simul->total_mb[source][ind_nmb]->mbspecies[layer][p][e][esub]; 
        } 
    } 
}

void Prose_duplication_total_mass_balance_annexe(int layer, int e,int esub, int source, int dest, FILE *fp) 
{
    int p, ind_nmb;

    for(ind_nmb = 0; ind_nmb < Simul->counter_bio->nmb; ind_nmb++)
    {
        for(p = HYDROLYSIS; p < NVAR_MB; p++)  
        {
            Simul->total_mb[dest][ind_nmb]->mbannex[layer][p][e][esub] = Simul->total_mb[source][ind_nmb]->mbannex[layer][p][e][esub]; 
        } 
    } 
}
/*
void Prose_extraction_all_mass_balance(s_species *pspecies, double *spmb, FILE *fp) 
{
    int p;

    for(p = HYDROLYSIS; p < NVAR_MB; p++)  
    {
        spmb[p] = pspecies->mb->deltamb[p];  
    } 
}

void Prose_assign_all_mass_balance(s_species *pspecies, double *spmb, FILE *fp) 
{
    int p;

    for(p = HYDROLYSIS; p < NVAR_MB; p++)  
    {
        pspecies->mb->deltamb[p] = spmb[p];  
    } 
}

void Prose_extraction_all_total_mass_balance(s_simul *psimulbio, int ind_nmb,int layer, int e,int esub, double *spmb_total, FILE *fp) 
{
    int p;

    for(p = HYDROLYSIS; p < NVAR_MB; p++)  
    {
        spmb_total[p] = psimulbio->total_mb[ind_nmb]->mbspecies[layer][p][e][esub];  
    } 
}

void Prose_assign_all_total_mass_balance(s_simul *psimulbio, int ind_nmb,int layer, int e,int esub, double *spmb_total, FILE *fp) 
{
    int p;

    for(p = HYDROLYSIS; p < NVAR_MB; p++)  
    {
        psimulbio->total_mb[ind_nmb]->mbspecies[layer][p][e][esub] = spmb_total[p];  
    } 
}*/

/* This function is used to calculate covariance matrix used in likelihood function*/
s_matrix_la * Prose_cov_matrix(s_carac_assim *passim, FILE *fp)
{
    int nparticules, num_t_obs;
    s_matrix_la *cov_matrix;
    
    int col = 0, row;
    double cov,sigma;
    //double *diff1, *diff2;
	//double Obs;
    nparticules = passim->N_particules;
    num_t_obs = passim->num_t_obs;
    cov_matrix = LA_matrix_calloc(num_t_obs,num_t_obs);

    //MH 03/05/2022: launching the weighted station function to produce weight for each station
    if (passim->weight_calc_assim->weighted_stations == YES_TS)
       {
	 Prose_create_obs_station_weight(passim, Simul->poutputs);
       }

    for(row = 0; row < passim->N_obs; row++)
    {   
        if(passim->pobs[row]->answer_obs == YES_TS) // there is a observation of row thd stations
		{

		  if (passim->weight_calc_assim->weighted_stations == YES_TS)
                  {

		    LP_printf(Simul->poutputs,"Norm_Weight @ obs_station_%d = %f \n",row+1, passim->pobs[row]->weight_station);

		    sigma = passim->error_obs_sigma*passim->pobs[row]->Obs / passim->pobs[row]->weight_station; // MH 03/05/2022 dividing sigma by the weight of each station which is based on the OM conc and the distance of the station from the OM source
                  }
		 else {
                      sigma = passim->error_obs_sigma*passim->pobs[row]->Obs ;
		      }
		cov = pow(sigma,2);
        //printf("row = %d val = %.10f\n",row,cov);
	    cov_matrix->data[col][col] = cov;
		col++;
		}
    }
   return cov_matrix;
}

s_matrix_la * Prose_inv_cov_matrix(s_matrix_la *cov_matrix, int num_t_obs, FILE *fp)
{
    s_matrix_la *inv_cov_matrix;
    int signum, row;
    double var;
    
    inv_cov_matrix = LA_matrix_calloc(num_t_obs,num_t_obs);
   
    /* in case of non diagonal matrix */
    //LA_LU_decomposition_L(cov_matrix, fp);
    //LA_inverse_matrix(cov_matrix, fp);

    for(row = 0; row < num_t_obs; row++)
    {
        var = cov_matrix->data[row][row]; // SW 08/11/2018 diagonal matrix
        var = 1./var;
	inv_cov_matrix->data[row][row] = var;
    }

    return inv_cov_matrix;
}


double Prose_det_cov_matrix(s_matrix_la *cov_matrix, FILE *fp)
{
    int signum,row;
    double det = 1.,var;
 
    for(row = 0; row < cov_matrix->row_size; row++)
    {
        var = cov_matrix->data[row][row];
        det *= var;
    }
	
    return det;
}

s_matrix_la *Prose_calc_matrix_produt(s_carac_assim *passim, s_matrix_la *inv_cov_matrix, int np, int num_t_obs, FILE *fp)
{
    int nobs,nobs_t = 0;
    s_matrix_la *diff_one_par, *trans_diff_one_par;
    s_matrix_la *diff_one_par_left; /* diff^T * inv_cov_matrix in likelihood function*/
    s_matrix_la *result;

    diff_one_par = LA_matrix_calloc(1,num_t_obs);
       
    for(nobs = 0; nobs < passim->N_obs; nobs++)
    {
        if(passim->pobs[nobs]->answer_obs == YES_TS)
        {
            diff_one_par->data[0][nobs_t] = passim->pobs[nobs]->difference[np]; // per particule
            nobs_t++;
        }
    }

    /* diff^T * inv_cov_matrix in likelihood function*/

    diff_one_par_left = LA_calc_mutiply_matrix(diff_one_par, inv_cov_matrix, fp); // memory allocation inside the function
    
    /*diff_one_par_left * diff  = (diff^T * inv_cov_matrix) * diff : diff = trans_diff_one_par */
    trans_diff_one_par = LA_transpose_matrix(diff_one_par, fp); // memory allocation inside the function
    result = LA_calc_mutiply_matrix(diff_one_par_left, trans_diff_one_par, fp); // memory allocation inside the function
        
    LA_free_matrix(diff_one_par, fp);
    LA_free_matrix(diff_one_par_left, fp);
    LA_free_matrix(trans_diff_one_par, fp);

    return result;
}

void Prose_calc_weight_all_particules(s_carac_assim *passim, FILE *fp)
{
    int np, num_t_obs, nparticules,num_threads,taille;
    double log_likelihood, det, produt_matrix,variance;
    s_matrix_la *inv_cov_matrix;
    s_matrix_la *cov_matrix;
    s_matrix_la *result;
    int nobs,nobs_t = 0; // SW debug 02/01/2023

    num_t_obs = passim->num_t_obs;
    nparticules = passim->N_particules;

    cov_matrix = Prose_cov_matrix(passim, fp);
    inv_cov_matrix = Prose_inv_cov_matrix(cov_matrix, num_t_obs, fp);
    //det = Prose_det_cov_matrix(cov_matrix, fp); // SW 16/09/2022 remove constant terms in likelihood

    //#ifdef OMP
    //num_threads = Simul->num_threads_par;
    //taille = PC_set_chunk_size_silent(Simul->poutputs,nparticules,num_threads);
    //omp_set_num_threads(num_threads);
    //#pragma omp parallel for schedule(dynamic,taille) shared(nparticules,inv_cov_matrix,passim,num_t_obs,det,fp) private(np,result,produt_matrix)
    //#endif	
    for(np = 0; np < nparticules; np++)
    {
        result = Prose_calc_matrix_produt(passim, inv_cov_matrix, np, num_t_obs, fp);

        produt_matrix = result->data[0][0];
        //LP_printf(fp,"logpi = %f\n",log(2*M_PI));
        
        log_likelihood = -0.5 * produt_matrix;

        // SW debug 02/01/2023
        /*if(np == 1 || np == 2)
        {
        LP_printf(fp,"np = %d, explog = %.10f, product_matrix = %f\n",np, exp(log_likelihood),produt_matrix);

        for(nobs = 0; nobs < passim->N_obs; nobs++)
        {
            if(passim->pobs[nobs]->answer_obs == YES_TS)
            {
                LP_printf(fp,"np = %d, nobs = %d, diff = %f\n", np, nobs, passim->pobs[nobs]->difference[np]); // per particule
                nobs_t++;
            }
        }
        
        if(np == 1)
        {
            LP_printf(fp, "print cov_matrix :\n");
            LA_print_matrix_type(cov_matrix, MATRIX_A, fp);
            LP_printf(fp, "print inv_cov_matrix :\n");
            LA_print_matrix_type(inv_cov_matrix, MATRIX_A, fp);
        }
        // SW end debug 02/01/2023
        }*/

        passim->omega[np] = passim->omega_prev[np] * exp(log_likelihood);
        LA_free_matrix(result, fp); // free memory of matrix result
    }

    LA_free_matrix(cov_matrix, fp); // free memory of covariance matrix
    LA_free_matrix(inv_cov_matrix, fp); // free memory of covariance matrix

}
 
void Prose_perturbation_parameters(s_simul ***psimul_bio, s_carac_assim *passim, int nele, int state, FILE *fp)
{
	switch(state){
		case HORS_BLOOM:{
			Prose_perturb_bact_param(psimul_bio, passim, nele, fp);
			break;
		}
		case BLOOM: {
			Prose_perturb_phy_param(psimul_bio, passim, nele, fp);
			Prose_perturb_bact_param(psimul_bio, passim, nele, fp);
			break;
		}
	}
}

void Prose_weights_to_weights_prev(s_carac_assim *passim, int np, FILE *fp)
{
	passim->omega_prev[np] = passim->omega_normlized[np];
}

void Prose_determine_trophic_state(s_simul ***psimul_bio, s_carac_assim *passim, FILE *fp)
{
	int np, id_ele, nobs;
	double conc_mean = 0.,c_chla;
	
	for(nobs = 0; nobs < passim->N_obs; nobs++)
	{
		id_ele = passim->pobs[nobs]->id_ele_obs;
	    for(np = 0; np < passim->N_particules; np++)
	    {
			c_chla = psimul_bio[np][id_ele]->section->compartments[WATER][0]->pspecies[PHY][0]->nut_C[NCOMP_RIVE];
		    conc_mean += psimul_bio[np][id_ele]->section->compartments[WATER][0]->pspecies[PHY][0]->C * passim->omega_normlized[np] * c_chla;
	    }
	    passim->pobs[nobs]->chla_mean = conc_mean;
		/* if all stations are in bloom, the state will be bloom; otherwise it is horsbloom*/
	    if(conc_mean < passim->seuil_chla)
		{
			passim->state = HORS_BLOOM;
			break;
		}
		else
			passim->state = BLOOM;
	}
}

/*used defore resampling*/
void PROSE_extraction_parameters(s_simul ***psimul_bio, int nparticules,int param, FILE *fp)
{
	double val;
	int np;
switch(param) {
   case MAINT_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->respiration->resp[MAINT20];
	       Simul->passim->param[MAINT_DA][np] = val; // before resampling
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(MAINT_DA));        
	   break;}
   case ALPHA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[A_RIVE];
	       Simul->passim->param[ALPHA_DA][np] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ALPHA_DA));		   
	   break;}
   case PMAX_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[PMAX20];
	       Simul->passim->param[PMAX_DA][np] = val;	
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(PMAX_DA));		   
	   break;}
   case ETA_CHLA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA_CHLA];
	       Simul->passim->param[ETA_CHLA_DA][np] = val;	
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ETA_CHLA_DA));		   
	   break;}
   case C_CHLA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->nut_C[NCOMP_RIVE]; //C/chla
	       Simul->passim->param[C_CHLA_DA][np] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(C_CHLA_DA));		   
	   break;}
   case ETA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA]; //eta_water
	       Simul->passim->param[ETA_DA][np] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ETA_DA));		   
	   break;}
   case TOPT_PHY_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->paraml[TOPT]; // Topt
	       Simul->passim->param[TOPT_PHY_DA][np] = val;	
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(TOPT_PHY_DA));		   
	   break;}
   case MU_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {

	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[MUMAX20]; // mumax
	       Simul->passim->param[MU_BACT_DA][np] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(MU_BACT_DA));	   
	   break;}
   case Y_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[YIELD];
	       Simul->passim->param[Y_BACT_DA][np] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(Y_BACT_DA));		   
	   break;}
   case MORT_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->mortality->mort[MORT20]; // mort
	       Simul->passim->param[MORT_BACT_DA][np] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(MORT_BACT_DA));		   
	   break;}
   case TOPT_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->paraml[TOPT]; // Topt
           Simul->passim->param[TOPT_BACT_DA][np] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(TOPT_BACT_DA));		   
	   break;}
   case KNAVIG_DA : {
        if(Simul->counter_bio->nsubspecies[O2] > 0)
        {
	   for(np = 0; np < nparticules; np++)
	   {
		   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[O2][0]->dissolved->gas->reaeration->rea[REA_NAVIG]; // Knavig
           Simul->passim->param[KNAVIG_DA][np] = val;
	   }	
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species O2 is not defined\n", PROSE_name_param(KNAVIG_DA));	   
	   break;}
       case B1_RIVER_DA : {
	  if ( (Simul->counter_bio->nmacrospecies > EPS_TS) && (Simul->passim->param_range[B1_RIVER_DA][PARAM_UP] > EPS_TS) ) {
	   for(np = 0; np < nparticules; np++)
	   {
	     val = Simul->passim->p_macro_da[TOC][np]->degradOrgMat[MACMOD][B]->val; // b1_river
	     //val = psimul_bio[np][0]->section->compartments[WATER][0]->pmacro[TOC]->degradOrgMat[MACMOD][B]->val;
	     Simul->passim->param[B1_RIVER_DA][np] = val;
	   }
	  }		   
	   break;}  
   default : {
	   LP_error(Simul->poutputs,"Unknown parameter for assimilation.\n");
	   break;}
   }
}

void Prose_perturb_phy_param(s_simul ***psimul_bio, s_carac_assim *passim, int nele, FILE *fp)
{
	double variance;
	double sd_sigma;
	int nparam;
	int np;
	double s, epsilon,val;
	for(nparam = 0; nparam < PARAM_PHY; nparam++)
	{
            // SW 24/01/2022 add check if parameter is assimilated and PHY species is defined
            if((Simul->counter_bio->nsubspecies[PHY] >0) && (passim->param_yesOrno[nparam] == YES_TS))
            {
		
		sd_sigma = (passim->param_range[nparam][PARAM_UP] - passim->param_range[nparam][PARAM_DOWN]) * passim->s_percent[nparam]; // MH 10/03/2022 : [NPARAMDA] added to facilitate different random walk for each param
	

		//parameters value after resampling
		if(passim->Neff < (passim->alpha*passim->N_particules))
		    PROSE_extraction_parameters(Simul->psimul_bio, passim->N_particules, nparam, fp);

		for(np = 0; np < passim->N_particules; np++)
		{
			epsilon = Prose_generate_random_norm_param(sd_sigma);

			/* SW 20/05/2022 calculation of the density for epsilon, used to recalculate the weight */
                        passim->param_perturbation_density[nparam][np] = PROSE_density_normal(epsilon, 0.0, sd_sigma, fp);

			
			val = passim->param[nparam][np] + epsilon;
                        if(passim->random_walk == LOOP) // SW 06/08/2020
                        {			
				if(val > passim->param_range[nparam][PARAM_UP]) // out of maximum
					val = passim->param_range[nparam][PARAM_DOWN] + (val - passim->param_range[nparam][PARAM_UP]);
				else if(val < passim->param_range[nparam][PARAM_DOWN]) // out of min
					val = passim->param_range[nparam][PARAM_UP] - (passim->param_range[nparam][PARAM_DOWN] - val);
                        }
                        else
                        {			
				if((val > passim->param_range[nparam][PARAM_UP]) || (val < passim->param_range[nparam][PARAM_DOWN]))// out of bound
					val = passim->param[nparam][np]; // SW 21/11/2018 dont change 
                        }

			//passim->param[nparam][np] = val;
			PROSE_assign_parameters_val(psimul_bio, nele, np, nparam, val, fp);

		}
            }
	}
		
}

void Prose_perturb_bact_param(s_simul ***psimul_bio, s_carac_assim *passim, int nele, FILE *fp)
{
	double variance;
	double sd_sigma;
	int nparam;
	int np;
	double s, epsilon,val;
	for(nparam = MU_BACT_DA; nparam < NPARAMDA; nparam++)
	{
            // SW 24/01/2022 add check if parameter is assimilated
            if((Simul->counter_bio->nsubspecies[BACT] >0) && (passim->param_yesOrno[nparam] == YES_TS))
            {

		//variance = gsl_stats_variance(passim->param[nparam],1,passim->N_particules);
		sd_sigma = (passim->param_range[nparam][PARAM_UP] - passim->param_range[nparam][PARAM_DOWN]) * passim->s_percent[nparam];     // MH 10/03/2022 : [NPARAMDA] added to facilitate different random walk for each param

		//LP_printf(Simul->poutputs,"random walk while perturbing of %s = %3.4f \n",PROSE_name_param(nparam),passim->s_percent[nparam]);
		
		//s = passim->s_percent;
		//if(variance == 0)
			//sd_sigma = (passim->param_range[nparam][PARAM_UP] - passim->param_range[nparam][PARAM_DOWN]) * passim->s_percent;
		//else
		   // sd_sigma = s * sqrt(variance);
		//parameters value after resampling
		if(passim->Neff < (passim->alpha*passim->N_particules))
		    PROSE_extraction_parameters(Simul->psimul_bio, passim->N_particules, nparam, fp);

		for(np = 0; np < passim->N_particules; np++)
		{
			epsilon = Prose_generate_random_norm_param(sd_sigma);

			/* SW 20/05/2022 calculation of the density for epsilon, used to recalculate the weight */
                        passim->param_perturbation_density[nparam][np] = PROSE_density_normal(epsilon, 0.0, sd_sigma, fp);

			
			val = passim->param[nparam][np] + epsilon;
                        if(passim->random_walk == LOOP) // SW 06/08/2020
                        {			
				if(val > passim->param_range[nparam][PARAM_UP]) // out of maximum
					val = passim->param_range[nparam][PARAM_DOWN] + (val - passim->param_range[nparam][PARAM_UP]);
				else if(val < passim->param_range[nparam][PARAM_DOWN]) // out of min
					val = passim->param_range[nparam][PARAM_UP] - (passim->param_range[nparam][PARAM_DOWN] - val);
                        }
                        else
                        {			
				if((val > passim->param_range[nparam][PARAM_UP]) || (val < passim->param_range[nparam][PARAM_DOWN]))// out of bound
					val = passim->param[nparam][np]; // SW 21/11/2018 dont change 
                        }

			//passim->param[nparam][np] = val;
			PROSE_assign_parameters_val(psimul_bio, nele, np, nparam, val, fp);

		}
            }
	}		
}

void PROSE_assign_epsilon_parameters(s_simul ***psimul_bio, int nele,int np, int param, double epsilon, FILE *fp)
{
	double val;
	int ne;
switch(param) {
   case MAINT_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->respiration->resp[MAINT20] + epsilon;
	   //val = val > Simul->passim->param_range[MAINT_DA][PARAM_UP] ? Simul->passim->param_range[MAINT_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[MAINT_DA][PARAM_DOWN] ? Simul->passim->param_range[MAINT_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[MAINT_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[MAINT_DA][PARAM_UP]))
	   {
	   for(ne = 0; ne < nele; ne++)
	   {
		   psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->respiration->resp[MAINT20] = val;
                   if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
		       psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->respiration->resp[MAINT20] = val;
	   }
	   
	   }
	   break;
   }
   case ALPHA_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[A_RIVE] + epsilon;
	   //val = val > Simul->passim->param_range[ALPHA_DA][PARAM_UP] ? Simul->passim->param_range[ALPHA_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[ALPHA_DA][PARAM_DOWN] ? Simul->passim->param_range[ALPHA_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[ALPHA_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[ALPHA_DA][PARAM_UP]))
	   {	   
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[A_RIVE] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
                   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[A_RIVE] = val;
	   }		   
	   }
	   break;
   }
   case PMAX_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[PMAX20] + epsilon;	   
	   //val = val > Simul->passim->param_range[PMAX_DA][PARAM_UP] ? Simul->passim->param_range[PMAX_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[PMAX_DA][PARAM_DOWN] ? Simul->passim->param_range[PMAX_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[PMAX_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[PMAX_DA][PARAM_UP]))
	   {
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[PMAX20] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
                  psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[PMAX20] = val;
	   }	   
	   }
	   break;
   }
   case ETA_CHLA_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA_CHLA] + epsilon;
	   //val = val > Simul->passim->param_range[ETA_CHLA_DA][PARAM_UP] ? Simul->passim->param_range[ETA_CHLA_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[ETA_CHLA_DA][PARAM_DOWN] ? Simul->passim->param_range[ETA_CHLA_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[ETA_CHLA_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[ETA_CHLA_DA][PARAM_UP]))
	   {	   
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA_CHLA] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
                   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA_CHLA] = val;
	   }	   
	   }
	   break;
   }
   case C_CHLA_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->nut_C[NCOMP_RIVE] + epsilon;
	   //val = val > Simul->passim->param_range[C_CHLA_DA][PARAM_UP] ? Simul->passim->param_range[C_CHLA_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[C_CHLA_DA][PARAM_DOWN] ? Simul->passim->param_range[C_CHLA_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[C_CHLA_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[C_CHLA_DA][PARAM_UP]))
	   {	   
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->nut_C[NCOMP_RIVE] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
		   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->nut_C[NCOMP_RIVE] = val;
	   }		   
	   }
	   break;
   }
   case ETA_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA] + epsilon;
	   //val = val > Simul->passim->param_range[ETA_DA][PARAM_UP] ? Simul->passim->param_range[ETA_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[ETA_DA][PARAM_DOWN] ? Simul->passim->param_range[ETA_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[ETA_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[ETA_DA][PARAM_UP]))
	   {	   
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
                   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA] = val;
	   }	   
	  }
	   break;
   }
   case TOPT_PHY_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->paraml[TOPT] + epsilon;
	   //val = val > Simul->passim->param_range[TOPT_PHY_DA][PARAM_UP] ? Simul->passim->param_range[TOPT_PHY_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[TOPT_PHY_DA][PARAM_DOWN] ? Simul->passim->param_range[TOPT_PHY_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[TOPT_PHY_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[TOPT_PHY_DA][PARAM_UP]))
	   {	   
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->paraml[TOPT] = val; // Topt
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
		   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->paraml[TOPT] = val; // Topt
	   }		   
	   }
	   break;
   }
   case MU_BACT_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[MUMAX20] + epsilon;
	   //val = val > Simul->passim->param_range[MU_BACT_DA][PARAM_UP] ? Simul->passim->param_range[MU_BACT_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[MU_BACT_DA][PARAM_DOWN] ? Simul->passim->param_range[MU_BACT_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[MU_BACT_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[MU_BACT_DA][PARAM_UP]))
	   {	   
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[MUMAX20] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
		   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[BACT][0]->particulate->living->growth->growth[MUMAX20] = val;
	   }
	   
	   }
	   break;
   }
   case Y_BACT_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[YIELD] + epsilon;
	   //val = val > Simul->passim->param_range[Y_BACT_DA][PARAM_UP] ? Simul->passim->param_range[Y_BACT_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[Y_BACT_DA][PARAM_DOWN] ? Simul->passim->param_range[Y_BACT_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[Y_BACT_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[Y_BACT_DA][PARAM_UP]))
	   {	   
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[YIELD] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
		   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[BACT][0]->particulate->living->growth->growth[YIELD]  = val;
	   }		   
	   }
	   break;
   }
   case MORT_BACT_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->mortality->mort[MORT20] + epsilon;
	   //val = val > Simul->passim->param_range[MORT_BACT_DA][PARAM_UP] ? Simul->passim->param_range[MORT_BACT_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[MORT_BACT_DA][PARAM_DOWN] ? Simul->passim->param_range[MORT_BACT_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[MORT_BACT_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[MORT_BACT_DA][PARAM_UP]))
	   {	   
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->mortality->mort[MORT20] = val; // mort
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)	
                   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[BACT][0]->particulate->living->mortality->mort[MORT20] = val; // mort
	   }		   
	   }
	   break;
   }
   case TOPT_BACT_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->paraml[TOPT] + epsilon;
	   //val = val > Simul->passim->param_range[TOPT_BACT_DA][PARAM_UP] ? Simul->passim->param_range[TOPT_BACT_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[TOPT_BACT_DA][PARAM_DOWN] ? Simul->passim->param_range[TOPT_BACT_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[TOPT_BACT_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[TOPT_BACT_DA][PARAM_UP]))
	   {	   
	   for(ne = 0; ne < nele; ne++)
	   {
           psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->paraml[TOPT] = val; // Topt
           if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
	       psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[BACT][0]->particulate->living->paraml[TOPT] = val; // Topt
	   }		   
	   }
	   break;
   }
   case KNAVIG_DA : {
	   val = psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[O2][0]->dissolved->gas->reaeration->rea[REA_NAVIG] + epsilon;
	   //val = val > Simul->passim->param_range[TOPT_BACT_DA][PARAM_UP] ? Simul->passim->param_range[TOPT_BACT_DA][PARAM_UP] : val;
	   //val = val < Simul->passim->param_range[TOPT_BACT_DA][PARAM_DOWN] ? Simul->passim->param_range[TOPT_BACT_DA][PARAM_DOWN] : val;
	   if((val >= Simul->passim->param_range[REA_NAVIG][PARAM_DOWN]) && (val <= Simul->passim->param_range[REA_NAVIG][PARAM_UP]))
	   {	   
	   for(ne = 0; ne < nele; ne++)
	   {
           psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[O2][0]->dissolved->gas->reaeration->rea[REA_NAVIG] = val; // Topt
           if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
	       psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[O2][0]->dissolved->gas->reaeration->rea[REA_NAVIG] = val; // Topt
	   }		   
	   }
	   break;
   }

   case B1_RIVER_DA : {
      if ( (Simul->counter_bio->nmacrospecies > EPS_TS) && (Simul->passim->param_range[B1_RIVER_DA][PARAM_UP] > EPS_TS ) ) {
     //val = psimul_bio[np][0]->section->compartments[WATER][0]->pmacro[TOC]->degradOrgMat[MACMOD][B]->val + epsilon;
       val = Simul->passim->p_macro_da[TOC][np]->degradOrgMat[MACMOD][B]->val + epsilon;
	   if((val >= Simul->passim->param_range[B1_RIVER_DA][PARAM_DOWN]) && (val <= Simul->passim->param_range[B1_RIVER_DA][PARAM_UP]))
	   {	   
	     /* for(ne = 0; ne < nele; ne++)
	      {
		psimul_bio[np][ne]->section->compartments[WATER][0]->pmacro[TOC]->degradOrgMat[MACMOD][B]->val = val;
		psimul_bio[np][ne]->section->compartments[VASE][0]->pmacro[TOC]->degradOrgMat[MACMOD][B]->val = val;		
		}*/
	    Simul->passim->p_macro_da[TOC][np]->degradOrgMat[MACMOD][B]->val = val;

	   }
      }
	   break;
   }     
   default : {
	   LP_error(Simul->poutputs,"Unknown parameter for assimilation. line 1001\n");
	   break;}
   }
}


void PROSE_assign_parameters_val(s_simul ***psimul_bio, int nele,int np, int param, double val, FILE *fp)
{
	//double val;
	int ne;
switch(param) {
   case MAINT_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
		   psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->respiration->resp[MAINT20] = val;
                   if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
		       psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->respiration->resp[MAINT20] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(MAINT_DA));
	break;}

   case ALPHA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[A_RIVE] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
                   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[A_RIVE] = val;
	   }	
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ALPHA_DA));	   
	break;}
   case PMAX_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[PMAX20] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
                   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[PMAX20] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(PMAX_DA));	   
	break;}
   case ETA_CHLA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA_CHLA] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
                   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA_CHLA] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ETA_CHLA_DA));	   
	   break;}
   case C_CHLA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->nut_C[NCOMP_RIVE] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
		   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->nut_C[NCOMP_RIVE] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(C_CHLA_DA));		   
	break;}
   case ETA_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
                   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->photosynthesis->phot[ETA] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(ETA_DA));	   
	break;}
   case TOPT_PHY_DA : {
        if(Simul->counter_bio->nsubspecies[PHY] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][0]->particulate->living->paraml[TOPT] = val; // Topt
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
		   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[PHY][0]->particulate->living->paraml[TOPT] = val; // Topt
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species PHY is not defined\n", PROSE_name_param(TOPT_PHY_DA));		   
	break;}
   case MU_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[MUMAX20] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
		   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[BACT][0]->particulate->living->growth->growth[MUMAX20] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(MU_BACT_DA));	   
	   break;}
   case Y_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->growth->growth[YIELD] = val;
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
		   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[BACT][0]->particulate->living->growth->growth[YIELD] = val;
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(Y_BACT_DA));		   
	   break;}
   case MORT_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
	       psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->mortality->mort[MORT20] = val; // mort
               if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)	
                   psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[BACT][0]->particulate->living->mortality->mort[MORT20] = val; // mort
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(MORT_BACT_DA));		   
	   break;}
   case TOPT_BACT_DA : {
        if(Simul->counter_bio->nsubspecies[BACT] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
           psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[BACT][0]->particulate->living->paraml[TOPT] = val; // Topt
           if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
	       psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[BACT][0]->particulate->living->paraml[TOPT] = val; // Topt
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species BACT is not defined\n", PROSE_name_param(TOPT_BACT_DA));		   
	   break;}
   case KNAVIG_DA : {
        if(Simul->counter_bio->nsubspecies[O2] > 0)
        {
	   for(ne = 0; ne < nele; ne++)
	   {
           psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[O2][0]->dissolved->gas->reaeration->rea[REA_NAVIG] = val; // Topt
           if(psimul_bio[np][ne]->section->nsublayers[VASE] > 0)
	       psimul_bio[np][ne]->section->compartments[VASE][0]->pspecies[O2][0]->dissolved->gas->reaeration->rea[REA_NAVIG] = val; // Topt
	   //LP_printf(Simul->poutputs,"np = %i) Knavig (passim) =  %3.4f =  %3.4f = Knavig (compart) \n",np, Simul->passim->param[KNAVIG_DA][np], val ); // dont why it gives 0  0
	   
	   }
        }
        else
            LP_warning(fp,"Parameter %s is defined for DA, but species O2 is not defined\n", PROSE_name_param(KNAVIG_DA));		   
	   break;}
    case B1_RIVER_DA : {
      if ( (Simul->counter_bio->nmacrospecies > EPS_TS) && (Simul->passim->param_range[B1_RIVER_DA][PARAM_UP] > EPS_TS) ) {
      Simul->passim->p_macro_da[TOC][np]->degradOrgMat[MACMOD][B]->val = val; // b1_river
      //LP_printf(Simul->poutputs,"Assign param func) np = %i) b1(passim_param) =  %3.4f , b1(passim_pmacro_degr) = %3.4f \n",np, Simul->passim->param[B1_RIVER_DA][np] , Simul->passim->p_macro_da[TOC][np]->degradOrgMat[MACMOD][B]->val );
      }
	  
	    break;}
   default : {
	   LP_error(Simul->poutputs,"Unknown parameter for assimilation. line 1091\n");
	   break;}
 }
}





/* SW 20/05/2022 calculation of density of a sample, drawn from a normal distribution N(mean, sd) */

double PROSE_density_normal(double sample, double mean, double sd, FILE *fp)
{
    double density;

    density = 1.0 / (sd * sqrt(2 * M_PI)) * exp(-1.0/2.0 * pow((sample - mean)/sd, 2));
    //LP_printf(fp, "density of %.12f with mean = %f sd = %.12f is %f\n", sample, mean, sd, density);

    return density;

}


/* SW 20/05/2022 production of the perturbation densities for each particle and return the sum*/

double PROSE_sum_perturbation_densties_productions(s_carac_assim *passim, FILE *fp)
{
    int np, nparam;
    double sum_density = 0., product_density;

    for(np = 0; np < passim->N_particules; np++)
    {
        for(nparam = 0; nparam < NPARAMDA; nparam++)
        {
            passim->param_perturbation_density_product[np] = 1.0;
            if(passim->param_yesOrno[nparam] == YES_TS)
                passim->param_perturbation_density_product[np] *= passim->param_perturbation_density[nparam][np];
        }

        sum_density += passim->param_perturbation_density_product[np];

    }

    return sum_density;
}

/* SW 20/05/2022 new weights obtained by normalizing the density production */

void PROSE_new_weights_after_perturbation(s_carac_assim *passim, FILE *fp)
{
    int np;
    double sum_density;

    sum_density = PROSE_sum_perturbation_densties_productions(passim, fp);

    for(np = 0; np < passim->N_particules; np++)
    {
        passim->omega[np] = passim->param_perturbation_density_product[np] / sum_density;
        passim->omega_normlized[np] = passim->param_perturbation_density_product[np] / sum_density;
        passim->omega_prev[np] = passim->omega_normlized[np];

        //LP_printf(fp, "After resampling, np = %d, old weight = %f, new weight = %f\n", np, 1./passim->N_particules, passim->omega_prev[np]);     
    }
}
