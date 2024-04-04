/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: manage_simulation.c
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
//#include <libprint.h>
#include <libprint.h>
#include "time_series.h"
#include "libpc.h"
#include "IO.h"
#include "GC.h"
#include "CHR.h"
//#include "spa.h"
#include "HYD.h"
#include "TTC.h"
#include "RIVE.h"
#include "SEB.h"
#include "TUB.h"
#include "MB.h"
#include "LA.h"
#include "PROSE.h"

s_simul_PROSE *PROSE_init_simulation()
{
  s_simul_PROSE *psimul;
  int i, nparam;
  psimul=new_simulation_PROSE();
  bzero((char *)psimul,sizeof(s_simul_PROSE));
  
  /*inialization of hydraulic part*/
  psimul->pchyd=HYD_create_chyd();
  psimul->pchyd->settings=HYD_create_sets();
  HYD_init_chyd(psimul->pchyd);
  
  psimul->regime = STEADY;
  /*initialization of transport part*/
  //psimul->pcarac_ttc=TTC_init_carac(); // SW 17/09/2018 mettre dans input.y
  /*initialization of bio part*/
  //psimul->section = init_section()  
  
  psimul->settings = init_settings();
  psimul->counter_bio = new_counter_rive();
  bzero((char *)psimul->counter_bio,sizeof(s_counter));
  psimul->counter_bio->nsubspecies[MOD] = 1;
  psimul->counter_bio->nsubspecies[MOP] = 1;
  psimul->counter_bio->nsubspecies[SI_D] = 0; //SW 25/02/2020 SI_D est null par défaut

  //for ( i = NH4; i < NSPECIES; i++) //SW 25/02/2020 SI_D est null par défaut
  //  psimul->counter_bio->nsubspecies[i] = 1;

  /*** SW 22/03/2023 SODA TA pH DIC CO2 not simulated by default ***/
  //psimul->counter_bio->nsubspecies[SODA] = 0;

  //for ( i = NH4; i < TA; i++)
  //  psimul->counter_bio->nsubspecies[i] = 1;
  
  //for ( i = O2; i < CO2; i++) 
   // psimul->counter_bio->nsubspecies[i] = 1;
  
  /*inialization of chorno part*/
  psimul->chronos = CHR_init_chronos();
  psimul->clock = CHR_init_clock();
  //psimul->clock->begin = time(NULL);

  /* initialization of assimilation part*/
  psimul->passim = new_assimilation();
  bzero((s_carac_assim *)psimul->passim,sizeof(s_carac_assim));
  psimul->passim->alpha = 0.3;
   // MH 10/03/2022 : [NPARAMDA] added to facilitate different random walk for each param   
  for(i=0; i < NPARAMDA; i++)
    {
      psimul->passim->s_percent[i] = 0.1;
      
    }
  //psimul->passim->s_percent = 0.1;
  psimul->passim->error_obs_sigma = 0.01;
  psimul->passim->seuil_chla = 0.;
  psimul->passim->random_walk = LOOP;
  psimul->passim->method = PF_PROSE;
  psimul->passim->nstep = 1; // SW 04/04/2024
  
  // SW 24/01/2022
  for(nparam = 0; nparam < NPARAMDA; nparam++)
      psimul->passim->param_yesOrno[nparam] = NO_TS;
  psimul->passim->num_Of_assimilated_param = 0;

  /*output part*/
  psimul->pinout=IO_create_inout_set(NOUTPUTS);
  psimul->pinout->init_from_file[IQ_IO] = NO_TS;
  psimul->pinout->init_from_file[IZ_IO] = NO_TS;
  psimul->pout = new_out_io();
  psimul->pout->init_from_file=NO_IO;
  //psimul->final_state = YES_TS;
  psimul->psmp = new_smp();
  psimul->psmp->nthreads = 1;
  psimul->passim->N_particules = 1;
  psimul->num_threads_par = 1;
  psimul->npk_mb_bio = 0; // SW 22/05/2019

  psimul->date_format = FR_TS;
  
  psimul->year0_meteo = INITIAL_YEAR_JULIAN_DAY_TS;

  psimul->default_t_inflows = 20; //31/05/2021
  //psimul->total_mb = new_several_total_mb();
  //psimul->total_conc = new_several_conc();
  for (i = 0; i < NOUTPUTS; i++)
    psimul->pinout->calc[i] = NO_TS;
  for(i = 0; i < NMODE; i++)
  {
	  if(i == HYD)
	      psimul->calc_mode[i] = YES_TS;
      else
		  psimul->calc_mode[i] = NO_TS;
  }
  
  for(i = 0; i < NANNEX_VAR; i++)
	  psimul->calc_bio_annexvar[i] = NO_TS; // SW 04/06/2019
  psimul->outputs[PRINT_PK] = (s_output_hyd **)calloc(1,sizeof(s_output_hyd *));
  psimul->outputs[FINAL_STATE] = (s_output_hyd **)calloc(1,sizeof(s_output_hyd *));
  psimul->outputs[ITERATIONS] = (s_output_hyd **)calloc(1,sizeof(s_output_hyd *));
  psimul->outputs[MB_ELE] = (s_output_hyd **)calloc(1,sizeof(s_output_hyd *));//LV 15/06/2012
  psimul->plec=(s_lec_tmp_hyd **)calloc(NB_LEC_MAX,sizeof(s_lec_tmp_hyd *));
 
  psimul->pout_tube = new_tube_out_put();
  bzero((s_output_tube *)psimul->pout_tube, sizeof(s_output_tube));
  psimul->pout_tube->calc_output = NO_TS;

  /* MH 12/12/2022: initialization of moving average variables   */ 
  psimul->passim->da_mov_ave = NO_TS;
  psimul->passim->mov_ave_theta = CODE_TS ;
  psimul->passim->lneigh = CODE_TS ;
  psimul->passim->nval_min = CODE_TS ;

   /* MH 31/05/2023: initialization of station weights: removal of bug */
  psimul->passim->weight_calc_assim = PROSE_init_calc_stat_weight();
  psimul->passim->weight_calc_assim->weighted_stations = NO_TS;
  
  for (i=0; i < psimul->passim->N_obs; i++)
    {
       psimul->passim->pobs[i]->weight_i = 0.;
       psimul->passim->pobs[i]->weight_station = 0.;
    }
 


  return psimul;
}
