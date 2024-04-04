/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: main_PROSE.c
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

//Standard libraries
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <malloc.h>
#include <math.h>
//#include <fenv.h>
//#include <signal.h>
#ifdef OMP
#include <omp.h>
#endif
//Dependency
#include <libprint.h>
#include <time_series.h>
#include "libpc.h"
#include "IO.h"
#include "GC.h"
#include "spmatrix.h"
#include "CHR.h"
//#include "spa.h"
#include "HYD.h"
#include "TTC.h"
#include "RIVE.h"
#include "SEB.h"
#include "TUB.h"
#include "MB.h"
#include "LA.h"
//prose-pa .h
#include "PROSE.h"
//#include "global_PROSE.h"
#include "ext_PROSE.h"

int main(int argc, char **argv)
{               
  /* Date of the beginning of the simulation */
  char *date_of_simulation;
  /* Structure from time.h containing date, time, etc */
  static struct timeval current_time;
  /* Structure containing computer clock time */
  struct tm *clock_time;
  /* Current time in seconds when programm starts */
  time_t clock;
  /* Error report */
  int timi;
  /* Current time in simulation */
  double t,dt,t_abs,dt_day,t_da;
  /* Time of end of the simulation in s */
  double tend;
  double *tempe;
  /*number of species for transport*/
  int nspecies, nele, nparticules, np, nparm, nreaches;
  int Neff,answer_obs = NO_TS;
  int num_threads,taille,ns2;
  double tempe1,Osat;
  //int ra,rd;
  /*** HT & SEB ***/
  int i, e, r;
  int index_nele = 0;
  s_met **p_met;
  s_rts **p_rts;
  s_rs_input **p_rside;
  int n_meteocell;
  int *id_meteo;
  s_carac_seb **p_surf_heat_flux;
  double H_tot;
  double *H_flux_ttc;
  double d;
  FILE *outfile1, *outfile2;
  int month_ind = 0, h_ind = 0, d_ind = 0;
  int num_met; double temp = 0., temp_air = 0., temp2 = 0., temp2_air = 0.;
  /****************/
  /* tube */
  int t_index_tube = 0, t_index2_tube = 0;
  /*test*/
  //FILE *fpno3 = NULL;
  //s_species_ttc *pspecies;
  //s_param_calc_ttc *pparam_calc_ttc;
  s_mb_hyd *pmb;// Mass balance structure containing the mass balance at t if asked for and chained with the previous mass balance pmb2
  s_mb_hyd *pmb2;// Mass balance structure containing the mass balance at t-1 if asked for
  pmb=NULL;
  pmb2=NULL;
  
  /* Creates the simulation structure */
  //fprintf(stdout,"argc = %d",argc);
  Simul=PROSE_init_simulation();
  
  /* Beginning of the simulation */
  Simul->clock->begin = time(NULL);
  
  /* Date */
  timi = gettimeofday(&current_time,NULL);
  clock = (time_t)current_time.tv_sec;
  ctime(&clock);
  clock_time = localtime(&clock);
  date_of_simulation = asctime(clock_time);
  //feenableexcept(FE_OVERFLOW); // | FE_DIVBYZERO);
  
  /*** Essentials for temperature modelling ***/
  Simul->init_from_file = NO_TS;
  Simul->calc_mode[H_T] = NO_TS;
  Simul->calc_mode[SEB] = NO_TS;
  /******************/
  
  /* Opening of command file */
  if (argc != 3) 
    LP_error(stderr,"WRONG COMMAND LINE : Start the code properly with 2 arguments\n > %s%4.2f CommandFileName DebugFileName\n",CODE_NAME,NVERSION_PROSE);
  
  if ((Simul->poutputs = fopen(argv[2],"w")) == 0) 
    LP_error(stderr,"Impossible to open the debug file %s\n",argv[2]);
  
    LP_log_file_header(Simul->poutputs,CODE_NAME,NVERSION_PROSE,date_of_simulation,argv[1]);
  
  /********************************************************/
  /*    Now input files are read  --> see file input.y    */
  /********************************************************/
  
  CHR_begin_timer();
  lecture(argv[1],Simul->poutputs);
  Simul->clock->time_spent[LEC_CHR] += CHR_end_timer();//LV 3/09/2012
  
  /********************************************************/
  /*           Starting HYDRO calculations                */
  /********************************************************/

  /*** Time parameters ***/
  t = Simul->chronos->t[BEGINNING_CHR];
  Simul->chronos->t[CUR_CHR] = Simul->chronos->t[BEGINNING_CHR];
  if(Simul->calc_mode[SEB] == YES_TS) // SW 04/12/2019
    t_abs = Simul->t_extrema[BEGINNING_CHR]; // Takes its origin in function of the day given by the user. SW 27/04/2021 julian day with origin = year0_meteo for t
  tend = Simul->chronos->t[END_CHR];
  dt = Simul->chronos->dt;
  dt_day = dt/(24.*3600.);
  CHR_calculate_simulation_date(Simul->chronos,t);

  /* Initialization of the hydraulic characteristics in the network *///LV nov2014 déplacé ici
  CHR_begin_timer();

  HYD_initialize_hydro(Simul->pchyd,Simul->pinout,Simul->pout,Simul->chronos,Simul->poutputs);
  HYD_assign_river(Simul->pchyd);//LV nov2014

  Simul->clock->time_spent[INIT_CHR] += CHR_end_timer();//LV 3/09/2012
 
  CHR_begin_timer();//LV 3/09/2012  
 
  Simul->pinout->fmb=HYD_create_header_mass_balance_file(Simul->poutputs);

  if (Simul->pinout->calc[PRINT_PK] == YES_TS)
    HYD_print_pk(Simul->outputs,Simul->pinout,Simul->pchyd,Simul->poutputs);
  
  if (Simul->pinout->calc[MB_ELE] == YES_TS)
    HYD_create_mb_at_elements_file(Simul->outputs,Simul->pinout,Simul->poutputs);

  HYD_print_outputs_formats(Simul->outputs,Simul->pinout,Simul->pchyd,Simul->poutputs);
  if(Simul->calc_mode[H_T] == YES_TS) //SW 26/01/2021 add if
      PROSE_print_outputs_formats_temp(Simul->outputs,Simul->pinout,Simul->pchyd,Simul->poutputs);//NF 12/10/2020
  
  Simul->clock->time_spent[OUTPUTS_CHR] += CHR_end_timer();//LV 3/09/2012
  //LV nov2014 : initialisation des coefs de la matrice de calcul de l'hydro
  HYD_initialize_GC(Simul->pchyd,Simul->chronos,Simul->poutputs);
  
  //printf("yes = %d \n",YES_TS);
  /*transport*/
  nparticules = Simul->passim->N_particules;

  Neff = nparticules; // SW 06/01/2023

  if(Simul->pcarac_ttc != NULL)
      nele = Simul->pcarac_ttc[0]->count[NELE_TTC];

  if((Simul->pcarac_ttc != NULL) && (Simul->solver == SP_PROSE)) // SW 28/02/2019
  {
    nspecies = Simul->pcarac_ttc[0]->count[NSPECIES_TTC];
    
    PROSE_allocate_spmatrix(nspecies,nele,nparticules,Simul->poutputs);
    PROSE_allocate_spmatrix_phy(nele, nparticules,Simul->poutputs);
  }
  /*** SW 05/12/2019 add initialization of sparse matrix***/
  if((Simul->pcarac_heat_ttc != NULL) && (Simul->calc_mode[H_T] == YES_TS) && (Simul->solver == SP_PROSE))
  {
	int error;
    if (Simul->calc_mode[TTC] == NO_TS)
	  nele = Simul->pcarac_heat_ttc->count[NELE_TTC];
	Simul->RHS_b_T = (spREAL *)calloc(nele+1,sizeof(double));
	Simul->mat_adv_T = spCreate(nele,0,&error);
    bzero((char *) Simul->RHS_b_T,(nele+1)*sizeof(double));	  
  }
  //tempe = (double *) malloc(nparticules*sizeof(double));
  
  if(Simul->calc_mode[RIVE] == YES_TS) // SW 28/02/2019
    //if(Simul->calc_mode[RIVE] == YES_TS || Simul->calc_mode[TTC] == YES_TS)  // SW 30/04/2019 // NFSW 21/01/2021
    {
      PROSE_init_output_simulation(nparticules, Simul->poutputs);
      
      /*#ifdef OMP
      num_threads = Simul->num_threads_par;
      taille = PC_set_chunk_size_silent(Simul->poutputs,nparticules,num_threads);	
      omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic,taille)  private(np)
#endif*/
      for(np = 0; np < nparticules; np++)
	{
	  //tempe[np] = 0.0;
	  PROSE_print_outputs_formats_bio(Simul->outputs,Simul->pinout,Simul->pchyd,np,Simul->poutputs);
	}
    }
  if(Simul->calc_mode[DA] == YES_TS)
    {
      Prose_init_assimilation(Simul->passim, nparticules, Simul->poutputs);
      Prose_init_rand_param(Simul->psimul_bio, Simul->passim, nele, Simul->poutputs);
      //LP_printf(Simul->poutputs,"dbug1\n");
      PROSE_create_files_weights(nparticules, Simul->poutputs);
      PROSE_create_files_parameters(nparticules, Simul->poutputs);
	  // SW 03/02/2020 print initial parameters and weights
	  for(nparm = 0; nparm < NPARAMDA; nparm++)
          {
              // SW 25/01/2022 check if parameter is assimilated
              if(Simul->passim->param_yesOrno[nparm] == YES_TS)
                  PROSE_print_parameters(Simul->psimul_bio, nparticules, nparm, t, Simul->poutputs);
          }
	  PROSE_print_weights(t, nparticules, Simul->poutputs);	
	  }
  //if(Simul->calc_mode[TTC] == YES_TS){
    if(Simul->calc_mode[RIVE] == YES_TS){
    /*#ifdef CDA     
#ifdef OMP
    num_threads = Simul->num_threads_par;
    taille = PC_set_chunk_size_silent(Simul->poutputs,nparticules,num_threads);
    omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic,taille) shared(nspecies,nele) private(np)
#endif
#endif*/
    for(np = 0; np < nparticules; np++) // possible to parallize
      {
        //if(np == 0)
	//{
	PROSE_init_reoxy_species(Simul->pcarac_ttc[np]->p_species, nspecies, nele, np,Simul->poutputs);
	PROSE_fill_param_base_all_species(Simul->pcarac_ttc[np]->p_species, nspecies, nele, np, Simul->poutputs);
	PROSE_fill_param_base_all_annex_species(Simul->p_phy_species[np], nele, np,Simul->poutputs);
	PROSE_alloc_u_dist_all_species(Simul->pcarac_ttc[np]->p_species, nspecies, nele, Simul->pchyd,Simul->poutputs);
	PROSE_alloc_u_dist_all_annex_species(Simul->p_phy_species[np], nele , Simul->pchyd,Simul->poutputs);
	PROSE_fill_neighbor_all_species(Simul->pcarac_ttc[np]->p_species, nspecies, nele, Simul->pchyd,Simul->poutputs);
	PROSE_fill_neighbor_all_annex_species(Simul->p_phy_species[np], nele, Simul->pchyd,Simul->poutputs);
	//PROSE_set_iapplic_all_species(Simul->pcarac_ttc->p_species, nspecies,nele,Simul->poutputs);
	//PROSE_set_iapplic_all_annex_species(Simul->p_phy_species, nele,Simul->poutputs);	 
	//}
      }
    LP_printf(Simul->poutputs,"number of species = %d, nele = %d\n",nspecies,nele);
  }
  
  /*heat transport + surface exchange budget*/
  if(Simul->calc_mode[H_T] == YES_TS)
    {
      /*** Essentials ***/
      //if (Simul->pcarac_heat_ttc->regime==STEADY_TTC)
      Simul->pcarac_heat_ttc->theta = 1.; 
      /******************/
      
      /*** Specie "Temp" initialization ***/
      if (Simul->calc_mode[TTC] == NO_TS)
	 nele = Simul->pcarac_heat_ttc->count[NELE_TTC];

      printf("nele: %d\n",nele);
      PROSE_init_temp_specie(Simul->pcarac_heat_ttc, nele);
      
      /*** Creatig lmat5,mat4 & mat5 ***/
      Simul->pcarac_heat_ttc->p_species[0]->pgc=GC_create_gc(TTC_GC);

      /*** Filling in the parameters for the heat transport (similar to the other species above) ***/ 
      PROSE_init_heat_transport(Simul->pcarac_heat_ttc->p_species[0], Simul->pchyd, nele, t, Simul->poutputs);

     
      /*** Meteorological parameters: filling, linking and calculating ***/
      if(Simul->calc_mode[SEB] == YES_TS)
	{
	  p_met = Simul->p_met;
	  p_rts = Simul->p_rts;
	  n_meteocell = Simul->n_meteocell;
	  
	  /****************************************************************************/
	  /*** Table of size nele, that links id_abs (of the elements) to id_meteo: ***/
	  /*** the indices correspond to id_abs and the values are int id_meteo. ******/
	  /*** This is an important function... ***************************************/
	  /****************************************************************************/
	  nreaches = Simul->pchyd->counter->nreaches;
	  id_meteo = PROSE_link_icell_to_imet_SEB(nele, n_meteocell, nreaches, p_rts, Simul->pchyd->p_reach, p_met); // SW 14/06/2021 add p_met
	  /*** p_rts est libérée à fin de la routine. ***/

	  /***************************/
	  /*** Initialization: t=0 ***/
	  /***************************/
	  
	  p_surf_heat_flux = SEB_alloc_shf(nele);
	  p_rside = SEB_alloc_rs_input(nreaches);
	  H_flux_ttc = (double *)malloc(nele*sizeof(double));
	  for (i=0; i<nele; i++)
	    {
	      p_rside[0]->val[VTS] = 1.0 ; //VIEWTOWHISKY; /*** !!! To be filled in through the input.y !!! --> parameter of vegetal cover ***/
	      p_surf_heat_flux[i]->pH_inputs->r_side = p_rside[0]->val;
	      p_surf_heat_flux[i]->pH_inputs->Tw = Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i];
	    }
          // SW 14/06/2021 add here allocation and initialization of meteo
          PROSE_alloc_inputs_meteo_all(nele, p_surf_heat_flux, NMETEO);

	}
    }
  
  /*biology*/
  if(Simul->calc_mode[RIVE] == YES_TS) // SW 25/01/2021
  //if(Simul->calc_mode[TTC] == YES_TS) 	  
    {
      /*#ifdef CDA
#ifdef OMP
      num_threads = Simul->num_threads_par;
      taille = PC_set_chunk_size_silent(Simul->poutputs,nparticules,num_threads);
      omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic,taille) shared(nele,t,nparticules) private(np)
#endif	
#endif	*/
      for(np = 0; np < nparticules; np++) // possible to parallize
	{
	  //if(np == 0)
	  //{
	  //tempe = 0.0;
	  PROSE_init_param_reaction_state_comp(Simul->psimul_bio[np], nele, Simul->poutputs);
	  //if(Simul->calc_mode[RIVE] == YES_TS){
	  /*** SW 05/12/2019 initialize temperature if heat transport is calulated***/
      if(Simul->calc_mode[H_T] == YES_TS)
	  {
		  for (i = 0; i < nele; i++)
			  Simul->psimul_bio[np][i]->section->meteo->tempe_value = Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i] - T_0;	 
	  }
	  PROSE_init_O2_conc_sat(Simul->psimul_bio[np], t, nele,Simul->poutputs);
	  PROSE_link_hyd_rive_geom(Simul->psimul_bio[np], Simul->pchyd, Simul->psmp->nthreads,Simul->poutputs);
	  // files bilan de masse a faire
	  PROSE_link_hyd_rive_hydro(Simul->psimul_bio[np],t, Simul->pchyd, Simul->psmp->nthreads,Simul->poutputs);
	  PROSE_calc_volume_water_all_sections(Simul->psimul_bio[np], nele, t, Simul->psmp->nthreads,Simul->poutputs);
	  PROSE_create_files_mb_bio(nparticules, np, Simul->poutputs);
	  //} 
	}	  
    }	  
  
  //if(Simul->calc_mode[TTC] == YES_TS) // for PHY species to PHYF PHYR PHYS
  //{
  //for(np = 0; np < nparticules; np++) // possible to parallize
  //{
  // PROSE_fill_var_all_annex_species(Simul->p_phy_species[np],nele , np,Simul->poutputs);
  /*fill init val*/
  //PROSE_fill_var_all_species(Simul->pcarac_ttc[np]->p_species, nspecies, nele, np,Simul->poutputs);
  //}
  //}
  
  if (Simul->pchyd->settings->calc_state == STEADY) {
    if(Simul->pchyd->settings->schem_type != MUSKINGUM)
      {
        //LP_printf(Simul->poutputs, "ok 111\n");
	HYD_steady_state(dt,pmb,pmb2,Simul->pinout->fmb,Simul->pinout,Simul->outputs,Simul->pchyd,Simul->clock,Simul->chronos,Simul->poutputs);//hydraulique_permanent(dt_calc) in ProSe
	//print_LHS_RHS(Simul->pchyd->calcul_hydro);
	CHR_begin_timer();//LV 3/09/2012
	//HYD_print_outputs_old(t,Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos,Simul->poutputs);
	PROSE_print_outputs_bio(t,Simul->psimul_bio[0],Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos, 0, Simul->poutputs); // SW 20/05/2019
	Simul->clock->time_spent[OUTPUTS_CHR] += CHR_end_timer();//LV 3/09/2012
      }
    else
      LP_error(Simul->poutputs,"No steady state for %s scheme \n check input files",HYD_name_scheme(Simul->pchyd->settings->schem_type));

    /* tube */
    if(Simul->calc_mode[TUB] == YES_TS)
      {
	Simul->p3_tube = build_tube_domain_for_all_reaches_TUB(Simul->pchyd,Simul->ntube_default);
	build_singularities_linkage_TUB(Simul->pchyd, Simul->p3_tube);
	write_wtk_file_for_sig_TUB(Simul->pinout, t, Simul->pchyd, Simul->chronos, Simul->p3_tube, Simul->poutputs);
      }
	
    /*transport*/
    if(Simul->calc_mode[TTC] == YES_TS) // SW 25/01/2021
    { 
      if(Simul->calc_mode[RIVE] == YES_TS) 
      { 
       np = 0;
	  
      tempe1 = calc_temperature(t,Simul->psimul_bio[0][0]);
      calc_O2_sat(tempe1,Simul->psimul_bio[0][0]);
      Osat = Simul->psimul_bio[0][0]->section->compartments[WATER][0]->pspecies[O2][0]->dissolved->gas->Csat;
      
      PROSE_link_hyd_rive_hydro(Simul->psimul_bio[np],t, Simul->pchyd, Simul->psmp->nthreads,Simul->poutputs);
      PROSE_link_hyd_rive_geom(Simul->psimul_bio[np], Simul->pchyd, Simul->psmp->nthreads,Simul->poutputs);
      //PROSE_calc_volume_water_all_sections(Simul->psimul_bio[np], nele, t, Simul->poutputs);	   
      PROSE_ttc_all_species(Simul->pcarac_ttc[np], Simul->pchyd, dt, t, np,tempe1,Osat,Simul->poutputs);
      PROSE_update_phy_ctot(Simul->psimul_bio[np], nele, np,Simul->poutputs);
      PROSE_update_conc_bio_all_species(Simul->pcarac_ttc[np]->p_species, nspecies, nele, Simul->psmp->nthreads,np,Simul->poutputs);
      //PROSE_update_conc_bio_all_annex_species(Simul->p_phy_species, nele, Simul->poutputs);		  
      //}
      PROSE_print_conc_init(nele, np, Simul->pinout, Simul->poutputs);
      }

    /* heat transport specifically */
    if(Simul->calc_mode[H_T] == YES_TS) // SW 25/01/2021
      {	
	if(Simul->calc_mode[SEB] == YES_TS)
	  PROSE_update_meteo_for_HT(nele, n_meteocell, t_abs, Simul->pchyd, p_surf_heat_flux, p_met, id_meteo, H_flux_ttc);
	PROSE_update_heat_transport(Simul->pcarac_heat_ttc->p_species[0], Simul->pchyd, H_flux_ttc, nele, t, Simul->poutputs); // ! Caution, time corresponds to the exact date and hour..!
	
        /*** SW 21/05/2021 test iteration for steady state when explicit uptake ***/
        double temp_var[nele], error;
        int niter;
        for(i=0; i<nele; i++)
             temp_var[i] = Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i];

	/*** SW 05/12/2019 add solver sparse***/
	if(Simul->solver == SP_PROSE)
	  PROSE_ttc_one_species_sparse(Simul->pcarac_heat_ttc,Simul->pcarac_heat_ttc->p_species[0]->plink->pparam_calc_ttc,Simul->pcarac_heat_ttc->p_species[0],dt,t,Simul->RHS_b_T, Simul->mat_adv_T,1,Simul->poutputs);
	else
	  PROSE_ttc_one_species(Simul->pcarac_heat_ttc,Simul->pcarac_heat_ttc->p_species[0]->plink->pparam_calc_ttc,Simul->pcarac_heat_ttc->p_species[0],dt,t,Simul->poutputs);
        
        /*** SW 21/05/2021 test iteration for steady state when explicit uptake ***/
            error = 0;
            for(i=0; i<nele; i++)
            {
                error = error > fabs(Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i] - temp_var[i]) ? error : fabs(Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i] - temp_var[i]);
                temp_var[i] = Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i];
            }
            niter = 1;
            while(error > Simul->settings->epsilon && niter < 1000)
            {
                
                PROSE_ttc_one_species_sparse(Simul->pcarac_heat_ttc,Simul->pcarac_heat_ttc->p_species[0]->plink->pparam_calc_ttc,Simul->pcarac_heat_ttc->p_species[0],dt,t,Simul->RHS_b_T, Simul->mat_adv_T,1,Simul->poutputs);
                niter++;
                for(i=0; i<nele; i++)
                {
                    error = error > fabs(Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i] - temp_var[i]) ? error : fabs(Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i] - temp_var[i]);
                    temp_var[i] = Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i];
                }
            }
            if(error > Simul->settings->epsilon)
                LP_printf(Simul->poutputs, "No convergence of the steady-state T after %d iterations.\n", niter);
            else
                LP_printf(Simul->poutputs, "Steady-state of T calculated with %d iterations\n", niter);
        /*** SW 21/05/2021 test iteration for steady state when explicit uptake ***/

	/*** Ouputs for heat transport ***/
	PROSE_print_outputs_heat(t, p_surf_heat_flux, Simul->outputs, Simul->pinout, Simul->pchyd, Simul->chronos, Simul->poutputs);

	//for (i=0; i<nele; i++) // SW 25/01/2021 no use for steady
	  //{
	   // fprintf(Simul->poutputs, "ind: %d T: %lf\n", i, Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i]);
	    //if(Simul->calc_mode[RIVE] == YES_TS)
	      //Simul->psimul_bio[np][i]->section->meteo->tempe_value = Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i] - T_0;	
	  //}
      }
    }
  } 
  
  if (Simul->pchyd->settings->calc_state == TRANSIENT) {
    
    CHR_begin_timer();//LV 3/09/2012
    
    CHR_update_specific_day_from_time(Simul->chronos,CUR_CHR,Simul->poutputs); /* SW 06/01/2021 update current date */
    //HYD_print_outputs_old(t,Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos,Simul->poutputs);//LV test 26/07/2012
    for(np = 0; np < nparticules; np++)
      {
	PROSE_print_outputs_bio(t,Simul->psimul_bio[np],Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos, np, Simul->poutputs); // SW 20/05/2019
      }
    // SW 10/01/2022 bug print
    PROSE_update_outputs_hyd_bio_heat_cur_io(t,Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos, 0, Simul->poutputs); // SW 26/01/2021

  
    Simul->clock->time_spent[OUTPUTS_CHR] += CHR_end_timer();//LV 3/09/2012
    //if(Simul->calc_mode[TTC] == YES_TS)
    //fpno3 = fopen("C_no3_simul_hyd_ttc.txt","w"); // SW test
    
    while (t <= tend) {
      
      LP_printf(Simul->poutputs,"t = %f days\n",t/3600/24);
      
      //t += dt; // SW 23/05/2018 elle est deja dans HYD_transient_hydraulics
      t_da += Simul->chronos->dt * Simul->passim->nstep;
      HYD_transient_hydraulics(&t,dt,pmb,pmb2,Simul->pinout->fmb,Simul->pchyd,Simul->clock,Simul->chronos,Simul->outputs,Simul->pinout,Simul->poutputs);//hydrau_moy(&t,&k,dt_calc) in ProSe 
      //CHR_refresh_t(Simul->chronos); /* SW 06/01/2021 update current time */ // SW 29/01/2021 bug, since we update t[BEGINNING_CHR] also which is used for initialization of mat4 and mat5
      Simul->chronos->t[CUR_CHR] += Simul->chronos->dt;
      CHR_update_specific_day_from_time(Simul->chronos,CUR_CHR,Simul->poutputs); /* SW 06/01/2021 update current date);*/
      
      if((Simul->calc_mode[RIVE] == NO_TS) && (Simul->calc_mode[HYD] == YES_TS) && (Simul->calc_mode[TUB] == NO_TS) ) // SW 20/05/2019 only hydraulic simulation
	{
	  CHR_begin_timer();//LV 3/09/2012
	  PROSE_print_outputs_bio(t,Simul->psimul_bio[0],Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos, 0, Simul->poutputs); // SW 20/05/2019
	  Simul->clock->time_spent[OUTPUTS_CHR] += CHR_end_timer();//LV 3/09/2012
	}
      //LP_printf(Simul->poutputs,"t3 = %f days\n",t/3600/24);
      CHR_calculate_simulation_date(Simul->chronos,dt);

      /*** AB 15/09/2020: TUBEs' generation ***/
      /*** SW 30/03/2021: print in PROSE_print_outputs_bio ***/
      if(Simul->calc_mode[TUB] == YES_TS)
	{
	  if (t - dt < Simul->chronos->t[BEGINNING_CHR] + 0.001)
	    {
	      Simul->p3_tube = build_tube_domain_for_all_reaches_TUB(Simul->pchyd,Simul->ntube_default);
	      build_singularities_linkage_TUB(Simul->pchyd, Simul->p3_tube);
	    }
	  //if (t_index_tube + EPS_HYD >= t_index2_tube)
            
           if(Simul->pout_tube->calc_output == YES_TS)
            {
                /* SW 02/04/2021 update tube only when need to print tube */
                double t0;
                int print_tube_flag = NO_TS, ntb_out;
           
                for(ntb_out = 0; ntb_out < Simul->pout_tube->nout; ntb_out++)
                { 
                    t0 = Simul->pout_tube->poutput_tube_type[ntb_out]->t_out[CUR_IO];
                    if((t >= Simul->pout_tube->poutput_tube_type[ntb_out]->t_out[BEGINNING]) && (t <= t0) && 
	               (t + Simul->chronos->dt > t0) && (t <= Simul->pout_tube->poutput_tube_type[ntb_out]->t_out[END]))
                    {
                        print_tube_flag = YES_TS;
                        break;
                    }
                }
               if(print_tube_flag == YES_TS)
               {
	          update_tube_domain_TUB(Simul->p3_tube, Simul->pchyd); // SW 30/03/2021 if we use tube for transport and bio, update_tube_domain_TUB must be done at each time
	          //time_geo_outputs_TUB(Simul->pinout, t_index_tube, Simul->pchyd, Simul->p3_tube, Simul->poutputs);
	          //t_index2_tube = t_index2_tube + 24*3600;
                  if(Simul->calc_mode[RIVE] == NO_TS) // SW 30/03/2021 if rive = YES, print after rive calculation, need to implement the code
                      PROSE_print_outputs_bio(t,Simul->psimul_bio[0],Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos, 0, Simul->poutputs);
               }

	    }
	 //t_index_tube = t_index_tube + dt;
	}
       
      
      /*** SW 04/12/2019 moved heat transport here, same for all particles ***/
      if(Simul->calc_mode[H_T] == YES_TS)
	{
	  /*****************************************/
          if(Simul->calc_mode[EB_HEAT] == YES_TS) // SW 04/05/2021 add energy balance, need SEB
              PROSE_calc_energy_init(t, Simul->chronos->dt, Simul->pcarac_heat_ttc->p_species[0], Simul->mb_heat, Simul->pchyd, Simul->poutputs);

	  if(Simul->calc_mode[SEB] == YES_TS)
          {
              CHR_begin_timer();
	      //LP_printf(Simul->poutputs, "debug, entering libseb\n");	      
	      PROSE_update_meteo_for_HT(nele, n_meteocell, t_abs, Simul->pchyd, p_surf_heat_flux, p_met, id_meteo, H_flux_ttc);
              Simul->clock->time_spent[SEB_CHR] += CHR_end_timer();
	      //LP_printf(Simul->poutputs, "debug, leaving libseb\n");

              if(Simul->calc_mode[EB_HEAT] == YES_TS) // SW 04/05/2021 add energy balance, need SEB
                  PROSE_calc_flux_from_seb(t, Simul->chronos->dt, Simul->pcarac_heat_ttc->p_species[0], Simul->mb_heat, Simul->pchyd, p_surf_heat_flux, Simul->poutputs);
          }
	  PROSE_update_heat_transport(Simul->pcarac_heat_ttc->p_species[0], Simul->pchyd, H_flux_ttc, nele, t, Simul->poutputs);
	  
	  /*** SW 05/12/2019 add solver sparse***/
          CHR_begin_timer();

	  if(Simul->solver == SP_PROSE)
	    PROSE_ttc_one_species_sparse(Simul->pcarac_heat_ttc,Simul->pcarac_heat_ttc->p_species[0]->plink->pparam_calc_ttc,Simul->pcarac_heat_ttc->p_species[0],dt,t,Simul->RHS_b_T, Simul->mat_adv_T,1,Simul->poutputs);
	  else
	    PROSE_ttc_one_species(Simul->pcarac_heat_ttc,Simul->pcarac_heat_ttc->p_species[0]->plink->pparam_calc_ttc,Simul->pcarac_heat_ttc->p_species[0],dt,t,Simul->poutputs);
	 
          Simul->clock->time_spent[SOLVE_TTC] += CHR_end_timer();
          /*****************************************/
          
          if(Simul->calc_mode[EB_HEAT] == YES_TS) // SW 04/05/2021 add energy balance, need SEB
              PROSE_calc_all_ebs_temperature_and_energy_end(Simul->pcarac_heat_ttc, t, dt, H_flux_ttc, Simul->poutputs);
          
          CHR_begin_timer();

	  PROSE_print_outputs_heat(t, p_surf_heat_flux, Simul->outputs, Simul->pinout, Simul->pchyd, Simul->chronos, Simul->poutputs);
	  
	  Simul->clock->time_spent[OUTPUTS_CHR] += CHR_end_timer();

	  /*** AB end of 2019 update water temp to libseb ***/
          CHR_begin_timer();
	  if(Simul->calc_mode[SEB] == YES_TS)
	    for (i=0; i<nele; i++)
	      p_surf_heat_flux[i]->pH_inputs->Tw = Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i];
          Simul->clock->time_spent[SEB_CHR] += CHR_end_timer();
	  
	  
	  /*** SW 04/12/2019 update temperature to rive***/
	  if(Simul->calc_mode[RIVE] == YES_TS)
	    {
	      for(np = 0; np < nparticules; np++)
		{
		  for (i = 0; i < nele; i++)
		    Simul->psimul_bio[np][i]->section->meteo->tempe_value = Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[i] - T_0;	 
		}
	    }
	  if(Simul->calc_mode[SEB] == YES_TS)
          {
           //   PROSE_free_inputs_meteo_all(nele, p_surf_heat_flux); /* SW 28/01/2021 need to free input meteo pointer */
              t_abs = t_abs + dt_day;
          }
	}
      //LP_printf(Simul->poutputs,"t4 = %f days\n",t/3600/24);
      /*transport*/
      if(Simul->calc_mode[TTC] == YES_TS && Simul->calc_mode[RIVE] == YES_TS)
	{
	  tempe1 = calc_temperature(t,Simul->psimul_bio[0][0]);
	  calc_O2_sat(tempe1,Simul->psimul_bio[0][0]);
	  Osat = Simul->psimul_bio[0][0]->section->compartments[WATER][0]->pspecies[O2][0]->dissolved->gas->Csat;
	  CHR_begin_timer();		
          //LP_printf(Simul->poutputs,"t5 = %f days\n",t/3600/24);

#ifdef CDA
#ifdef OMP
	  num_threads = Simul->num_threads_par;
	  taille = PC_set_chunk_size_silent(Simul->poutputs,nparticules,num_threads);
	  omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic,taille) shared(nele,dt,t,tempe1,Osat) private(np)
#endif	
#endif		
	  for(np = 0; np < nparticules; np++)
	    {    
	      PROSE_ttc_all_species(Simul->pcarac_ttc[np], Simul->pchyd, dt, t, np,tempe1,Osat,Simul->poutputs);     
	      //On met la solution de gc dans le tableau var a linstant t qui servira d initialisation pour literation temporelle suivante
	      PROSE_update_phy_ctot(Simul->psimul_bio[np], nele, np,Simul->poutputs);
	      //}
	      //debug
	      //SW 17/10/2018 add mass balance, only mb in total domaine validate
	      if(Simul->calc_mode[MB_BIO] == YES_TS)
		PROSE_calc_mb_all_species(Simul->pcarac_ttc[np], nspecies, t, dt, np, H_flux_ttc, Simul->poutputs); // SW 05/05/2021 add H_flux_ttc for heat
	    }
	  Simul->clock->time_spent[SOLVE_TTC] += CHR_end_timer();
  
	}
      //LP_printf(Simul->poutputs,"t6 = %f days\n",t/3600/24);

      if(Simul->calc_mode[RIVE] == YES_TS)
	{
	  //fprintf(Simul->poutputs,"yes1 \n");
          //t -= dt; // SW 24/05/2018 dans libhyd t += dt, mais pour crive c'est apres
          CHR_begin_timer();
#ifdef CDA
#ifdef OMP
          num_threads = Simul->num_threads_par;
          taille = PC_set_chunk_size_silent(Simul->poutputs,nparticules,num_threads);
          omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic,taille) shared(nele,nspecies,dt,t) private(np)
#endif
#endif		  
          for(np = 0; np < nparticules; np++)  
	    {	
	      //if(np == 0)
	      //{	  
	      PROSE_update_conc_bio_all_species(Simul->pcarac_ttc[np]->p_species, nspecies, nele, Simul->psmp->nthreads,np,Simul->poutputs);
	      PROSE_link_hyd_rive_hydro(Simul->psimul_bio[np], t, Simul->pchyd, Simul->psmp->nthreads,Simul->poutputs);
	      PROSE_link_hyd_rive_geom(Simul->psimul_bio[np], Simul->pchyd, Simul->psmp->nthreads,Simul->poutputs);
	      PROSE_calc_volume_water_all_sections(Simul->psimul_bio[np], nele, t,Simul->psmp->nthreads, Simul->poutputs);
	      PROSE_calc_interface_water_sed(Simul->psimul_bio[np], Simul->pchyd, Simul->psmp->nthreads,Simul->poutputs);
              //LP_printf(Simul->poutputs,"t7 = %f days\n",t/3600/24);
		  PROSE_all_sections(Simul->psimul_bio[np], t, dt, nele, np, Simul->poutputs); // SW 04/12/2019
	      //t += dt;
              //LP_printf(Simul->poutputs,"t8 = %f days\n",t/3600/24);
	      if(Simul->calc_mode[MB_BIO] == YES_TS)
		{		
		  
		  //for(ns2 = 0; ns2 < nele; ns2++)
		    //{
		      //s_simul *psimulbio2;
		      //psimulbio = psimul_bio[][ns2];
		      PROSE_calc_mb_sections(t, Simul->psimul_bio[np], nele, np, Simul->poutputs);
		    //}
		  /*calulation of final mass for all sections and the sum for */	
		  PROSE_calc_mb_mend_all_sections(Simul->psimul_bio[np], t, dt, nele,np, Simul->poutputs);
		}		  
	      //CHR_begin_timer();
	      //fprintf(Simul->poutputs,"yes1 np = %d\n",np);
              //LP_printf(Simul->poutputs,"t9 = %f days\n",t/3600/24);
	      if(Simul->calc_mode[DA] == NO_TS)
		{
		  //PROSE_print_outputs_bio(t,Simul->psimul_bio[np],Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos,np,Simul->poutputs);
		  Simul->clock->time_spent[SOLVE_RIVE] += CHR_end_timer();
		  CHR_begin_timer();
		  PROSE_print_outputs_bio(t,Simul->psimul_bio[np],Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos,np,Simul->poutputs); // SW 20/05/2019
		  Simul->clock->time_spent[OUTPUTS_CHR] += CHR_end_timer();
		}
	      //fprintf(Simul->poutputs,"yes2 np = %d\n",np);
              //LP_printf(Simul->poutputs,"t10 = %f days\n",t/3600/24);
	      //Simul->clock->time_spent[OUTPUTS_CHR] += CHR_end_timer();
	      if(Simul->calc_mode[MB_BIO] == YES_TS && Simul->calc_mode[DA] == NO_TS)
	      {
		   CHR_begin_timer();
		   PROSE_print_mb_bio_domaine(Simul->psimul_bio[np], t, dt, nele, np, Simul->poutputs);
		   Simul->clock->time_spent[OUTPUTS_CHR] += CHR_end_timer();
	      }
	    }
	    //if(Simul->calc_mode[DA] == YES_TS)
			//Simul->clock->time_spent[SOLVE_RIVE] += CHR_end_timer();
	}
      //LP_printf(Simul->poutputs,"t7 = %f days\n",t/3600/24);
      /**/
      if(Simul->calc_mode[DA] == YES_TS)
	{
          //Particle filter SW 16/11/2021 to use a function
         //LP_printf(Simul->poutputs, "enkf md = %d\n",Simul->passim->method);
        if(Simul->passim->method == PF_PROSE)
        {
	  // extraction des simuls
	  if(fabs(t_da - t) < EPS_TS) // SW 04/04/2024 data assimilation time step
	      answer_obs = Prose_calc_difference_obs_simul(Simul->passim, Simul->psimul_bio, nparticules, t, Simul->poutputs);
	  // check if we have observations at time t
	  
	  if(answer_obs == YES_TS)
	    {
	      
	      // calc weights
	      Prose_calc_weight_all_particules(Simul->passim, Simul->poutputs);
	      // calc normalized weights and effective sample size
	      Neff = Prose_calc_normalized_weights_and_sample_size(Simul->passim, Simul->poutputs); 
	      //in order to calculate variance of all particules for each parameter
	      
              PROSE_print_resampling_size(t, Neff,Simul->passim->pout_sampling_size);
	      // may be after resampling SW 14/11/2019
	      //print param	
#ifdef OMP
              //num_threads = Simul->num_threads_par;
	      num_threads = NPARAMDA;
              taille = PC_set_chunk_size_silent(Simul->poutputs,NPARAMDA,num_threads);
              omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic,taille) shared(nparticules) private(nparm)
#endif
	      for(nparm = 0; nparm < NPARAMDA; nparm++)
	      {
                  // SW 25/01/2022 check if parameter is assimilated
                  if(Simul->passim->param_yesOrno[nparm] == YES_TS)
		      PROSE_print_extract_parameters(Simul->psimul_bio, nparticules, nparm, t, Simul->poutputs);
	      }
	      //determine trophic state // SW 22/11/2018 before resampling
              //Prose_determine_trophic_state(Simul->psimul_bio, Simul->passim, Simul->poutputs);
	      
	      //print weights
	      PROSE_print_weights(t, nparticules, Simul->poutputs);			  
	      
              //#ifdef OMP
              //num_threads = Simul->num_threads_par;
	      //num_threads = NPARAMDA;
              //taille = PC_set_chunk_size_silent(Simul->poutputs,NPARAMDA,num_threads);
              //omp_set_num_threads(num_threads);
              //#pragma omp parallel for schedule(dynamic,taille) shared(nparticules,t) private(nparm)
              //#endif
	      //for(nparm = 0; nparm < NPARAMDA; nparm++)
	      //PROSE_print_parameters(Simul->psimul_bio, nparticules, nparm, t, Simul->poutputs);
	      //openmp
              //#ifdef OMP
              //num_threads = Simul->num_threads_par;
              //taille = PC_set_chunk_size_silent(Simul->poutputs,nparticules,num_threads);
              //omp_set_num_threads(num_threads);
              //#pragma omp parallel for schedule(dynamic,taille) shared(nparticules,t,dt,nele) private(np)
	      // #endif
	      
	      // SW 03/05/2019 print concentrations and MB_BIO before resample
	      for(np = 0; np < nparticules; np++)
		{
		  //print concentrations
		  PROSE_print_outputs_bio(t,Simul->psimul_bio[np],Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos,np,Simul->poutputs);
		  if(Simul->calc_mode[MB_BIO] == YES_TS)
		    //print mass_balance
		    PROSE_print_mb_bio_domaine(Simul->psimul_bio[np], t, dt, nele, np, Simul->poutputs);				  
		}			 
	      for(np = 0; np < nparticules; np++)
		Prose_weights_to_weights_prev(Simul->passim, np, Simul->poutputs);			 		  
	      
	      if(Neff < (Simul->passim->alpha*Simul->passim->N_particules))
		{
		  //LP_printf(Simul->poutputs,"begin resampling neff = %d\n",Neff);
		  Prose_fill_re_sample_elim_dupli(Simul->passim, Simul->poutputs);
		  // resample
		  //LP_printf(Simul->poutputs,"after fill resample\n");
		  Prose_resample_particules(Simul->passim, Simul->poutputs);
		  Prose_reinitialization_weights(Simul->passim, Simul->poutputs);
		  Prose_perturbation_parameters(Simul->psimul_bio, Simul->passim, nele, BLOOM, Simul->poutputs); // SW 31/05/2019 all parameters

		  // SW 20/05/2022 recalculate the particle weights after perturbation of parameters
                  PROSE_new_weights_after_perturbation(Simul->passim, Simul->poutputs);

		}
	      // SW 23/11/2018 after calculation of parameter variance
	      //for(nparm = 0; nparm < NPARAMDA; nparm++)
	      //PROSE_extraction_parameters(Simul->psimul_bio, nparticules, nparm, Simul->poutputs);
	      
	      //output
	      
	      //if(Simul->passim->state == BLOOM)
	      //LP_printf(Simul->poutputs,"algal blooms\n");
	      //perturbation
	      //Prose_perturbation_parameters(Simul->psimul_bio, Simul->passim, nele, Simul->passim->state, Simul->poutputs);
	      
	      PROSE_reset_no_answer_obs(Simul->passim,Simul->poutputs);
	    }
	  
	  else // SW 05/06/2019
	    {
              PROSE_print_resampling_size(t, Neff,Simul->passim->pout_sampling_size);
#ifdef OMP
              //num_threads = Simul->num_threads_par;
	      num_threads = NPARAMDA;
              taille = PC_set_chunk_size_silent(Simul->poutputs,NPARAMDA,num_threads);
              omp_set_num_threads(num_threads);
#pragma omp parallel for schedule(dynamic,taille) shared(nparticules) private(nparm)
#endif
	      for(nparm = 0; nparm < NPARAMDA; nparm++)
	      {
                  // SW 25/01/2022 check if parameter is assimilated
                  if(Simul->passim->param_yesOrno[nparm] == YES_TS)
		      PROSE_print_extract_parameters(Simul->psimul_bio, nparticules, nparm, t, Simul->poutputs);
	      }
	      //determine trophic state // SW 22/11/2018 before resampling
              //Prose_determine_trophic_state(Simul->psimul_bio, Simul->passim, Simul->poutputs);
	      
	      //print weights
	      PROSE_print_weights(t, nparticules, Simul->poutputs);			  
	      
	      // SW 03/05/2019 print concentrations and MB_BIO at the end
	      for(np = 0; np < nparticules; np++)
		{
		  //print concentrations
		  PROSE_print_outputs_bio(t,Simul->psimul_bio[np],Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos,np,Simul->poutputs);
		  if(Simul->calc_mode[MB_BIO] == YES_TS)
		    //print mass_balance
		    PROSE_print_mb_bio_domaine(Simul->psimul_bio[np], t, dt, nele, np, Simul->poutputs);				  
		}
	      for(np = 0; np < nparticules; np++)
		Prose_weights_to_weights_prev(Simul->passim, np, Simul->poutputs);			 		  
	    }
		Simul->clock->time_spent[SOLVE_RIVE] += CHR_end_timer();
	}
	else // SW 16/11/2021
    {
      Prose_processus_enkf(Simul->psimul_bio, nele, Simul->passim, nparticules, t, t_da, Simul->poutputs);
        //LP_printf(Simul->poutputs, "enkf t = %f\n", t);
    }
    }
    

    PROSE_update_outputs_hyd_bio_heat_cur_io(t,Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos, 0, Simul->poutputs); // SW 26/01/2021  update t_out[CUR_IO]
    } //end while t	
    // masse balance bio a faire
    //if(Simul->calc_mode[TTC] == YES_TS)
    //fclose(fpno3);

  }
  
  CHR_begin_timer();//LV 3/09/2012
  pmb=HYD_print_mass_balance(pmb,t,Simul->pinout->fmb,Simul->chronos,Simul->poutputs);//NF 11/9/06
  fclose(Simul->pinout->fmb);//NF 11/9/06

  if (Simul->pinout->calc[MB_ELE] == YES_TS)//LV 15/06/2012
    HYD_print_mb_at_elements(t,Simul->outputs,Simul->pchyd,Simul->chronos,Simul->poutputs);
  
  if (Simul->pinout->calc[FINAL_STATE] == YES_TS)
    {
      if(Simul->calc_mode[RIVE] ==YES_TS)
	PROSE_print_conc_volume_final(nele,Simul->pinout, Simul->poutputs);
      HYD_print_final_state(Simul->outputs,Simul->pinout,Simul->pchyd,Simul->poutputs);
      if(Simul->calc_mode[H_T] == YES_TS) // AB 04/06/20
	PROSE_print_final_temp(Simul->outputs, Simul->pinout, Simul->pchyd, Simul->poutputs);
    }
  
  Simul->clock->time_spent[OUTPUTS_CHR] += CHR_end_timer();//LV 3/09/2012
  
  printf("date of end of simulation : %d %s %d %f h\n",Simul->chronos->day_d,TS_name_month(Simul->chronos->month,Simul->poutputs),Simul->chronos->year[BEGINNING],Simul->chronos->hour_h);

  /********************************************************/
  /*        End, summary of the computation length        */
  /********************************************************/
  
  CHR_calculate_calculation_length(Simul->clock,Simul->poutputs);
  printf("Time of calculation : %f s\n",Simul->clock->time_length);
  #ifdef CHR_CHECK
  PROSE_calculate_check_time(nele, Simul->poutputs);
  #endif 
  
  CHR_print_times(Simul->poutputs,Simul->clock);
  
  //LP_printf(Simul->poutputs,"data read ! dt = %f tend = %f\n",dt,tend);
  
  fclose(Simul->poutputs);

  return 0;
}


