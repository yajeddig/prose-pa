/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: struct_simul_PROSE.h
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

typedef struct simulation_PROSE s_simul_PROSE;
typedef struct carac_assimilation s_carac_assim;
typedef struct carac_obs_assimilation s_carac_obs_assim;
typedef struct weight_calc_for_assim s_weight_calc_assim;
/* Simulation structure */
struct simulation_PROSE {
  /* Name of the simulation */
  char *name;
  /*regime of the simulation*/
  int regime;
  int calc_mode[NMODE];
  /* parallel caculation of n particules*/
  int num_threads_par;
  /*omp struct libpc.h*/
  s_smp *psmp;
  int solver;
  /* Time characteristics */
  s_chronos_CHR *chronos;
  /* Timer */
  s_clock_CHR *clock;
  /* Date format in input data*/
  int date_format;
  /*structur of hydraulic characteristics*/
  s_chyd *pchyd;
  /*structur of transport N particule*/
  s_carac_ttc **pcarac_ttc;
  double ***RHS_b;
  char ***mat_adv;
  double ***RHS_b_phy;
  char ***mat_adv_phy;
  double *RHS_b_T;
  char *mat_adv_T;
  /*structur of species of PHY, N particules*/
  s_species_ttc ****p_phy_species;
  /*structur of bio N particules*/
  s_simul ***psimul_bio;
  /*structur counter for bio*/
  s_counter *counter_bio;
  /*setting biology*/
  s_settings *settings;
   /* Table containing the pointers towards the macrospecies that are used for data assimilation of boundary conditions such as TOC  */ // MH: 10/09/2021
  //s_macrospecies_RIV **p_macro_da[NMACROSPECIES]; // in case of [TOC][np]
  s_macrospecies_RIV *p_macrospecies[NMACROSPECIES]; 
  
  

  //NF 12/10/2020 It is required to create a meta structur for the energy balance, and then include all of those forcings in it
  /* Structure for the heat transport */
  s_carac_ttc *pcarac_heat_ttc;
  /* For temperature initialization */
  s_ft *init_T;
  int init_from_file;
  double default_t_inflows; // SW 10/03/2021
  
  // SW 04/02/2021 add a meta structur for libseb
  //s_meteodata_seb *pmeteodata_seb;

  /* Corresp id_meteo to id_reach */
  s_rts **p_rts; // links reaches to safran's cell
  /* Safran data */
  s_met **p_met; // Size corresponds to the number of the different meteo data sets being used -> usually, number of safran's cells
  s_input_seb **p_safran;
  long n_meteocell;
  /* Time delineation for gathering the right time frame of meteo data */
  double *t_extrema; // This variable has been added, and is important!
  int year0_meteo; //SW 04/02/2021

  s_mbheat_mb **mb_heat; // SW 03/05/2021
  
  /* Structure for the tubes */
  s_def_tub ****p3_tube;
  int ntube_default;
  s_output_tube *pout_tube;

  /* Debug file */
  FILE *poutputs;
  FILE *fphy;
  /*output structur for hydraulics for N simulations */
  s_output_hyd **outputs[NOUTPUTS];
  s_total_mb ***total_mb; // mass balance for biology
  s_lp_pk_hyd **lp_pk;  // lp_pk for biology mass balance
  int npk_mb_bio;
  //s_conc **total_conc;
  /*output structur for biology Output files (mass balances) for N simulations*/
  FILE ****mass_balances[NLAYERS+1][NSPECIES+NANNEX_VAR];
  FILE ****ads_mass_balances[NLAYERS+1][NDISS];
  FILE ****concentrations[NLAYERS+1];
  FILE ****ads_concentrations[NLAYERS+1][NDISS];  
  int calc_bio_annexvar[NANNEX_VAR]; // SW 04/06/2018
  double unit_bio_annexvar[NANNEX_VAR]; // SW 04/06/2018
  FILE ***transv_profil_temp;//NF 12/10/2020 pointers to transv profiles (time series) of temperature. For now only for the water compartment and independent of the particule filter
  /*settings for input and outputs */
  s_inout_set_io *pinout;
  s_lec_tmp_hyd **plec;
  /*Simulation output */ // SW 24/12/2018
  s_out_io *pout;
  
  /*struct of data assimilation*/
  s_carac_assim *passim;
   
};

/*Structure for configuring the calculation of the weight related to an assimilation station. For now the parameters are related to oxygen data assimilation related to organic matter content*/ 
struct weight_calc_for_assim {
  /*    Potentially user defined parameters for weight calculation*/
  int weighted_stations; // yes or no to using weighteds stations in data assimilation
  double vel_low_flow ;
  double bact_growth_rate ;
  double bact_yield ;
  double fast_bdom_low_flow ;
  double bact_biomass;
  double OM_flux_threshold;
  double OM_decay_rate;

  int N_flows_significant; // total number of inflows including the upstream in the system that passes the threshold
  double weight_max;
  
};


/*struct of data assimilation*/
struct carac_assimilation {
	/*the number of particules*/
        int method; // SW 16/11/2011 PF method and EnKF method
	int N_particules;
	int N_obs;
	int Neff;
	int num_t_obs;
	int state;
	int nstep;
	double alpha;
	double seuil_chla;
        int random_walk;

        /* SW 20/05/2022 add perturbation density to recalculate the particle weights */
        double *param_perturbation_density[NPARAMDA];
        double *param_perturbation_density_product;

	/*parameters of all particules at t - 1*/
	double *param[NPARAMDA];
        /* Corresponding gaussian parameter values which follow a gaussian distribution N(0, 1) SW 13/01/2022 */
        double *paramGaussian[NPARAMDA];
	double units_param[NPARAMDA];
	double param_range[NPARAMDA][PARAMR]; // UP and DWON
        // SW 24/01/2022 add param_yesOrno for assimilable parameters' implementation
        int param_yesOrno[NPARAMDA];
        int num_Of_assimilated_param;

	/*the number of noeud used for MPI*/
	//int N_noeud;
	/*the number of particules sur un noeud N_particules/N_noeud*/
	//int nparticules;
	/* the error percent associated for observation defined by user for exemple 5% */
	double error_obs_sigma;
	/* check if perturb the obeservations */
	//int obs_perturb_answer;
	/* s_percent, percent of parameter range for perturbation*/ 
  double s_percent[NPARAMDA];  // MH 10/03/2022 : [NPARAMDA] added to facilitate different random walk for each param
	/*weights of particules, length = N_particules*/
	double *omega_prev;
	double *omega;
	/*normlized weights of particules, length = N_particules*/
	double *omega_normlized;
	/*Cumulative distribution function, length = N_particules, used for main thread*/
	double *cdf;
	/*ensemble weights of all particules for main thread length = N_particules*/
	//double *omega_all;
	/*rank of processeurs, length = N_particules*/
	//int *rank;
	/*for re-sample, check if the particule is saved, 1 = eliminated 0 = saved*/
	int *eliminated;
	/*the number of duplications of a particule */
	int *ndupli;

         /* Table containing the pointers towards the macrospecies that are used for data assimilation of boundary conditions such as TOC  */ // MH: 10/09/2021
        s_macrospecies_RIV **p_macro_da[NMACROSPECIES]; // in case of [TOC][np]

     /* Moving average parameters*/ // MH: 26/04/2022
        int da_mov_ave;  //yes or no if we want to use moving average in our data assimilation
        double mov_ave_theta; // moving average alignment: centered, left or right
        double lneigh; //length of neighborhood window or the jump for moving average
        int do_time_step; // oxygen data frequency
        int nval_min; // minimum number of values to be used for moving average at the two ends in case of using the centered method

  /* Weighted stations for data assimilation  */  // MH: 03/05/2022
  s_weight_calc_assim *weight_calc_assim;
  

    /*observations in M stations*/
	s_carac_obs_assim **pobs;
	/*simulation for M stations, N particules*/
	//s_ft ***sim;
	FILE *pout_weight;
	FILE **pout_param;
	FILE *pout_sampling_size;
};


struct obs_station_weight {
  int  nupstream;
  int nsources;
  //s_inflow **p_inflows;
  double *pk_sources;
  double *mass;
  double *dist;
};


struct carac_obs_assimilation {

	/*observations one stations*/
    s_ft *obs;	
	/*simulation of N particules*/
	//s_ft **y;
	double *difference;
        /* SW 19/01/2022 store epsilon Y_true = Y^* + epsilon for EnKF*/
        double *epsilon;
	/*observation at time t for one station*/
	double Obs;
	/*check if we have the observation at current time t YES_TS or NO_TS*/
	int answer_obs;
	int var;
	int num;
	double chla_mean;
	/* time step for assimilation dt_da = nstep*dt*/
	int nstep;
	int pk_type;
	/* grid identity of one station*/
	int id_ele_obs;
    char *river;
    int reach_nb;
    int branch_nb;
    double pk;
  /*weighted stations */
  double weight_station; //normalized
    double weight_i;

  
	s_carac_obs_assim *prev;
	s_carac_obs_assim *next;
	
};

#define new_simulation_PROSE() ((s_simul_PROSE *) malloc(sizeof(s_simul_PROSE)))
#define new_assimilation() ((s_carac_assim *) malloc(sizeof(s_carac_assim)))
#define new_obs_assimilation() ((s_carac_obs_assim *) malloc(sizeof(s_carac_obs_assim)))
#define new_weight_calc_assim() ((s_weight_calc_assim *) malloc(sizeof(s_weight_calc_assim)))
