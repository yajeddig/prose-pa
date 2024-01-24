/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: input.y
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

%{



#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <errno.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <signal.h>
#include <time.h>
#include <sys/time.h>
#include <libprint.h>
#include <time_series.h>
#include "libpc.h"
#include "IO.h"
#include "GC.h"
#include "CHR.h"
#include "TTC.h"
#include "SEB.h"
#include "HYD.h"
#include "RIVE.h"
#include "TUB.h"
#include "MB.h"
#include "LA.h"
#include "PROSE.h"
//#include "global_PROSE.h"
#include "ext_PROSE.h"


  
  /* Previous time read in a time series */
  double told = 0.;
  /* Time */
  double tnew = 0.;
  /* Unit of the time or of ft->t in a s_ft structure */
  double unit_t = 1.;
  /* Unit of ft->f in a s_ft structure */
  double unit_f = 1.;
  double unit_x,unit_y,unit_z,unit_l,unit_kappa;
  /* Pointer towards the currently read reach and the total chain of reaches */
  s_reach_hyd *preach, *preachtot;
  /* Reach in which the bathymetry is currently being defined */
  s_reach_hyd *preach_bat;
  /* Reach in which the inflow arrives */
  s_reach_hyd *preach_infl;
  /* Pointer towards a transversal face */
  s_face_hyd *pface;
  /* Pointer towards a calculation element */
  s_element_hyd *pele;
  /* Pointer towards an inflow */
  s_inflow_hyd *pinflow;
  /* Pointer towards the currently read singularity and the total chain of singularities */
  s_singularity_hyd *psing, *psingtot;
  /* Current BC characteristics */
  s_BC_char_hyd *pwork, *pworktot;
  /* Nb of hydraulic works constituting the BC */
  int nworks;
  /* Type of hydraulic work*/
  int worktype;
  /* Geometry of the face described in Lambert coordinates */
  s_pointAbscZ_hyd *pptAbscZ;
  /* Geometry of the face described in XYZ coordinates */
  s_pointXYZ_hyd *pptXYZ;
  /* Number of steps in a time series s_ft */
  int nbsteps;
  /* Loop indexes */
  int s,r,i,p,f,e;
  /* Outputs */
  s_output_hyd *pout_hyd;
  s_output_hyd *pts=NULL,*plp=NULL,*pmb=NULL;
  s_ts_pk_hyd *pts_pk,*pts_pktot;
  s_carac_obs_assim *pts_obs_pktot, *pts_pk_obs;
  s_lp_pk_hyd *plp_pk,*plp_pktot,*pmb_pk,*pmb_pktot, *plp_pk_tube, *plp_pk_tube_tot, *pmb_pk_heat;
  int output_type;
  int nout_tube = 0;
  s_output_tube_type *tb_out;
  char *river;
  char**FnodeTnode;
  int nb_outlet;
  int nb_reach=0;
  int dim_tmp_lec = NB_LEC_MAX;
  int nb_realloc = 1;
  s_date_ts *pd;

  /*macrospecies*/  //  09/09/2021
  
  /* Pointer towards a macrospecies*/
  s_macrospecies_RIV *pmacspe;//NF 21/7/2021
  s_stoch_param_RIV *pstoch; // MH 19/08/2021
  int rev_threshold = NO_TS;//NF 19/8/2021 for macrospecies TOC
  int imacspe, idegval;//NF 19/8/2021 for macrospecies TOC
  double t_param, b1_param, b2_param, s1_param, s2_param; // MH 25/08/2021 for calculating the share of MOD1 .... MOP3 
  double sharemod1, sharemod2, sharemod3, sharemop1, sharemop2, sharemop3;
  s_c_flows *pconstflux;

  s_weight_calc_assim *pstat_weight; //16/05/2022 struct for station weight parameters

  
  //moving average
  s_ft *pft_obs;
  s_ft *pft_mov_ave;
  double t0,tf;
  double lneigh,theta_mave;
  int nval_min;

  /*libttc*/
  
  /* Loop indexes */
  int s,r,i,p,f,e;
  s_ft *pft;
  int nbsteps;
  /*Variables for the transport*/
  s_species_ttc *pspec_ttc;
  s_species_ttc **pspecies[NSPECIES];
  s_link_ttc *plink;
  int nspecies;
  int itype_t=1;
  int idir_t=1;
  /*Variables for the mesh*/
  int nlayer;
  int nele;
  int nele_tot;
  int intern_id;
  int coordo;
  int verif_layer=0;
  int verif_specie=0;
  int verif_ele=0;
  int kindof_transm,kindof_BC,kindof_SOURCE;
  int  kindof_attribut,kindof;
  int idir=1;
  double surf_riv=0;
  double ep_riv_bed = 0;
  double TZ_riv_bed = 0;
  double temp_var;
  double stot;
  double ltot;
  int nelebound = 0;
  int general_out = CODE;
  
  /*crive*/
  
  /* Last and new times in ft chains */  
  double told;
  /* Time unit variable */
  double unit_t;
  /* Variable unit variable */
  double unit_f;
  /* Number of the subspecies */
  int num;
  /* Number of mass balances asked by the user */
  int num_mb = 0;
  int num_conc = 0;
  /* Number of time steps in a time series */
  int nbsteps;
  /* Kind of species */
  int var,class;
  /* Type of layer */
  int lay,layspe,layfl;
  /* Number of the sublayer when defining flows at interfaces */
  int nlfl;
  /* Table of integers defining for each layer whether or not its species were defined */
  int def_spe[NLAYERS];
  /* Type of the 1st layer in which species were defined */
  int type_1st_layer = CODE;
  /* Type of flow */
  int fl;
  /* Loop index */
  int i,e,j,nl,k;
  /* Biogeochemical process */
  int proc;
  /* Intermediary species variables */
  s_species *spec2;
  s_species **pspec2[NLAYERS][NSPECIES];
  /* Intermediary section variable */
  s_section *sect2;
  /* Intermediary output variables */
  s_total_mb *mb2;
  s_conc *conc2;
  /* Intermediary compartment */
  s_compartment *comp2;
  /* Intermediary interface */
  s_interface *int2;
  /* Intermediary reaction variables */
  s_reac *reaction2;
  char **other_reactors_old;
  char **other_reactors_new;
  double *stoechio_other_reactors_old;
  double *stoechio_other_reactors_new;
  char **products_old;
  char **products_new;
  double *stoechio_products_old;
  double *stoechio_products_new;
  int nother_reactors;
  int nproducts;
  int nsections;
  int ns;

  /*libseb*/

  int i;
  int check = 0;
  int nout_mb_heat = 0;
  
  s_carac_ttc *pcarac_heat_ttc;
  double theta_T;

  //libmb
  s_mbheat_mb *mb_heat_temp;
%}


%union {
  double real;
  int integer;
  char *string;
  s_ft *function;
  s_date_ts *date_ts;
  s_singularity_hyd *sing;
  s_reach_hyd *reach;
  s_pointXYZ_hyd *pointXYZ;
  s_pointAbscZ_hyd *pointAbscZ;
  s_face_hyd *face;
  s_inflow_hyd *inflow;
  s_BC_char_hyd *boundary;
  s_ts_pk_hyd *ts_pk;
  s_lp_pk_hyd *lp_pk;
  
  
  //libttc
  s_species_ttc *ttc_species;
  // int struct_id;
  s_id_io *struct_id;
  s_link_ttc *plink;
  
  //crive
  s_species *species;
  s_total_mb *mass_balance;
  s_conc *conc_output;
  s_compartment *comp;
  s_interface *interface;
  s_reac *reaction;

  //macrospecies
  s_macrospecies_RIV *macrospecies;//NF 21/7/2021 introducing the reading of macrospecies for DA of boundary conditions
  s_stoch_param_RIV *stoch_param;//NF 19/8/2021 introducing the reading of stochastic param

  
  //DA
  s_carac_obs_assim *ts_pk_obs;
  //s_weight_calc_assim *pstat_weight; // MH 16/05/2022: introducing the reading of station weight calc params
  

 //libtube
 s_output_tube_type *out_tube;

 //libmb
 s_mbheat_mb *mb_heat;
}

%token <real> LEX_DOUBLE LEX_A_UNIT LEX_VAR_UNIT LEX_A_UNIT_MOL
%token <integer> LEX_INT
%token LEX_INIT LEX_INIT_Z LEX_INIT_Q
%token LEX_OPENING_BRACE LEX_CLOSING_BRACE LEX_OPENING_BRACKET LEX_CLOSING_BRACKET LEX_REVERSE_ARROW
%token LEX_HARD_BRACKET_OPEN LEX_HARD_BRACKET_CLOSE
%token LEX_INV LEX_POW LEX_COMA LEX_EQUAL LEX_ARROW LEX_COLON LEX_SEMI_COLON
%token <string> LEX_NAME LEX_DATE_DAY LEX_DATE_HH_MM LEX_DATE_HH_MM_SS
%token LEX_INPUT_FOLDER LEX_OUTPUT_FOLDER LEX_FOLDER
%token LEX_COMM LEX_LOG LEX_SEMICOLON
%token LEX_SIMUL LEX_TIMESTEP LEX_TIMESTEP_DA LEX_CHRONOS LEX_STRICKLER LEX_BEGIN LEX_DIVDT LEX_DATEFORMAT 
%token LEX_NUM LEX_NUM_THREAD LEX_NUM_THREAD_PAR LEX_NUM_PARTICULES LEX_NUM_TUBE LEX_ALPHA_DA LEX_SEUIL_CHLA 
%token LEX_PARAM_RANGE LEX_ERROR_OBS LEX_S_PERCENT LEX_LIMIT_FACTOR LEX_RANDOM_WALK LEX_DA_METHOD
%token <integer> LEX_PHY2  LEX_LIMIT_FACTOR_VAL LEX_RANDOM_WALK_VAL LEX_DA_METHOD_VAL
%token LEX_SPECIES LEX_ANNEXVAR LEX_ADS_SPECIES LEX_EXCH
%token <integer> LEX_ONESPECIES LEX_TEMPSPECIES LEX_COMP LEX_ONEANNEXVAR LEX_PARAM_DA LEX_SED_VAR
%token <integer> LEX_PROCESS
%token <integer> LEX_PARAMP LEX_PARAML LEX_PARAMM LEX_PARAMG LEX_PARAMD
%token <integer> LEX_PHOT LEX_GROWTH LEX_MORT LEX_RESP LEX_GRAZ LEX_EXCR 
%token <integer> LEX_REA LEX_ADS_SENEQUE LEX_ADS_FR LEX_HYDR LEX_RAD_DECAY 
%token LEX_KMICH LEX_KLIM LEX_NUTC LEX_ADS_ON LEX_REACTIONS LEX_REACTION LEX_DEFAULT_C_INFLOWS
%token LEX_OTHERS LEX_PRODUCTS LEX_CONDITION
%token <integer> LEX_INF_SUP
%token LEX_SECTION LEX_GEOMETRY LEX_ELEVATION LEX_WIDTH LEX_LENGTH
%token LEX_METEO LEX_TEMPERATURE LEX_WIND LEX_RADIATION LEX_PHOTOPERIOD
%token LEX_MEAN LEX_AMPLITUDE LEX_DELAY LEX_ATTENUATION
%token LEX_HYD LEX_TUNIT LEX_PERIMETER LEX_SURFACE
%token LEX_DISCHARGE LEX_VELOCITY LEX_HYDRAD
%token LEX_HEIGHT LEX_SCOURING LEX_RET
%token <integer> LEX_LAYERS LEX_LAYER
%token LEX_MASS LEX_VOLUME LEX_PHI LEX_RHO LEX_IC LEX_FC LEX_SAT LEX_RATIO_WATER LEX_SURFLAYER
%token <integer> LEX_FLOWS
%token LEX_INTERFACES

%token <integer> LEX_MONTH LEX_YEAR0 LEX_YEAR0_METEO
%token <integer> LEX_ANSWER LEX_TIME LEX_METHOD LEX_SETYPE LEX_PARAM_EROS LEX_TCTYPE LEX_TC_PARAM LEX_DATEFORMAT_FR_TS 
%token LEX_SET LEX_EPS LEX_COMPPHY LEX_DZ_RIVE LEX_DBO LEX_CALC_SE LEX_CALC_TC
%token LEX_DIM LEX_TYPE LEX_CALC_CURVATURE 
%token <integer> LEX_CALC_STATE LEX_GENERALP LEX_CALC_MODE
%token LEX_NETWORK LEX_PK LEX_HYDWORKS LEX_POSITION LEX_HOLLER
%token <integer> LEX_FION LEX_DIR
%token LEX_PROFILES LEX_INFLOWS LEX_SING LEX_REACH LEX_TRAJECTORY
%token <integer> LEX_INFLOW_TYPE
%token LEX_ABSC LEX_X LEX_Y LEX_HYD_VAR LEX_SHAPE
%token <integer> LEX_BC_TYPE LEX_CROSS_SECTION_TYPE
%token <integer> LEX_WORK LEX_WORK_PARAM
%token LEX_FILE_NAME LEX_POINTS LEX_EXTENT LEX_POINTS_OBS
%token LEX_OUTPUTS LEX_VAR LEX_GRAPHICS
%token LEX_MB LEX_CONC LEX_NSTEPS
%token <integer> LEX_ONE_VAR LEX_OUTPUT_TYPE LEX_GRAPH_TYPE
%token <integer> LEX_MUSK_PAR LEX_CALC_MUSK_PAR LEX_DEF_SCHEM_TYPE
%token LEX_AREA LEX_SLOPE LEX_SCHEM_TYPE LEX_NETWORK_MUSK LEX_REACH_MUSK LEX_ELE_MUSK

%token <integer>  LEX_FORMAT LEX_PREF
%token LEX_NEWSAM  LEX_AQUIFER LEX_AVIEW LEX_SETUP_LAYER LEX_SET_UP 
%token LEX_SETUP_SPECIES
%token LEX_BOUND LEX_SOURCE LEX_TIME_S LEX_PARAM LEX_PARAM_RIV LEX_COND LEX_INI_STATE
%token <integer> LEX_MESH_ATT LEX_TYPE_BC LEX_TYPE_MEAN LEX_TYPE_SOURCE LEX_SIM  LEX_TYPE_COND
%token <integer> LEX_INI_FLOW_TYPE LEX_INI_TRANS_TYPE
%token LEX_SETTINGS_TTC LEX_TRANSPORT_EQ LEX_SETUP_TRANSPORT LEX_INI_FLOW LEX_INI_TRANS
%token <integer> LEX_TRANSPORT_ATT LEX_TRANSPORT_BASE_ATT LEX_TRANSPORT_U_ATT
%token <integer> LEX_INIT_VAR_ATT
%token <integer> LEX_BOUND_TYPE LEX_BOUND_DIR
%token <integer> LEX_KINDOF_BOUND
%token LEX_HYDRO LEX_BIOLOGY LEX_TRANSPORT LEX_REGIME LEX_MED_TYPE LEX_THETA LEX_POROSITY
%token LEX_DIFF_MOL LEX_SOLVER
%token LEX_WATER
%token <integer> LEX_MEDIA_TYPE LEX_PROC_TRANSP LEX_TRANSPORT_MODE LEX_SOLVER_TYPE

%token LEX_HEAT_TRANSPORT LEX_INIT_T LEX_ID LEX_TEMP_INFLOW LEX_DEFAULT_T_INFLOWS LEX_THETA_T
%token LEX_TUBE_OUT LEX_EB_HEAT LEX_HEAT_UNIT
%token <integer> LEX_MET LEX_CODE_TS LEX_TUBE_OUT_TYPE
%token LEX_MACROSPECIES LEX_RELATED_TO LEX_THRESHOLD LEX_VAL LEX_RANGES
%token <integer> LEX_TYPEOF_MACROSPECIES LEX_BIODEG LEX_SHARE_MO
%token LEX_DA_MOV_AVE LEX_MOV_AVE_THETA LEX_LNEIGH LEX_DO_TIME_STEP LEX_NVAL_MIN
%token LEX_WEIGHTED_STATIONS LEX_CONFIGURE_WEIGHTS LEX_VEL_LOW_FLOW LEX_BACT_YIELD LEX_BACT_GROWTH_RATE LEX_FAST_BDOM LEX_BACT_BIOMASS LEX_OM_FLUX_THRESHOLD LEX_OM_DECAY_RATE

%type <real> flottant mesure 
%type <real> a_unit units one_unit a_unit_value all_units read_units a_unit_f
%type <function> f_ts f_t date_f_ts date_f_t
%type <mass_balance> def_mbs def_one_mb
%type <species> subspecies one_subspecies
//%type <comp> sub_layers sub_layer
%type <reaction> reacs reac
%type <boundary> fion hydworks hydwork hydstructures
%type <sing> singularity singularities
%type <reach> reach reaches 
%type <pointAbscZ> sectionAbscZ sectionAbscZs 
%type <pointXYZ> sectionXYZ sectionXYZs
%type <face> cross_sections cross_section bathymetry trajectory traj_point traj_points 
%type <inflow> inflow inflows
%type <lp_pk> extents
%type <ts_pk> points one_pk pk_list
%type <ts_pk_obs> points_obs one_pk_obs pk_list_obs
%type <out_tube> def_outs_tube def_out_tube
%type <mb_heat> def_out_mb_heat def_outs_mb_heat
%type <macrospecies> macrospecs macrospec


%start beginning
%%


beginning : begin
            outputs;


begin : paths model_settings
| model_settings 
;

/* This part enables to define specific input and output folders */
paths : path paths
| path
;

path : input_folders
| output_folder
;

input_folders : LEX_INPUT_FOLDER folders
;

/* List of input folders' names */
folders : folder folders
| folder 
;

folder : LEX_EQUAL LEX_NAME
{
  if (folder_nb >= NPILE) {
    fprintf(Simul->poutputs,"Only %d folder names are available as input folders\n",NPILE);
    printf("Only %d folder names are available as input folders\n",NPILE);
  }
  name_out_folder[folder_nb++] = strdup($2);
  fprintf(Simul->poutputs,"path : %s\n",name_out_folder[folder_nb-1]);
  printf("path : %s\n",name_out_folder[folder_nb-1]);
  free($2);
  pd = new_date_ts(); // SW 05/01/2021
} 
;

/* Definition of the folder where outputs are stored */
output_folder : LEX_OUTPUT_FOLDER LEX_EQUAL LEX_NAME
{
  int noname,str_len;
  char *new_name;
  char cmd[MAXCHAR_PROSE];
  FILE *fp;

  fp=Simul->poutputs;

  Simul->pinout->name_outputs = $3;
  str_len=strlen($3)+8;
  new_name = (char *)calloc(str_len,sizeof(char));
  sprintf(new_name,"RESULT=%s",$3);
  fprintf(fp,"%s\n",new_name);
  printf("%s\n",new_name);
  noname = putenv(new_name); 
  if (noname == -1)
    LP_error(fp,"File %s, line %d : undefined variable RESULT\n",
		current_read_files[pile].name,line_nb); 
  sprintf(cmd,"%s",getenv("RESULT")); //NF 25/3/2019 mkdir recurssive to avoid issues for dummies
  //system(cmd);
  IO_mkdir(cmd); // SW 19/10/2022 mkdir recurssive to avoid issues for dummies, system may has problems sometimes

} 
;

/* This part defines global parameters of the model and stores them in the simulation structure */
model_settings : LEX_SIMUL LEX_EQUAL LEX_OPENING_BRACE simul_name atts_simul LEX_CLOSING_BRACE
{
  LP_printf(Simul->poutputs,"Global parameters of the simulation have been defined\n\n");
  
}
;

simul_name : LEX_NAME
{
  Simul->name = $1;
  LP_printf(Simul->poutputs,"\nSimulation name = %s\n",Simul->name);
  
}
;

atts_simul : att_simul atts_simul
| att_simul
; 

att_simul : def_chronos
| def_settings
;

def_chronos : LEX_CHRONOS LEX_EQUAL LEX_OPENING_BRACE def_times LEX_CLOSING_BRACE

def_times : def_time def_times
| def_time
;

def_time : LEX_DATEFORMAT LEX_EQUAL LEX_DATEFORMAT_FR_TS
{
   Simul->date_format = $3;
}
| LEX_YEAR0 LEX_EQUAL LEX_INT
{
   Simul->chronos->yr0 = $3;
}
| LEX_TIMESTEP LEX_EQUAL mesure
{
  Simul->chronos->dt = $3;
  
LP_printf(Simul->poutputs,"time step = %f s\n",Simul->chronos->dt);
  
} 
| LEX_TIME LEX_EQUAL mesure
{
  Simul->chronos->t[$1] = $3;
  LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]);
  Simul->chronos->pd[$1] = TS_convert_julian2date_decimal(TS_convert_seconds2days(Simul->chronos->t[$1],Simul->poutputs),Simul->chronos->yr0,Simul->poutputs);
  //LP_error(Simul->poutputs,"month = %d\n",Simul->chronos->pd[$1]->mm);
  Simul->chronos->pd[$1]->mm += 1;
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_DATE_HH_MM_SS
{
  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &Simul->chronos->pd[$1]->dd, &Simul->chronos->pd[$1]->mm,&Simul->chronos->pd[$1]->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &Simul->chronos->pd[$1]->mm, &Simul->chronos->pd[$1]->dd,&Simul->chronos->pd[$1]->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    sscanf($4,"%d:%d:%d", &Simul->chronos->pd[$1]->hh, &Simul->chronos->pd[$1]->min,&Simul->chronos->pd[$1]->ss);
    
    //Simul->chronos->pd[$1]->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  Simul->chronos->t[$1] = TS_date2julian_dd_hm(Simul->chronos->pd[$1],Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]); 
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_DATE_HH_MM
{
  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &Simul->chronos->pd[$1]->dd, &Simul->chronos->pd[$1]->mm,&Simul->chronos->pd[$1]->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &Simul->chronos->pd[$1]->mm, &Simul->chronos->pd[$1]->dd,&Simul->chronos->pd[$1]->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    sscanf($4,"%d:%d", &Simul->chronos->pd[$1]->hh, &Simul->chronos->pd[$1]->min);
    Simul->chronos->pd[$1]->ss = 0.;
    
    //Simul->chronos->pd[$1]->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  Simul->chronos->t[$1] = TS_date2julian_dd_hm(Simul->chronos->pd[$1],Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]); 
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_INT
{
  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &Simul->chronos->pd[$1]->dd, &Simul->chronos->pd[$1]->mm,&Simul->chronos->pd[$1]->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &Simul->chronos->pd[$1]->mm, &Simul->chronos->pd[$1]->dd,&Simul->chronos->pd[$1]->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    Simul->chronos->pd[$1]->hh = $4;
    Simul->chronos->pd[$1]->min = 0.;
    Simul->chronos->pd[$1]->ss = 0.;
    
    //Simul->chronos->pd[$1]->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  Simul->chronos->t[$1] = TS_date2julian_dd_hm(Simul->chronos->pd[$1],Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]); 
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY
{
  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &Simul->chronos->pd[$1]->dd, &Simul->chronos->pd[$1]->mm,&Simul->chronos->pd[$1]->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &Simul->chronos->pd[$1]->mm, &Simul->chronos->pd[$1]->dd,&Simul->chronos->pd[$1]->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    Simul->chronos->pd[$1]->hh = 0.;
    Simul->chronos->pd[$1]->min = 0.;
    Simul->chronos->pd[$1]->ss = 0.;
    
    //Simul->chronos->pd[$1]->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  Simul->chronos->t[$1] = TS_date2julian_dd_hm(Simul->chronos->pd[$1],Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]); 
}
| LEX_BEGIN LEX_EQUAL LEX_INT LEX_MONTH LEX_INT flottant flottant
{
  Simul->chronos->day_d = $3;
  Simul->chronos->month = $4;
  Simul->chronos->year[BEGINNING] = $5;
  Simul->chronos->hour_h = $6 + $7 / 60;
  LP_printf(Simul->poutputs,"simulation starts on %d %s %d at %f h\n",
	    Simul->chronos->day_d,TS_name_month(Simul->chronos->month,Simul->poutputs),
	  Simul->chronos->year[BEGINNING],Simul->chronos->hour_h);
  //Simul->t_extrema = time_extrema(Simul->chronos); // !!! Must be added into ProSe !!! // SW 04/02/2021 do it after
  //LP_printf(Simul->poutputs,"First julian day is %lf and last julian day is %lf\n\n",Simul->t_extrema[BEGINNING_TS],Simul->t_extrema[END_TS]); 
}
| LEX_TIMESTEP_DA LEX_EQUAL LEX_INT
{
  Simul->passim->nstep = $3;
}
;

def_settings : LEX_SET LEX_EQUAL LEX_OPENING_BRACE settings LEX_CLOSING_BRACE
;

settings : setting settings
| setting
;

setting : set_hydro
| set_transport
| set_biology
| set_param_da
| def_species_and_macrospecies
| set_heat_transport
| set_safran_access
| regime
| calc_modes
| set_bios
| set_params_general
;

regime : LEX_REGIME LEX_EQUAL LEX_CALC_STATE
{
    Simul->regime = $3;
	Simul->pchyd->settings->calc_state = $3;
	LP_printf(Simul->poutputs,"calculation in %s state\n",HYD_name_calc_state(Simul->pchyd->settings->calc_state));
 
}
;
calc_modes : calc_mode calc_modes
| calc_mode
;
calc_mode : LEX_CALC_MODE LEX_EQUAL LEX_ANSWER
{
   Simul->calc_mode[$1] = $3;
}
;
//set_obs_pertub : LEX_OBS_PERTUB LEX_EQUAL LEX_ANSWER
//{
//   Simul->passim->obs_perturb_answer = $3;
//}
//;

set_params_general : LEX_NUM_PARTICULES LEX_EQUAL LEX_INT
{
   //Simul->num_threads = $3;
   Simul->passim->N_particules = $3;
}
|
LEX_NUM_THREAD LEX_EQUAL LEX_INT
{
   //Simul->num_threads = $3;
   Simul->psmp->nthreads = $3;
}
| LEX_NUM_THREAD_PAR LEX_EQUAL LEX_INT
{
   Simul->num_threads_par = $3;
}
| LEX_NUM_TUBE LEX_EQUAL LEX_INT
{
   Simul->ntube_default = $3;
}
| LEX_SEUIL_CHLA LEX_EQUAL mesure
{
    Simul->passim->seuil_chla = $3;
}
| LEX_ERROR_OBS LEX_EQUAL flottant
{
    Simul->passim->error_obs_sigma = $3;
}
| LEX_S_PERCENT LEX_EQUAL flottant
{
  for(i=0; i < NPARAMDA; i++) // MH 10/03/2022 : [NPARAMDA] added to facilitate different random walk for each param   
    {
      Simul->passim->s_percent[i] = $3;
      LP_printf(Simul->poutputs,"Default random walk of  %s = %3.4f \n",PROSE_name_param(i),Simul->passim->s_percent[i]);
      
    }

}
| LEX_RANDOM_WALK LEX_EQUAL LEX_RANDOM_WALK_VAL
{
    Simul->passim->random_walk = $3;
}

| LEX_ALPHA_DA LEX_EQUAL LEX_DOUBLE
{
    Simul->passim->alpha = $3;
}
| LEX_DA_METHOD LEX_EQUAL LEX_DA_METHOD_VAL
{
    // SW 16/11/2021 add da method setting
    Simul->passim->method = $3;
}
| LEX_SOLVER LEX_EQUAL LEX_SOLVER_TYPE
{
   Simul->solver = $3;
}
| LEX_DEFAULT_T_INFLOWS LEX_EQUAL mesure // SW 10/03/2021
{
   
   Simul->default_t_inflows = $3; // SW 10/03/2021 c'est mieux de ne pas mettre dans la structure Simul 
   LP_printf(Simul->poutputs, "default values T = %f\n", Simul->default_t_inflows);
}
// MH 26/04/2022: the reading of moving average parameters
| LEX_DA_MOV_AVE LEX_EQUAL LEX_ANSWER
{
  Simul->passim->da_mov_ave = $3;
}
| LEX_MOV_AVE_THETA LEX_EQUAL flottant
{
  Simul->passim->mov_ave_theta = $3;
}
| LEX_LNEIGH LEX_EQUAL mesure
{
  Simul->passim->lneigh = $3;
}
| LEX_NVAL_MIN LEX_EQUAL flottant
{
  Simul->passim->nval_min = $3;
}


;

set_bios : set_bio set_bios
| set_bio
;
set_bio : LEX_EPS LEX_EQUAL flottant
{
    temp_var = $3;
    Simul->settings->epsilon = temp_var;
  LP_printf(Simul->poutputs,"epsilon = %f\n",temp_var);
}
| LEX_COMPPHY LEX_EQUAL flottant
{
    temp_var = $3;
    Simul->settings->nb_comp_phy = temp_var;
        if (temp_var == 3) {
          Simul->settings->phy2[PHYF] = PHY2PHYF;
          Simul->settings->phy2[PHYS] = PHY2PHYS;
          Simul->settings->phy2[PHYR] = PHY2PHYR;
  }

  LP_printf(Simul->poutputs,"nb_comp_phy = %d\n",Simul->settings->nb_comp_phy);

}
| LEX_PHY2 LEX_EQUAL flottant
{
  temp_var = $3;
  if (Simul->settings->nb_comp_phy == 3) 
      Simul->settings->phy2[$1] = temp_var;

}
| LEX_DBO LEX_EQUAL mesure
{
  temp_var = $3;
  Simul->settings->dbo_oxy = temp_var;
  LP_printf(Simul->poutputs,"dbo_oxy = %f\n",temp_var);
}
| LEX_LIMIT_FACTOR LEX_EQUAL LEX_LIMIT_FACTOR_VAL
{
  Simul->settings->limit_factor = $3;
  LP_printf(Simul->poutputs,"limitting factor for growth of phytoplankton is set as %s \n",name_limit_factor($3));
}

;


/*** Safran accessing: (i) Directory that contains all the Safran's files, 
(ii) the file that makes reaches and safran cells correspond to each other, 
(iii) the optional file containing for each reach the vegetal cover or some 
other zone parameters. ***/
set_safran_access : LEX_METEO LEX_EQUAL brace saccess brace
{
  Simul->calc_mode[SEB] = YES_TS; // Activation du module de transfert thermique depuis l'atmo.
}
;

saccess : corresp meteo_access // "riverside_param" to possibly add
{
}
;

corresp : LEX_ID LEX_EQUAL LEX_NAME // (i)
{
  int endoffile ;
  FILE *fp;
  int icheck=0, n_reach;
  char r_amont[1000] ="", r_aval[1000]="";
  int r_voie;
  long id_safran;

  n_reach = Simul->pchyd->counter->nreaches;
  Simul->p_rts = SEB_alloc_rts(n_reach);
  printf("n_reach: %d\n",n_reach);

  fp=LP_openr_file(file_path($3),Simul->poutputs,ERR_LP);
  fscanf(fp, "%*[^\n]\n", NULL); // Skip first "en-tête" line
  endoffile=fscanf(fp,"%s %s %d %lu\n", r_amont, r_aval, &r_voie, &id_safran);
  
  printf("Reach ID <-> Safran ID\n");
  while(endoffile!=EOF)
    {
	  strcpy(Simul->p_rts[icheck]->pid_reach->amont,r_amont);
	  strcpy(Simul->p_rts[icheck]->pid_reach->aval,r_aval);
	  Simul->p_rts[icheck]->pid_reach->voie=r_voie;
	  Simul->p_rts[icheck]->id_meteo=id_safran;
	  printf("%s %s %d <-> %lu\n", Simul->p_rts[icheck]->pid_reach->amont, Simul->p_rts[icheck]->pid_reach->aval, Simul->p_rts[icheck]->pid_reach->voie, Simul->p_rts[icheck]->id_meteo);
      icheck = icheck + 1;
      endoffile=fscanf(fp,"%s %s %d %lu\n", r_amont, r_aval, &r_voie, &id_safran);
    }

  fclose(fp);
  printf("%s closed\n",$3);

  if (icheck == n_reach)
    printf("Successful operation: the number of extracted safran's correspondencies exactly matches the number of reaches\n");
  else
    {
      printf("Failure: the number of extracted safran's correspondencies does not match the number of reaches...\n");
	  return 0;
	}

  Simul->n_meteocell = safran_cell_count(n_reach, Simul->p_rts);
  printf("Number of safran's cells: %lu\n",Simul->n_meteocell);
}
;

meteo_access : LEX_ID LEX_EQUAL brace file_access brace
;

file_access : met_dirctory_year_init met_files // SW 04/02/2021 add year_init_meteo
;

met_dirctory_year_init : met_directory year_init_meteo
{
  /* SW 04/02/2021 add here the calculation of julian day of begin and end dates, with origin = LEX_YEAR0_METEO.
     These julian day are used to extract the meteo data. it is nice to save this year0_meteo. Need a meta structur*/  
  Simul->t_extrema = (double*)malloc(NEXTREMA_TS*sizeof(double));
  
  Simul->t_extrema[BEGINNING_TS] = TS_date2julian_dd_hm(Simul->chronos->pd[BEGINNING_TS],Simul->year0_meteo,Simul->poutputs);
  Simul->t_extrema[END] = TS_date2julian_dd_hm(Simul->chronos->pd[END_TS],Simul->year0_meteo,Simul->poutputs);

}
;
met_directory : LEX_FOLDER LEX_EQUAL LEX_NAME // (ii)
{
  int noname;
  char *new_name;

  new_name = (char *)calloc(strlen($3) + 17,sizeof(char));
  sprintf(new_name,"FOLDER_SAFRAN=%s",$3);
  printf("%s\n",new_name);
  noname = putenv(new_name);
  if (noname == -1)
    LP_error(Simul->poutputs,"File %s, line %d : undefined variable FOLDER_SAFRAN\n",current_read_files[pile].name,line_nb);
}
;

year_init_meteo : LEX_YEAR0_METEO LEX_EQUAL LEX_INT
{
  Simul->year0_meteo = $3; // default 1850
}
;
met_files : LEX_ID LEX_EQUAL brace sfiles brace
;

sfiles : sfile sfiles
| sfile
;

sfile : LEX_MET LEX_EQUAL LEX_NAME
{
  char name[STRING_LENGTH], name2[STRING_LENGTH];
  FILE *fp;
  long n_reach, n_safran;
  double t_beg, t_end;
  
  n_reach = Simul->pchyd->counter->nreaches; //n_reach = N_REACH;
  n_safran = Simul->n_meteocell;
  t_beg = Simul->t_extrema[BEGINNING_TS];
  t_end = Simul->t_extrema[END_TS];
  sprintf(name,"%s",$3);
  sprintf(name2,"%s/%s",getenv("FOLDER_SAFRAN"),name);
  fp=LP_openr_file(name2,Simul->poutputs,ERR_LP);

  printf("Reading meteo variable: %s\n",SEB_name_meteo($1));

  if (check == 0)
	Simul->p_met = SEB_alloc_met(n_safran);
  if ($1 == P_ATM)
	PROSE_Patm_reading_in_inputy(check, fp, n_reach, n_safran, $1, t_beg, t_end, Simul->p_rts, Simul->p_met);
  else
	PROSE_safransreading_in_inputy(check, fp, n_reach, n_safran, $1, t_beg, t_end, Simul->p_rts, Simul->p_met);
  fclose(fp);
  printf("%s closed\n\n",$3);
  check++; 
}
;

/*** To trigger the heat transport ***/
set_heat_transport : heat_transport_intro brace temperature_objects brace //SW 10/03/2021
;

heat_transport_intro : LEX_HEAT_TRANSPORT LEX_EQUAL // AB 02.10.2019 for heat transport, and libseb
{
  pcarac_heat_ttc = TTC_init_carac(); // Initialize the whole structure !!!
  //pcarac_heat_ttc->theta = theta_T; // SW 10/03/2021
  pcarac_heat_ttc->regime = Simul->regime;
  pcarac_heat_ttc->count[NELE_TTC] = Simul->pchyd->counter->nele_tot;
  Simul->pcarac_heat_ttc = pcarac_heat_ttc;
  //  Simul->calc_mode[H_T] = YES_TS;//NF 12/10/2020 controlled only in the settings as well as rive, ttc, etc
}
;
// SW 10/03/2021 modifying

temperature_objects : temperature_object temperature_objects
| temperature_object
;

temperature_object : LEX_THETA_T LEX_EQUAL flottant
{
  
  pcarac_heat_ttc->theta = $3;
  LP_printf(Simul->poutputs, "theta for T = %f\n", $3);
  
}
| LEX_INIT_T LEX_EQUAL a_unit_f f_ts
{
  Simul->init_T = $4; // SW 06/10/2020 add a_unit_f for controlling unit_t and unit_f in f_ts.
  Simul->init_from_file = YES_TS;
  printf("Temperature init. file has been read!\n");
}
| LEX_DEFAULT_T_INFLOWS LEX_EQUAL mesure // SW 10/03/2021
{
   LP_printf(Simul->poutputs, "default values T = %f\n", Simul->default_t_inflows);
   Simul->default_t_inflows = $3; // SW 10/03/2021 c'est mieux de ne pas mettre dans la structure Simul 
}
;

/* This part describes the modelled species and their parameters */
def_species_and_macrospecies : def_species
| def_species def_macrospecies
;

def_species : intro_species layers_species LEX_CLOSING_BRACE //SW 01/03/2017 if we define several layers in section species
{
  for (i = 0; i < NLAYERS; i++) {
    if (def_spe[i] == NO_TS) {
      for (e = 0; e < NSPECIES; e++) {
	  if(pspec2[WATER][e] != NULL) // SW 22/02/2019 if we define a species in water, allocate for sediment
	     pspec2[i][e] = (s_species **)calloc(Simul->counter_bio->nsubspecies[e],sizeof(s_species *));

	  for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++) {
	  //if (pspec2[type_1st_layer][e][j] != NULL){
	  if (pspec2[type_1st_layer][e] != NULL){ // SW 22/02/2019 remove [j]
	    pspec2[i][e][j] = copy_species(pspec2[i][e][j],pspec2[type_1st_layer][e][j],Simul->counter_bio, Simul->settings->nb_comp_phy,Simul->settings->dbo_oxy);
        pspec2[i][e][j]->diff_mol_ttc = pspec2[type_1st_layer][e][j]->diff_mol_ttc;
	  }
	  }
    }
  }
}
 LP_printf(Simul->poutputs,"\nspecies defined : \n");
}
;

intro_species :  LEX_SPECIES LEX_EQUAL LEX_OPENING_BRACE
{
  def_spe[WATER] = def_spe[VASE] = def_spe[PERIPHYTON] = NO_TS;
  LP_printf(Simul->poutputs,"\nEntering species definition...\n");
  /*for (i = 0; i < NLAYERS; i++) {
    for (e = 0; e < NSPECIES; e++){  SW 22/02/2019 remove two for
      pspec2[i][e] = new_species();

	  }
  }*/
}
;

layers_species : layer_species layers_species //SW 01/03/2017 if we define several layeys in section species
|layer_species
;

layer_species : intro_layer_species species LEX_CLOSING_BRACE
{
  //NF 11/10/2020 setting up the nb of transported variables
  int nb_var=nspecies;//NF 11/10/2020
  
  if(layspe==WATER){//NF 11/10/2020
    //if(Simul->calc_mode[H_T]==YES_TS)//NF 11/10/2020//NF 12/10/2020 would have been important if pcarac_heat_ttc wouldn't exist and only pcarac_ttc would !
    //nb_var++;//NF 11/10/2020
  //Simul->counter_bio->nspecies=nb_var;//NF 11/10/2020
  }//NF 11/10/2020
}
;

intro_layer_species : LEX_LAYER LEX_EQUAL LEX_OPENING_BRACE
{
  char *name_lay;
  layspe = $1;
  def_spe[layspe] = YES_TS;

  if (type_1st_layer == CODE) 
    type_1st_layer = layspe;

  //else for (e = 0; e < NSPECIES; e++) { // SW 16/09/2020 allocation is below, when we declare different physiological propertiesn, this is a bug
    //pspec2[layspe][e] = (s_species **)calloc(Simul->counter_bio->nsubspecies[e],sizeof(s_species *));
	//for (j = 0; j < Simul->counter->nsubspecies[e]; j++) {
      //if (pspec2[type_1st_layer][e][j] != NULL)
	//pspec2[layspe][e][j] = copy_species(pspec2[layspe][e][j],pspec2[type_1st_layer][e][j],Simul->settings->nb_comp_phy);
    //}
  //}
  name_lay = name_layer(layspe);
  LP_printf(Simul->poutputs,"\nspecies in %s compartment : \n",name_layer(layspe));
  free(name_lay);
}
;
species : one_species species
| one_species
;

one_species : intro_onespecies subspecies
{
    char *name_spe;
	name_spe = name_species(var);
	LP_printf(Simul->poutputs,"%d subspecies of %s read\n",num,name_spe); //NF 22/2/2011
    free(name_spe);
    //int j;
    pspec2[layspe][var] = (s_species **)calloc(num,sizeof(s_species *));
	pspecies[var] = (s_species_ttc **)calloc(num,sizeof(s_species_ttc *));
	for(i = 0; i < num; i++)
	{
		pspecies[var][i] = TTC_create_species();
		pspecies[var][i]->pchronos=CHR_copy_chronos(Simul->chronos);//Each species will have its own chronometer
    }
  pspec2[layspe][var][0] = $2;
  pspecies[var][0]->name = pspec2[layspe][var][0]->name;
  pspecies[var][0]->type = pspec2[layspe][var][0]->type_ttc;
  //pspecies[var][0]->plink->pbase_ttc->psolute->diff_mol = pspec2[layspe][var][0]->diff_mol_ttc;
  pspecies[var][0]->media_type = pspec2[layspe][var][0]->media_type_ttc;
  pspecies[var][0]->calc_process[DIFFUSION_TTC] = pspec2[layspe][var][0]->calc_process_ttc[DIFFUSION_TTC];
  pspecies[var][0]->calc_process[CONVECTION_TTC] = pspec2[layspe][var][0]->calc_process_ttc[CONVECTION_TTC];
  pspecies[var][0]->oxygen = NO_TS;
  LP_printf(Simul->poutputs,"ra = %d,rd = %d layspe = %d\n",pspec2[layspe][var][0]->calc_process_ttc[CONVECTION_TTC],pspec2[layspe][var][0]->calc_process_ttc[DIFFUSION_TTC],layspe);
  for (i = 1; i < Simul->counter_bio->nsubspecies[var]; i++){
    pspec2[layspe][var][i] = pspec2[layspe][var][i-1]->next;
    pspecies[var][i]->name = pspec2[layspe][var][i]->name;
    pspecies[var][i]->type = pspec2[layspe][var][i]->type_ttc;
    pspecies[var][i]->media_type = pspec2[layspe][var][i]->media_type_ttc;
    pspecies[var][i]->calc_process[DIFFUSION_TTC] = pspec2[layspe][var][i]->calc_process_ttc[DIFFUSION_TTC];
    pspecies[var][i]->calc_process[CONVECTION_TTC] = pspec2[layspe][var][i]->calc_process_ttc[CONVECTION_TTC];	
    pspecies[var][i]->oxygen = NO_TS;
	}
}
;

intro_onespecies : LEX_ONESPECIES
{
    char *name_spe;
	var = $1;
    num = 0;
    Simul->counter_bio->nsubspecies[var] = 0;
    name_spe = name_species(var);
    LP_printf(Simul->poutputs,"\nsub-species of %s :\n",name_spe);
	free(name_spe);

}
;

subspecies : one_subspecies subspecies
{
  $$ = chain_subspecies($1,$2);
}
| one_subspecies
{
  $$ = $1;
}
;

one_subspecies : intro_subspecies atts_species LEX_CLOSING_BRACE
{
  $$ = spec2;
  spec2 = NULL; 
}
| subspecies_name 
{
  $$ = spec2;
  spec2 = NULL; 
}
;

intro_subspecies : subspecies_name LEX_EQUAL LEX_OPENING_BRACE
;

subspecies_name : LEX_NAME
{
  spec2 = init_subspecies(var,Simul->settings->dbo_oxy); 
  spec2->num = ++num;
  Simul->counter_bio->nsubspecies[var]++;
  spec2->mb = new_mass_balance_species();
  
  sprintf(spec2->name,"%s",$1);//SW 21/02/2017
  LP_printf(Simul->poutputs,"\n  %s\n",spec2->name);
  free($1);
  
  if(layspe==WATER)//NF 11/10/2020 Otherwise the definition of an existing specie in another compartment will increment the counter generating an error
    nspecies++; // pour libttc
  
}
;

atts_species : att_species atts_species
| att_species
;

att_species : LEX_KMICH LEX_ONESPECIES LEX_EQUAL mesure
{
  spec2->kmich_species[$2] = $4;
  LP_printf(Simul->poutputs,"kmich %s = %f g/m^3\n",name_species($2),spec2->kmich_species[$2]);
  //printf("kmich %s = %f g/m^3\n",name_species($2),spec2->kmich_species[$2]);
}
| LEX_KMICH LEX_COMP LEX_EQUAL mesure
{
  spec2->kmich_comp[$2] = $4;
  LP_printf(Simul->poutputs,"kmich %s = %f g/m^3\n",name_component($2),spec2->kmich_comp[$2]);
  //printf("kmich %s = %f g/m^3\n",name_component($2),spec2->kmich_comp[$2]);
}
| LEX_KMICH LEX_ONEANNEXVAR LEX_EQUAL mesure
{
  spec2->kmich_annex[$2] = $4;
  LP_printf(Simul->poutputs,"kmich %s = %f g/m^3\n",name_annex_var($2),spec2->kmich_annex[$2]);
  //printf("kmich %s = %f g/m^3\n",name_annex_var($2),spec2->kmich_annex[$2]);
}
| LEX_KLIM LEX_ONESPECIES LEX_EQUAL mesure
{
  spec2->klim_species[$2] = $4;
  LP_printf(Simul->poutputs,"klim %s = %f g/m^3\n",name_species($2),spec2->klim_species[$2]);
  //printf("klim %s = %f g/m^3\n",name_species($2),spec2->klim_species[$2]);
}
| LEX_NUTC LEX_COMP LEX_EQUAL mesure 
{
  if ($4 > 0.) spec2->nut_C[$2] = 1. / $4;
  else spec2->nut_C[$2] = 0.;
  LP_printf(Simul->poutputs,"%s/C = %f g/g\n",name_component($2),spec2->nut_C[$2]);
  //printf("%s/C = %f g/g\n",name_component($2),spec2->nut_C[$2]);
}
| LEX_NUTC LEX_ONEANNEXVAR LEX_EQUAL mesure 
{
  if ($4 > 0.) spec2->nut_C[NCOMP_RIVE] = 1. / $4; // SW 02/03/2017 $2 ou NCOMP_RIVE ?, c'est pas le même nombre, remplace NCOMP_RIVE par $2
  else spec2->nut_C[NCOMP_RIVE] = 0.;
  LP_printf(Simul->poutputs,"chla/C = %f g/g\n",spec2->nut_C[NCOMP_RIVE]);
  //printf("chla/C = %f g/g\n",spec2->nut_C[NCOMP_RIVE]);
}
| LEX_NUTC LEX_ONEANNEXVAR LEX_EQUAL read_units f_ts
{
   s_ft *pft;
   //pft = $5;
   pft = TS_browse_ft($5,BEGINNING_TS);
   pft->ft = 1./pft->ft; // SW inversion chla/C   
   if(pft->next != NULL)
   {   
   pft = pft->next;
   while(pft->next != NULL)
   {
   pft->ft = 1./pft->ft; // SW inversion chla/C
   pft = pft->next;   
   }
   pft->ft = 1./pft->ft; // SW inversion chla/C
   pft = TS_browse_ft(pft,BEGINNING_TS);
   }
   spec2->nut_C_variable[NCOMP_RIVE] = pft;
   unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
   unit_t = 1.;
}
| LEX_NUTC LEX_ONEANNEXVAR LEX_EQUAL read_units date_f_ts
{
   s_ft *pft;
   //pft = $5;
   pft = TS_browse_ft($5,BEGINNING_TS);
   pft->ft = 1./pft->ft; // SW inversion chla/C   
   if(pft->next != NULL)
   {   
   pft = pft->next;
   while(pft->next != NULL)
   {
   pft->ft = 1./pft->ft; // SW inversion chla/C
   pft = pft->next;   
   }
   pft->ft = 1./pft->ft; // SW inversion chla/C
   pft = TS_browse_ft(pft,BEGINNING_TS);
   }
   spec2->nut_C_variable[NCOMP_RIVE] = pft;
   unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
   unit_t = 1.;
}
| LEX_PARAMD LEX_EQUAL mesure
{
  spec2->dissolved->paramd[$1] = $3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_paramd($1),spec2->dissolved->paramd[$1],unit_paramd($1));
  //printf("%s = %f %s\n",name_paramd($1),spec2->dissolved->paramd[$1],unit_paramd($1));
}
| LEX_PARAMP LEX_EQUAL mesure
{
  spec2->particulate->paramp[$1] = $3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_paramp($1),spec2->particulate->paramp[$1],unit_paramp($1));
  //printf("%s = %f %s\n",name_paramp($1),spec2->particulate->paramp[$1],unit_paramp($1));
}
| LEX_PARAML LEX_EQUAL mesure
{
  spec2->particulate->living->paraml[$1] = $3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_paraml($1),spec2->particulate->living->paraml[$1],unit_paraml($1));
  //printf("%s = %f %s\n",name_paraml($1),spec2->particulate->living->paraml[$1],unit_paraml($1));
}
| LEX_PARAML LEX_EQUAL read_units f_ts // SW 03/02/2020 time variant parameters
{
  spec2->particulate->living->paraml_variable[$1] = $4;
  LP_printf(Simul->poutputs,"time variant parameter %s defined\n",name_paraml($1));
  unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
  unit_t = 1.;
}
| LEX_PARAML LEX_EQUAL read_units date_f_ts // SW 03/02/2020 time variant parameters
{
  spec2->particulate->living->paraml_variable[$1] = $4;
  LP_printf(Simul->poutputs,"time variant parameter %s defined\n",name_paraml($1));
  unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
  unit_t = 1.;
}
| LEX_PARAMM LEX_EQUAL mesure
{
  if (var < NPART)
    spec2->particulate->mineral->paramm[$1] = $3;
  else spec2->dissolved->mineral->paramm[$1] = $3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_paramm($1),$3,unit_paramm($1));
  //printf("%s = %f %s\n",name_paramm($1),$3,unit_paramm($1));
}
| LEX_DEFAULT_C_INFLOWS LEX_EQUAL mesure // SW 11/03/2021
{
  spec2->default_C_inflows = $3;
  LP_printf(Simul->poutputs, "default c value = %f \n for species %s\n", spec2->default_C_inflows, spec2->name);
}
| process
| reactions
{
  if (spec2->type == PARTICULATE) {
    spec2->particulate->mineral->reactions = (s_reac **)calloc(spec2->particulate->mineral->nreactions,sizeof(s_reac *));
    for (i = 0; i < spec2->particulate->mineral->nreactions; i++) {
      spec2->particulate->mineral->reactions[i] = reaction2;
      reaction2 = reaction2->next;
    }
  }
  else {
    spec2->dissolved->mineral->reactions = (s_reac **)calloc(spec2->dissolved->mineral->nreactions,sizeof(s_reac *));
    for (i = 0; i < spec2->dissolved->mineral->nreactions; i++) {
      spec2->dissolved->mineral->reactions[i] = reaction2;
      reaction2 = reaction2->next;
    }
  }
}
| LEX_TRANSPORT_EQ LEX_EQUAL LEX_TRANSPORT_MODE
{
  spec2->type_ttc=$3;
  LP_printf(Simul->poutputs,"\ttype of problem %s transport\n",TTC_name_equa_type($3));
}
| LEX_MED_TYPE LEX_EQUAL LEX_MEDIA_TYPE
{

  spec2->media_type_ttc=$3;
  LP_printf(Simul->poutputs,"\tType of media %s\n",TTC_environnement($3));
}
| LEX_PROC_TRANSP LEX_EQUAL LEX_ANSWER
{
  spec2->calc_process_ttc[$1]=$3;
  LP_printf(Simul->poutputs,"\tAccounting for %s = %s\n",TTC_name_process($1),LP_answer(spec2->calc_process_ttc[$1],Simul->poutputs));
}
| LEX_DIFF_MOL LEX_EQUAL mesure
{
  spec2->diff_mol_ttc = $3;
  LP_printf(Simul->poutputs,"\t %s diff_mol= %f\n",spec2->name,spec2->diff_mol_ttc);
}
;

process : intro_process atts_process LEX_CLOSING_BRACE
{
  assign_process_function(spec2,proc);
}
| intro_process LEX_CLOSING_BRACE//LV 01/06/2011
{
  assign_process_function(spec2,proc);
}
;

intro_process : LEX_PROCESS LEX_EQUAL LEX_OPENING_BRACE
{
  proc = $1;
  spec2->calc_process[proc] = YES_RIVE;
  LP_printf(Simul->poutputs,"--%s-- %d\n",name_process(proc),YES_RIVE); // SW 30/05/2018 ici YES_RIVE = 1
  //printf("--%s--\n",name_process(proc));
  spec2 = init_proc(proc,spec2,3);

  if (proc == ADS_DESORPTION) {
    for (e = 0; e < NPART; e++) 
      spec2->dissolved->mineral->ads_desorption->adsorbs_on[e] = (int *)calloc(Simul->counter_bio->nsubspecies[e],sizeof(int));
  }
 }
;

atts_process : att_process atts_process
| att_process
;

att_process : LEX_PHOT LEX_EQUAL mesure
{
  spec2->particulate->living->photosynthesis->phot[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_phot($1),$3,unit_param_phot($1));
  //printf("%s = %f %s\n",name_param_phot($1),$3,unit_param_phot($1));
}
| LEX_PHOT LEX_EQUAL read_units f_ts
{
  spec2->particulate->living->photosynthesis->phot_variable[$1] = $4;
  LP_printf(Simul->poutputs,"%s \n",name_param_phot($1));
  unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
  unit_t = 1.;
}
| LEX_PHOT LEX_EQUAL read_units date_f_ts
{
  spec2->particulate->living->photosynthesis->phot_variable[$1] = $4;
  LP_printf(Simul->poutputs,"%s \n",name_param_phot($1));
  unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
  unit_t = 1.;
}
| LEX_MORT LEX_EQUAL mesure
{
  spec2->particulate->living->mortality->mort[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_mort($1),$3,unit_param_mort($1));
  //printf("%s = %f %s\n",name_param_mort($1),$3,unit_param_mort($1));
}
| LEX_MORT LEX_EQUAL read_units f_ts
{
  spec2->particulate->living->mortality->mort_variable[$1] = $4;
  unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
  unit_t = 1.;
}
| LEX_MORT LEX_EQUAL read_units date_f_ts
{
  spec2->particulate->living->mortality->mort_variable[$1] = $4;
  unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
  unit_t = 1.;
}
| LEX_RESP LEX_EQUAL mesure
{
  spec2->particulate->living->respiration->resp[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_resp($1),$3,unit_param_resp($1));
  //printf("%s = %f %s\n",name_param_resp($1),$3,unit_param_resp($1));
} 
| LEX_RESP LEX_EQUAL read_units f_ts
{
  spec2->particulate->living->respiration->resp_variable[$1] = $4;
  LP_printf(Simul->poutputs,"%s \n",name_param_resp($1));
  unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
  unit_t = 1.;
}
| LEX_RESP LEX_EQUAL read_units date_f_ts
{
  spec2->particulate->living->respiration->resp_variable[$1] = $4;
  LP_printf(Simul->poutputs,"%s \n",name_param_resp($1));
  unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
  unit_t = 1.;
}
| LEX_GROWTH LEX_EQUAL mesure
{
  spec2->particulate->living->growth->growth[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_growth($1),$3,unit_param_growth($1));
  //printf("%s = %f %s\n",name_param_growth($1),$3,unit_param_growth($1));
}
| LEX_GROWTH LEX_EQUAL read_units f_ts
{
  spec2->particulate->living->growth->growth_variable[$1] = $4;
  LP_printf(Simul->poutputs,"%s \n",name_param_growth($1));
  unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
  unit_t = 1.;
}
| LEX_GROWTH LEX_EQUAL read_units date_f_ts
{
  spec2->particulate->living->growth->growth_variable[$1] = $4;
  LP_printf(Simul->poutputs,"%s \n",name_param_growth($1));
  unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
  unit_t = 1.;
}
| LEX_GRAZ LEX_EQUAL mesure
{
  spec2->particulate->living->grazing->graz[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_graz($1),$3,unit_param_graz($1));
  //printf("%s = %f %s\n",name_param_graz($1),$3,unit_param_graz($1));
}
| LEX_EXCR LEX_EQUAL mesure
{
  spec2->particulate->living->excretion->excr[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_excr($1),$3,unit_param_excr($1));
  //printf("%s = %f %s\n",name_param_excr($1),$3,unit_param_excr($1));
}
| LEX_REA LEX_EQUAL mesure
{
  spec2->dissolved->gas->reaeration->rea[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_rea($1),$3,unit_param_rea($1));
  //printf("%s = %f %s\n",name_param_rea($1),$3,unit_param_rea($1));
}
| LEX_REA LEX_EQUAL read_units f_ts
{
  spec2->dissolved->gas->reaeration->rea_variable[$1] = $4;
  //LP_printf(Simul->poutputs,"%s \n",name_param_rea($1));
  LP_printf(Simul->poutputs,"%s t = %f val = %f\n",name_param_rea($1),spec2->dissolved->gas->reaeration->rea_variable[$1]->t,spec2->dissolved->gas->reaeration->rea_variable[$1]->ft);
   unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
   unit_t = 1.;
}
| LEX_REA LEX_EQUAL read_units date_f_ts
{
  spec2->dissolved->gas->reaeration->rea_variable[$1] = $4;
  LP_printf(Simul->poutputs,"%s \n",name_param_rea($1));
   unit_f = 1.; // SW 04/02/2020 need to reset to 1. for after usage
   unit_t = 1.;
}
| LEX_HYDR LEX_EQUAL mesure
{
  if (spec2->type == DISSOLVED)
    spec2->dissolved->mineral->hydrolysis->hydr[$1]=$3;
  else
    spec2->particulate->mineral->hydrolysis->hydr[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_hydr($1),$3,unit_param_hydr($1));
  //printf("%s = %f %s\n",name_param_hydr($1),$3,unit_param_hydr($1));
}
| LEX_ADS_SENEQUE LEX_EQUAL mesure
{
  spec2->dissolved->mineral->ads_desorption->ads_des[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_ads_des_seneque($1),$3,unit_param_ads_des_seneque($1));
  //printf("%s = %f %s\n",name_param_ads_des_seneque($1),$3,unit_param_ads_des_seneque($1));
}
| LEX_ADS_FR LEX_EQUAL mesure
{
  spec2->dissolved->mineral->ads_desorption->ads_des[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_ads_des_freundlich($1),$3,unit_param_ads_des_freundlich($1));
  //printf("%s = %f %s\n",name_param_ads_des_freundlich($1),$3,unit_param_ads_des_freundlich($1));
  spec2->dissolved->mineral->ads_desorption->adsorption_desorptionfunction = &ads_desorption_freundlich;
}
| adsorbs_on
;

adsorbs_on : LEX_ADS_ON LEX_EQUAL LEX_OPENING_BRACE list_adsorbs_on_species LEX_CLOSING_BRACE
;

list_adsorbs_on_species : LEX_ONESPECIES LEX_INT list_adsorbs_on_species
{
  spec2->dissolved->mineral->ads_desorption->adsorbs_on[$1][$2-1] = YES_RIVE; // SW ic YES_RIVE == 1
  LP_printf(Simul->poutputs,"adsorbs on %s\n",pspec2[layspe][$1][$2-1]->name);
  //printf("adsorbs on %s\n",pspec2[layspe][$1][$2-1]->name);
}
| LEX_ONESPECIES LEX_INT 
{
  spec2->dissolved->mineral->ads_desorption->adsorbs_on[$1][$2-1] = YES_RIVE; // SW ic YES_RIVE == 1
  LP_printf(Simul->poutputs,"adsorbs on %s\n",pspec2[layspe][$1][$2-1]->name);
  //printf("adsorbs on %s\n",pspec2[layspe][$1][$2-1]->name);
}
;

reactions : intro_reactions reacs LEX_CLOSING_BRACE
;

intro_reactions : LEX_REACTIONS LEX_EQUAL LEX_OPENING_BRACE
{
  spec2->calc_process[RAD_DECAY] = YES_RIVE; // SW 31/05/2018 ici YES_RIVE = 1
  spec2->proc[RAD_DECAY] = reaction;//NF 7/10/2020 warning  that the operation is of incompatible types. It is true that in the species structur in c-rive there is only a pointer toward void (*proc[PHOT])(double,int,int,int,s_species *,s_simul *); // SW 26/04/2018//NF 7/10/2020 why proc[PHOT] only, what abour the other processes ?
  //printf("Definition of the species' reactions...\n");
  LP_printf(Simul->poutputs,"Definition of the species' reactions...\n");
}
;

reacs : reac reacs
{
  $$ = chain_reactions($1,$2);
}
| reac
{
  $$ = $1;
}
;

reac : intro_reaction atts_reac LEX_CLOSING_BRACE
{
  $$ = reaction2;
}
;

intro_reaction : LEX_REACTION LEX_EQUAL LEX_OPENING_BRACE
{
  if (spec2->type == PARTICULATE) {
    spec2->particulate->mineral->nreactions++;
  }
  else {
    spec2->dissolved->mineral->nreactions++;
  }
  reaction2 = init_reaction();
  //printf("reaction1\n");
  LP_printf(Simul->poutputs,"reaction1\n");
}
;

atts_reac : att_reac atts_reac
| att_reac
;

att_reac : LEX_RAD_DECAY LEX_EQUAL mesure
{
  reaction2->reac[$1]=$3;
  LP_printf(Simul->poutputs,"%s = %f %s\n",name_param_rad_decay($1),$3,unit_param_rad_decay($1));
  //printf("%s = %f %s\n",name_param_rad_decay($1),$3,unit_param_rad_decay($1));
}
| other_reactors
| products
| control
;

other_reactors : intro_others others LEX_CLOSING_BRACE
{
  reaction2->nother_reactors = nother_reactors;
  reaction2->name_other_reactors = (char **)calloc(nother_reactors,sizeof(char *));
  reaction2->name_other_reactors = other_reactors_new;
  
  reaction2->stoechio_other_reactors = (double *)calloc(nother_reactors,sizeof(double));
  reaction2->stoechio_other_reactors = stoechio_other_reactors_new;
}
;

intro_others : LEX_OTHERS LEX_EQUAL LEX_OPENING_BRACE
{
  nother_reactors = 0;
  //printf("\tother reactors :\n");
  LP_printf(Simul->poutputs,"\tother reactors :\n");
}
;

others : other others
| other
;

other : LEX_NAME mesure
{
  nother_reactors++;
  
  other_reactors_new = (char **)calloc(nother_reactors,sizeof(char *));
  stoechio_other_reactors_new = (double *)calloc(nother_reactors,sizeof(double));
  
  for (i = 0; i < nother_reactors - 1; i++) {
    other_reactors_new[i] = other_reactors_old[i];
    stoechio_other_reactors_new[i] = stoechio_other_reactors_old[i];
  }
  other_reactors_new[nother_reactors - 1] = $1;
  stoechio_other_reactors_new[nother_reactors - 1] = $2;
  other_reactors_old = other_reactors_new;
  stoechio_other_reactors_old = stoechio_other_reactors_new;
  //printf("\t\t%s\n",$1);
  LP_printf(Simul->poutputs,"\t\t%s\n",$1);
}
;

products : intro_products prods LEX_CLOSING_BRACE
{
  reaction2->nproducts = nproducts;
  reaction2->name_products = (char **)calloc(nproducts,sizeof(char *));
  reaction2->name_products = products_new;
  
  reaction2->stoechio_products = (double *)calloc(nproducts,sizeof(double));
  reaction2->stoechio_products = stoechio_products_new;
}
;

intro_products : LEX_PRODUCTS LEX_EQUAL LEX_OPENING_BRACE
{
  nproducts = 0;
  //printf("\tproducts :\n");
  LP_printf(Simul->poutputs,"\tproducts :\n");
}
;

prods : prod prods
| prod
;

prod : LEX_NAME mesure
{
  nproducts++;
  products_new = (char **)calloc(nproducts,sizeof(char *));
  stoechio_products_new = (double *)calloc(nproducts,sizeof(double));
  for (i = 0; i < nproducts - 1; i++) {
    products_new[i] = products_old[i];
    stoechio_products_new[i] = stoechio_products_old[i];
  }
  products_new[nproducts - 1] = $1;
  stoechio_products_new[nproducts - 1] = $2;
  products_old = products_new;
  stoechio_products_old = stoechio_products_new;
  //printf("\t\t%s\n",$1);
  LP_printf(Simul->poutputs,"\t\t%s\n",$1);
}
;

control : intro_control condition LEX_CLOSING_BRACE
;

intro_control : LEX_CONDITION LEX_EQUAL LEX_OPENING_BRACE
{
  reaction2->reaction_control = YES_RIVE; // SW 31/05/2018 ici YES_RIVE = 1
  //printf("\tcontrol species :\n");
  LP_printf(Simul->poutputs,"\tcontrol species :\n");
}
;

condition : LEX_NAME LEX_INF_SUP mesure
{
  //reaction2->name_control_species = $1;
  reaction2->name_control_species = (char *) malloc(MAXCHAR_PROSE*sizeof(char));
  sprintf(reaction2->name_control_species,"%s",$1);//SW 21/02/2017
  reaction2->inf_sup = $2;
  reaction2->limit_conc = $3;
  //printf("\t\t%s\n",$1);
  LP_printf(Simul->poutputs,"\t\t%s\n",$1);
  free($1);
}
;


/* This part describes the macrospecies if any */

def_macrospecies : intro_macrospec macrospecs LEX_CLOSING_BRACE
{
  int ii;
  int i;
  pmacspe=$2;
   for (ii=0;ii<Simul->counter_bio->nmacrospecies;ii++)
    {
      RIV_print_macrospecies(pmacspe,Simul->poutputs);
      Simul->p_macrospecies[pmacspe->type]=pmacspe;
      LP_printf(Simul->poutputs,"b1_input.y =  %3.2f \n ", Simul->p_macrospecies[pmacspe->type]->degradOrgMat[MACMOD][B]->val);

      // MH 25/11/2021: to store macrospecies for DA purpose
      if (Simul->counter_bio->nmacrospecies > EPS_TS) {
	Simul->passim->p_macro_da[pmacspe->type] = (s_macrospecies_RIV **)calloc(Simul->passim->N_particules,sizeof(s_macrospecies_RIV *));
	
	for(i = 0; i < Simul->passim->N_particules; i++)
	{   	   
	  Simul->passim->p_macro_da[pmacspe->type][i] = RIV_init_macrospecies(pmacspe->type,Simul->poutputs);
	  Simul->passim->p_macro_da[pmacspe->type][i]->degradOrgMat[MACMOD][B] =  RIV_init_stoch_param(Simul->poutputs);	  
	  //LP_printf(Simul->poutputs,"b1_1_input.y =  %3.2f \n ", Simul->passim->p_macro_da[pmacspe->type][i]->degradOrgMat[MACMOD][B]->val);
	}
      }// if
      
      pmacspe=pmacspe->prev;
    }
}
;

intro_macrospec : LEX_MACROSPECIES LEX_EQUAL LEX_OPENING_BRACE
{
  pmacspe=NULL;
  pstoch = NULL;
}

macrospecs : macrospec macrospecs
{ $$ = RIV_chain_macrospecies($1,$2); }
| macrospec
{ $$ = $1;}
;

macrospec : intro_onemacro atts_onemacro  LEX_CLOSING_BRACE
{
  $$=pmacspe;
}   
;

intro_onemacro : LEX_TYPEOF_MACROSPECIES LEX_EQUAL  LEX_OPENING_BRACE
{
  Simul->counter_bio->nmacrospecies++;
  pmacspe = RIV_init_macrospecies($1,Simul->poutputs);
}
;

atts_onemacro : att_onemacro atts_onemacro
| att_onemacro
;

att_onemacro : LEX_NAME
{
  pmacspe->name=$1;
  //printf("Macrospecies name %s\n",$1);
}
|  relateto LEX_COMA thresh
|  share_speciess
   //|  share_species_mop
   //|  share_species_mod
;

relateto : LEX_RELATED_TO LEX_ONESPECIES LEX_ONESPECIES
{
  int idtoc=CODE_TS;
  int id1,id2;
  id1 = $2;
  id2= $3;

  /*  if (id1 == MOD)
    LP_printf(Simul->poutputs,"id1 = %s",name_species(MOD));
  else
  LP_printf(Simul->poutputs,"id1 = %s",name_species(0));*/
    
  
  LP_printf(Simul->poutputs,"Reading Macrospecies %s related to %s (int val %d) and %s (int val %d)\n",RIV_name_macrospecies(pmacspe->type,Simul->poutputs),name_species(id1),id1,name_species(id2),id2);
  
  if (pmacspe->type==TOC){
    
    if(id1 == MOD) idtoc=MACMOD;
    else {
      idtoc = MACMOP;
      rev_threshold = YES_TS;
    }
    pmacspe->irelatedspec[idtoc]=id1;
    
    if (id2==MOP){
      idtoc=MACMOP;
    }
    else{
      idtoc = MACMOD;
      rev_threshold = YES_TS;
    }
    pmacspe->irelatedspec[idtoc]=id2;
  }
  else {
    LP_warning(Simul->poutputs,"The macrospecies number %d is not fully integrated yet in librive, please contact the developpers\n",pmacspe->type);
  }
}
;

thresh : LEX_THRESHOLD LEX_EQUAL flottant
{
  //pmacspe->threshold = $3;
  if (rev_threshold == NO_TS) {
    pmacspe->threshold->val = $3;
  }
  else
    pmacspe->threshold->val = 1-$3;
    
  //printf("Macrospecies threshold %3.2f\n",$3);

}
;
share_speciess : share_species share_speciess
| share_species
;

share_species : intro_share_mo atts_share LEX_CLOSING_BRACE
;


//share_species_mop : intro_share_mop atts_share_mop LEX_CLOSE_CURVED_BRACKET
//{
  // printf("The share of MOP was read");
//  $$ = pstoch;
//}
//;


intro_share_mo : LEX_SHARE_MO LEX_EQUAL LEX_OPENING_BRACE
{
  imacspe=$1;
  //Simul->counter->nsharespecies++;  //MH : add to counter
  LP_printf(Simul->poutputs,"The share of %s is read now \n",name_species(pmacspe->irelatedspec[imacspe]));
}
;
atts_share : att_share atts_share
| att_share
;

att_share : intro_att_share atts_stoch_param LEX_CLOSING_BRACE
{
  pmacspe->degradOrgMat[imacspe][idegval]=pstoch;
  //LP_printf(Simul->poutputs,"%s = %3.2f \n",RIV_name_stoch(idegval,Simul->poutputs),pstoch->val);
  //LP_printf(Simul->poutputs," %s =  %3.2f \n",RIV_name_stoch(idegval,Simul->poutputs),pmacspe->degradOrgMat[imacspe][idegval]->val);

   if (pstoch->range[MIN_RIV] ==0 && pstoch->range[MAX_RIV]==0)
    {
    LP_printf(Simul->poutputs,"with no defined range \n");
    } else {
     //LP_printf(Simul->poutputs,"and it ranges between %3.2f and %3.2f \n",pstoch->range[MIN_RIV],pstoch->range[MAX_RIV]);
    LP_printf(Simul->poutputs,"Min = %3.2f , Max = %3.2f \n",pmacspe->degradOrgMat[imacspe][idegval]->range[MIN_RIV],pmacspe->degradOrgMat[imacspe][idegval]->range[MAX_RIV]); 

    }
}
;

intro_att_share : LEX_BIODEG LEX_EQUAL LEX_OPENING_BRACE
{
  idegval=$1;
  pstoch = RIV_init_stoch_param(Simul->poutputs);
  //pmacspe->degradOrgMat[imacspe][idegval] = RIV_init_stoch_param(Simul->poutputs);
}
;

 
atts_stoch_param : att_stoch_param atts_stoch_param
| att_stoch_param
;


att_stoch_param : LEX_VAL LEX_EQUAL flottant 
{
  pstoch->val = $3;
  //pstoch->varying_YesorNo = NO_RIVE; // SW 08/01/2023 it is initialized in RIV_init_stoch_param with NO_RIVE
  //pmacspe->degradOrgMat[imacspe][idegval]->val = $3;
}
| LEX_VAL LEX_EQUAL read_units f_ts //MH 23/05/2022: for time varying (dynamic) b1
{
  pstoch->val_variable = $4;
  pstoch->varying_YesorNo = YES_RIVE; // SW 08/01/2023 update YES_RIVE
  LP_printf(Simul->poutputs,"time varying b1 created: t= %f, b1= %f \n", pstoch->val_variable->t/NSEC_DAY_TS, pstoch->val_variable->ft);
  unit_f = 1.;
  unit_t = 1.;
}
| LEX_VAL LEX_EQUAL read_units date_f_ts
{
  pstoch->val_variable = $4;
  pstoch->varying_YesorNo = YES_RIVE; // SW 08/01/2023 update YES_RIVE
  LP_printf(Simul->poutputs,"time varying b1 created: t= %f, b1= %f \n", pstoch->val_variable->t/NSEC_DAY_TS, pstoch->val_variable->ft);
  unit_f = 1.;
  unit_t = 1.;
}


| LEX_RANGES LEX_EQUAL LEX_OPENING_BRACE flottant LEX_COMA flottant LEX_CLOSING_BRACE
{
  
  pstoch->range[MIN_RIV]=$4;
  pstoch->range[MAX_RIV]=$6;
  /*
  int a;
  a = MIN_RIV;
  pmacspe->degradOrgMat[imacspe][idegval]->range[MIN_RIV] = $4;
  pmacspe->degradOrgMat[imacspe][idegval]->range[MAX_RIV] = $6;
  */
  // LP_printf(Simul->poutputs,"min = %3.2f , max = %3.2f \n ", pstoch->range[MIN_RIV], pstoch->range[MAX_RIV]); 
  
}
;


set_hydro : LEX_HYDRO LEX_EQUAL LEX_OPENING_BRACE hydro_atts LEX_CLOSING_BRACE
{
  Simul->pchyd->settings->solver_type = Simul->solver;
  Simul->pchyd->settings->psmp = Simul->psmp;
  Simul->pcarac_ttc = (s_carac_ttc **) malloc(Simul->passim->N_particules * sizeof(s_carac_ttc *));
  Simul->psimul_bio = (s_simul ***)malloc(Simul->passim->N_particules*sizeof(s_simul **));
  for(i = 0; i < Simul->passim->N_particules; i++)
       Simul->pcarac_ttc[i] = TTC_init_carac();
  LP_printf(Simul->poutputs,"Hydraulics simulation configured solver = %d\n",Simul->pchyd->settings->solver_type);
}
;

hydro_atts : hydro_att hydro_atts
| hydro_att
;

hydro_att : LEX_DIM LEX_EQUAL LEX_INT
{
  Simul->pchyd->settings->ndim = $3; // BL ici !!!
  LP_printf(Simul->poutputs,"ndim = %d\n",Simul->pchyd->settings->ndim);
 
}
| LEX_TYPE LEX_EQUAL LEX_CALC_STATE
{
  Simul->pchyd->settings->calc_state = $3; 
}
| LEX_CALC_CURVATURE LEX_EQUAL LEX_ANSWER
{
  Simul->pchyd->settings->calc_curvature = $3;
  LP_printf(Simul->poutputs,"calculation of the curvature radius = %s\n",
	  HYD_name_answer(Simul->pchyd->settings->calc_curvature));
  
}
| LEX_GENERALP LEX_EQUAL mesure
{
  Simul->pchyd->settings->general_param[$1] = $3;
  LP_printf(Simul->poutputs,"!!! %s = %f \n",HYD_name_general_param($1),Simul->pchyd->settings->general_param[$1]);
  
}
| LEX_SCHEM_TYPE LEX_EQUAL LEX_DEF_SCHEM_TYPE
{
  Simul->pchyd->settings->schem_type=$3;
  LP_printf(Simul->poutputs,"schem_type = %s\n",HYD_name_scheme($3));
}
| def_initialization
| datas_hydro
;

def_initialization : LEX_INIT LEX_EQUAL brace def_inits brace

def_inits : def_init def_inits
| def_init
;

def_init : init_Z
| init_Q
;

init_Z : LEX_INIT_Z LEX_EQUAL f_ts
{
  Simul->pchyd->settings->init_Z = $3;
  Simul->pinout->init_from_file[0] = YES_TS;
}
;

init_Q : LEX_INIT_Q LEX_EQUAL f_ts
{
  Simul->pchyd->settings->init_Q = $3;
  Simul->pinout->init_from_file[1] = YES_TS;
}
;

/* Data giving information on the river network : singularities, reaches, inflows, cross_sections... */
datas_hydro : data_hydro datas_hydro
| data_hydro 
;

data_hydro : LEX_NETWORK LEX_EQUAL brace network_attributes brace
{
  
  HYD_BC_faces(Simul->pchyd,Simul->poutputs);
  //HYD_print_sing(Simul->pchyd,Simul->poutputs);
  HYD_verif_singularities(Simul->pchyd,Simul->poutputs);
  HYD_calculate_pk_sing(Simul->pchyd,Simul->poutputs);
  HYD_calculate_pk_faces_centers(Simul->pchyd,Simul->poutputs);
  if (Simul->pchyd->settings->calc_curvature == YES_TS)
    HYD_calculate_curvature_radius(Simul->pchyd);
  else if (Simul->pchyd->settings->calc_curvature == NO_TS)
    HYD_zero_curvature(Simul->pchyd,Simul->poutputs);
  //assign_river();//LV nov2014 déplacé après initialisation hydro
}
;

network_attributes : network_att network_attributes
| network_att
;

network_att : def_singularities
| def_reaches
| bathymetry //type face
{
  Simul->clock->time_spent[LEC_CHR] += CHR_end_timer();//LV 3/09/2012
  CHR_begin_timer();//LV 3/09/2012
  HYD_table_hydro_faces($1,Simul->pchyd,Simul->chronos,Simul->poutputs);
  HYD_create_faces($1,Simul->pchyd,Simul->chronos,Simul->poutputs);
  LP_printf(Simul->poutputs,"ntot_faces = %d\n",Simul->pchyd->counter->nfaces_tot);
  HYD_create_elements(Simul->pchyd,Simul->poutputs);
  HYD_table_hydro_elements(Simul->pchyd,Simul->chronos,Simul->poutputs);
  //HYD_print_hydro_table(Simul->pchyd->p_reach[35]->p_ele[0]->center->hydro);
  Simul->clock->time_spent[TAB_CHR] += CHR_end_timer();//LV 3/09/2012
  CHR_begin_timer();//LV 3/09/2012
}
| def_inflows
;

/* Singularities attributes
 * - name
 * - pk 
 * - weir characteristics for dams
 *
 * Syntax : upstream_point;
 *          downstream_dam {
 * 	                     pk = [km] 100. 
 *	                     z(t) = { [d]  [m] 
 *		                      0.0 17.40
 *	                            }
 *                         };
 */
def_singularities : LEX_SING LEX_EQUAL brace singularities brace
{
  psingtot = $4;
  HYD_create_singularities(psingtot,Simul->pchyd,Simul->poutputs);
}
;

singularities : singularity singularities
{
  $$ = HYD_chain_singularities($1,$2);
}
| singularity
{
  $$ = $1;
}
;

singularity : read_sing_name sing_attributes 
{
  $$ = new_singularity(); //curieux non ?
  $$ = psing;
  LP_printf(Simul->poutputs,"The singularity %s has been read\n",psing->name);
  
}
;

read_sing_name : LEX_NAME 
{
  psing = HYD_create_singularity($1,Simul->pchyd);
  nworks = 0;
  free($1);
}
;

sing_attributes : brace sequence_sing_attributes brace LEX_SEMI_COLON 
| LEX_SEMI_COLON 
;

sequence_sing_attributes : sing_attribute sequence_sing_attributes
| sing_attribute 
;

sing_attribute : LEX_PK LEX_EQUAL mesure
{
  psing->pk = $3;
  psing->passage = YES_TS;
}
| LEX_TYPE LEX_EQUAL LEX_BC_TYPE
{
  psing->type = $3;
}
| LEX_POSITION LEX_EQUAL LEX_DIR LEX_INT
{
  if (psing->BC_char == NULL)
    psing->BC_char = ((s_BC_char_hyd **) calloc(1,sizeof(s_BC_char_hyd *)));
  psing->BC_char[0]->position[0] = $3;
  psing->BC_char[0]->position[1] = $4 - 1;
}
| LEX_HOLLER LEX_EQUAL flottant
{
   psing->holler = $3;
   psing->formule_rea = HOLLER;
}
| fion
{
  if (psing->BC_char == NULL)
    psing->BC_char = ((s_BC_char_hyd **) calloc(1,sizeof(s_BC_char_hyd *)));
  psing->BC_char[0] = $1;
  psing->nworks = 1;
}
| hydstructures
{
  pworktot = $1;
  psing->nworks = nworks;
  psing->BC_char = ((s_BC_char_hyd **) calloc(nworks,sizeof(s_BC_char_hyd *)));
  psing->BC_char[0] = pworktot;

  for (i = 1; i < nworks; i++) {
    pworktot = pworktot->next;
    psing->BC_char[i] = pworktot;
  }
}
;

hydstructures : LEX_HYDWORKS LEX_EQUAL brace hydworks brace
{
  $$ = $4;
}
;

hydworks : hydwork hydworks
{
  $$ = HYD_chain_hydworks($1,$2);
}
| hydwork
{
  $$ = $1;
}
;

hydwork : intro_work work_atts brace
{
  $$ = pwork;
}
;


intro_work : LEX_WORK LEX_EQUAL brace
{
  worktype = $1;
  pwork = HYD_create_hydwork(worktype);
  nworks++;
}
;

work_atts : work_att work_atts
| work_att
;

work_att : LEX_WORK_PARAM LEX_EQUAL mesure
{
  pwork->hydwork_param[$1] = $3;
}
| LEX_POSITION LEX_EQUAL LEX_DIR LEX_INT
{
  pwork->position[0] = $3;
  pwork->position[1] = $4 - 1;
}
| fion
{
  //LP_printf(Simul->poutputs,"Dans le fion !! \n");//NF 25/3/2019 For transfer, not a very clear message !
  pwork->fion_type = $1->fion_type;
  pwork->fion = $1->fion;
  pwork->fion_t = TS_function_value_t(Simul->chronos->t[INI],pwork->fion,Simul->poutputs);//LV nov2014
}
;

fion : LEX_FION LEX_EQUAL brace read_units f_ts brace
{
  $$ = HYD_create_hydwork(NONE);
  $$->fion_type = $1;
  $$->fion = $5;
}
| LEX_FION LEX_EQUAL brace read_units date_f_ts brace
{
  $$ = HYD_create_hydwork(NONE);
  $$->fion_type = $1;
  $$->fion = $5;
}
;


/* Reaches' characteristics
 * - name of the upsteam limit
 * - name of the downstream limit
 * - Strickler coefficients depending on the discharge
 * 
 * Syntax : upstream_name -> downstream_name 1 {
 *                                               min_radius = [m] 1000 
 *                                               strickler = 40.
 *                                             };	  
 */
def_reaches : LEX_REACH LEX_EQUAL brace reaches brace
{
  preachtot = $4;
  HYD_create_reaches(preachtot,Simul->pchyd);
}
;

reaches : reach reaches
{
  $$ = HYD_chain_reaches($1,$2);
}
| reach
{
  $$ =$1;
}
;

reach : read_reach_name reach_attributes
{
  $$ = new_reach();
  $$ = preach;
  LP_printf(Simul->poutputs,"The reach %s %s %d has been read\n",
	  preach->limits[ONE]->name,preach->limits[TWO]->name,preach->branch_nb);
  
}
;

read_reach_name : LEX_NAME LEX_ARROW LEX_NAME LEX_INT 
{
  preach = HYD_create_reach($1,$3,$4,Simul->pchyd,Simul->poutputs); // doit envoyer pschd car cherche singularity à reach
  free($1);
  free($3);
}
;

reach_attributes : brace sequence_reach_attributes brace LEX_SEMI_COLON 
| LEX_SEMI_COLON 
;

sequence_reach_attributes : reach_attribute sequence_reach_attributes
| reach_attribute 
;

reach_attribute : strickler
;

/*strickler : LEX_STRICKLER LEX_EQUAL mesure 
{
  preach->strickler = $3;
}
;*/

//LV 05/06/2012 : pour que ce soit décrit comme dans ProSe
strickler : intro_strickler f_ts brace
{
  preach->strickler = $2;
  //preach->strickler->prev = NULL;//LV 05/07/2012
}
;

intro_strickler : LEX_STRICKLER LEX_EQUAL brace
{
  nbsteps = 0;//Il n'y a pas d'unités à lire...
  unit_t = 1.0; // SW 29/01/2018 strickler en fonction du Q m3/s
  unit_f = 1.0;
}
;

/* Bathymetry
 * An element is created for each defined cross_section.
 * Cross_Sections' characteristics :
 * - name
 * - reach (upstream limit and branch nb)
 * - type of description : Absc Z or X Y Z
 *   - X Y Z description (X Y in Lambert coordinates)
 *     calculation of the transversal cross_section with the values of X and Y
 *     the first point is taken as the referance
 *     (function transversal_abscisse()).
 *   - ABSC Z description 
 *     ABSC = curvilinear distance from a reference point
 * - Length of the cross_section. A length of zero is only used for the 
 *   interpolation of the sections that are not used as calculation 
 *   elements (usually at the end of a reach)
 *  
 * Synthaxe : cross_section_name <- limit_name i {
 *                         dx = 700 [m]
 *                         type  = ABSC	 Z
 *  		                   0.0	 20.97
 *  		                   10.7	 15.64
 *  		                   60.8	 15.64
 *  		                   71.5	 20.97
 *                         } ;
 */
bathymetry : LEX_PROFILES LEX_EQUAL brace cross_sections brace
{
  $$ = $4;
}
| LEX_PROFILES LEX_EQUAL brace cross_sections shape brace
{
  $$ = $4;
}
| LEX_TRAJECTORY LEX_EQUAL brace trajectory brace
{
  FILE *fp;
  fp=Simul->poutputs;
  $$ = $4;
  if (Simul->pchyd->counter->nreaches > 1) {
    LP_error(fp,"There cannot be more than one reach if a trajectory shape is given\n");
  }
  if (Simul->pchyd->p_reach[0]->nfaces > 1) {
    HYD_calculate_faces_traj($4,Simul->pchyd);
  }
  else {
    LP_error(fp,"The trajectory shape must be defined by more than one point\n");
  }
}
;

cross_sections : cross_section cross_sections
{
  $$ = HYD_chain_faces($1,$2);
}
| cross_section
{
  $$ = $1;
}
;

cross_section : read_geometry_name brace cross_section_options brace  LEX_SEMI_COLON
{
  $$ = new_face();
  $$ = pface;
  pface = NULL;
}
;

read_geometry_name : LEX_NAME LEX_REVERSE_ARROW LEX_NAME LEX_INT 
{
  preach_bat = HYD_find_reach($3,$4,Simul->pchyd);
  preach_bat->nfaces++;
  free($3);
  pface = HYD_create_face(X_HYD,Simul->pchyd);
  pface->name = $1;
  pface->def = RAW_SECTION;
  pface->reach = preach_bat;
  pface->geometry = HYD_create_geometry();
  LP_printf(Simul->poutputs,"The cross-section %s in reach %s %s %d has been read\n",
	  $1,preach_bat->limits[ONE]->name,preach_bat->limits[TWO]->name,preach_bat->branch_nb);
  
}
; 

cross_section_options : cross_section_option cross_section_options
| cross_section_option 
;

cross_section_option : LEX_GENERALP LEX_EQUAL mesure 
{  
  if ($1 == DX) {
    pface->geometry->dx = $3;
    preach_bat->length += $3;
  }

  else if ($1 == CURVATURE) {
    pface->description->curvature = $3;
  }
}
| intro_sectionAbscZ sectionAbscZs 
{
  pptAbscZ = $2;
  HYD_create_ptsAbscZ(pface,pptAbscZ);
  // !!!!! il faut libérer pptAbscZ car recrée tableau dans HYD_create_ptsAbscZ(pface,pptAbscZ); --> fuite memoire !!!
}
| intro_sectionXYZ sectionXYZs 
{
  pptXYZ = $2;
  HYD_create_ptsXYZ(pface,pptXYZ);
  HYD_calculate_AbscZ_pts(pface->geometry);
// !!!!! J'imagine idem plus haut il faut libérer pptAbscZ car recrée tableau dans HYD_create_ptsAbscZ(pface,pptAbscZ); --> fuite memoire !!!
}
;

intro_sectionAbscZ : LEX_TYPE LEX_EQUAL LEX_ABSC LEX_HYD_VAR
{
  pface->geometry->type = ABSC_Z;
}
;

intro_sectionXYZ :LEX_TYPE LEX_EQUAL LEX_X LEX_Y LEX_HYD_VAR
{
  pface->geometry->type = X_Y_Z;
}
;

sectionAbscZs : sectionAbscZ sectionAbscZs 
{
  $$ = HYD_chain_ptsAbscZ($1,$2);
}
| sectionAbscZ 
{
  $$ = $1; 
} 
;
		

sectionXYZs : sectionXYZ sectionXYZs 
{
  $$ = HYD_chain_ptsXYZ($1,$2);
}
| sectionXYZ
{
  $$ = $1; 
} 
;
		
sectionAbscZ : flottant flottant 
{ 
   $$ = HYD_create_ptAbscZ($1,$2);
   pface->geometry->npts++;
} 
;

sectionXYZ : flottant flottant flottant 
{ 
  $$ = HYD_create_ptXYZ($1,$2,$3);
  pface->geometry->npts++;
} 
;

shape : LEX_SHAPE LEX_EQUAL LEX_CROSS_SECTION_TYPE
{
  Simul->pchyd->settings->cs_shape = $3;
}
;

trajectory : traj_units traj_points
{
  $$ = $2;
}
;

traj_units : a_unit a_unit a_unit a_unit a_unit
{
  unit_x = $1;
  unit_y = $2;
  unit_z = $3;
  unit_l = $4;
  unit_kappa = $5;
}
;

traj_points : traj_point traj_points
{
  $$ = HYD_chain_faces($1,$2);
}
| traj_point
{
  $$ = $1;
}
;

traj_point : flottant flottant flottant flottant flottant
{
  preach_bat = Simul->pchyd->p_reach[0];
  preach_bat->nfaces++;

  pface = HYD_create_face(X_HYD,Simul->pchyd);
  pface->def = RAW_SECTION;
  pface->reach = preach_bat;
  pface->geometry = HYD_create_geometry();
  pface->description->xc = $1 * unit_x;
  pface->description->yc = $2 * unit_y;
  pface->description->Zbottom = $3 * unit_z;
  pface->description->l = $4 * unit_l;
  pface->description->curvature = $5 * unit_kappa;
  $$ = pface;
}
;


/* Inflows' characteristics : 
 * - TYPE : point or diffuse
 * - name
 * - upstream reach limit
 * - number of the river branche i
 * - distance x from upstream limit
 * - length dx on which the inflow occurs (diffuse inflow)
 * - transversal localization y (0=left bank to 1=right bank)
 * - discharge
 *
 * Syntax : name_inf1 <- name_singularity1 i1 {
 *                           type = UPSTREAM_INFLOW
 *                           q = { [d] [m^3/s] 0. 60.0 }
 *                   };
 *          name_inf2 <- name_singularity2 i2 {
 *                           type = PT_INFLOW
 *                           x = [km] 10.
 *                           q = { [d] [m^3/s] 0. 60.0 }
 *                   };                    
 *          name_inf3 <- name_singularity3 i3 {
 *                           type = DIFFUSE_INFLOW
 *                           x = [km] 10.
 *                           dx = [km] 30.
 *                           q = { [d] [m^3/s] 0. 60.0 }
 *                   };                    
 */
def_inflows : LEX_INFLOWS LEX_EQUAL brace  inflows brace
{
  // HYD_create_inflows($4,Simul->pchyd,nspecies,Simul->poutputs);//NF 7/10/2020 nspecies is required in the last version of libhyd. Watch how T is handled. Included in nspecies or a special status. If it is the case how is it handled in the inflow structure if any ?
  HYD_create_inflows($4,Simul->pchyd,0,Simul->poutputs);//NF 7/10/2020 nspecies is required in the last version of libhyd. Watch how T is handled. Included in nspecies or a special status. If it is the case how is it handled in the inflow structure if any ?//NF 12/10/2020 discussion with SW, NG said to put at 0, this is only for cawaqs
  // loop for all sections (Simul->pchyd->p_reach[nr]->face[Y_HYD][side]->pt_inflows[nt]->app_bio[nspecies][subs]) in all particles
};

inflows : inflow inflows
{
  $$ = HYD_chain_inflows($1,$2);
}
| inflow
{
  $$ = $1;
}
;

inflow : read_name_inflow brace options_inflow brace LEX_SEMI_COLON
{
  //$$ = new_inflow();//NF 7/10/2020 This line creates a memory leak because the object is lost next line
  $$ = pinflow;
  pinflow = NULL;
}
;

read_name_inflow : LEX_NAME LEX_COLON LEX_NAME LEX_INT
{
  preach_infl = HYD_find_reach($3,$4,Simul->pchyd);
  preach_infl->ninflows++;
  free($3);
  pinflow = HYD_create_inflow($1,Simul->pchyd);
  pinflow->reach = preach_infl;
  
	
  LP_printf(Simul->poutputs,"The inflow %s in reach %s %s %d has been read\n",
	  $1,preach_infl->limits[ONE]->name,preach_infl->limits[TWO]->name,preach_infl->branch_nb);
  
  
 free($1);
} 
;

options_inflow : option_inflow options_inflow
| option_inflow 
;

option_inflow : LEX_TYPE LEX_EQUAL LEX_INFLOW_TYPE 
{
  pinflow->type = $3;
  if (pinflow->type == DIFFUSE_INFLOW)
    pinflow->diff_inflow = new_diffuse_inflow();
  else {
    pinflow->pt_inflow = new_pt_inflow();
    pinflow->pt_inflow->name = pinflow->name;
     pinflow->pt_inflow->temperature = NULL; // SW 06/10/2020 a voir, pourquoi pinflow->pt_inflow->temperature != NULL parfois apr�s new_pt_inflow()
    //if(strcmp(pinflow->name,"fresnes_choisy") == 0)
	//printf("debug");
	for (e = 0; e < NSPECIES; e++) 
    {
	 pinflow->pt_inflow->app_bio[e] = (s_ft **)malloc(Simul->counter_bio->nsubspecies[e]*sizeof(s_ft*));
        for(i = 0; i < Simul->counter_bio->nsubspecies[e]; i++)
		  //pinflow->pt_inflow->app_bio[e][i] = (s_ft *)malloc(sizeof(s_ft));
                    pinflow->pt_inflow->app_bio[e][i] = TS_init_ft(); // SW 06/10/2020 maybe memory leak
                  
		  
    }
	/* MH : 10/09/2021 memory allocation for macrospecies and subspecies */
	// allocate NMACROSPECIES flow_in_marcospe pinflow->pt_inflow->flow_in_macrospeciescies

	 pinflow->pt_inflow->flow_in_macrospecies = (s_ft **)malloc(NMACROSPECIES*sizeof(s_ft*));
	 for(e = 0; e < NMACROSPECIES; e++)
	   pinflow->pt_inflow->flow_in_macrospecies[e] = NULL;
  }
}
| LEX_X LEX_EQUAL mesure 
{
  pinflow->x = $3;
}
| LEX_GENERALP LEX_EQUAL mesure
{
  pinflow->diff_inflow->dx = $3;
}
| LEX_Y LEX_EQUAL flottant 
{
  pinflow->transversal_position = $3;
}
| LEX_HYD_VAR LEX_EQUAL brace read_units f_ts brace 
{
  pinflow->discharge = $5;
  if (pinflow->type != DIFFUSE_INFLOW)
    pinflow->pt_inflow->discharge = $5;
}
| LEX_HYD_VAR LEX_EQUAL brace read_units date_f_ts brace 
{
  pinflow->discharge = $5;
  if (pinflow->type != DIFFUSE_INFLOW)
    pinflow->pt_inflow->discharge = $5;
}
| LEX_TEMP_INFLOW LEX_EQUAL brace read_units f_ts brace 
{
  pinflow->temperature = $5;
  if (pinflow->type != DIFFUSE_INFLOW)
  {
    pinflow->pt_inflow->temperature = $5;

      /* SW 15/03/2021 to find CODE_TS and replace it by default value */
      if((int)pinflow->pt_inflow->temperature->ft == CODE_TS)
      {
         pinflow->pt_inflow->temperature->ft = Simul->default_t_inflows;
         //LP_printf(Simul->poutputs, "inflowsssss %s T = %f\n", pinflow->pt_inflow->name, pinflow->pt_inflow->temperature->ft);
      } 
     while(pinflow->pt_inflow->temperature->next != NULL)
      {
         pinflow->pt_inflow->temperature = pinflow->pt_inflow->temperature->next;
         if((int)pinflow->pt_inflow->temperature->ft == CODE_TS)
             pinflow->pt_inflow->temperature->ft = Simul->default_t_inflows;
      }
      
      pinflow->pt_inflow->temperature = TS_browse_ft(pinflow->pt_inflow->temperature, BEGINNING_TS);
  }
}
| LEX_TEMP_INFLOW LEX_EQUAL brace read_units date_f_ts brace 
{
  pinflow->temperature = $5;
  if (pinflow->type != DIFFUSE_INFLOW)
  {
    pinflow->pt_inflow->temperature = $5;

      /* SW 15/03/2021 to find CODE_TS and replace it by default value */
      if((int)pinflow->pt_inflow->temperature->ft == CODE_TS)
      {
         pinflow->pt_inflow->temperature->ft = Simul->default_t_inflows;
         //LP_printf(Simul->poutputs, "inflowsssss %s T = %f\n", pinflow->pt_inflow->name, pinflow->pt_inflow->temperature->ft);
      }
     
      while(pinflow->pt_inflow->temperature->next != NULL)
      {
         pinflow->pt_inflow->temperature = pinflow->pt_inflow->temperature->next;
         if((int)pinflow->pt_inflow->temperature->ft == CODE_TS)
             pinflow->pt_inflow->temperature->ft = Simul->default_t_inflows;
      }
      
      pinflow->pt_inflow->temperature = TS_browse_ft(pinflow->pt_inflow->temperature, BEGINNING_TS);
  }
}
| LEX_ONESPECIES LEX_INT LEX_EQUAL brace read_units f_ts brace
{
   e = $1;
   j = $2;
   LP_printf(Simul->poutputs,"apport %s %d\n",name_species(e),j);
   if(j <= Simul->counter_bio->nsubspecies[e])
   {
      pinflow->pt_inflow->app_bio[e][j-1] = $6;
      
      /* SW 15/03/2021 to find CODE_TS and replace it by default value */
      if((int)pinflow->pt_inflow->app_bio[e][j-1]->ft == CODE_TS)
         pinflow->pt_inflow->app_bio[e][j-1]->ft = pspec2[WATER][e][j-1]->default_C_inflows;
      while(pinflow->pt_inflow->app_bio[e][j-1] ->next != NULL)
      { 
         pinflow->pt_inflow->app_bio[e][j-1] = pinflow->pt_inflow->app_bio[e][j-1]->next;
         if((int)pinflow->pt_inflow->app_bio[e][j-1]->ft == CODE_TS)
             pinflow->pt_inflow->app_bio[e][j-1]->ft = pspec2[WATER][e][j-1]->default_C_inflows;
      }
      
      pinflow->pt_inflow->app_bio[e][j-1] = TS_browse_ft(pinflow->pt_inflow->app_bio[e][j-1], BEGINNING_TS);	 
   }
	  
   else
     LP_error(Simul->poutputs, "No species defined for %s %d",name_species(e),j);
   	 
}
| LEX_ONESPECIES LEX_INT LEX_EQUAL brace read_units date_f_ts brace
{
   e = $1;
   j = $2;
   LP_printf(Simul->poutputs,"apport %s %d\n",name_species(e),j);
   if(j <= Simul->counter_bio->nsubspecies[e])
   {
      pinflow->pt_inflow->app_bio[e][j-1] = $6;

      /* SW 15/03/2021 to find CODE_TS and replace it by default value */
      if((int)pinflow->pt_inflow->app_bio[e][j-1]->ft == CODE_TS)
         pinflow->pt_inflow->app_bio[e][j-1]->ft = pspec2[WATER][e][j-1]->default_C_inflows;
      while(pinflow->pt_inflow->app_bio[e][j-1] ->next != NULL)
      {
         pinflow->pt_inflow->app_bio[e][j-1] = pinflow->pt_inflow->app_bio[e][j-1]->next;
         if((int)pinflow->pt_inflow->app_bio[e][j-1]->ft == CODE_TS)
             pinflow->pt_inflow->app_bio[e][j-1]->ft = pspec2[WATER][e][j-1]->default_C_inflows;
      }
      
      pinflow->pt_inflow->app_bio[e][j-1] = TS_browse_ft(pinflow->pt_inflow->app_bio[e][j-1], BEGINNING_TS);		 
   }
	  
   else
     LP_error(Simul->poutputs, "No species defined for %s %d",name_species(e),j);
   	 
}

| LEX_TYPEOF_MACROSPECIES LEX_EQUAL brace read_units f_ts brace
//MH : 13/09/2021 reading of TOC and storing it at pt_inflow in struct_hydro.h
{
   e = $1;
   // allocate flow_in_macrospecies[$1]

   pinflow->pt_inflow->flow_in_macrospecies[$1] = TS_init_ft(); 

   LP_printf(Simul->poutputs,"READING OF %s \n",RIV_name_macrospecies($1,Simul->poutputs));
   
   if (e <= NMACROSPECIES)
   {
      pinflow->pt_inflow->flow_in_macrospecies[e] = $5;
      LP_printf(Simul->poutputs,"TOC:  %3.2f %3.2f \n",pinflow->pt_inflow->flow_in_macrospecies[$1]->t/86400,pinflow->pt_inflow->flow_in_macrospecies[$1]->ft/83.3);
      // LP_printf(Simul->poutputs,"TOC:  %3.2f %3.2f \n",pinflow->pt_inflow->flow_in_macrospecies[$1]->next->next->t/86400,pinflow->pt_inflow->flow_in_macrospecies[$1]->next->next->ft/83.3);
 
      pinflow->pt_inflow->flow_in_macrospecies[e] = TS_browse_ft(pinflow->pt_inflow->flow_in_macrospecies[e], BEGINNING_TS);	 
   
   LP_printf(Simul->poutputs,"END OF reading %s \n",RIV_name_macrospecies($1,Simul->poutputs));
   }
   else
    LP_error(Simul->poutputs, "No species defined for %s",RIV_name_macrospecies($1,Simul->poutputs));

   // function to calculate MOD1 ...MOP3 and store in pt_inflow struct of pinflow using val and params

   if (Simul->p_macrospecies[TOC]->degradOrgMat[MACMOD][B]->varying_YesorNo == YES_RIVE) // SW 08/01/2023 replace YES by YES_RIVE, varying_YesorNo is initialized with NO_RIVE
     {
       PROSE_toc_to_mod_mop_fract_var_b1($1, Simul, pinflow, Simul->poutputs); 
       
       LP_printf(Simul->poutputs,"varying b1 function implemented");
     }
   else {
     PROSE_toc_to_mod_mop_fract($1, Simul, pinflow, Simul->poutputs);
     LP_printf(Simul->poutputs,"const b1 function implemented");

        }
     

   
}

;

set_transport : LEX_TRANSPORT LEX_EQUAL brace transport_atts brace
{
  int ind, var, ns_ttc,nele, ind_phy,nb_var;
  nele = Simul->pchyd->counter->nele_tot;
  Simul->p_phy_species = (s_species_ttc ****)malloc(Simul->passim->N_particules * sizeof(s_species_ttc ***));

  nb_var=nspecies;//NF 11/10/2020
  //if(Simul->calc_mode[H_T]==YES_TS)//NF 11/10/2020//NF 12/10/2020 would be useful if we would have only pcarac_ttc
      //nb_var++;//NF 11/10/2020
  
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  ns_ttc = 0;
  Simul->pcarac_ttc[i]->regime = Simul->regime; 
  Simul->pcarac_ttc[i]->count[NSPECIES_TTC] = nspecies;
  Simul->pcarac_ttc[i]->count[NELE_TTC]=nele;
  if(Simul->pcarac_ttc[i]->theta == 0)
     Simul->pcarac_ttc[i]->theta = THET_TTC;
  
  Simul->pcarac_ttc[i]->p_species = (s_species_ttc **)malloc(nb_var*sizeof(s_species_ttc*));//NF 11/10/2020
  Simul->p_phy_species[i] = (s_species_ttc ***)malloc(Simul->counter_bio->nsubspecies[PHY]*sizeof(s_species_ttc **));
  
  for(ind = 0; ind < Simul->counter_bio->nsubspecies[PHY]; ind++)
     Simul->p_phy_species[i][ind] = ((s_species_ttc **)calloc(3,sizeof(s_species_ttc*))); // 3 compartiments

  for(var = 0; var < NSPECIES; var++){
    if(pspecies[var] != NULL){
    for(ind = 0; ind < Simul->counter_bio->nsubspecies[var]; ind++){
	  Simul->pcarac_ttc[i]->p_species[ns_ttc] = TTC_init_species(nele);
	  Simul->pcarac_ttc[i]->p_species[ns_ttc]->name = pspecies[var][ind]->name;
	  Simul->pcarac_ttc[i]->p_species[ns_ttc]->type = pspecies[var][ind]->type;
	  Simul->pcarac_ttc[i]->p_species[ns_ttc]->media_type = pspecies[var][ind]->media_type;
	  Simul->pcarac_ttc[i]->p_species[ns_ttc]->calc_process[DIFFUSION_TTC] = pspecies[var][ind]->calc_process[DIFFUSION_TTC];
	  Simul->pcarac_ttc[i]->p_species[ns_ttc]->theta=Simul->pcarac_ttc[i]->theta;
      Simul->pcarac_ttc[i]->p_species[ns_ttc]->calc_process[CONVECTION_TTC] = pspecies[var][ind]->calc_process[CONVECTION_TTC];
      Simul->pcarac_ttc[i]->p_species[ns_ttc]->pchronos=CHR_copy_chronos(Simul->chronos);
      ns_ttc++;
	  
	  if(var == PHY){
	  for(ind_phy = 0; ind_phy < 3; ind_phy++)
	  {
	     Simul->p_phy_species[i][ind][ind_phy] = TTC_init_species(nele);
		 Simul->p_phy_species[i][ind][ind_phy]->name = pspecies[var][ind]->name;
		 Simul->p_phy_species[i][ind][ind_phy]->type = pspecies[var][ind]->type;
		 Simul->p_phy_species[i][ind][ind_phy]->media_type = pspecies[var][ind]->media_type;
		 Simul->p_phy_species[i][ind][ind_phy]->calc_process[DIFFUSION_TTC] = pspecies[var][ind]->calc_process[DIFFUSION_TTC];
		 Simul->p_phy_species[i][ind][ind_phy]->calc_process[CONVECTION_TTC] = pspecies[var][ind]->calc_process[CONVECTION_TTC];
         Simul->p_phy_species[i][ind][ind_phy]->theta = Simul->pcarac_ttc[i]->theta;
		 Simul->p_phy_species[i][ind][ind_phy]->pchronos=CHR_copy_chronos(Simul->chronos);
		 }
		 }
      //free(pspecies[var][ind]->pgc->b);
	  //free(pspecies[var][ind]->pgc);
	  //pspecies[var][ind]->pgc = NULL;
	  
	}
	}
  }	 
  }

  LP_printf(Simul->poutputs,"Transport simulation configured: %d species\n",ns_ttc);
}
; 
transport_atts : transport_att transport_atts
| transport_att
;

transport_att : LEX_THETA LEX_EQUAL flottant
{
  for(i = 0; i < Simul->passim->N_particules; i++)
      Simul->pcarac_ttc[i]->theta=$3;
  LP_printf(Simul->poutputs,"Theta = %f \n",Simul->pcarac_ttc[0]->theta);
};

set_biology : intro_set_biology biology_atts LEX_CLOSING_BRACE
{
  
   LP_printf(Simul->poutputs,"Biology simulation configured\n");
}
;

intro_set_biology : LEX_BIOLOGY LEX_EQUAL LEX_OPENING_BRACE
{
    
	nsections = Simul->pchyd->counter->nele_tot;
	for(i = 0; i < Simul->passim->N_particules; i++)
	{
    Simul->psimul_bio[i] = (s_simul **)calloc(nsections,sizeof(s_simul *));
	LP_printf(Simul->poutputs,"%d sections defined for biology\n", nsections);
    for(ns = 0; ns < nsections; ns++){
        Simul->psimul_bio[i][ns] = init_simulation();
		free(Simul->psimul_bio[i][ns]->counter);
		Simul->psimul_bio[i][ns]->poutputs = Simul->poutputs;
		Simul->psimul_bio[i][ns]->counter = Simul->counter_bio;
		if(Simul->psimul_bio[i][ns]->settings != NULL)
		  free(Simul->psimul_bio[i][ns]->settings);
		Simul->psimul_bio[i][ns]->settings = Simul->settings;
		Simul->psimul_bio[i][ns]->section = init_section();		
		Simul->psimul_bio[i][ns]->section->id_section = ns + 1; //corresponding to id_abs_ele pchyd
		init_annex_var(Simul->psimul_bio[i][ns]->section, Simul->psimul_bio[i][ns]);
		//Simul->psimul_bio[ns]->section->nsublayers[WATER] = 1;
		//Simul->psimul_bio[ns]->section->nsublayers[VASE] = 1;
		//Simul->psimul_bio[ns]->section->compartments[WATER] = (s_compartment **)calloc(1,sizeof(s_compartment *));
		//Simul->psimul_bio[ns]->section->compartments[VASE] = (s_compartment **)calloc(1,sizeof(s_compartment *));
		//Simul->psimul_bio[ns]->section->compartments[WATER][0] = init_compartment(WATER,1,Simul->psimul_bio[ns]->section, Simul->settings->dbo_oxy);
		//Simul->psimul_bio[ns]->section->compartments[WATER][0]->num = 1;
		//Simul->psimul_bio[ns]->section->compartments[WATER][0]->down = (s_interface **)calloc(1,sizeof(s_interface *));
		//Simul->psimul_bio[ns]->section->compartments[WATER][0]->down[0] = init_interface();
		//Simul->psimul_bio[ns]->section->compartments[WATER][0]->n_int_down = 1;
		//Simul->psimul_bio[ns]->section->compartments[WATER][0]->down[0]->upper_comp = Simul->psimul_bio[ns]->section->compartments[WATER][0];
		//Simul->psimul_bio[ns]->section->compartments[VASE][0] = init_compartment(VASE,1,Simul->psimul_bio[ns]->section,Simul->settings->dbo_oxy);
		//Simul->psimul_bio[ns]->section->compartments[VASE][0]->num = 1;
		//Simul->psimul_bio[ns]->section->compartments[WATER][0]->down[0]->lower_comp = Simul->psimul_bio[ns]->section->compartments[VASE][0];        
		}
		}
  // LP_printf(Simul->poutputs,"%d sections defined for biology\n", nsections);
}
;

biology_atts : biology_att biology_atts
| biology_att
;

biology_att : intro_layers layers LEX_CLOSING_BRACE
| meteo
| exchanges
| LEX_NUM LEX_EQUAL LEX_METHOD
{
  char *name_m;
  name_m = name_method($3);
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  if ($3 == EXPLICIT) {
   for(ns = 0; ns < nsections; ns++){
    Simul->psimul_bio[i][ns]->numerical_method->name = name_m;
    Simul->psimul_bio[i][ns]->numerical_method->kmax = $3;
    Simul->psimul_bio[i][ns]->numerical_method->coef_RK[0] = 1.;
    Simul->psimul_bio[i][ns]->numerical_method->coef_RK[1] = Simul->psimul_bio[i][ns]->numerical_method->coef_RK[2] = Simul->psimul_bio[i][ns]->numerical_method->coef_RK[3] = 0.;
    Simul->psimul_bio[i][ns]->numerical_method->dt_RK[0] = 1.;
    Simul->psimul_bio[i][ns]->numerical_method->dt_RK[1] = Simul->psimul_bio[i][ns]->numerical_method->dt_RK[2] = Simul->psimul_bio[i][ns]->numerical_method->dt_RK[3] = 0.;
  }
  }
  }
  free(name_m);
  LP_printf(Simul->poutputs,"numerical method = %s\n",Simul->psimul_bio[0][0]->numerical_method->name);
  //printf("numerical method = %s\n",Simul->numerical_method->name);
  LP_printf(Simul->poutputs,"kmax = %d\n",Simul->psimul_bio[0][0]->numerical_method->kmax);
  //printf("kmax = %d\n",Simul->numerical_method->kmax);
}
| LEX_DIVDT LEX_EQUAL flottant
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++)
   Simul->psimul_bio[i][ns]->numerical_method->max_div_dt = $3;
   }
  LP_printf(Simul->poutputs,"max_div_dt = %f\n",Simul->psimul_bio[0][0]->numerical_method->max_div_dt);
  
}
| LEX_DZ_RIVE LEX_EQUAL mesure
{
  temp_var = $3;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {  
  for(ns = 0; ns < nsections; ns++)
     Simul->psimul_bio[i][ns]->settings->dz = temp_var;
	 }
  LP_printf(Simul->poutputs,"dz = %f\n",temp_var);
} 
;

intro_layers : LEX_LAYERS LEX_EQUAL LEX_OPENING_BRACE
;

layers : layer layers
| layer
;

layer : intro_layer sub_layers
{
  //LP_printf(Simul->poutputs,"\n\t%d sublayers of %s read\n",sect2->nsublayers[lay],name_layer(lay));
  //printf("\n\t%d sublayers of %s read\n",sect2->nsublayers[lay],name_layer(lay));
  
  //Simul->psimul_bio[ns]->section->compartments[lay] = (s_compartment **)calloc(sect2->nsublayers[lay],sizeof(s_compartment *));
  //sect2->compartments[lay][0] = $2;
  //for (i = 1; i < sect2->nsublayers[lay]; i++)
    //sect2->compartments[lay][i] = sect2->compartments[lay][i-1]->next;
 for(i = 0; i < Simul->passim->N_particules; i++)
  {
	for(ns = 0; ns < nsections; ns++){
  sect2 = Simul->psimul_bio[i][ns]->section;
  if (sect2->compartments[WATER] != NULL) {  
    init_O2(sect2->compartments[WATER][0]->pspecies[O2][0]->dissolved->gas);
    init_N2O(sect2->compartments[WATER][0]->pspecies[N2O][0]->dissolved->gas);
  }//LV 05/05/2011 rajouté sinon si l'utilisateur définit l'oxygène les variable pour calcul de la saturation ne sont pas définies
  }
  }
}
;

intro_layer : LEX_LAYER
{
  lay = $1;
  //sect2->nsublayers[lay] = 0.;
  LP_printf(Simul->poutputs,"\n\tlayers of %s\n",name_layer(lay));
  //printf("\n\tlayers of %s\n",name_layer(lay));
  
	nsections = Simul->pchyd->counter->nele_tot;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
    for(ns = 0; ns < nsections; ns++){
        Simul->psimul_bio[i][ns]->section->nsublayers[lay] = 1;
		Simul->psimul_bio[i][ns]->section->compartments[lay] = (s_compartment **)calloc(1,sizeof(s_compartment *));
        Simul->psimul_bio[i][ns]->section->compartments[lay][0] = init_compartment(lay,1,Simul->psimul_bio[i][ns]->section, Simul->settings->dbo_oxy);		
        Simul->psimul_bio[i][ns]->section->compartments[lay][0]->num = 1;
		if(lay == VASE){
		Simul->psimul_bio[i][ns]->section->compartments[WATER][0]->down = (s_interface **)calloc(1,sizeof(s_interface *));
		Simul->psimul_bio[i][ns]->section->compartments[WATER][0]->down[0] = init_interface();
		Simul->psimul_bio[i][ns]->section->compartments[WATER][0]->down[0]->lower_comp = Simul->psimul_bio[i][ns]->section->compartments[VASE][0];
		Simul->psimul_bio[i][ns]->section->compartments[WATER][0]->n_int_down = 1;
		Simul->psimul_bio[i][ns]->section->compartments[WATER][0]->down[0]->upper_comp = Simul->psimul_bio[i][ns]->section->compartments[WATER][0];

		/* MH: 25/11/2021  memory for the values of random_params   */
		/* 
		Simul->psimul_bio[i][ns]->section->compartments[WATER][0]->pmacro[TOC] =  RIV_init_macrospecies(TOC,Simul->poutputs);
		Simul->psimul_bio[i][ns]->section->compartments[WATER][0]->pmacro[TOC]->degradOrgMat[MACMOD][B]=  RIV_init_stoch_param(Simul->poutputs);
		Simul->psimul_bio[i][ns]->section->compartments[VASE][0]->pmacro[TOC] =  RIV_init_macrospecies(TOC,Simul->poutputs);
		Simul->psimul_bio[i][ns]->section->compartments[VASE][0]->pmacro[TOC]->degradOrgMat[MACMOD][B]=  RIV_init_stoch_param(Simul->poutputs);


*/


		}
		}
  }
}
;

sub_layers : sub_layer sub_layers
| sub_layer
;
//{
//  $$ = chain_compartments($1,$2);
//}

//{
 // $$ = $1;
//}


sub_layer : intro_sub_layer atts_compartment LEX_CLOSING_BRACE
;
//{
  /* Initialization of the compartment's volume */ //SW 03/03/2017 why ??????
  //if (lay != WATER) {
   // if (comp2->state->mass > EPS)  // SW 03/03/2017 attention ne pas reinitialiser a 0 pour volume ou mass defini en bas
   // {
	//  if(comp2->state->volume < EPS){}
	//	//comp2->state->volume = comp2->state->mass / (comp2->state->rho *(1-comp2->state->phi)); // SW 03/03/2017 pas de fonction calc_volume_vase
	//}
   // else if(comp2->state->volume > EPS){}
	//   //comp2->state->mass = comp2->state->rho * (1-comp2->state->phi) * comp2->state->volume;
   // else 
	//{
	//  LP_printf(Simul->poutputs,"Compartiment %s 's volume and mass are zero !!! ",comp2->name);
	//}
 // }

  //$$ = comp2;
//}


intro_sub_layer : LEX_NAME LEX_EQUAL LEX_OPENING_BRACE
{
  s_compartment *pcomp;
  int np;
  //sect2->nsublayers[lay]++;
  //comp2 = init_compartment(lay,sect2->nsublayers[lay],sect2);
  //comp2->name = $1;
 
  for(np = 0; np < Simul->passim->N_particules; np++)
  { 
  for(ns = 0; ns < nsections; ns++){
    pcomp = Simul->psimul_bio[np][ns]->section->compartments[lay][0];
    sprintf(pcomp->name,"%s",$1);
	//if(lay == VASE){
	  //LP_printf(Simul->poutputs,"%s\n",pcomp->name);
	  //}
    for (e = 0; e < NSPECIES; e++) {
      //if (pspec2[lay][e][0] != NULL) {
	  if (pspec2[lay][e] != NULL) { //SW 22/02/2019 remove [j]
	  if(e < NH4) // SW 25/02/2020
		pcomp->pspecies[e] = (s_species **)calloc(Simul->counter_bio->nsubspecies[e],sizeof(s_species *));
        for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++)
       {
	   pcomp->pspecies[e][j] = copy_species(pcomp->pspecies[e][j],pspec2[lay][e][j],Simul->counter_bio,Simul->settings->nb_comp_phy,Simul->settings->dbo_oxy);
		//pspec2[lay][e][j] = NULL;
	   }		
	
    }	
  }
  for(e = 0; e < NPART; e++){
    if (pcomp->pspecies[e] != NULL) {
    for (j = 0; j < Simul->counter_bio->nsubspecies[e]; j++) {
	  	for (i = 0; i < NDISS; i++) {
	    pcomp->pspecies[e][j]->particulate->adsorbedC[i] = (double *)calloc(Simul->counter_bio->nsubspecies[NPART+i],sizeof(double));
	    for (k = 0; k < 4; k++)
	      pcomp->pspecies[e][j]->particulate->intermediate_adsorbedC[k][i] = (double *)calloc(Simul->counter_bio->nsubspecies[NPART+i],sizeof(double));
	  }
}
}
}	  

    //if ((lay == WATER) && (sect2->nsublayers[lay] == 1))
  if ((lay == WATER))
    pcomp->pspecies[O2][0]->C = CODE;
  }
  }

  free($1);
}
;

/* This part describes the compartments' attributes */

atts_compartment : att_compartment atts_compartment
| att_compartment
;

att_compartment : LEX_MASS LEX_EQUAL a_unit_f f_ts
{
  s_compartment *pcomp;
  s_ft *pf_ts;
  if(TS_length_ts($4) == nsections){
  
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  pf_ts = TS_browse_ft($4,BEGINNING_TS);
  while(pf_ts != NULL)
  {
  ns = (int)pf_ts->t;
  pcomp = Simul->psimul_bio[i][ns-1]->section->compartments[lay][0];
  pcomp->state->mass = pf_ts->ft;
  pf_ts = pf_ts->next;
  }
  }
  }
  else if(TS_length_ts($4) == 1){
    pf_ts = $4;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {	
    for(ns = 0; ns < nsections; ns++){
	   pcomp = Simul->psimul_bio[i][ns]->section->compartments[lay][0];
	   pcomp->state->mass = pf_ts->ft;
	}
  }
   LP_printf(Simul->poutputs,"All sections with a same mass %f\n",pf_ts->ft);
  }
}
| LEX_MASS LEX_EQUAL a_unit_f date_f_ts
{
  s_compartment *pcomp;
  s_ft *pf_ts;
  if(TS_length_ts($4) == nsections){
  
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  pf_ts = TS_browse_ft($4,BEGINNING_TS);
  while(pf_ts != NULL)
  {
  ns = (int)pf_ts->t;
  pcomp = Simul->psimul_bio[i][ns-1]->section->compartments[lay][0];
  pcomp->state->mass = pf_ts->ft;
  pf_ts = pf_ts->next;
  }
  }
  }
  else if(TS_length_ts($4) == 1){
    pf_ts = $4;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {	
    for(ns = 0; ns < nsections; ns++){
	   pcomp = Simul->psimul_bio[i][ns]->section->compartments[lay][0];
	   pcomp->state->mass = pf_ts->ft;
	}
  }
   LP_printf(Simul->poutputs,"All sections with a same mass %f\n",pf_ts->ft);
  }
}
| LEX_VOLUME LEX_EQUAL a_unit_f f_ts
{
  s_compartment *pcomp;
  s_ft *pf_ts;
  if(TS_length_ts($4) == nsections){
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  pf_ts = TS_browse_ft($4,BEGINNING_TS);
  while(pf_ts != NULL)
  {
  ns = (int)pf_ts->t;
  pcomp = Simul->psimul_bio[i][ns-1]->section->compartments[lay][0];
  pcomp->state->volume = pf_ts->ft;
  pf_ts = pf_ts->next;
  }
  }
  }
  else if(TS_length_ts($4) == 1){
    pf_ts = $4;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
    for(ns = 0; ns < nsections; ns++){
	   pcomp = Simul->psimul_bio[i][ns]->section->compartments[lay][0];
	   pcomp->state->volume = pf_ts->ft;
	}
  }
   LP_printf(Simul->poutputs,"All sections with a same volume %f\n",pf_ts->ft);
  }
} 
| LEX_VOLUME LEX_EQUAL a_unit_f date_f_ts
{
  s_compartment *pcomp;
  s_ft *pf_ts;
  if(TS_length_ts($4) == nsections){
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  pf_ts = TS_browse_ft($4,BEGINNING_TS);
  while(pf_ts != NULL)
  {
  ns = (int)pf_ts->t;
  pcomp = Simul->psimul_bio[i][ns-1]->section->compartments[lay][0];
  pcomp->state->volume = pf_ts->ft;
  pf_ts = pf_ts->next;
  }
  }
  }
  else if(TS_length_ts($4) == 1){
    pf_ts = $4;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
    for(ns = 0; ns < nsections; ns++){
	   pcomp = Simul->psimul_bio[i][ns]->section->compartments[lay][0];
	   pcomp->state->volume = pf_ts->ft;
	}
  }
   LP_printf(Simul->poutputs,"All sections with a same volume %f\n",pf_ts->ft);
  }
}
| LEX_RHO LEX_EQUAL mesure
{
  s_compartment *pcomp;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
    pcomp = Simul->psimul_bio[i][ns]->section->compartments[lay][0];
	pcomp->state->rho = $3;
  }
  } 
  LP_printf(Simul->poutputs,"\trho = %f g/m^3\n",$3);
  //printf("\trho = %f g/m^3\n",comp2->state->rho);
}
| LEX_PHI LEX_EQUAL flottant
{
    s_compartment *pcomp;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
    pcomp = Simul->psimul_bio[i][ns]->section->compartments[lay][0];
	pcomp->state->phi = $3;
  }
  }
  LP_printf(Simul->poutputs,"\tporosity = %f\n",$3);
  //printf("\tporosity = %f\n",comp2->state->phi);
}
| initial_concentrations
;

a_unit_f : a_unit
{
 unit_f = $1;
 unit_t = 1;
}
;

initial_concentrations : LEX_IC LEX_EQUAL LEX_OPENING_BRACE ics LEX_CLOSING_BRACE
;

ics : ic ics
|ic
;

ic : LEX_ONESPECIES LEX_INT LEX_EQUAL a_unit_f f_ts
{
  s_compartment *pcomp;
  s_ft *pf_ts;
  //LP_printf(Simul->poutputs,"%d\n",$1);
  //if($1==12)
   //LP_printf(Simul->poutputs,"%d\n",$1);
  if(TS_length_ts($5) == nsections){
     
     if (Simul->counter_bio->nsubspecies[$1] >= $2) {
  for(i = 0; i < Simul->passim->N_particules; i++)
  {	 
        pf_ts = TS_browse_ft($5,BEGINNING_TS);
        while(pf_ts != NULL){
		ns = (int)pf_ts->t;
		pcomp = Simul->psimul_bio[i][ns-1]->section->compartments[lay][0];
		if(pcomp->pspecies[$1][$2-1] != NULL)
		{
		    pcomp->pspecies[$1][$2-1]->C = pf_ts->ft;
		    pf_ts = pf_ts->next;
		}
		else
		{
          LP_warning(Simul->poutputs,"\tAn initial concentration of %s %d is given in %s but the species is not defined.\n",name_species($1),$2,name_layer(pcomp->type));
         break;
		}
		}
  }
      }  
    else
      LP_warning(Simul->poutputs,"\tAn initial concentration of %s %d is given in %s but the species is not defined.\n",name_species($1),$2,name_layer(pcomp->type));
  }
  else if(TS_length_ts($5) == 1){
  if (Simul->counter_bio->nsubspecies[$1] >= $2) {
    pf_ts = $5;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
	for(ns = 0; ns < nsections; ns++){
	pcomp = Simul->psimul_bio[i][ns]->section->compartments[lay][0];
	pcomp->pspecies[$1][$2-1]->C = pf_ts->ft;	
	}
  }
  }
  else
      LP_warning(Simul->poutputs,"\tAn initial concentration of %s %d is given in %s but the species is not defined.\n",name_species($1),$2,name_layer(pcomp->type));  
  }
  else
    LP_error(Simul->poutputs,"Length of initial concentrations is not good for species %s\n",name_species($1));
}
| LEX_ONESPECIES LEX_INT LEX_EQUAL a_unit_f date_f_ts
{
  s_compartment *pcomp;
  s_ft *pf_ts;
  //LP_printf(Simul->poutputs,"%d\n",$1);
  //if($1==12)
   //LP_printf(Simul->poutputs,"%d\n",$1);
  if(TS_length_ts($5) == nsections){
     
     if (Simul->counter_bio->nsubspecies[$1] >= $2) {
  for(i = 0; i < Simul->passim->N_particules; i++)
  {	 
        pf_ts = TS_browse_ft($5,BEGINNING_TS);
        while(pf_ts != NULL){
		ns = (int)pf_ts->t;
		pcomp = Simul->psimul_bio[i][ns-1]->section->compartments[lay][0];
		if(pcomp->pspecies[$1][$2-1] != NULL)
		{
		    pcomp->pspecies[$1][$2-1]->C = pf_ts->ft;
		    pf_ts = pf_ts->next;
		}
		else
		{
          LP_warning(Simul->poutputs,"\tAn initial concentration of %s %d is given in %s but the species is not defined.\n",name_species($1),$2,name_layer(pcomp->type));
         break;
		}
		}
  }
      }  
    else
      LP_warning(Simul->poutputs,"\tAn initial concentration of %s %d is given in %s but the species is not defined.\n",name_species($1),$2,name_layer(pcomp->type));
  }
  else if(TS_length_ts($5) == 1){
  if (Simul->counter_bio->nsubspecies[$1] >= $2) {
    pf_ts = $5;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
	for(ns = 0; ns < nsections; ns++){
	pcomp = Simul->psimul_bio[i][ns]->section->compartments[lay][0];
	pcomp->pspecies[$1][$2-1]->C = pf_ts->ft;	
	}
  }
  }
  else
      LP_warning(Simul->poutputs,"\tAn initial concentration of %s %d is given in %s but the species is not defined.\n",name_species($1),$2,name_layer(pcomp->type));  
  }
  else
    LP_error(Simul->poutputs,"Length of initial concentrations is not good for species %s\n",name_species($1));
}
;

meteo : LEX_METEO LEX_EQUAL LEX_OPENING_BRACE atts_meteo LEX_CLOSING_BRACE
;

atts_meteo : att_meteo atts_meteo
| att_meteo
;

att_meteo : LEX_TEMPERATURE LEX_EQUAL LEX_OPENING_BRACE temperature LEX_CLOSING_BRACE
| LEX_WIND LEX_EQUAL LEX_OPENING_BRACE wind LEX_CLOSING_BRACE
| LEX_RADIATION LEX_EQUAL LEX_OPENING_BRACE radiation LEX_CLOSING_BRACE
| LEX_PHOTOPERIOD LEX_EQUAL LEX_OPENING_BRACE atts_photoperiod LEX_CLOSING_BRACE
;

temperature : read_units f_ts
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->temperature->function_meteo = $2;
  Simul->psimul_bio[i][ns]->section->meteo->temperature->calc = NO_RIVE;
  }
  }
  //LP_printf(Simul->poutputs,"temperature time series defined\n");
  //printf("temperature time series defined\n");
}
| read_units date_f_ts
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->temperature->function_meteo = $2;
  Simul->psimul_bio[i][ns]->section->meteo->temperature->calc = NO_RIVE;
  }
  }
  //LP_printf(Simul->poutputs,"temperature time series defined\n");
  //printf("temperature time series defined\n");
}
| atts_temperature
;

atts_temperature : att_temperature atts_temperature
| att_temperature
;

att_temperature : LEX_MEAN LEX_EQUAL mesure
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->temperature->mean = $3;
  }
  }
  //LP_printf(Simul->poutputs,"temperature : mean = %f °C\n",$3);
  
}
| LEX_AMPLITUDE LEX_EQUAL mesure
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->temperature->amplitude = $3;
  }
  }
  //LP_printf(Simul->poutputs,"temperature : amplitude = %f °C\n",$3);

}
| LEX_DELAY LEX_EQUAL mesure
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->temperature->delay = $3;
  }
  }
  //LP_printf(Simul->poutputs,"temperature : delay = %f °C\n",$3);
  //printf("temperature : delay = %f s\n",sect2->meteo->temperature->delay);
}
;

wind : read_units f_ts
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->wind->function_meteo = $2;
  Simul->psimul_bio[i][ns]->section->meteo->wind->calc = NO_RIVE;
  }
  }
  //LP_printf(Simul->poutputs,"wind time series defined length = %d\n",TS_length_ts($2));
  //printf("wind time series defined\n");
}
| read_units date_f_ts
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->wind->function_meteo = $2;
  Simul->psimul_bio[i][ns]->section->meteo->wind->calc = NO_RIVE;
  }
  }
  //LP_printf(Simul->poutputs,"wind time series defined length = %d\n",TS_length_ts($2));
  //printf("wind time series defined\n");
}

;

radiation : read_units f_ts
{
  s_ft *pft_st;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {  
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->radiation->function_meteo = $2;
  Simul->psimul_bio[i][ns]->section->meteo->radiation->calc = NO_RIVE;
  }
  }
  pft_st = TS_browse_ft($2,END_TS);
  //LP_printf(Simul->poutputs,"radiation time series defined length = %d\n",TS_length_ts($2));
  //printf("radiation time series defined\n");
}
| read_units date_f_ts
{
  s_ft *pft_st;
  for(i = 0; i < Simul->passim->N_particules; i++)
  {  
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->radiation->function_meteo = $2;
  Simul->psimul_bio[i][ns]->section->meteo->radiation->calc = NO_RIVE;
  }
  }
  pft_st = TS_browse_ft($2,END_TS);
  //LP_printf(Simul->poutputs,"radiation time series defined length = %d\n",TS_length_ts($2));
  //printf("radiation time series defined\n");
}

| atts_radiation
;

atts_radiation : att_radiation atts_radiation
| att_radiation
;

att_radiation : LEX_MEAN LEX_EQUAL mesure
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->radiation->mean = $3;
  }
  }
  //LP_printf(Simul->poutputs,"radiation : mean = %f µE/m^2/s\n",$3);
  //printf("radiation : mean = %f µE/m^2/s\n",sect2->meteo->radiation->mean);
}
| LEX_AMPLITUDE LEX_EQUAL mesure
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->radiation->amplitude = $3;
  }
  }
  //LP_printf(Simul->poutputs,"radiation : amplitude = %f µE/m^2/s\n",$3);
  //printf("radiation : amplitude = %f µE/m^2/s\n",sect2->meteo->radiation->amplitude);
}
| LEX_ATTENUATION LEX_EQUAL flottant
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->radiation->attenuation = $3;
  }
  }
  //LP_printf(Simul->poutputs,"attenuation factor : %f\n",$3);
  //printf("attenuation factor : %f\n",sect2->meteo->radiation->attenuation);
}
;

atts_photoperiod : att_photoperiod atts_photoperiod
| att_photoperiod
;

att_photoperiod : LEX_MEAN LEX_EQUAL mesure
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->photoperiod->mean = $3;
  }
  }
  //LP_printf(Simul->poutputs,"photoperiod : mean = %f s\n",$3);
  //printf("photoperiod : mean = %f s\n",sect2->meteo->photoperiod->mean);
}
| LEX_AMPLITUDE LEX_EQUAL mesure
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->meteo->photoperiod->amplitude = $3;
  } 
  }
  //LP_printf(Simul->poutputs,"photoperiod : amplitude = %f s\n",$3);
  //printf("photoperiod : amplitude = %f s\n",sect2->meteo->photoperiod->amplitude);
}
;

exchanges : intro_exchanges atts_exchanges LEX_CLOSING_BRACE
;

intro_exchanges : LEX_EXCH LEX_EQUAL LEX_OPENING_BRACE
;

atts_exchanges : att_exchanges atts_exchanges
| att_exchanges
;

att_exchanges : LEX_CALC_SE LEX_EQUAL LEX_SETYPE
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->exchange_settings->calc_sed_eros = $3;
  }
  }
}
| LEX_PARAM_EROS LEX_EQUAL mesure 
{
  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  for(ns = 0; ns < nsections; ns++){
  Simul->psimul_bio[i][ns]->section->exchange_settings->param_ero[$1] = $3;
  }
  }
}
;

set_param_da : LEX_PARAM_RANGE LEX_EQUAL LEX_OPENING_BRACE atts_param_da LEX_CLOSING_BRACE
{
    LP_printf(Simul->poutputs, "parameters' ranges read\n");
}
;
atts_param_da : att_param_da atts_param_da
| att_param_da
;

att_param_da : LEX_PARAM_DA a_unit flottant flottant
{
    if($1 == C_CHLA_DA)
	{
            Simul->passim->units_param[$1] = $2;
            Simul->passim->param_range[$1][PARAM_DOWN] = 1 / ($2 * $4);
	    Simul->passim->param_range[$1][PARAM_UP] = 1 / ($2 * $3);
            // SW 24/01/2022
            Simul->passim->param_yesOrno[$1] = YES_TS;
            Simul->passim->num_Of_assimilated_param++;
            //chla /C = g / g 		
	}
	else
	{
	    Simul->passim->units_param[$1] = $2;
            Simul->passim->param_range[$1][PARAM_UP] = $2 * $4;
	    Simul->passim->param_range[$1][PARAM_DOWN] = $2 * $3;

            // SW 24/01/2022
            Simul->passim->param_yesOrno[$1] = YES_TS;
            Simul->passim->num_Of_assimilated_param++;
		
	}
}
| LEX_PARAM_DA flottant flottant
{
    Simul->passim->units_param[$1] = 1.0;
    Simul->passim->param_range[$1][PARAM_UP] = $3;
    Simul->passim->param_range[$1][PARAM_DOWN] = $2;
    // SW 24/01/2022
    Simul->passim->param_yesOrno[$1] = YES_TS;
    Simul->passim->num_Of_assimilated_param++;

}

/* MH 10/03/2022 : added the 3rd flottant to facilitate different random walk for each param */
| LEX_PARAM_DA a_unit flottant flottant flottant
{
    if($1 == C_CHLA_DA)
	{
            Simul->passim->units_param[$1] = $2;
            Simul->passim->param_range[$1][PARAM_DOWN] = 1 / ($2 * $4);
	    Simul->passim->param_range[$1][PARAM_UP] = 1 / ($2 * $3);
            // SW 24/01/2022
            Simul->passim->param_yesOrno[$1] = YES_TS;
            Simul->passim->num_Of_assimilated_param++;
            //chla /C = g / g

	    // MH 10/03/2022
	    Simul->passim->s_percent[$1] = $5;
	}
	else
	{
	    Simul->passim->units_param[$1] = $2;
            Simul->passim->param_range[$1][PARAM_UP] = $2 * $4;
	    Simul->passim->param_range[$1][PARAM_DOWN] = $2 * $3;

            // SW 24/01/2022
            Simul->passim->param_yesOrno[$1] = YES_TS;
            Simul->passim->num_Of_assimilated_param++;
	     // MH 10/03/2022
	    Simul->passim->s_percent[$1] = $5;
	    LP_printf(Simul->poutputs,"Target random walk of %s = %3.4f \n",PROSE_name_param($1),Simul->passim->s_percent[$1]);
	}
}
| LEX_PARAM_DA flottant flottant flottant
{
    Simul->passim->units_param[$1] = 1.0;
    Simul->passim->param_range[$1][PARAM_UP] = $3;
    Simul->passim->param_range[$1][PARAM_DOWN] = $2;
    // SW 24/01/2022
    Simul->passim->param_yesOrno[$1] = YES_TS;
    Simul->passim->num_Of_assimilated_param++;

    // MH 10/03/2022
    Simul->passim->s_percent[$1] = $4;
    LP_printf(Simul->poutputs,"Target random walk of %s = %3.4f \n",PROSE_name_param($1),Simul->passim->s_percent[$1]);

}

// MH 29/10/2021 : assigning the range of b1_river from the macrospecies properties to DA param table
/*
| LEX_PARAM_DA
{
    Simul->passim->units_param[$1] = 1.0;
    Simul->passim->param_range[$1][PARAM_UP] = Simul->psimul_bio->p_macrospecies[TOC]->degradOrgMat[MACMOD][B]->range[MAX_RIV];
    Simul->passim->param_range[$1][PARAM_DOWN] =  Simul->psimul_bio->p_macrospecies[TOC]->degradOrgMat[MACMOD][B]->range[MIN_RIV];
    LP_printf(Simul->poutputs,"%s range added to the DA param table \n ",PROSE_name_param(B1_RIVER_DA));
    }*/
;

/* Definition of the outputs */
outputs : intro_outputs def_outputs LEX_CLOSING_BRACE
{
  char *cmd;
  
  cmd = (char *) malloc(MAXCHAR_PROSE * sizeof(char));
  bzero((char *)cmd, sizeof(char));

  if (pts != NULL) {
    HYD_create_time_series(pts,Simul->outputs,Simul->pchyd->counter->nts);
    sprintf(cmd,"%s/time_series",getenv("RESULT")); 
    IO_mkdir(cmd); // SW 19/10/2022 replace system(cmd) by mkdir d function

    LP_printf(Simul->poutputs,"time_series directory %s\n",cmd);
    bzero((char *)cmd, sizeof(char));

  }
  if (plp != NULL) {
    HYD_create_long_profiles(plp,Simul->outputs,Simul->pchyd->counter->nlp);
    sprintf(cmd,"%s/longitudinal_profiles",getenv("RESULT")); 
    IO_mkdir(cmd); // SW 19/10/2022 replace system(cmd) by mkdir c function which is called in IO_mkdir function
    
    LP_printf(Simul->poutputs,"longitudinal_profiles directory %s Information: %s \n",cmd, strerror(errno));
    bzero((char *)cmd, sizeof(char));
  }
  if (pmb != NULL) {
    HYD_create_mass_balances(pmb,Simul->outputs,Simul->pchyd->counter->nmb);
    sprintf(cmd,"%s/mass_balances",getenv("RESULT")); 
    IO_mkdir(cmd); // SW 19/10/2022 replace system(cmd) by mkdir c function which is called in IO_mkdir function
    
    LP_printf(Simul->poutputs,"/mass_balances directory %s Information: %s \n",cmd, strerror(errno));
    bzero((char *)cmd, sizeof(char));
  }
  LP_printf(Simul->poutputs,"End of outputs formats reading\n"); 
}
;

intro_outputs : LEX_OUTPUTS LEX_EQUAL LEX_OPENING_BRACE
{LP_printf(Simul->poutputs,"Begin outputs formats reading\n"); }
;

def_outputs : def_output def_outputs
| def_output
;

def_output : intro_output output_options brace
{
  if (output_type == TRANSV_PROFILE) {
    //Simul->outputs[TRANSV_PROFILE] = pout;
    if (pts == NULL)
      pts = pout_hyd;
    else
      pts = HYD_chain_outputs(pout_hyd,pts);
    HYD_create_output_points(pts,pts_pktot);
    Simul->pchyd->counter->nts++;
  }

  else if (output_type == LONG_PROFILE) {
    //Simul->outputs[LONG_PROFILE] = pout;
    if (plp == NULL)
      plp = pout_hyd;
    else
      plp = HYD_chain_outputs(pout_hyd,plp);
    HYD_create_output_extents(plp,plp_pktot);
    Simul->pchyd->counter->nlp++;
  }

  else if (output_type == MASS_BALANCE) {
    if (pmb == NULL)
      pmb = pout_hyd;
    else
      pmb = HYD_chain_outputs(pout_hyd,pmb);
    HYD_create_output_extents(pmb,pmb_pktot);
    Simul->pchyd->counter->nmb++;
  }

  else Simul->outputs[output_type][0] = pout_hyd;

  if(Simul->calc_mode[DA] == YES_TS)
      Prose_create_obs_points(Simul->passim, pts_obs_pktot);
  plp_pktot = NULL;
  pmb_pktot = NULL;
  pts_pktot = NULL;
  pts_obs_pktot = NULL;
}
| def_mass_balances
| def_tube_output
| def_mass_balances_heat
;

intro_output : LEX_OUTPUT_TYPE LEX_EQUAL brace
{
  Simul->pinout->calc[$1] = YES_TS;
  pout_hyd = HYD_initialize_output(Simul->chronos);
  // SW 30/05/2018 pour soties des bios

  for(i = 0; i < NSPECIES_IO; i++)
  {
    //pout_hyd->pout->biovar[i] = 1;
    pout_hyd->pout->biovar[i] = NO_TS;//NF 11/10/2020 watch this is from c-rive where YES and NO are inversed compared to libts, meaning NO_TS = YES and NO_TS = YES
    pout_hyd->pout->biovar_unit[i] = 1.;	  
  }
  for(i = 0; i < NTEMP_IO; i++)//NF 11/10/2020
  {
    pout_hyd->pout->tempvar[i] = NO_TS;//NF 11/10/2020 watch this is from c-rive where YES and NO are inversed compared to libts, meaning NO_TS = YES and NO_TS = YES
    pout_hyd->pout->tempvar_unit[i] = 1.;	  
  }
  for(i = 0; i < NSEDVAR_IO; i++) // SW 23/01/2024
  {
    pout_hyd->pout->biosedvar[i] = NO_TS;//SW 23/01/2024 watch this is from c-rive where YES and NO are inversed compared to libts, meaning NO_TS = YES and NO_TS = YES
    pout_hyd->pout->biosedvar_unit[i] = 1.;	  
  }

  output_type = $1;
}
;

output_options : output_option output_options
| output_option
;

output_option : LEX_ANSWER
{
  Simul->pinout->calc[output_type] = $1;
}
| LEX_FILE_NAME LEX_EQUAL LEX_NAME
{
  /*Simul->outputs[output_type] = new_output();
  Simul->outputs[output_type]->fic = new_file();
  Simul->outputs[output_type]->fic->name = $3;*/
  pout_hyd->fic = new_file_io();
  pout_hyd->fic->name = $3;
}
| time_unit
| out_chronos
| variables
| points
{
  if (pts_pktot != NULL)
    pts_pktot = HYD_chain_ts_pk(HYD_browse_ts_pk(pts_pktot,END),$1);
  else 
    pts_pktot = $1; 
}
| points_obs
{
  if (pts_obs_pktot != NULL)
    pts_obs_pktot = Prose_chain_ts_obs(Prose_browse_ts_obs(pts_obs_pktot,END),$1);
  else 
    pts_obs_pktot = $1;
}
| extents
{
  //plp_pktot = new_lp_pk(); // SW 25/01/2018
  if (output_type == LONG_PROFILE) {
    if (plp_pktot != NULL)
      plp_pktot = HYD_chain_lp_pk(HYD_browse_lp_pk(plp_pktot,END),$1);
    else
      plp_pktot = $1;
  }

  else if (output_type == MASS_BALANCE) {
    if (pmb_pktot != NULL)
      pmb_pktot = HYD_chain_lp_pk(HYD_browse_lp_pk(pmb_pktot,END),$1);
    else
      pmb_pktot = $1;
  }
}
| graphics
| configure_weight

;

time_unit : LEX_TUNIT LEX_EQUAL a_unit
{
  pout_hyd->pout->time_unit = 1 / $3;
}
;

out_chronos : LEX_CHRONOS LEX_EQUAL brace out_times brace
{
  //pout_hyd->pout->t_out[INI_IO] = Simul->chronos->t[BEGINNING];
  pout_hyd->pout->t_out[CUR_IO] = pout_hyd->pout->t_out[INI_IO]; // SW 20/05/2019 output extent
}
;

out_times : out_time out_times
| out_time
;

out_time : LEX_TIME LEX_EQUAL mesure
{
  pout_hyd->pout->t_out[$1] = $3;
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_DATE_HH_MM_SS
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();
  
  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    sscanf($4,"%d:%d:%d", &pd_tmp->hh, &pd_tmp->min,&pd_tmp->ss);
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  pout_hyd->pout->t_out[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),pout_hyd->pout->t_out[$1]); 
  
  free(pd_tmp);
  pd_tmp = NULL;
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_DATE_HH_MM
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();

  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    sscanf($4,"%d:%d", &pd_tmp->hh, &pd_tmp->min);
    pd_tmp->ss = 0.;
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  pout_hyd->pout->t_out[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  //LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]); 

  free(pd_tmp);
  pd_tmp = NULL;
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_INT
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();

  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    pd_tmp->hh = $4;
    pd_tmp->min = 0.;
    pd_tmp->ss = 0.;
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  pout_hyd->pout->t_out[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  //LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]);
  free(pd_tmp);
  pd_tmp = NULL;

}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();

  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    pd_tmp->hh = 0.;
    pd_tmp->min = 0.;
    pd_tmp->ss = 0.;
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  pout_hyd->pout->t_out[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  //LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]);
  free(pd_tmp);
  pd_tmp = NULL;
 
}
| LEX_TIMESTEP LEX_EQUAL mesure
{
  pout_hyd->pout->deltat = $3 < Simul->chronos->dt ? Simul->chronos->dt : $3;
}
;

variables : LEX_VAR LEX_EQUAL brace var_list brace
| LEX_ANNEXVAR LEX_EQUAL LEX_OPENING_BRACE list_annex_var LEX_CLOSING_BRACE
;

list_annex_var : one_annex_var_bio list_annex_var
| one_annex_var_bio
;

one_annex_var_bio : LEX_ONEANNEXVAR a_unit
{
    Simul->calc_bio_annexvar[$1] = YES_TS;
	Simul->unit_bio_annexvar[$1] = $2;
}

var_list : one_var var_list
| one_var
| one_var_bio var_list
| one_var_bio
| one_var_temp var_list //NF 12/10/2020
| one_var_temp //NF 12/10/2020
| one_var_sediment var_list // SW 23/01/2024
| one_var_sediment // SW 23/01/2024
;

one_var : LEX_ONE_VAR a_unit
{
  pout_hyd->pout->hydvar[$1] = YES_TS;
  pout_hyd->pout->hydvar_unit[$1] = $2;
}
;

one_var_bio : LEX_ONESPECIES a_unit
{
  pout_hyd->pout->biovar[$1] = YES_TS;
  pout_hyd->pout->biovar_unit[$1] = $2;
}
;

//SW 23/01/2024
one_var_sediment : LEX_SED_VAR a_unit
{
  pout_hyd->pout->biosedvar[$1] = YES_TS;
  pout_hyd->pout->biosedvar_unit[$1] = $2;
}
;

//NF 11/10/2020
one_var_temp : LEX_TEMPSPECIES a_unit
{
  pout_hyd->pout->tempvar[$1] = YES_TS;
  pout_hyd->pout->tempvar_unit[$1] = $2;
  if (Simul->calc_mode[SEB] != YES_TS && $1 == TA_IO)
    {
      LP_warning(Simul->poutputs,"\tInconsistancy in output Tair request: Tthere is no T_air variable accounted into the calculations\nPlease check that you didn't forget the atmospheric forcings.\n So far, only SAFRAN forcings are implemented.\nRequest ignored NO TAIR OUTPUTS\n");
      pout_hyd->pout->tempvar[$1] = NO_TS;
    }
}
;

points : LEX_POINTS LEX_EQUAL brace river pk_list brace
{
  $$ = $5;
}
;

points_obs : LEX_POINTS_OBS LEX_EQUAL brace river pk_list_obs brace
{
  $$ = $5;
}
;

river : LEX_NAME
{
  river = $1;
  LP_lowercase(river);
}
;

pk_list : one_pk pk_list
{
  $$ = HYD_chain_ts_pk($1,$2);
}
| one_pk
{
  $$ = $1;
}
;

pk_list_obs : one_pk_obs pk_list_obs
{
  $$ = Prose_chain_ts_obs($1,$2);
}
| one_pk_obs
{
  $$ = $1;
}
;

one_pk : mesure LEX_INT 
{
  pts_pk = new_ts_pk();
  bzero((s_ts_pk_hyd *)pts_pk,sizeof(s_ts_pk_hyd)); // SW 26/01/2018 il faut initialiser a NULL
  pts_pk->river = river;
  pts_pk->pk = $1;
  pts_pk->pk_type = PK_HYD; // SW 26/01/2018 
  pts_pk->branch_nb = $2;
  //find_element_reach(ts_pk);
  $$ = pts_pk;
  pout_hyd->npk++;
}
;

one_pk_obs : mesure LEX_INT LEX_ONESPECIES LEX_INT read_units f_ts
{
  pts_pk_obs = new_obs_assimilation();
  bzero((s_carac_obs_assim *)pts_pk_obs,sizeof(s_carac_obs_assim)); // SW 26/01/2018 il faut initialiser a NULL
  pts_pk_obs->river = river;
  pts_pk_obs->pk = $1;
  pts_pk_obs->answer_obs = NO_TS;
  pts_pk_obs->pk_type = PK_HYD; // SW 26/01/2018 
  pts_pk_obs->branch_nb = $2;
  pts_pk_obs->obs = $6;

  /* MH 20/12/2022 moving average implemented on the obs data */
  LP_printf(Simul->poutputs, "t = %f  o2 = %f , length = %d \n ", pts_pk_obs->obs->t/NSEC_DAY_TS, pts_pk_obs->obs->ft, TS_length_ts(pts_pk_obs->obs));
  
  if (Simul->passim->da_mov_ave == YES_TS) {

    //getting the input arguments of the TS_moving_average function
    lneigh = Simul->passim->lneigh;
    theta_mave =  Simul->passim->mov_ave_theta;
    nval_min =  Simul->passim->nval_min;
    //nval_min = ceil(( (lneigh/Simul->passim->do_time_step) + 1)/2);
    pft_obs = pts_pk_obs->obs;
    pft_obs =TS_browse_ft(pft_obs,END_TS);
    tf=pft_obs->t;
    pft_obs=TS_browse_ft(pft_obs,BEGINNING_TS);
    t0=pft_obs->t;
        
    pft_mov_ave = TS_moving_average(pft_obs,(double) lneigh, (double) theta_mave,t0,tf, nval_min,Simul->poutputs);

    pft_mov_ave =TS_browse_ft(pft_mov_ave,BEGINNING_TS);

    /*
    
    LP_printf(Simul->poutputs, " t_obs = %f  o2_obs = %f , length_obs = %d \n ", pft_obs->t/NSEC_DAY_TS, pft_obs->ft/31.25, TS_length_ts(pft_obs));

    LP_printf(Simul->poutputs, " t_mva = %f  o2_mva = %f , length_mva = %d \n ", pft_mov_ave->t/NSEC_DAY_TS, pft_mov_ave->ft/31.25, TS_length_ts(pft_mov_ave));

    
    //printing t and o2 at every time step
    pft_mov_ave = TS_browse_ft(pft_mov_ave,BEGINNING_TS);
    pft_obs = TS_browse_ft(pft_obs,BEGINNING_TS);
    
    LP_printf(Simul->poutputs, "obs_sdv = %f \n  mov_ave_sdv = %f \n" , TS_stdev(pft_obs,(TS_average(pft_obs))),  TS_stdev(pft_mov_ave,(TS_average(pft_mov_ave))));
    
    double t = pft_obs->t;
      while( (t <=tf) && (pft_obs != NULL) )
	{
	  LP_printf(Simul->poutputs, "t_obs = %f, t_mov = %f, o2_obs = %f, o2_mov = %f \n", t/NSEC_DAY_TS, pft_mov_ave->t/NSEC_DAY_TS, pft_obs->ft/31.25, pft_mov_ave->ft/31.25);
	 
	  if (pft_obs->next != NULL)
	    {
	      pft_mov_ave = pft_mov_ave->next;
	      pft_obs = pft_obs->next;
	      t = pft_obs->t;
	    }
	  else {

	    LP_printf(Simul->poutputs, " End of printing \n ");

	    break;
	  }
	} 
    */
    TS_free_ts(pft_obs,Simul->poutputs);
    pts_pk_obs->obs = pft_mov_ave;
    //LP_printf(Simul->poutputs, "After attaching pft_mov_ave to obs: t = %f  o2_new = %f , length = %d \n ", pts_pk_obs->obs->t/NSEC_DAY_TS, pts_pk_obs->obs->ft/31.25, TS_length_ts(pts_pk_obs->obs));

    // now to print out the moving average data in the output folder

    char cmd[MAXCHAR_PROSE];
    char filename[MAXCHAR_PROSE];
    
    sprintf(cmd,"mkdir %s/moving_average",getenv("RESULT"));
    system(cmd);
    
    sprintf(filename,"%s/moving_average/mov_ave_pk%4.2f",getenv("RESULT"), pts_pk_obs->pk/1000);

    FILE *filename_2 = fopen(filename,"w");
    TS_printf_ts(pts_pk_obs->obs, filename_2);
    fclose(filename_2);

    }
 
  pts_pk_obs->var = $3;
  pts_pk_obs->num = $4 - 1;
  $$ = pts_pk_obs;
  Simul->passim->N_obs++;
}
| mesure LEX_INT LEX_ONESPECIES LEX_INT read_units date_f_ts
{
  pts_pk_obs = new_obs_assimilation();
  bzero((s_carac_obs_assim *)pts_pk_obs,sizeof(s_carac_obs_assim)); // SW 26/01/2018 il faut initialiser a NULL
  pts_pk_obs->river = river;
  pts_pk_obs->pk = $1;
  pts_pk_obs->answer_obs = NO_TS;
  pts_pk_obs->pk_type = PK_HYD; // SW 26/01/2018 
  pts_pk_obs->branch_nb = $2;
  pts_pk_obs->obs = $6;


  /* MH 20/12/2022 moving average implemented on the obs data */
  LP_printf(Simul->poutputs, "t = %f  o2 = %f , length = %d \n ", pts_pk_obs->obs->t/NSEC_DAY_TS, pts_pk_obs->obs->ft, TS_length_ts(pts_pk_obs->obs));
  
  if (Simul->passim->da_mov_ave == YES_TS) {

    //getting the input arguments of the TS_moving_average function
    lneigh = Simul->passim->lneigh;
    theta_mave =  Simul->passim->mov_ave_theta;
    nval_min =  Simul->passim->nval_min;
    //nval_min = ceil(( (lneigh/Simul->passim->do_time_step) + 1)/2);
    pft_obs = pts_pk_obs->obs;
    pft_obs =TS_browse_ft(pft_obs,END_TS);
    tf=pft_obs->t;
    pft_obs=TS_browse_ft(pft_obs,BEGINNING_TS);
    t0=pft_obs->t;
     
    pft_mov_ave = TS_moving_average(pft_obs,(double) lneigh, (double) theta_mave,t0,tf, nval_min,Simul->poutputs);

    pft_mov_ave =TS_browse_ft(pft_mov_ave,BEGINNING_TS);

    TS_free_ts(pft_obs,Simul->poutputs);
    pts_pk_obs->obs = pft_mov_ave;
    //LP_printf(Simul->poutputs, "After attaching pft_mov_ave to obs: t = %f  o2_new = %f , length = %d \n ", pts_pk_obs->obs->t/NSEC_DAY_TS, pts_pk_obs->obs->ft/31.25, TS_length_ts(pts_pk_obs->obs));

    // now to print out the moving average data in the output folder

    char cmd[MAXCHAR_PROSE];
    char filename[MAXCHAR_PROSE];
    
    sprintf(cmd,"mkdir %s/moving_average",getenv("RESULT"));
    system(cmd);
    
    sprintf(filename,"%s/moving_average/mov_ave_pk%4.2f",getenv("RESULT"), pts_pk_obs->pk/1000);

    FILE *filename_2 = fopen(filename,"w");
    TS_printf_ts(pts_pk_obs->obs, filename_2);
    fclose(filename_2);
    
    }
  
  pts_pk_obs->var = $3;
  pts_pk_obs->num = $4 - 1;
  $$ = pts_pk_obs;
  Simul->passim->N_obs++;
}
/* MH 16/05/2022 implementation of user defined parameters for calculation of station weights */ 
// configure_weight as an output option
;
configure_weight : intro_configure_weight atts_configure_weight LEX_CLOSING_BRACE
{
  LP_printf(Simul->poutputs,"All station weights parameters read \n");
  Simul->passim->weight_calc_assim = pstat_weight;

}
;
intro_configure_weight : LEX_CONFIGURE_WEIGHTS LEX_EQUAL LEX_OPENING_BRACE
{
  
  LP_printf(Simul->poutputs,"reading of station weight configuration params started \n");
  pstat_weight = PROSE_init_calc_stat_weight() ;

}
;
atts_configure_weight : att_configure_weight atts_configure_weight
| att_configure_weight
;
att_configure_weight : LEX_WEIGHTED_STATIONS LEX_EQUAL LEX_ANSWER
{
  if ($3 == YES_TS)
    {
      LP_printf(Simul->poutputs, "Let's weight the assimilation stations \n");
      //pstat_weight = PROSE_init_calc_stat_weight() ;
    }
  pstat_weight->weighted_stations = $3;

}
| LEX_VEL_LOW_FLOW LEX_EQUAL flottant
{
  //psimul->passim->weight_calc_assim->vel_low_flow
  LP_printf(Simul->poutputs,"reading of velocity at low flow \n");
  pstat_weight->vel_low_flow = $3;

}
| LEX_BACT_YIELD LEX_EQUAL flottant
{
  LP_printf(Simul->poutputs,"reading of yield \n");
  pstat_weight->bact_yield = $3;

}
|  LEX_BACT_GROWTH_RATE LEX_EQUAL flottant
{
  LP_printf(Simul->poutputs,"reading of mu \n");
  pstat_weight->bact_growth_rate = $3;

}
|  LEX_FAST_BDOM LEX_EQUAL flottant
{
  LP_printf(Simul->poutputs,"reading of mod1 \n");
  pstat_weight->fast_bdom_low_flow = $3 ;

}
| LEX_BACT_BIOMASS LEX_EQUAL flottant
{
  LP_printf(Simul->poutputs,"reading of bact \n");
  pstat_weight->bact_biomass = $3 ;

}
| LEX_OM_FLUX_THRESHOLD LEX_EQUAL flottant
{
  LP_printf(Simul->poutputs,"reading of mod1_flux_threshold \n");
  pstat_weight->OM_flux_threshold = $3;

}
| LEX_OM_DECAY_RATE LEX_EQUAL flottant
{
  LP_printf(Simul->poutputs,"reading of decay rate \n");
  pstat_weight->OM_decay_rate = $3;
  
}
;



extents : LEX_EXTENT LEX_EQUAL brace river two_pk brace
{
  if (output_type == LONG_PROFILE) {
    $$ = plp_pk;
  }

  else if ((output_type == MASS_BALANCE) || (output_type == MASS_BALANCE_BIO_OUT)) {
    $$ = pmb_pk;
  }
  else if ((output_type == TUBE_MESH) || (output_type == TUBE_HYD)) { // SW 30/03/2021
    $$ = plp_pk_tube;
  }
    else if(output_type == ENERGY_BALANCE_HEAT_OUT)
    $$ = pmb_pk_heat; 
}
;

two_pk : mesure LEX_INT mesure LEX_INT
{
  if (output_type == LONG_PROFILE) {
    plp_pk = HYD_init_lp_pk(river,$1,$2,$3,$4,Simul->pchyd);
    pout_hyd->npk++;
	//LP_printf(Simul->poutputs,"long file lu \n");
  }

  else if (output_type == MASS_BALANCE) {
    pmb_pk = HYD_init_lp_pk(river,$1,$2,$3,$4,Simul->pchyd);
    pout_hyd->npk++;
  }
  else if (output_type == MASS_BALANCE_BIO_OUT) {
    pmb_pk = HYD_init_lp_pk(river,$1,$2,$3,$4,Simul->pchyd);
    Simul->npk_mb_bio++;
  }
  else if ((output_type == TUBE_MESH) || (output_type == TUBE_HYD)) { // SW 30/03/2021
    plp_pk_tube = HYD_init_lp_pk(river,$1,$2,$3,$4,Simul->pchyd);
  }
  else if(output_type = ENERGY_BALANCE_HEAT_OUT) {
   pmb_pk_heat = HYD_init_lp_pk(river,$1,$2,$3,$4,Simul->pchyd); // SW 03/05/2021
  }

  //if(Simul->calc_mode[MB_BIO] == YES_TS) // SW 30/03/2021 maybe a bug when bio mass balance defined before LONG_PROFILE or MASS_BALANCE
    // Simul->npk_mb_bio++;
  //else
     //pout_hyd->npk++;
}
;

graphics : LEX_GRAPHICS LEX_EQUAL LEX_ANSWER
{
  pout_hyd->graphics = NO_TS;
}
| LEX_GRAPHICS LEX_EQUAL LEX_GRAPH_TYPE
{
  pout_hyd->graphics = $3;

  if ($3 == GNUPLOT) {
    char cmd[MAXCHAR_PROSE];
    sprintf(cmd,"%s/gnuplot",getenv("RESULT")); 
    //system(cmd);
    IO_mkdir(cmd); // SW 19/10/2022 mkdir recurssive to avoid issues for dummies, system may has problems sometimes
  }
}
;

def_mass_balances : intro_mb_bio def_mbs LEX_CLOSING_BRACE
{
  s_total_mb *mb_temp;
  mb_temp = $2;
  Simul->total_mb = (s_total_mb ***) malloc(Simul->passim->N_particules * sizeof(s_total_mb **));

  for(i = 0; i < Simul->passim->N_particules; i++)
  {
  Simul->total_mb[i] = (s_total_mb **) malloc(num_mb * sizeof(s_total_mb *)); 
  //Simul->total_mb[i][0] = Prose_copy_total_mb(mb_temp);
  }
  //for(j = 1; j < num_mb; j++)
     //Simul->total_mb[i][j] = Prose_copy_total_mb(mb_temp->next);
  j = 0;
  while(mb_temp != NULL)
  {
     for(i = 0; i < Simul->passim->N_particules; i++)
     {
         Simul->total_mb[i][j] = Prose_copy_total_mb(mb_temp);
     }
	 mb_temp = mb_temp->next;
	 j++;
  }
  char cmd[MAXCHAR_PROSE];
  sprintf(cmd,"%s/mass_balance_bio",getenv("RESULT")); 
  //system(cmd);

  IO_mkdir(cmd); // SW 19/10/2022 mkdir recurssive to avoid issues for dummies, system may has problems sometimes
  //SW 06/12/2019 ATTENTION : coherence entre mb_temp et pmb_pktot, à modifier un jour plus proprement
  if(output_type == MASS_BALANCE_BIO_OUT) {
	PROSE_create_output_extents_bio(Simul->npk_mb_bio,pmb_pktot);
  }

}
;

intro_mb_bio : LEX_MB LEX_EQUAL LEX_OPENING_BRACE
{
   //output_type = MASS_BALANCE;
   output_type = MASS_BALANCE_BIO_OUT;
   Simul->calc_mode[MB_BIO] = YES_TS;
}
;

def_mbs : def_one_mb def_mbs
{
  $$ = chain_mbs($1,$2);
}
| def_one_mb
{
  $$ = $1;
}
;

def_one_mb : intro_mb atts_mb LEX_CLOSING_BRACE
{
  mb2->t0 = mb2->t[BEGINNING] + mb2->ndt;
  $$ = mb2;
  mb2 = NULL;
}
;

intro_mb : LEX_NAME LEX_EQUAL LEX_OPENING_BRACE
{
  num_mb++;
  //for(ns = 0; ns < nsections; ns++) 
     Simul->counter_bio->nmb++;
  mb2 = init_mass_balance(Simul->psimul_bio[0][0],Simul->chronos->t[BEGINNING], Simul->chronos->t[END], Simul->chronos->dt);
  //mb2->name = $1;
  sprintf(mb2->name,"%s",$1); //SW 21/02/2017
  LP_printf(Simul->poutputs,"\ncharacteristics of %s\n",mb2->name);
  //printf("\ncharacteristics of %s\n",mb2->name);
}
;

atts_mb : att_mb atts_mb
| att_mb
;

att_mb : LEX_TIME LEX_EQUAL mesure
{
  mb2->t[$1] = $3;
  LP_printf(Simul->poutputs,"t %s = %f s\n",name_extremum($1),mb2->t[$1]);
  //printf("t %s = %f s\n",name_extremum($1),mb2->t[$1]);
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_DATE_HH_MM_SS
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();
  
  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    sscanf($4,"%d:%d:%d", &pd_tmp->hh, &pd_tmp->min,&pd_tmp->ss);
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  mb2->t[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),pout_hyd->pout->t_out[$1]); 
  
  free(pd_tmp);
  pd_tmp = NULL;
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_DATE_HH_MM
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();

  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    sscanf($4,"%d:%d", &pd_tmp->hh, &pd_tmp->min);
    pd_tmp->ss = 0.;
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  mb2->t[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  //LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]); 

  free(pd_tmp);
  pd_tmp = NULL;
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_INT
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();

  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    pd_tmp->hh = $4;
    pd_tmp->min = 0.;
    pd_tmp->ss = 0.;
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  mb2->t[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  //LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]);
  free(pd_tmp);
  pd_tmp = NULL;

}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();

  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    pd_tmp->hh = 0.;
    pd_tmp->min = 0.;
    pd_tmp->ss = 0.;
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  mb2->t[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  //LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]);
  free(pd_tmp);
  pd_tmp = NULL;
 
}
| LEX_TUNIT LEX_EQUAL a_unit
{
  mb2->time_unit = $3;
}
| LEX_NSTEPS LEX_EQUAL flottant
{
  mb2->ndt = $3 * Simul->chronos->dt;
  LP_printf(Simul->poutputs,"ndt = %f s\n",mb2->ndt);
  //printf("ndt = %f s\n",mb2->ndt);
}
| LEX_SPECIES LEX_EQUAL LEX_OPENING_BRACE list_species_mb LEX_CLOSING_BRACE
| LEX_ANNEXVAR LEX_EQUAL LEX_OPENING_BRACE list_annex_var_mb LEX_CLOSING_BRACE
| LEX_ADS_SPECIES LEX_EQUAL LEX_OPENING_BRACE list_adsorbed_var_mb LEX_CLOSING_BRACE
| extents
{
   if(output_type == MASS_BALANCE_BIO_OUT) {
     if(pmb_pktot != NULL)
       pmb_pktot = HYD_chain_lp_pk(HYD_browse_lp_pk(pmb_pktot,END),$1);
     else
      pmb_pktot = $1;
  }
}
;

list_species_mb : LEX_ONESPECIES a_unit list_species_mb
{
  if (Simul->counter_bio->nsubspecies[$1] == 0) {
    LP_warning(Simul->poutputs,"%s mass balance asked for but the species wasn't defined\n",name_species($1));
  }
  else {
    for (j = 0; j < Simul->counter_bio->nsubspecies[$1]; j++) {
      mb2->calc_mb_species[$1][j] = YES_RIVE; // SW 31/05/2018 ici No == 1
      mb2->unit_mb_species[$1][j] = $2;
    }
    LP_printf(Simul->poutputs,"%s mass balance asked for\n",name_species($1));
    //printf("%s mass balance asked for\n",name_species($1));
  }
}
| LEX_ONESPECIES a_unit
{
  for (j = 0; j < Simul->counter_bio->nsubspecies[$1]; j++) {
    mb2->calc_mb_species[$1][j] = YES_RIVE; // SW 31/05/2018 ici No == 1
    mb2->unit_mb_species[$1][j] = $2;
  }
  LP_printf(Simul->poutputs,"%s mass balance asked for\n",name_species($1));
  //printf("%s mass balance asked for\n",name_species($1));
}
;

list_annex_var_mb : LEX_ONEANNEXVAR a_unit list_annex_var_mb
{
  for (j = 0; j < Simul->counter_bio->nsubannex_var[$1]; j++) {
    mb2->calc_mb_annex_var[$1][j] = YES_RIVE; // SW 31/05/2018 ici No == 1
    mb2->unit_mb_annex_var[$1][j] = $2;
  }
  LP_printf(Simul->poutputs,"%s mass balance asked for\n",name_annex_var($1));
  //printf("%s mass balance asked for\n",name_annex_var($1));
}
| LEX_ONEANNEXVAR a_unit
{
  for (j = 0; j < Simul->counter_bio->nsubannex_var[$1]; j++) {
    mb2->calc_mb_annex_var[$1][j] = YES_RIVE; // SW 31/05/2018 ici No == 1
    mb2->unit_mb_annex_var[$1][j] = $2;
  }
  LP_printf(Simul->poutputs,"%s mass balance asked for\n",name_annex_var($1));
  //printf("%s mass balance asked for\n",name_annex_var($1));
}
;

list_adsorbed_var_mb : LEX_ONESPECIES a_unit list_adsorbed_var_mb
{
  for (j = 0; j < Simul->counter_bio->nsubspecies[$1]; j++) {
    mb2->calc_mb_adsorbed_species[$1-NPART][j] = YES_RIVE; // SW 31/05/2018 ici No == 1
    mb2->unit_mb_adsorbed_species[$1-NPART][j] = $2;
  }
  LP_printf(Simul->poutputs,"adsorbed %s mass balance asked for\n",name_species($1));
  //printf("adsorbed %s mass balance asked for\n",name_species($1));
}
| LEX_ONESPECIES a_unit
{
  for (j = 0; j < Simul->counter_bio->nsubspecies[$1]; j++) {
    mb2->calc_mb_adsorbed_species[$1-NPART][j] = YES_RIVE; // SW 31/05/2018 ici No == 1
    mb2->unit_mb_adsorbed_species[$1-NPART][j] = $2;
  }
  LP_printf(Simul->poutputs,"adsorbed %s mass balance asked for\n",name_species($1));
  //printf("adsorbed %s mass balance asked for\n",name_species($1));
}
;

def_mass_balances_heat : intro_mb_heat def_outs_mb_heat LEX_CLOSING_BRACE
{
  s_mbheat_mb *mb_heat;
  mb_heat = $2;
  
  Simul->mb_heat = (s_mbheat_mb **) malloc(nout_mb_heat * sizeof(s_mbheat_mb *)); 

  j = 0;
  while(mb_heat != NULL)
  {
         Simul->mb_heat[j] = mb_heat;
         Simul->mb_heat[j]->nmb = nout_mb_heat;
	 mb_heat = mb_heat->next;
	 j++;
  }
  
}
;

intro_mb_heat : LEX_EB_HEAT LEX_EQUAL LEX_OPENING_BRACE
{
   char cmd[MAXCHAR_PROSE];
   Simul->pout_tube->calc_output = YES_TS;
   sprintf(cmd,"%s/energy_balance_heat", getenv("RESULT"));
   IO_mkdir(cmd); // SW 19/10/2022 mkdir recurssive to avoid issues for dummies, system may has problems sometimes
   
   Simul->calc_mode[EB_HEAT] = YES_TS;
   output_type = ENERGY_BALANCE_HEAT_OUT;
}
;

def_outs_mb_heat : def_out_mb_heat def_outs_mb_heat
{
   $$ = MB_chain_heat_mbs($1, $2);
}
| def_out_mb_heat
{
  $$ = $1;
}
;

def_out_mb_heat : intro_mb_heat_out atts_heat_out LEX_CLOSING_BRACE
{
   
   mb_heat_temp->t0 = mb_heat_temp->t[BEGINNING] + mb_heat_temp->ndt;

   /*** SW 16/06/2023 check if we define t_ini and t_end. if not, set as the ones of simulation. ***/
   if((mb_heat_temp->t[BEGINNING] < EPS_TS) && fabs(mb_heat_temp->t[END] < EPS_TS)) // check if they are both 0
   {
       mb_heat_temp->t[BEGINNING] = Simul->chronos->t[BEGINNING];
       mb_heat_temp->t[END] = Simul->chronos->t[END];
   }

   $$ = mb_heat_temp;
   mb_heat_temp = NULL;   
}
;

intro_mb_heat_out : LEX_NAME LEX_EQUAL LEX_OPENING_BRACE
{
   nout_mb_heat++;
   mb_heat_temp = MB_init_mb_heat();
   sprintf(mb_heat_temp->name,"%s",$1);
   
}
;

atts_heat_out : att_heat_out atts_heat_out
| att_heat_out
;

att_heat_out : LEX_TIME LEX_EQUAL mesure
{
  mb_heat_temp->t[$1] = $3;
  LP_printf(Simul->poutputs, "t mab_heat = %f\n", $3); 
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_DATE_HH_MM_SS
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();
  
  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    sscanf($4,"%d:%d:%d", &pd_tmp->hh, &pd_tmp->min,&pd_tmp->ss);
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  mb_heat_temp->t[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),pout_hyd->pout->t_out[$1]); 
  
  free(pd_tmp);
  pd_tmp = NULL;
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_DATE_HH_MM
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();

  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    sscanf($4,"%d:%d", &pd_tmp->hh, &pd_tmp->min);
    pd_tmp->ss = 0.;
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  mb_heat_temp->t[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  //LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]); 

  free(pd_tmp);
  pd_tmp = NULL;
}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY LEX_INT
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();

  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    pd_tmp->hh = $4;
    pd_tmp->min = 0.;
    pd_tmp->ss = 0.;
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  mb_heat_temp->t[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  //LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]);
  free(pd_tmp);
  pd_tmp = NULL;

}
| LEX_TIME LEX_EQUAL LEX_DATE_DAY
{
  s_date_ts *pd_tmp;
  pd_tmp = new_date_ts();

  switch(Simul->date_format){
  case FR_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->dd, &pd_tmp->mm,&pd_tmp->yyyy);
    break;
  case US_TS:
    sscanf($3,"%d/%d/%d", &pd_tmp->mm, &pd_tmp->dd,&pd_tmp->yyyy);
    
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    pd_tmp->hh = 0.;
    pd_tmp->min = 0.;
    pd_tmp->ss = 0.;
    
    //pd_tmp->mm -= 1; // the list of months in pd (s_date_ts) startes at 0 and not 1, it is done in TS_date2julian_dd_hm

  mb_heat_temp->t[$1] = TS_date2julian_dd_hm(pd_tmp,Simul->chronos->yr0,Simul->poutputs) * NSEC_DAY_TS;
  //LP_printf(Simul->poutputs,"t %s = %f s\n",HYD_name_extremum($1),Simul->chronos->t[$1]);
  free(pd_tmp);
  pd_tmp = NULL;
 
}
| LEX_TUNIT LEX_EQUAL a_unit
{
  mb_heat_temp->time_unit = 1 / $3;
}
| LEX_NSTEPS LEX_EQUAL flottant
{
  mb_heat_temp->ndt = $3 * Simul->chronos->dt;
}
| LEX_HEAT_UNIT LEX_EQUAL a_unit
{
  mb_heat_temp->unit_mb_heat = 1 / $3;
}
| extents
{
      if(output_type == ENERGY_BALANCE_HEAT_OUT)
          mb_heat_temp->pk_mb_heat = $1;
      LP_printf(Simul->poutputs, "mab_heat %s pks defined\n", mb_heat_temp->name); 

}
;


/* SW 30/03/2021 define tube output MESH or HYD*/
def_tube_output : intro_tube_output def_outs_tube LEX_CLOSING_BRACE
{
   
   int ntb;
   s_output_tube_type *pout_tube_temp;
   
   pout_tube_temp = $2;
   Simul->pout_tube->nout = nout_tube;
   LP_printf(Simul->poutputs, "number of outputs for tube = %d\n", Simul->pout_tube->nout);
   Simul->pout_tube->poutput_tube_type = (s_output_tube_type **) malloc(nout_tube * sizeof(s_output_tube_type *));
  
   for(ntb = 0; ntb < nout_tube; ntb++)
   {
      Simul->pout_tube->poutput_tube_type[ntb] = pout_tube_temp;
      pout_tube_temp = pout_tube_temp->next;
   }
}
;

intro_tube_output : LEX_TUBE_OUT LEX_EQUAL LEX_OPENING_BRACE
{
   char cmd[MAXCHAR_PROSE];
   Simul->pout_tube->calc_output = YES_TS;
   sprintf(cmd,"%s/tube", getenv("RESULT"));
   IO_mkdir(cmd); // SW 19/10/2022 mkdir recurssive to avoid issues for dummies, system may has problems sometimes
}
;

def_outs_tube : def_out_tube def_outs_tube
{
  // to chain 
  $$ = TUB_chain_output_bck_tube($1, $2);
}
| def_out_tube
{
  $$ = $1;
}
;

def_out_tube : intro_tube_out atts_tube LEX_CLOSING_BRACE
{
   $$ = tb_out;
   nout_tube++;
   tb_out->t_out[CUR_IO] = tb_out->t_out[INI_IO];
   tb_out = NULL;
}
;

intro_tube_out : LEX_TUBE_OUT_TYPE LEX_EQUAL LEX_OPENING_BRACE
{
   output_type = $1;
   tb_out = new_tube_output_type();
   bzero((s_output_tube_type *)tb_out, sizeof(s_output_tube_type));
   tb_out->output_type = $1;
   

}
;
atts_tube : att_tube atts_tube
| att_tube
;

att_tube : LEX_TIME LEX_EQUAL mesure
{
  tb_out->t_out[$1] = $3;
}
| LEX_TIMESTEP LEX_EQUAL mesure
{
  tb_out->dt = $3 < Simul->chronos->dt ? Simul->chronos->dt : $3;
}
| LEX_TUNIT LEX_EQUAL a_unit
{
  tb_out->time_unit = 1. / $3;
}
| extents
{
      // SW 30/03/2021 if we want to define several extents, we need to chain with HYD_chain_lp_pk
      tb_out->pk_tube = $1;
}
;


/* This part describes how to deal with numbers and units*/
read_units : a_unit a_unit 
{
  unit_t = $1;
  unit_f = $2;
  nbsteps = 0;
}
;

a_unit : LEX_OPENING_BRACKET units LEX_CLOSING_BRACKET {$$ = $2;}
;


/* This part allows a_unit conversion in the programm units (SI), understanding power, division... */
units : one_unit units 
{
  $$ = $1 * $2;
}
| one_unit 
{
  $$ = $1;
}
;

one_unit : a_unit_value 
{
  $$ = $1;
}
| LEX_INV a_unit_value 
{
  $$ = 1.0/$2;
}
; /* inversion */


a_unit_value : all_units 
{
  $$ = $1;
} /* simple a_unit */
| all_units LEX_POW flottant 
{
  $$ = pow($1,$3);
}
; /* power a_unit */

all_units : LEX_A_UNIT  {$$ = $1;}
| LEX_A_UNIT_MOL {$$ = $1;}
| LEX_A_UNIT LEX_VAR_UNIT {$$ = $1*$2;}
| LEX_A_UNIT_MOL LEX_VAR_UNIT {$$ = $1;}
;

mesure : a_unit flottant 
{
  $$ = $1 * $2;
} 
| flottant a_unit//LV 06/06/2012 pour ne pas a avoir à modifier tous les fichiers de profils ProSe !
{
  $$ = $1 * $2;
}
|flottant 
{
  $$ = $1;
}
;

flottant : LEX_DOUBLE
{
  $$ = $1;
}
| LEX_INT
{
  $$ = (double)$1;
}
;


/* Examples of chaining function of t s_ft */
f_ts : f_t f_ts 
{
  $$ = TS_chain_bck_ts($1,$2);
}
|  f_t 
{ 
  $$ = $1;
};

f_t : flottant flottant 
{
  nbsteps++;
  told = tnew;
  tnew = $1 * unit_t;
  if ((line_nb >= 1) && (tnew < told) && (nbsteps >=2)) 
    LP_warning(Simul->poutputs,"%s%4.2f -> Inversion of the time steps in %s line %d, tnew = %f told = %f\n",CODE_NAME,NVERSION_PROSE,current_read_files[pile].name,line_nb+1,tnew,told);
  
  if((int)$2 == CODE_TS) // SW 15/03/2021 no unit for CODE_TS (-9999)
    $$ = TS_create_function($1*unit_t,$2);
  else
    $$ = TS_create_function($1*unit_t,$2*unit_f);
}
| flottant LEX_CODE_TS // SW 11/03/2021 define 0 CODE_TS in inflows when data are not available. CODE_TS will be replaced by default_T_value or default_C_value defined by user
{
  nbsteps++;
  told = tnew;
  tnew = $1 * unit_t;
  if ((line_nb >= 1) && (tnew < told) && (nbsteps >=2)) 
    LP_warning(Simul->poutputs,"%s%4.2f -> Inversion of the time steps in %s line %d, tnew = %f told = %f\n",CODE_NAME,NVERSION_PROSE,current_read_files[pile].name,line_nb+1,tnew,told);
	$$ = TS_create_function($1*unit_t,$2); // pas d'unit� pour CODE_TS = -9999
}
;

/* Examples of chaining function of t s_ft from date */
date_f_ts : date_f_t date_f_ts
{
  $$ = TS_chain_bck_ts($1,$2);
}
| date_f_t
{
  $$ = $1;
};

date_f_t : LEX_DATE_DAY LEX_DATE_HH_MM_SS flottant
{
  double julian_day;

  switch(Simul->date_format){
  case FR_TS:
    sscanf($1,"%d/%d/%d", &pd->dd, &pd->mm,&pd->yyyy);
    break;
  case US_TS:
    sscanf($1,"%d/%d/%d", &pd->mm, &pd->dd,&pd->yyyy);
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    sscanf($2,"%d:%d:%d", &pd->hh, &pd->min,&pd->ss);

  julian_day = TS_date2julian_dd_hm(pd,Simul->chronos->yr0,Simul->poutputs); // pd->mm commence par 1
  julian_day *= NSEC_DAY_TS; // in second
  $$ = TS_create_function(julian_day,$3*unit_f);
}
| LEX_DATE_DAY LEX_DATE_HH_MM flottant
{
  double julian_day;

  switch(Simul->date_format){
  case FR_TS:
    sscanf($1,"%d/%d/%d", &pd->dd, &pd->mm,&pd->yyyy);
    break;
  case US_TS:
    sscanf($1,"%d/%d/%d", &pd->mm, &pd->dd,&pd->yyyy);
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    sscanf($2,"%d:%d", &pd->hh, &pd->min);
    pd->ss = 0.;
  julian_day = TS_date2julian_dd_hm(pd,Simul->chronos->yr0,Simul->poutputs); // pd->mm commence par 1
  julian_day *= NSEC_DAY_TS; // in second
  $$ = TS_create_function(julian_day,$3*unit_f);
}
| LEX_DATE_DAY LEX_INT flottant
{
  double julian_day;

  switch(Simul->date_format){
  case FR_TS:
    sscanf($1,"%d/%d/%d", &pd->dd, &pd->mm,&pd->yyyy);
    break;
  case US_TS:
    sscanf($1,"%d/%d/%d", &pd->mm, &pd->dd,&pd->yyyy);
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    pd->hh = $2;
    pd->min = 0.;
    pd->ss = 0.;
  julian_day = TS_date2julian_dd_hm(pd,Simul->chronos->yr0,Simul->poutputs); // pd->mm commence par 1
  julian_day *= NSEC_DAY_TS; // in second
  $$ = TS_create_function(julian_day,$3*unit_f);
}
| LEX_DATE_DAY flottant
{
  double julian_day;

  switch(Simul->date_format){
  case FR_TS:
    sscanf($1,"%d/%d/%d", &pd->dd, &pd->mm,&pd->yyyy);
    break;
  case US_TS:
    sscanf($1,"%d/%d/%d", &pd->mm, &pd->dd,&pd->yyyy);
    break;
  default: LP_error(Simul->poutputs,"In %s%4.2f -> Encoded date unknown, impossible to read the line %d in %s. Use US_TS or FR_TS\n",CODE_NAME,NVERSION_PROSE,line_nb+1,current_read_files[pile].name);
    break;
   }
    pd->hh = 0.;
    pd->min = 0.;
    pd->ss = 0.;
  julian_day = TS_date2julian_dd_hm(pd,Simul->chronos->yr0,Simul->poutputs); // pd->mm commence par 1
  julian_day *= NSEC_DAY_TS; // in second
  $$ = TS_create_function(julian_day,$2*unit_f);
}
;


brace : LEX_OPENING_BRACE
| LEX_CLOSING_BRACE
;
%%


/* Procedure used to display errors
 * automatically called in the case of synthax errors
 * Specifies the file name and the line where the error was found
 */
	
#if defined GCC482 || defined GCC472 ||  defined GCC471
void yyerror(char const *s)
#else 
void yyerror(char *s)
#endif
{
  if (pile >= 0) {
    //LP_printf(Simul->poutputs,"File %s, line %d : %s\n",current_read_files[pile].name,line_nb+1,s);
    //exit -1;
    LP_error(Simul->poutputs,"File %s, line %d : %s\n",current_read_files[pile].name,line_nb+1,s);
  }
}


void lecture(char *name,FILE *fp) 
{
  pile = 0;
  current_read_files[pile].name = strdup(name);
  if ((yyin = fopen(name,"r")) == NULL)
    LP_error(fp,"File %s doesn't exist\n",name);
  current_read_files[pile].address = yyin;
  current_read_files[pile].line_nb = line_nb;

   LP_printf(fp,"\n****Reading input data****\n");
  yyparse();
  LP_printf(fp,"\n****Input data read****\n");
}

