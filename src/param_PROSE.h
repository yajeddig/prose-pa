/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: param_PROSE.h
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

/* Parameter for error messages with LEX and YACC */
#define YYERROR_VERBOSE 1 
/* Maximum number of lines that can be read in a file */
#define YYMAXDEPTH_PROSE 1500000
/* Maximum number of input files which can be read with include() */
#define NPILE 25
/* Maximum number of characters in a (char *) structure */
#define MAXCHAR_PROSE  STRING_LENGTH_LP
enum random_work_da {LOOP, NOT_LOOP};
enum calc_module {HYD, TTC, RIVE, MB_BIO, DA, H_T, SEB, TUB, EB_HEAT, NMODE}; // AB 25.05.20 -> TUB for libtube. AB 03.10.19 -> H_T: heat transport (temperature variable is activated) & SEB: account for atmospheric conditions in the temperature calculation.
enum calc_solver {GC_PROSE, SP_PROSE, NSOLVER};
enum param_type_assimilation{MAINT_DA, ALPHA_DA, PMAX_DA, ETA_CHLA_DA, C_CHLA_DA, ETA_DA, TOPT_PHY_DA, MU_BACT_DA, Y_BACT_DA, MORT_BACT_DA, TOPT_BACT_DA, KNAVIG_DA, B1_RIVER_DA , NPARAMDA};
enum trophic_state {HORS_BLOOM, BLOOM, NTROPHIC_STATE};
enum param_up_down {PARAM_UP, PARAM_DOWN, PARAMR};
enum type_formule_rea {HOLLER,AVERY};
#define PARAM_PHY MU_BACT_DA
#define NVERSION_PROSE 0.76

#define CODE_NAME "ProSe-PA"


enum output_types {MASS_BALANCE_BIO_OUT = 10, TUBE_MESH, TUBE_HYD, ENERGY_BALANCE_HEAT_OUT, NOUTPUT_TYPES_PROSE}; // SW 30/03/2021 need a bigger value for MASS_BALANCE_BIO_OUT
/* Very important: number of Safran's cells for Seine model,
correspond to the safran's cells stored into the safran's met' files */
#define N_SEINE_SAFRAN 26 // Number of Safran's cells for the "Seine model". May be modified in the future..

enum da_method {PF_PROSE, ENKF_PROSE, DA_METHOD_PROSE};

#define EPS_PROSE 1.e-12
