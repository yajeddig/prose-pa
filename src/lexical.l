/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: lexical.l
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

/********************************************************/
/*             Lexical words used in HyCu               */
/********************************************************/


%x incl str incl_str variable incl_variable unit

%{

#include <stdio.h>
#include <strings.h>
#include <string.h>
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
#include "input.h"




/* Local functions declaration */
void include_function(char *);
/*void yyerror(char *);*/
int yylex();

char *pname;



%}

%a 10000
%p 10000
%o 15000
%e 5000
%n 2000

number	[0-9]
to_power [DdEe][-+]?{number}+
alnum [0-9A-Za-z]

%%
#[^\n]*\n 	{ line_nb++;}	/* To not read the comments */
![^\n]*\n 	{ line_nb++;}	/* To not read the comments */
\n	{ line_nb++;}
[ \t]+ ;			/* To ignore spaces and tabulations */

\$ {BEGIN(variable);} 
\[ {BEGIN(unit); return LEX_OPENING_BRACKET;}
[Ii][Nn][Cc][Ll][Uu][DdRr][Ee] {BEGIN(incl);}
[Nn][oO] {yylval.integer = NO_TS; return LEX_ANSWER;}
[Yy][Ee][Ss] {yylval.integer = YES_TS; return LEX_ANSWER;}
[Gg][Ii][Vv][Ee][Nn] {yylval.integer = GIVEN; return LEX_ANSWER;}
"=" return LEX_EQUAL;
"," return LEX_COMA;
\{ return LEX_OPENING_BRACE;
\} return LEX_CLOSING_BRACE;
":" return LEX_COLON;
";" return LEX_SEMI_COLON;
"->" return LEX_ARROW;
">" return LEX_ARROW; 
"<-" return LEX_REVERSE_ARROW; 
"<" return LEX_REVERSE_ARROW; 

[Ii][Nn][Pp][Uu][Tt]"_"[Ff][Oo][Ll][Dd][Ee][Rr][Ss] return LEX_INPUT_FOLDER;
[Oo][Uu][Tt][Pp][Uu][Tt]"_"[Ff][Oo][Ll][Dd][Ee][Rr] return LEX_OUTPUT_FOLDER;

[Ss][Ii][Mm][Uu][Ll][Aa][Tt][Ii][Oo][Nn] return LEX_SIMUL;

[Ii][Nn][Ii][Tt][Ii][Aa][Ll][Ii][Zz][Aa][Tt][Ii][Oo][Nn] return LEX_INIT;
[Ii][Nn][Ii][Tt]"_"[Zz]"_"[Ff][Ii][Ll][Ee] return LEX_INIT_Z;
[Ii][Nn][Ii][Tt]"_"[Qq]"_"[Ff][Ii][Ll][Ee] return LEX_INIT_Q;

[Tt][Ii][Mm][Ee] return LEX_CHRONOS;
[Dd][Tt] return LEX_TIMESTEP;
[Nn][Ss][Tt][Ee][Pp]"_"[Dd][Aa] return LEX_TIMESTEP_DA;
[Ss][Ee][Uu][Ii][Ll]"_"[Cc][Hh][Ll][Aa] return LEX_SEUIL_CHLA;
[Aa][Ll][Pp][Hh][Aa]"_"[Dd][Aa] return LEX_ALPHA_DA;
[Ee][Rr][Rr][Oo][Rr]"_"[Oo][Bb][Ss] return LEX_ERROR_OBS;
[Ss]"_"[Pp][Ee][Rr][Cc][Ee][Nn][Tt] return LEX_S_PERCENT;
[Rr][Aa][Nn][Dd][Oo][Mm]"_"[Ww][Aa][Ll][Kk] return LEX_RANDOM_WALK;
[Ll][Oo][Oo][Pp] {yylval.integer = LOOP; return LEX_RANDOM_WALK_VAL;}
[Nn][Oo][Tt]"_"[Ll][Oo][Oo][Pp] {yylval.integer = NOT_LOOP; return LEX_RANDOM_WALK_VAL;}
[Pp][Ff] {yylval.integer = PF_PROSE; return LEX_DA_METHOD_VAL;}
[En][Nn][Kk][Ff] {yylval.integer = ENKF_PROSE; return LEX_DA_METHOD_VAL;}
[Dd][Aa]"_"[Mm][Ee][Tt][Hh][Oo][Dd] return LEX_DA_METHOD;
[Tt]"_"[Ii][Nn][Ii] {yylval.integer = BEGINNING; return LEX_TIME;}
[Tt]"_"[Ee][Nn][Dd] {yylval.integer = END; return LEX_TIME;}
[Dd][Aa][Tt][Ee]"_"[Ff][Oo][Rr][Mm][Aa][Tt] return LEX_DATEFORMAT;
[Ff][Rr]"_"[Tt][Ss] {yylval.integer = FR_TS; return LEX_DATEFORMAT_FR_TS;}
[Uu][Ss]"_"[Tt][Ss] {yylval.integer = US_TS; return LEX_DATEFORMAT_FR_TS;}
[Yy][Ee][Aa][Rr][0]"_"[Jj][Uu][Ll][Ii][Aa][Nn] return LEX_YEAR0;
[Yy][Ee][Aa][Rr][0]"_"[Mm][Ee][Tt][Ee][Oo] return LEX_YEAR0_METEO;
[0-9]+"/"[0-9]+"/"[0-9]+ {yylval.string = strdup(yytext); return LEX_DATE_DAY;}
[0-9]+":"[0-9]+":"[0-9]+ {yylval.string = strdup(yytext); return LEX_DATE_HH_MM_SS;}
[0-9]+":"[0-9]+ {yylval.string = strdup(yytext); return LEX_DATE_HH_MM;}


[Bb][Ee][Gg][Ii][Nn][Nn][Ii][Nn][Gg] return LEX_BEGIN;
[Jj][Aa][Nn][Uu][Aa][Rr][Yy] {yylval.integer = JANUARY; return LEX_MONTH;}
[Ff][Ee][Bb][Rr][Uu][Aa][Rr][Yy] {yylval.integer = FEBRUARY; return LEX_MONTH;}
[Mm][Aa][Rr][Cc][Hh] {yylval.integer = MARCH; return LEX_MONTH;}
[Aa][Pp][Rr][Ii][Ll] {yylval.integer = APRIL; return LEX_MONTH;}
[Mm][Aa][Yy] {yylval.integer = MAY; return LEX_MONTH;}
[Jj][Uu][Nn][Ee] {yylval.integer = JUNE; return LEX_MONTH;}
[Jj][Uu][Ll][Yy] {yylval.integer = JULY; return LEX_MONTH;}
[Aa][Uu][Gg][Uu][Ss][Tt] {yylval.integer = AUGUST; return LEX_MONTH;}
[Ss][Ee][Pp][Tt][Ee][Mm][Bb][Ee][Rr] {yylval.integer = SEPTEMBER; return LEX_MONTH;}
[Oo][Cc][Tt][Oo][Bb][Ee][Rr] {yylval.integer = OCTOBER; return LEX_MONTH;}
[Nn][Oo][Vv][Ee][Mm][Bb][Ee][Rr] {yylval.integer = NOVEMBER; return LEX_MONTH;}
[Dd][Ee][Cc][Ee][Mm][Bb][Ee][Rr] {yylval.integer = DECEMBER; return LEX_MONTH;}

[Ss][Pp][Ee][Cc][Ii][Ee][Ss] return LEX_SPECIES;
[Tt][Rr][Aa][Nn][Ss][Pp][Oo][Rr][Tt]"_"[Tt][Yy][Pp][Ee] return LEX_TRANSPORT_EQ;
[Hh][Ee][Aa][Tt] {yylval.integer = HEAT_TTC; return LEX_TRANSPORT_MODE;}
[Ss][Oo][Ll][Uu][Tt][Ee] {yylval.integer = SOLUTE_TTC; return LEX_TRANSPORT_MODE;}
[Pp][Oo][Rr][Oo][Uu][Ss] {yylval.integer = POROUS_TTC; return LEX_MEDIA_TYPE;}
[Ff][Rr][Ee][Ee]"_"[Ww][Aa][Tt][Ee][Rr] {yylval.integer = WATER_TTC; return LEX_MEDIA_TYPE;}
[Mm][Ee][Dd][Ii][Aa]"_"[Tt][Yy][Pp][Ee] return LEX_MED_TYPE;
[Aa][Dd][Vv][Ee][Cc][Tt][Ii][Oo][Nn] {yylval.integer = CONVECTION_TTC ; return LEX_PROC_TRANSP;}
[Dd][Ii][Ff][Ff][Uu][Ss][Ii][Oo][Nn] {yylval.integer = DIFFUSION_TTC ; return LEX_PROC_TRANSP;}
[Dd][Ii][Ff][Ff]_[Mm][Oo][Ll] return LEX_DIFF_MOL;

[Ss][Ee][Tt][Tt][Ii][Nn][Gg][Ss] return LEX_SET;
[Nn][Dd][Ii][Mm] return LEX_DIM;
[Tt][Yy][Pp][Ee] return LEX_TYPE;
[Rr][Ee][Gg][Ii][Mm][Ee] return LEX_REGIME;
[Tt][Rr][Aa][Nn][Ss][Ii][Ee][Nn][Tt] {yylval.integer = TRANSIENT; return LEX_CALC_STATE;}
[Ss][Tt][Ee][Aa][Dd][Yy] {yylval.integer = STEADY; return LEX_CALC_STATE;}
[Nn][Uu][Mm]"_"[Pp][Aa][Rr][Tt][Ii][Cc][Uu][Ll][Ee][Ss] return LEX_NUM_PARTICULES;

[Hh][Yy][Dd] {yylval.integer = HYD; return LEX_CALC_MODE;}
[Tt][Tt][Cc] {yylval.integer = TTC; return LEX_CALC_MODE;}
[Rr][Ii][Vv][Ee] {yylval.integer = RIVE; return LEX_CALC_MODE;}
[Dd][Aa] {yylval.integer = DA; return LEX_CALC_MODE;}
[Tt][Uu][Bb][Ee] {yylval.integer = TUB; return LEX_CALC_MODE;}
[Dd][Yy][Nn][Aa][Mm][Ii][Cc][Aa][Ll]"_"[Tt] {yylval.integer = H_T; return LEX_CALC_MODE;}
[Dd][Yy][Nn][Aa][Mm][Ii][Cc][Aa][Ll]"_"[Tt][Ee][Mm][Pp][Ee][Rr][Aa][Tt][Uu][Rr][Ee] {yylval.integer = H_T; return LEX_CALC_MODE;}
[Dd][Yy][Nn][Aa][Mm][Ii][Cc]"_"[Tt][Ee][Mm][Pp][Ee][Rr][Aa][Tt][Uu][Rr][Ee] {yylval.integer = H_T; return LEX_CALC_MODE;}
[Dd][Yy][Nn][Aa][Mm][Ii][Cc]"_"[Tt] {yylval.integer = H_T; return LEX_CALC_MODE;}
[Gg][Cc]"_"[Pp][Rr][Oo][Ss][Ee] {yylval.integer = GC_PROSE; return LEX_SOLVER_TYPE;}
[Ss][Pp]"_"[Pp][Rr][Oo][Ss][Ee] {yylval.integer = SP_PROSE; return LEX_SOLVER_TYPE;}
[Ss][Oo][Ll][Vv][Ee][Rr] return LEX_SOLVER;
[Nn][Uu][Mm]"_"[Tt][Hh][Rr][Ee][Aa][Dd] return LEX_NUM_THREAD;
[Nn][Uu][Mm]"_"[Tt][Hh][Rr][Ee][Aa][Dd]"_"[Pp][Aa][Rr] return LEX_NUM_THREAD_PAR;
[Nn][Uu][Mm]"_"[Tt][Uu][Bb][Ee]"_"[Dd][Ee][Ff] return LEX_NUM_TUBE;
[Cc][Aa][Ll][Cc][Uu][Ll][Aa][Tt][Ee]"_"[Cc][Uu][Rr][Vv][Aa][Tt][Uu][Rr][Ee] return LEX_CALC_CURVATURE;
[Uu][Pp][Ss][Tt][Rr][Ee][Aa][Mm]"_"[Hh][Mm][Ii][Nn] {yylval.integer = UP_HMIN;return LEX_GENERALP;}
[Dd][Oo][Ww][Nn][Ss][Tt][Rr][Ee][Aa][Mm]"_"[Hh][Mm][Aa][Xx] {yylval.integer = DOWN_HMAX;return LEX_GENERALP;}
[Dd][Xx] {yylval.integer = DX;return LEX_GENERALP;}
[Dd][Yy] {yylval.integer = DY;return LEX_GENERALP;}
[Dd][Zz] {yylval.integer = DZ;return LEX_GENERALP;}
[Cc][Uu][Rr][Vv][Aa][Tt][Uu][Rr][Ee] {yylval.integer = CURVATURE;return LEX_GENERALP;}
[Gg][Ll][Oo][Bb][Aa][Ll]"_"[Ss][Tt][Rr][Ii][Cc][Kk][Ll][Ee][Rr] {yylval.integer = STRICKLER;return LEX_GENERALP;}
[Tt][Hh][Ee][Tt][Aa] {yylval.integer = THETA;return LEX_GENERALP;}
[Ee][Pp][Ss]"_"[Qq] {yylval.integer = EPS_Q;return LEX_GENERALP;}
[Ee][Pp][Ss]"_"[Zz] {yylval.integer = EPS_Z;return LEX_GENERALP;}

[Hh][Yy][Dd][Rr][Aa][Uu][Ll][Ii][Cc] return LEX_HYDRO;
[Tt][Rr][Aa][Nn][Ss][Pp][Oo][Rr][Tt] return LEX_TRANSPORT;
[Tt][Hh][Ee][Tt][Aa]_[Tt][Tt][Cc] {return LEX_THETA;}
[Bb][Ii][Oo][Ll][Oo][Gg][Yy] return LEX_BIOLOGY;

[Pp][Aa][Rr][Aa][Mm]"_"[Rr][Aa][Nn][Gg][Ee]"_"[Dd][Aa] return LEX_PARAM_RANGE;
[Rr][Mm]"_"[Pp][Hh][Yy] {yylval.integer = MAINT_DA; return LEX_PARAM_DA;}
[Aa]"_"[Pp][Hh][Yy] {yylval.integer = ALPHA_DA; return LEX_PARAM_DA;}
[Pp][Mm][Aa][Xx]"_"[Pp][Hh][Yy] {yylval.integer = PMAX_DA; return LEX_PARAM_DA;}
[Ee][Tt][Aa]"_"[Cc][Hh][Ll][Aa]"_"[Pp][Hh][Yy] {yylval.integer = ETA_CHLA_DA; return LEX_PARAM_DA;}
[Cc]"_"[Cc][Hh][Ll][Aa]"_"[Pp][Hh][Yy] {yylval.integer = C_CHLA_DA; return LEX_PARAM_DA;}
[Ee][Tt][Aa]"_"[Ww][Aa][Tt][Ee][Rr] {yylval.integer = ETA_DA; return LEX_PARAM_DA;}
[Tt][Oo][Pp][Tt]"_"[Pp][Hh][Yy] {yylval.integer = TOPT_PHY_DA; return LEX_PARAM_DA;}
[Mm][Uu][Mm][Aa][Xx]"_"[Bb][Aa][Cc][Tt] {yylval.integer = MU_BACT_DA; return LEX_PARAM_DA;}
[Yy][Ii][Ee][Ll][Dd]"_"[Bb][Aa][Cc][Tt] {yylval.integer = Y_BACT_DA; return LEX_PARAM_DA;}
[Mm][Oo][Rr][Tt]"_"[Bb][Aa][Cc][Tt] {yylval.integer = MORT_BACT_DA; return LEX_PARAM_DA;}
[Tt][Oo][Pp][Tt]"_"[Bb][Aa][Cc][Tt] {yylval.integer = TOPT_BACT_DA; return LEX_PARAM_DA;}
[Kk][Rr][Ee][Aa]"_"[Nn][Aa][Vv][Ii][Gg] {yylval.integer = KNAVIG_DA; return LEX_PARAM_DA;}
[Bb][1]"_"[Rr][Ii][Vv][Ee][Rr] {yylval.integer = B1_RIVER_DA ; return LEX_PARAM_DA;}


[Nn][Uu][Mm][Ee][Rr][Ii][Cc][Aa][Ll]"_"[Mm][Ee][Tt][Hh][Oo][Dd] return LEX_NUM;
[Ee][Xx][Pp][Ll][Ii][Cc][Ii][Tt] {yylval.integer = EXPLICIT; return LEX_METHOD;}
[Rr][Uu][Nn][Gg][Ee]"_"[Kk][Uu][Tt][Tt][Aa] {yylval.integer = RUNGE_KUTTA; return LEX_METHOD;}
[Mm][Aa][Xx]"_"[Dd][Ii][Vv]"_"[Dd][Tt] return LEX_DIVDT;
[Ee][Pp][Ss][Ii][Ll][Oo][Nn] return LEX_EPS;
[Nn][Bb]"_"[Cc][Oo][Mm][Pp]"_"[Pp][Hh][Yy] return LEX_COMPPHY;
[Pp][Hh][Yy][2][Pp][Hh][Yy][Ff] {yylval.integer = PHYF; return LEX_PHY2;}
[Pp][Hh][Yy][2][Pp][Hh][Yy][Rr] {yylval.integer = PHYR; return LEX_PHY2;}
[Pp][Hh][Yy][2][Pp][Hh][Yy][Ss] {yylval.integer = PHYS; return LEX_PHY2;}
[Dd][Zz]_[Rr][Ii][Vv][Ee] return LEX_DZ_RIVE;
[Dd][Bb][Oo]"_"[Oo][Xx][Yy] return LEX_DBO;
[Ll][Ii][Mm][Ii][Tt]"_"[Ff][Aa][Cc][Tt][Oo][Rr] return LEX_LIMIT_FACTOR;
[Ll][Ii][Mm]"_"[Mm][Ii][Nn] {yylval.integer = LIM_MIN; return LEX_LIMIT_FACTOR_VAL;}
[Mm][Uu][Ll][Tt][Ii][Pp][Ll][Ii][Cc][Aa][Tt][Ii][Oo][Nn] {yylval.integer = MULTIPLICATION; return LEX_LIMIT_FACTOR_VAL;}

[Mm][Ee][Tt]"_"[Cc][Oo][Rr][Rr][Ee][Ss][Pp] return LEX_ID;
[Mm][Ee][Tt]"_"[Dd][Aa][Tt][Aa] return LEX_ID;
[Ff][Oo][Ll][Dd][Ee][Rr] return LEX_FOLDER;
[Ff][Oo][Rr][Cc][Ii][Nn][Gg][Ss] return LEX_ID;
[Tt]"_"[Aa][Ii][Rr] {yylval.integer = T_AIR; return LEX_MET;}
[Hh]"_"[Ss][Ww] {yylval.integer = H_SW; return LEX_MET;}
[Hh]"_"[Ll][Ww] {yylval.integer = H_LW; return LEX_MET;}
[Hh][Uu][Mm][Ii][Dd] {yylval.integer = HUMID; return LEX_MET;}
[P]"_"[Aa][Tt][Mm] {yylval.integer = P_ATM; return LEX_MET;}
[Uu]"_"[Ww][Ii][Nn][Dd] {yylval.integer = U_WIND; return LEX_MET;}
[Rr][Ii][Vv][Ee][Rr][Ss][Ii][Dd][Ee]"_"[Pp][Aa][Rr][Aa][Mm][Ee][Tt][Ee][Rr] return LEX_ID;

[Ii][Nn][Ii][Tt]"_"[Tt]"_"[Ff][Ii][Ll][Ee] return LEX_INIT_T;
[Hh][Ee][Aa][Tt]"_"[Tt][Rr][Aa][Nn][Ss][Pp][Oo][Rr][Tt] return LEX_HEAT_TRANSPORT;
[Tt][Hh][Ee][Tt][Aa]"_"[Tt] return LEX_THETA_T;
[Dd][Ee][Ff][Aa][Uu][Ll][Tt]"_"[Tt]"_"[Ii][Nn][Ff][Ll][Oo][Ww][Ss] return LEX_DEFAULT_T_INFLOWS;
[Cc][Oo][Dd][Ee]"_"[Tt][Ss] {yylval.integer = CODE_TS; return LEX_CODE_TS;}

[Ee][Xx][Cc][Hh][Aa][Nn][Gg][Ee][Ss] return LEX_EXCH;
[Cc][Aa][Ll][Cc]"_"[Ss][Ee][Dd][Ii][Mm]"_"[Ee][Rr][Oo][Ss] return LEX_CALC_SE;
[Cc][Oo][Nn][Tt][Ii][Nn][Uu][Oo][Uu][Ss]"_"[Mm][Oo][Nn][Oo][Pp][Rr][Oo][Cc][Ee][Ss][Ss] {yylval.integer = MONOPROCESS_C; return LEX_SETYPE;}
[Dd][Ii][Ss][Cc][Oo][Nn][Tt][Ii][Nn][Uu][Oo][Uu][Ss]"_"[Mm][Oo][Nn][Oo][Pp][Rr][Oo][Cc][Ee][Ss][Ss] {yylval.integer = MONOPROCESS_D; return LEX_SETYPE;}
[Ss][Ii][Mm][Uu][Ll][Tt][Aa][Nn][Ee][Oo][Uu][Ss] {yylval.integer = SIMULTANEOUS; return LEX_SETYPE;}
[Ss][Ii][Mm][Uu][Ll][Tt][Aa][Nn][Ee][Oo][Uu][Ss]"_"[Pp][Rr][Oo][Ss][Ee] {yylval.integer = SIMULTANEOUS_PROSE; return LEX_SETYPE;}
[Nn][Uu] {yylval.integer = NU; return LEX_PARAM_EROS;}
[Ee][Tt][Aa]"_"[Hh][Yy][Dd] {yylval.integer = ETA_HYD; return LEX_PARAM_EROS;}
[Pp][Nn][Aa][Vv][Ii][Gg] {yylval.integer = PNAVIG; return LEX_PARAM_EROS;}
[Kk]"_"[Ee][Rr][Oo] {yylval.integer = K_ERO; return LEX_PARAM_EROS;}
[Ll][Aa][Mm][Bb][Dd][Aa] {yylval.integer = LAMBDA; return LEX_PARAM_EROS;}
[Cc][Aa][Ll][Cc]"_"[Tt][Rr][Aa][Nn][Ss][Pp][Oo][Rr][Tt]"_"[Cc][Aa][Pp][Aa][Cc][Ii][Tt][Yy] return LEX_CALC_TC;
[Ss][Hh][Ee][Aa][Rr]"_"[Ss][Tt][Rr][Ee][Ss][Ss] {yylval.integer = SHEAR_STRESS; return LEX_TCTYPE;}
[Ss][Tt][Rr][Ee][Aa][Mm]"_"[Pp][Oo][Ww][Ee][Rr] {yylval.integer = STREAM_POWER; return LEX_TCTYPE;}
[Uu][Nn][Ii][Tt]"_"[Ss][Tt][Rr][Ee][Aa][Mm]"_"[Pp][Oo][Ww][Ee][Rr] {yylval.integer = UNIT_STREAM_POWER; return LEX_TCTYPE;}
[Ss][Ii][Mm][Pp][Ll][Ee]  {yylval.integer = SIMPLE; return LEX_TCTYPE;}
[Kk]"_"[Tt][Cc] {yylval.integer = K_TC; return LEX_TC_PARAM;}
[Bb][Ee][Tt][Aa]"_"[Tt][Cc] {yylval.integer = BETA_TC; return LEX_TC_PARAM;}
[Gg][Aa][Mm][Mm][Aa]"_"[Tt][Cc] {yylval.integer = GAMMA_TC; return LEX_TC_PARAM;}
[Aa][Nn][Nn][Ee][Xx]"_"[Vv][Aa][Rr] return LEX_ANNEXVAR;
[Aa][Dd][Ss][Oo][Rr][Bb][Ee][Dd]"_"[Ss][Pp][Ee][Cc][Ii][Ee][Ss] return LEX_ADS_SPECIES;

[Pp][Hh][Yy] {yylval.integer = PHY; return LEX_ONESPECIES;}
[Zz][Oo][Oo] {yylval.integer = ZOO; return LEX_ONESPECIES;}
[Bb][Aa][Cc][Tt] {yylval.integer = BACT; return LEX_ONESPECIES;}
[Bb][Aa][Cc][Tt][Nn] {yylval.integer = BACTN; return LEX_ONESPECIES;}
[Bb][Ii][Ff] {yylval.integer = BIF; return LEX_ONESPECIES;}
[Mm][Ee][Ss] {yylval.integer = MES; return LEX_ONESPECIES;}
[Mm][Oo][Pp] {yylval.integer = MOP; return LEX_ONESPECIES;}
[Mm][Oo][Dd] {yylval.integer = MOD; return LEX_ONESPECIES;}
[Ss][Ii]"_"[Dd] {yylval.integer = SI_D; return LEX_ONESPECIES;}
[Nn][Hh][4] {yylval.integer = NH4; return LEX_ONESPECIES;}
[Nn][Oo][3] {yylval.integer = NO3; return LEX_ONESPECIES;}
[Pp][Oo][4] {yylval.integer = PO4; return LEX_ONESPECIES;}
[Oo][2] {yylval.integer = O2; return LEX_ONESPECIES;}
[Nn][Oo][2] {yylval.integer = NO2; return LEX_ONESPECIES;}
[Nn][2][Oo] {yylval.integer = N2O; return LEX_ONESPECIES;}

[Mm][Aa][Cc][Rr][Oo][Ss][Pp][Ee][Cc][Ii][Ee][Ss] return LEX_MACROSPECIES; // NF 21/07/2021
[Tt][Oo][Cc] {yylval.integer = TOC; return LEX_TYPEOF_MACROSPECIES;} // MH 21/07/2021
"APOUB" {yylval.integer = APOUB; return LEX_TYPEOF_MACROSPECIES;} // MH 21/07/2021
[Rr][Ee][Ll][Aa][Tt][Ee][Dd]"_"[Tt][Oo] return LEX_RELATED_TO;//NF 21/7/2021
[Tt][Hh][Rr][Ee][Ss][Hh][Oo][Ll][Dd] return LEX_THRESHOLD;//NF 21/7/2021
[Ss][Hh][Aa][Rr][Ee]"_"[Mm][Oo][Pp] {yylval.integer = MACMOP; return LEX_SHARE_MO;} //MH 06/08/2021//NF 19/8/2021
[Ss][Hh][Aa][Rr][Ee]"_"[Mm][Oo][Dd] {yylval.integer = MACMOD; return LEX_SHARE_MO;} //MH 06/08/2021//NF 19/8/2021
[Bb][Ii][Oo][Dd][Ee][Gg][Rr][Aa][Dd][Aa][Bb][Ll][Ee] {yylval.integer = B; return LEX_BIODEG;} // MH 06/08/2021
[Ff][Aa][Ss][Tt]"_"[Bb][Ii][Oo][Dd][Ee][Gg][Rr][Aa][Dd][Aa][Bb][Ll][Ee] {yylval.integer = S; return LEX_BIODEG;} //MH 06/08/2021
[Vv][Aa][Ll] return LEX_VAL; // MH 06/08/2021
[Rr][Aa][Nn][Gg][Ee] return LEX_RANGES; // MH 06/08/2021

[Dd][Aa]"_"[Mm][Oo][Vv]"_"[Aa][Vv][Ee] return LEX_DA_MOV_AVE; // MH 26/04/2022 moving average
[Mm][Oo][Vv]"_"[Aa][Vv][Ee]"_"[Tt][Hh][Ee][Tt][Aa] return LEX_MOV_AVE_THETA; // MH 26/04/2022
[Ll][Nn][Ee][Ii][Gg][Hh] return LEX_LNEIGH;
[Dd][Oo]"_"[Tt][Ii][Mm][Ee]"_"[Ss][Tt][Ee][Pp] return LEX_DO_TIME_STEP; // MH 26/04/2022
[Nn][Vv][Aa][Ll]"_"[Mm][Ii][Nn] return LEX_NVAL_MIN; // MH 14/12/2022



[Ww][Ee][Ii][Gg][Hh][Tt][Ee][Dd]"_"[Ss][Tt][Aa][Tt][Ii][Oo][Nn][Ss] return LEX_WEIGHTED_STATIONS; // MH 02/05/2022
[Cc][Oo][Nn][Ff][Ii][Gg][Uu][Rr][Ee]"_"[Ww][Ee][Ii][Gg][Hh][Tt][Ss] return LEX_CONFIGURE_WEIGHTS; // MH 16/05/2022: station weight parameters
[Vv][Ee][Ll]"_"[Ll][Oo][Ww]"_"[Ff][Ll][Oo][Ww] return LEX_VEL_LOW_FLOW;
[Bb][Aa][Cc][Tt]"_"[Yy][Ii][Ee][Ll][Dd] return LEX_BACT_YIELD;
[Bb][Aa][Cc][Tt]"_"[Gg][Rr][Oo][Ww][Tt][Hh]"_"[Rr][Aa][Tt][Ee] return LEX_BACT_GROWTH_RATE;
[Ff][Aa][Ss][Tt]"_"[Bb][Dd][Oo][Mm] return LEX_FAST_BDOM;
[Bb][Aa][Cc][Tt]"_"[Bb][Ii][Oo][Mm][Aa][Ss][Ss] return LEX_BACT_BIOMASS;
[Oo][Mm]"_"[Ff][Ll][Uu][Xx]"_"[Tt][Hh][Rr][Ee][Ss][Hh][Oo][Ll][Dd] return LEX_OM_FLUX_THRESHOLD;
[Oo][Mm]"_"[Dd][Ee][Cc][Aa][Yy]"_"[Rr][Aa][Tt][Ee] return LEX_OM_DECAY_RATE;


[Ss][Ee][Dd][Vv][Oo][Ll] {yylval.integer = SED_VOL_IO; return LEX_SED_VAR;}
[Ss][Ee][Dd][Hh] {yylval.integer = SED_H_IO; return LEX_SED_VAR;}

[Tt][Ww] {yylval.integer = TW_IO; return LEX_TEMPSPECIES;}
[Tt][Aa] {yylval.integer = TA_IO; return LEX_TEMPSPECIES;}

[Cc] {yylval.integer = C; return LEX_COMP;}
[Oo] {yylval.integer = O; return LEX_COMP;}
[Nn] {yylval.integer = N_RIVE; return LEX_COMP;}
[Pp] {yylval.integer = P; return LEX_COMP;}
[Ss][Ii] {yylval.integer = SI; return LEX_COMP;}
[Hh] {yylval.integer = H_RIVE; return LEX_COMP;}

[Pp][Hh][Yy][Ff] {yylval.integer = PHYF; return LEX_ONEANNEXVAR;}
[Pp][Hh][Yy][Rr] {yylval.integer = PHYR; return LEX_ONEANNEXVAR;}
[Pp][Hh][Yy][Ss] {yylval.integer = PHYS; return LEX_ONEANNEXVAR;}
[Cc][Hh][Ll][Aa] {yylval.integer = CHLA; return LEX_ONEANNEXVAR;}
[Cc][Aa][Rr][Bb][Oo][Nn][Ee] {yylval.integer = CARBONE; return LEX_ONEANNEXVAR;}
[Aa][Cc][Tt][Bb][Aa][Cc][Tt] {yylval.integer = ACTBACT; return LEX_ONEANNEXVAR;}
[Aa][Cc][Tt][Bb][Nn] {yylval.integer = ACTBN; return LEX_ONEANNEXVAR;}
[Nn][Tt][Oo][Tt] {yylval.integer = NTOT; return LEX_ONEANNEXVAR;}
[Pp][Tt][Oo][Tt] {yylval.integer = PTOT; return LEX_ONEANNEXVAR;}

[Pp][Hh][Oo][Tt][Oo][Ss][Yy][Nn][Tt][Hh][Ee][Ss][Ii][Ss] {yylval.integer = PHOT; return LEX_PROCESS;}
[Rr][Ee][Ss][Pp][Ii][Rr][Aa][Tt][Ii][Oo][Nn] {yylval.integer = RESP; return LEX_PROCESS;}
[Gg][Rr][Aa][Zz][Ii][Nn][Gg] {yylval.integer = GRAZING; return LEX_PROCESS;}
[Mm][Oo][Rr][Tt][Aa][Ll][Ii][Tt][Yy] {yylval.integer = MORTALITY; return LEX_PROCESS;}
[Rr][Ee][Aa][Ee][Rr][Aa][Tt][Ii][Oo][Nn] {yylval.integer = REA; return LEX_PROCESS;}
[Gg][Rr][Oo][Ww][Tt][Hh] {yylval.integer = GROWTH; return LEX_PROCESS;}
[Ee][Xx][Cc][Rr][Ee][Tt][Ii][Oo][Nn] {yylval.integer = EXCR; return LEX_PROCESS;}
[Hh][Yy][Dd][Rr][Oo][Ll][Yy][Ss][Ii][Ss] {yylval.integer = HYDROLYSIS; return LEX_PROCESS;}
[Aa][Dd][Ss]"_"[Dd][Ee][Ss][Oo][Rr][Pp][Tt][Ii][Oo][Nn] {yylval.integer = ADS_DESORPTION; return LEX_PROCESS;}
[Rr][Ee][Aa][Cc][Tt][Ii][Oo][Nn][Ss] return LEX_REACTIONS;
[Rr][Ee][Aa][Cc][Tt][Ii][Oo][Nn] return LEX_REACTION;

[Dd][Ee][Ff][Aa][Uu][Ll][Tt]"_"[Cc]"_"[Ii][Nn][Ff][Ll][Oo][Ww][Ss] return LEX_DEFAULT_C_INFLOWS;

[Dd][Ss] {yylval.integer = DS; return LEX_PARAMD;}

[Rr][Hh][Oo] {yylval.integer = RHO; return LEX_PARAMP;}
[Pp][Hh][Ii] {yylval.integer = PHI; return LEX_PARAMP;}
[Vv][Ss][Ee][Dd] {yylval.integer = VSED; return LEX_PARAMP;}
[Ll][Ii][Ff][Tt][Ii][Nn][Gg] {yylval.integer = ARRACHE; return LEX_PARAMP;}

[Tt][Oo][Pp][Tt] {yylval.integer = TOPT; return LEX_PARAML;}
[Ss][Ii][Gg][Mm][Aa] {yylval.integer = SIGMA; return LEX_PARAML;}
[Dd][Bb][Oo] {yylval.integer = DBO; return LEX_PARAML;}
[Pp][Rr][Oo][Dd]"_"[Nn][2][Oo] {yylval.integer = PROD_N2O; return LEX_PARAML;}

[Cc][Oo][Mm][Pp][Oo][Ss][Ee]"_"[Mm][Oo] {yylval.integer = COMPOSE_MO; return LEX_PARAMM;}

[Ee][Tt][Aa] {yylval.integer = ETA; return LEX_PHOT;}
[Ee][Tt][Aa]"_"[Mm][Ee][Ss] {yylval.integer = ETA_MES; return LEX_PHOT;}
[Ee][Tt][Aa]"_"[Cc][Hh][Ll][Aa] {yylval.integer = ETA_CHLA; return LEX_PHOT;}
[Aa][Ll][Pp][Hh][Aa] {yylval.integer = A_RIVE; return LEX_PHOT;}
[Bb][Ee][Tt][Aa] {yylval.integer = B_RIVE; return LEX_PHOT;}
[Pp][Mm][Aa][Xx] {yylval.integer = PMAX20; return LEX_PHOT;}

[Ss][Rr] {yylval.integer = SR20; return LEX_GROWTH;}
[Cc][Rr] {yylval.integer = CR20; return LEX_GROWTH;}
[Mm][Uu][Mm][Aa][Xx] {yylval.integer = MUMAX20; return LEX_GROWTH;}
[Yy][Ii][Ee][Ll][Dd] {yylval.integer = YIELD; return LEX_GROWTH;}
[Pp][Hh][Ii]"_"[Ll][Ii][Mm]"_"[Gg][Rr] {yylval.integer = PHI_LIM_GR; return LEX_GROWTH;}
[Yy][Ii][Ee][Ll][Dd]"_"[Nn][Ii][Tt] {yylval.integer = YIELD_NIT; return LEX_GROWTH;}
[Ss][Tt][Oo][Ee][Cc][Hh][Ii][Oo]"_"[Nn][Ii][Tt] {yylval.integer = STOECHIO_NIT; return LEX_GROWTH;}
[Kk][Mm][Ii][Cc][Hh]"_"[Nn][Oo][2] {yylval.integer = KMICH_NO2; return LEX_GROWTH;}
[Cc][Mm][Ii][Nn]"_"[Oo][2] {yylval.integer = CMIN_O2; return LEX_GROWTH;}
[Ff]"_"[Bb][Ii][Ss] {yylval.integer = F_BIS; return LEX_GROWTH;}
[Cc][Ll][Ii][Mm]"_"[Bb][Aa][Cc][Tt] {yylval.integer = CLIM_BACT; return LEX_GROWTH;}
[Dd][Ee][Ll][Tt][Aa]"_"[Bb][Aa][Cc][Tt] {yylval.integer = DELTA_BACT; return LEX_GROWTH;}
[Kk]"_"[Bb][Aa][Cc][Tt] {yylval.integer = K_BACT; return LEX_GROWTH;}

[Mm][Aa][Ii][Nn][Tt] {yylval.integer = MAINT20; return LEX_RESP;}
[Ee][Nn][Ee][Rr][Gg] {yylval.integer = ENERG; return LEX_RESP;}

[Mm][Oo][Rr][Tt] {yylval.integer = MORT20; return LEX_MORT;}
[Pp][Hh][Ii]"_"[Ll][Ii][Mm]"_"[Mm][Oo][Rr][Tt] {yylval.integer = PHI_LIM_MORT; return LEX_MORT;}
[Dd][Ee][Ll][Tt][Aa] {yylval.integer = DELTA; return LEX_MORT;}
[Kk][Dd] {yylval.integer = KD; return LEX_MORT;}
[Kk][Pp] {yylval.integer = KP; return LEX_MORT;}
[Aa][Ll][Pp][Hh][Aa]"_"[Bb][Ii][Ff] {yylval.integer = ABIF; return LEX_MORT;}
[Mm]"_"[Bb][Ii][Ff] {yylval.integer = MBIF; return LEX_MORT;}
[Cc]"_"[Rr][Ee][Ff] {yylval.integer = CREF; return LEX_MORT;}

[Gg][Rr][Aa][Zz] {yylval.integer = GRAZ20; return LEX_GRAZ;}

[Ee][Xx][Cc][Rr]"_"[Cc][Ss][Tt] {yylval.integer = EXCR_CST; return LEX_EXCR;}
[Ee][Xx][Cc][Rr]"_"[Pp][Hh][Oo][Tt] {yylval.integer = EXCR_PHOT; return LEX_EXCR;}

[Rr][Ee][Aa]"_"[Nn][Aa][Vv][Ii][Gg] {yylval.integer = REA_NAVIG; return LEX_REA;}
[Rr][Ee][Aa]"_"[Ww][Ii][Nn][Dd] {yylval.integer = REA_WIND; return LEX_REA;}
[Ee][Mm] {yylval.integer = EM; return LEX_REA;}
[Dd]"_"[Rr][Ii][Cc][Hh][Ee][Yy] {yylval.integer = D_RICHEY; return LEX_REA;}
[Dd][Ee][Gg][Aa][Ss] {yylval.integer = KDEGAS; return LEX_REA;}


[Kk][Pp][Ss] {yylval.integer = KPS; return LEX_ADS_SENEQUE;}
[Pp][Aa][Cc] {yylval.integer = PAC; return LEX_ADS_SENEQUE;}
[Dd][Aa][Mm][Pp][Ii][Nn][Gg] {yylval.integer = DAMPING; return LEX_ADS_SENEQUE;}
[Aa]"_"[Ff][Rr] {yylval.integer = A_FR; return LEX_ADS_FR;}
[Bb]"_"[Ff][Rr] {yylval.integer = B_FR; return LEX_ADS_FR;}
[Dd]"_"[Ff][Rr] {yylval.integer = D_FR; return LEX_ADS_FR;}
[Tt][Cc]"_"[Ff][Rr] {yylval.integer = TC_FR; return LEX_ADS_FR;}
[Aa][Dd][Ss][Oo][Rr][Bb][Ss]"_"[Oo][Nn] return LEX_ADS_ON;

[Kk][Cc]"_"[Hh][Yy][Dd][Rr] {yylval.integer = KC_HYDR20; return LEX_HYDR;}

[Dd][Ee][Cc][Aa][Yy] {yylval.integer = DECAY; return LEX_RAD_DECAY;}
[Oo][Tt][Hh][Ee][Rr]"_"[Rr][Ee][Aa][Cc][Tt][Oo][Rr][Ss] return LEX_OTHERS;
[Pp][Rr][Oo][Dd][Uu][Cc][Tt][Ss] return LEX_PRODUCTS;
[Cc][Oo][Nn][Dd][Ii][Tt][Ii][Oo][Nn] return LEX_CONDITION;
[Ii][Nn][Ff] {yylval.integer = INF; return LEX_INF_SUP;}
[Ss][Uu][Pp] {yylval.integer = SUP; return LEX_INF_SUP;}

[Kk][Mm][Ii][Cc][Hh] return LEX_KMICH;
[Kk][Ll][Ii][Mm] return LEX_KLIM;
[C][/] return LEX_NUTC;

[Mm][Ee][Tt][Ee][Oo] return LEX_METEO;
[Tt][Ee][Mm][Pp][Ee][Rr][Aa][Tt][Uu][Rr][Ee] return LEX_TEMPERATURE;
[Ww][Ii][Nn][Dd] return LEX_WIND;
[Rr][Aa][Dd][Ii][Aa][Tt][Ii][Oo][Nn] return LEX_RADIATION;
[Pp][Hh][Oo][Tt][Oo][Pp][Ee][Rr][Ii][Oo][Dd] return LEX_PHOTOPERIOD;
[Mm][Ee][Aa][Nn] return LEX_MEAN;
[Aa][Mm][Pp][Ll][Ii][Tt][Uu][Dd][Ee] return LEX_AMPLITUDE;
[Dd][Ee][Ll][Aa][Yy] return LEX_DELAY;
[Aa][Tt][Tt][Ee][Nn][Uu][Aa][Tt][Ii][Oo][Nn]"_"[Ff][Aa][Cc][Tt][Oo][Rr] return LEX_ATTENUATION;


[Ll][Aa][Yy][Ee][Rr][Ss] return LEX_LAYERS;
[Ww][Aa][Tt][Ee][Rr] {yylval.integer = WATER; return LEX_LAYER;}
[Vv][Aa][Ss][Ee] {yylval.integer = VASE; return LEX_LAYER;}
[Pp][Ee][Rr][Ii][Pp][Hh][Yy][Tt][Oo][Nn] {yylval.integer = PERIPHYTON; return LEX_LAYER;}

[Vv]"_"[Vv][Ww][Aa][Tt][Ee][Rr] return LEX_RATIO_WATER;
[Rr][Ee][Tt][Ee][Nn][Tt][Ii][Oo][Nn] return LEX_RET;
[Ss][Cc][Oo][Uu][Rr][Ii][Nn][Gg] return LEX_SCOURING;
[Mm][Aa][Ss][Ss]"_"[Ll][Aa][Yy][Ee][Rr] return LEX_MASS;
[Vv][Oo][Ll][Uu][Mm][Ee]"_"[Ll][Aa][Yy][Ee][Rr] return LEX_VOLUME;
[Ss][Uu][Rr][Ff][Aa][Cc][Ee]"_"[Ll][Aa][Yy][Ee][Rr] return LEX_SURFLAYER;
[Pp][Oo][Rr][Oo][Ss][Ii][Tt][Yy] return LEX_PHI;
[Rr][Hh][Oo]"_"[Ll][Aa][Yy][Ee][Rr] return LEX_RHO;
[Ii][Nn][Ii][Tt][Ii][Aa][Ll]"_"[Cc][Oo][Nn][Cc][Ee][Nn][Tt][Rr][Aa][Tt][Ii][Oo][Nn][Ss] return LEX_IC;
[Ff][Ll][Oo][Ww]"_"[Ii][Nn] {yylval.integer = FLOWIN; return LEX_FLOWS;}
[Ff][Ll][Oo][Ww]"_"[Oo][Uu][Tt] {yylval.integer = FLOWOUT; return LEX_FLOWS;}
[Ff][Ll][Oo][Ww]"_"[Ii][Nn][Tt][Ee][Rr][Ff][Aa][Cc][Ee] {yylval.integer = FLOWINT; return LEX_FLOWS;}
[Ff][Oo][Rr][Cc][Ee][Dd]"_"[Cc][Oo][Nn][Cc][Ee][Nn][Tt][Rr][Aa][Tt][Ii][Oo][Nn][Ss] return LEX_FC;
[Ss][Aa][Tt][Uu][Rr][Aa][Tt][Ii][Oo][Nn] return LEX_SAT;

[Nn][Ee][Tt][Ww][Oo][Rr][Kk]"_"[Dd][Ee][Ss][Cc][Rr][Ii][Pp][Tt][Ii][Oo][Nn] return LEX_NETWORK;

[Ss][Ii][Nn][Gg][Uu][Ll][Aa][Rr][Ii][Tt][Ii][Ee][Ss] return LEX_SING;
[Oo][Pp][Ee][Nn] {yylval.integer = OPEN; return LEX_BC_TYPE;}
[Cc][Oo][Nn][Tt][Ii][Nn][Uu] {yylval.integer = CONTINU; return LEX_BC_TYPE;}
[Dd][Ii][Ss][Cc][Hh][Aa][Rr][Gg][Ee] {yylval.integer = DISCHARGE; return LEX_BC_TYPE;}
[Ww][Aa][Tt][Ee][Rr]"_"[Ll][Ee][Vv][Ee][Ll] {yylval.integer = WATER_LEVEL; return LEX_BC_TYPE;}
[Hh][Yy][Dd][Ww][Oo][Rr][Kk] {yylval.integer = HYDWORK; return LEX_BC_TYPE;}
[Hh][Oo][Ll][Ll][Ee][Rr] return LEX_HOLLER;
[Rr][Aa][Tt][Ii][Nn][Gg]"_"[Cc][Uu][Rr][Vv][Ee] {yylval.integer = RATING_CURVE; return LEX_BC_TYPE;}
[Pp][Kk] return LEX_PK;
[Pp][Oo][Ss][Ii][Tt][Ii][Oo][Nn] return LEX_POSITION;
[Uu][Pp][Ss][Tt][Rr][Ee][Aa][Mm] {yylval.integer = UPSTREAM; return LEX_DIR;}
[Dd][Oo][Ww][Nn][Ss][Tt][Rr][Ee][Aa][Mm] {yylval.integer = DOWNSTREAM; return LEX_DIR;}
[Hh][Yy][Dd]"_"[Ss][Tt][Rr][Uu][Cc][Tt][Uu][Rr][Ee][Ss] return LEX_HYDWORKS;
[Ww][Ee][Ii][Rr] {yylval.integer = WEIR; return LEX_WORK;}
[Gg][Aa][Tt][Ee] {yylval.integer = GATE; return LEX_WORK;}
[Zz]"("[Tt]")" {yylval.integer = ZT;return LEX_FION;}
[Zz][Ww]"("[Tt]")" {yylval.integer = ZWT;return LEX_FION;}
[Hh]"("[Tt]")" {yylval.integer = HT;return LEX_FION;}
[Mm][Uu] {yylval.integer = MU;return LEX_WORK_PARAM;}
[Ww][Ii][Dd][Tt][Hh] {yylval.integer = WIDTH;return LEX_WORK_PARAM;}

[Rr][Ee][Aa][Cc][Hh][Ee][Ss] return LEX_REACH;
[Rr][Ee][Aa][Cc][Hh][Ee][Ss]"_"[Mm][Uu][Ss][Kk] return LEX_REACH_MUSK;
[Ss][Tt][Rr][Ii][Cc][Kk][Ll][Ee][Rr] return LEX_STRICKLER;

[Pp][Rr][Oo][Ff][Ii][Ll][Ee][Ss] return LEX_PROFILES;
[Xx] return LEX_X;
[Yy] return LEX_Y;
[Aa][Bb][Ss][Cc] return LEX_ABSC;
[Tt][Rr][Aa][Jj][Ee][Cc][Tt][Oo][Rr][Yy] return LEX_TRAJECTORY;
[Ss][Hh][Aa][Pp][Ee] return LEX_SHAPE;
[Rr][Ee][Cc][Tt][Aa][Nn][Gg][Uu][Ll][Aa][Rr] {yylval.integer = RECTANGULAR; return LEX_CROSS_SECTION_TYPE;}
[Cc][Oo][Mm][Pp][Ll][Ee][Xx] {yylval.integer = COMPLEX; return LEX_CROSS_SECTION_TYPE;}

[Zz] return LEX_HYD_VAR;
[Qq] return LEX_HYD_VAR;

[Ii][Nn][Ff][Ll][Oo][Ww][Ss] return LEX_INFLOWS;
[Uu][Pp][Ss][Tt][Rr][Ee][Aa][Mm]_[Ii][Nn][Ff][Ll][Oo][Ww] {yylval.integer = UPSTREAM_INFLOW; return LEX_INFLOW_TYPE;}
[Ii][Nn][Ff][Ll][Uu][Ee][Nn][Tt] {yylval.integer = INFLUENT; return LEX_INFLOW_TYPE;}
[Ee][Ff][Ff][Ll][Uu][Ee][Nn][Tt] {yylval.integer = EFFLUENT; return LEX_INFLOW_TYPE;}
[Dd][Ii][Ff][Ff][Uu][Ss][Ee]_[Ii][Nn][Ff][Ll][Oo][Ww] {yylval.integer = DIFFUSE_INFLOW; return LEX_INFLOW_TYPE;}
[Tt] return LEX_TEMP_INFLOW; 

[Oo][Uu][Tt][Pp][Uu][Tt][Ss] return LEX_OUTPUTS;
[Tt][Ii][Mm][Ee]"_"[Uu][Nn][Ii][Tt] return LEX_TUNIT;
[Ii][Tt][Ee][rr][Aa][Tt][Ii][Oo][Nn][Ss] {yylval.integer = ITERATIONS; return LEX_OUTPUT_TYPE;}
[Pp][Rr][Ii][Nn][Tt]"_"[Pp][Kk] {yylval.integer = PRINT_PK; return LEX_OUTPUT_TYPE;}
[Ff][Ii][Nn][Aa][Ll]"_"[Ss][Tt][Aa][Tt][Ee] {yylval.integer = FINAL_STATE; return LEX_OUTPUT_TYPE;}
[Mm][Bb]"_"[Ee][Ll][Ee][Mm][Ee][Nn][Tt][Ss] {yylval.integer = MB_ELE; return LEX_OUTPUT_TYPE;}
[Ff][Ii][Ll][Ee]"_"[Nn][Aa][Mm][Ee] return LEX_FILE_NAME;
[Tt][Ii][Mm][Ee]"_"[Ss][Ee][Rr][Ii][Ee][Ss] {yylval.integer = TRANSV_PROFILE; return LEX_OUTPUT_TYPE;}
[Ll][Oo][Nn][Gg][Ii][Tt][Uu][Dd][Ii][Nn][Aa][Ll]"_"[Pp][Rr][Oo][Ff][Ii][Ll][Ee][Ss] {yylval.integer = LONG_PROFILE; return LEX_OUTPUT_TYPE;}
[Mm][Aa][Ss][Ss]"_"[Bb][Aa][Ll][Aa][Nn][Cc][Ee][Ss] {yylval.integer = MASS_BALANCE; return LEX_OUTPUT_TYPE;}

[Vv][Aa][Rr] return LEX_VAR;
[Pp][Oo][Ii][Nn][Tt][Ss] return LEX_POINTS;
[Pp][Oo][Ii][Nn][Tt][Ss]"_"[Oo][Bb][Ss] return LEX_POINTS_OBS;
[Ee][Xx][Tt][Ee][Nn][Tt] return LEX_EXTENT;
[Gg][Rr][Aa][Pp][Hh][Ii][Cc][Ss] return LEX_GRAPHICS;
[Gg][Nn][Uu][Pp][Ll][Oo][Tt] {yylval.integer = GNUPLOT; return LEX_GRAPH_TYPE;}

[Zz][Ff][Ss] {yylval.integer = ZFS_IO; return LEX_ONE_VAR;}
[Hh][Hh] {yylval.integer = H_IO; return LEX_ONE_VAR;}
[Ss][Uu][Rr][Ff] {yylval.integer = SURF_IO; return LEX_ONE_VAR;}
[Pp][Ee][Rr][Ii] {yylval.integer = PERI_IO; return LEX_ONE_VAR;}
[Rr][Hh] {yylval.integer = RH_IO; return LEX_ONE_VAR;}
[Vv][Ee][Ll] {yylval.integer = VEL_IO; return LEX_ONE_VAR;}
[Ww][Ff][Ss] {yylval.integer = WFS_IO; return LEX_ONE_VAR;}
[Dd][Ii][Ss][Cc][Hh] {yylval.integer = Q_IO; return LEX_ONE_VAR;}

[Tt][Uu][Bb][Ee]"_"[Oo][Uu][Tt][Pp][Uu][Tt] return LEX_TUBE_OUT;
[Mm][Ee][Ss][Hh]"_"[Tt][Uu][Bb][Ee] {yylval.integer = TUBE_MESH; return LEX_TUBE_OUT_TYPE;}
[Hh][Yy][Dd]"_"[Tt][Uu][Bb][Ee] {yylval.integer = TUBE_HYD; return LEX_TUBE_OUT_TYPE;}



[Mm][Aa][Ss][Ss]"_"[Bb][Aa][Ll][Aa][Nn][Cc][Ee][Ss]"_"[Bb][Ii][Oo] return LEX_MB;
[Ee][Nn][Ee][Rr][Gg][Yy]"_"[Bb][Aa][Ll][Aa][Nn][Cc][Ee][Ss]"_"[Hh][Ee][Aa][Tt] return LEX_EB_HEAT;
[Hh][Ee][Aa][Tt]"_"[Uu][Nn][Ii][Tt] return LEX_HEAT_UNIT;

[Nn]"_"[Ss][Tt][Ee][Pp][Ss] return LEX_NSTEPS;

[Kk]"_"[Mm][Uu][Ss][Kk] {yylval.integer = K; return LEX_MUSK_PAR;}
[Xx]"_"[Mm][Uu][Ss][Kk] {yylval.integer = ALPHA; return LEX_MUSK_PAR;}
[Aa][Rr][Ee][Aa] {yylval.integer = AREA; return LEX_CALC_MUSK_PAR;}
[Ss][Ll][Oo][Pp][Ee] {yylval.integer = SLOPE; return LEX_CALC_MUSK_PAR;}
[Mm][Uu][Ss][Kk][Ii][Nn][Gg][Uu][Mm] {yylval.integer = MUSKINGUM; return LEX_DEF_SCHEM_TYPE;}
[Ss][Tt]"_"[Vv][Ee][Nn][Aa][Nn][Tt] {yylval.integer = ST_VENANT; return LEX_DEF_SCHEM_TYPE;}
[Ss][Cc][Hh][Ee][Mm]"_"[Tt][Yy][Pp][Ee] return LEX_SCHEM_TYPE;
[Nn][Ee][Tt][Ww][Oo][Rr][Kk]"_"[Mm][Uu][Ss][Kk] return LEX_NETWORK_MUSK;
[Ee][Ll][Ee][Mm][Ee][Nn][Tt][Ss]"_"[Mm][Uu][Ss][Kk] return LEX_ELE_MUSK;



{number}+ 	|
"-"{number}+ {yylval.integer  = atoi(yytext); return LEX_INT;}
{number}+"."{number}*({to_power})?	| 
"-"{number}+"."{number}*({to_power})?	| 
{number}*"."{number}+({to_power})?	 |
"-"{number}*"."{number}+({to_power})?	{ yylval.real = atof(yytext); return LEX_DOUBLE;}
<incl>[ \t]*    /* skip spaces and tabulations */
<incl>[ \n]*    /* skip newlines */
<incl>\$ {BEGIN(incl_variable);} 
<incl>[^$ \t\n:;,\"]+ {
               char *name;
               name = strdup(yytext);
               include_function(name);
               BEGIN(INITIAL);
 }
\" {/* Allows to read string chains with spaces */ BEGIN(str);} 
<incl>\" {BEGIN(incl_str); /* Beginning of a string chain */ } 
<incl_str>\n { /* A string chain ends with \n */
         LP_error(Simul->poutputs,"HyCu %4.2f -> String chain not ended\n",NVERSION_HYCU);}
<incl_str>[^\\\n\"]+\" {
             int i,k;
             i = strlen(yytext);
             pname = (char *)malloc((i - 1) * sizeof(char));
             for (k = 0; k < i - 1; k++)
                   pname[k] = yytext[k];
             include_function(pname);
             BEGIN(INITIAL);
         }
<str>\" {BEGIN(INITIAL);} 
<str>\n { /* A string chain ends with \n */
         LP_error(Simul->poutputs,"HyCu %4.2f -> String chain not ended\n",NVERSION_HYCU);}
<str>[^\\\n\"]+ {
             yylval.string = strdup(yytext);
             return LEX_NAME;
}

<unit>"/" {return LEX_INV;}
<unit>"^" {return LEX_POW;}

<unit>[Â°][C] {yylval.real = 1.0; return LEX_A_UNIT;}
<unit>"K" {yylval.real = 1.0; return LEX_A_UNIT;}

<unit>"nb" {yylval.real = 1.0; return LEX_A_UNIT;}

<unit>"min" {yylval.real = 60.0; return LEX_A_UNIT;}
<unit>"s" {yylval.real = 1.0; return LEX_A_UNIT;}
<unit>"h" {yylval.real = 3600.0; return LEX_A_UNIT;}
<unit>"d" {yylval.real = 86400.0; return LEX_A_UNIT;}
<unit>"yr" {yylval.real = 31536000; return LEX_A_UNIT;} /* years on the bases of 365 days per year */


<unit>"m" {yylval.real = 1.0; return LEX_A_UNIT;}
<unit>"km" {yylval.real = 1000.0; return LEX_A_UNIT;}
<unit>"cm" {yylval.real = 0.01; return LEX_A_UNIT;}
<unit>"ha" {yylval.real = 10000; return LEX_A_UNIT;}
<unit>"mm" {yylval.real = 0.001; return LEX_A_UNIT;}

<unit>"l" {yylval.real = 0.001; return LEX_A_UNIT;}

<unit>"uE" {yylval.real = 1.0; return LEX_A_UNIT;}
<unit>"J" {yylval.real = J_UE * 0.36; return LEX_A_UNIT;}
<unit>"cal" {yylval.real = J_UE * 0.36 * 4.184; return LEX_A_UNIT;} /* 1 cal = 4.1868 J */

<unit>"ug" {yylval.real = 0.000001; return LEX_A_UNIT;}
<unit>"mg" {yylval.real = 0.001; return LEX_A_UNIT;}
<unit>"g" {yylval.real = 1.; return LEX_A_UNIT;}
<unit>"kg" {yylval.real = 1000.0; return LEX_A_UNIT;}
<unit>"t" {yylval.real = 1000000.0; return LEX_A_UNIT;}

<unit>"umol" {yylval.real = 0.001; return LEX_A_UNIT_MOL;} 
<unit>"mmol" {yylval.real = 1.0; return LEX_A_UNIT_MOL;}
<unit>"mol" {yylval.real = 1000.0; return LEX_A_UNIT_MOL;}

<unit>"N" {yylval.real=1000.0/14.0;return LEX_VAR_UNIT;}
<unit>"NH4" {yylval.real=1000.0/18.0;return LEX_VAR_UNIT;}
<unit>"NO2" {yylval.real=1000.0/46.0;return LEX_VAR_UNIT;}
<unit>"NO3" {yylval.real=1000.0/62.0;return LEX_VAR_UNIT;}
<unit>"P" {yylval.real=1000.0/31.0;return LEX_VAR_UNIT;}
<unit>"PO4" {yylval.real=1000.0/95.0;return LEX_VAR_UNIT;}
<unit>"P2O5" {yylval.real=1000.0/71.0;return LEX_VAR_UNIT;}
<unit>"O2" {yylval.real=1000.0/32.0;return LEX_VAR_UNIT;}
<unit>[Cc][Hh][Ll][Aa] {yylval.real=1.0;return LEX_VAR_UNIT;}
<unit>S[Ii][Oo]2 {yylval.real=1000.0/60.0;return LEX_VAR_UNIT;}
<unit>S[Ii] {yylval.real=1000.0/28.0;return LEX_VAR_UNIT;}
<unit>"C" {yylval.real=1000./12.;return LEX_VAR_UNIT;} 

<unit>"%" {yylval.real = 1.0; return LEX_A_UNIT;}
<unit>"-" {yylval.real = 1.0; return LEX_A_UNIT;}
<unit>{number}+ 	|
<unit>"-"{number}+ {yylval.integer  = atoi(yytext); return LEX_INT;}
<unit>{number}+"."{number}*({to_power})?	| 
<unit>"-"{number}+"."{number}*({to_power})?	| 
<unit>{number}*"."{number}+({to_power})?	 |
<unit>"-"{number}*"."{number}+({to_power})?	{ yylval.real = atof(yytext); return LEX_DOUBLE;}
<unit>\] {
  BEGIN(INITIAL);
  return LEX_CLOSING_BRACKET;
}
<variable>[[:alnum:]]+ {
  /* Inserts the contents of the given variable */
  char *valeur = getenv(yytext);
  if (valeur == NULL) {
    LP_error(Simul->poutputs,"HyCu %4.2f -> File \"%s\", variable \"%s\" undefined, line %d\n",NVERSION_HYCU,current_read_files[pile].name,yytext,line_nb);
  }
  {
    int i;
    for (i = strlen(valeur) - 1; i >= 0; i--)
      unput(valeur[i]);
  }
  BEGIN(INITIAL);
}
<incl_variable>[[:alnum:]]+ {
  /* Inserts the contents of the given variable */
  char *valeur = getenv(yytext);
  if (valeur == NULL) {
    LP_error(Simul->poutputs,"HYCU %4.2f --> File \"%s\", variable \"%s\" undefined, line %d\n",NVERSION_HYCU,current_read_files[pile].name,yytext,line_nb);
  }
  {
    int i;
    for (i = strlen(valeur) - 1; i >= 0; i--)
      unput(valeur[i]);
  }
  BEGIN(incl);
}
[^$ \t\n:;,\"\[\]=]+ {
yylval.string = strdup(yytext);
return LEX_NAME;
}

%%

/* Procedure used to find a file which name is given
* The folders defined by the user in the command file are searched in
* yyin points towards the file read by the parser
*/

void find_file(char *name)
{
   int i;
   char *new_name;

   yyin = fopen(name,"r");
   if (yyin == NULL) {
     for (i = 0; i <folder_nb; i++) {
       new_name = (char *)calloc(strlen(name_out_folder[i]) + strlen(name) + 2,sizeof(char));
       sprintf(new_name,"%s/%s",name_out_folder[i],name);
       /*printf("%s\n",new_name);*/
       yyin = fopen(new_name,"r");
       free(new_name);
       if (yyin != NULL) 
         break;
     }
   }
   if (yyin == NULL)
      LP_error(Simul->poutputs,"HyCu %4.2f -> No file %s\n",NVERSION_HYCU,name);
   
}

/* Procedure used to include a file
*
* - opening of the file which name is in argument
* - storage of the variables related to the file which was read (line_nb, file name, buffer's content)
* - incrementation of the variables of the file pile
*/

void include_file(char *name)
{
int i;
   char *new_name;

   current_read_files[pile].line_nb = line_nb;
   current_read_files[pile].buffer = YY_CURRENT_BUFFER;
   line_nb = 0;
   pile++;
   if (pile >= NPILE)
       LP_error(Simul->poutputs,"HyCu %4.2f -> File %s : Too many files included the ones in the others, maximum = %d\n",NVERSION_HYCU,name,NPILE);
   yyin = fopen(name,"r");
   if (yyin != NULL)
       current_read_files[pile].name = strdup(name);
   else {
     for (i = 0; i <folder_nb; i++) {
       new_name = (char *)calloc(strlen(name_out_folder[i])+strlen(name)+2,sizeof(char));
       sprintf(new_name,"%s/%s",name_out_folder[i],name);
       yyin = fopen(new_name,"r");
       current_read_files[pile].name = strdup(new_name);
       free(new_name);
       if (yyin != NULL) 
         break;
     }
   }
   if (yyin == NULL)
      LP_error(Simul->poutputs,"HyCu %4.2f -> Not possible to find the file %s\n",NVERSION_HYCU,name);
   
   current_read_files[pile].address = yyin;
   yy_switch_to_buffer(yy_create_buffer(yyin,YY_BUF_SIZE));
   LP_printf(Simul->poutputs,"File %s opened\n",current_read_files[pile].name);
   free(name);

}

void include_function(char *name)
{
int i;
   char *new_name;

   current_read_files[pile].line_nb = line_nb;
   current_read_files[pile].buffer = YY_CURRENT_BUFFER;
   line_nb = 0;
   pile++;
   if (pile >= NPILE)
       LP_error(Simul->poutputs,"HyCu %4.2f -> File %s : Too many files included the ones in the others, maximum = %d\n",NVERSION_HYCU,name,NPILE);
   yyin = fopen(name,"r");
   if (yyin != NULL)
       current_read_files[pile].name = strdup(name);
   else {
     for (i = 0; i <folder_nb; i++) {
       new_name = (char *)calloc(strlen(name_out_folder[i])+strlen(name)+2,sizeof(char));
       sprintf(new_name,"%s/%s",name_out_folder[i],name);
       yyin = fopen(new_name,"r");
       current_read_files[pile].name = strdup(new_name);

       if (yyin != NULL) 
         break;
     }
   }
   if (yyin == NULL)
      LP_error(Simul->poutputs,"HyCu %4.2f -> Not possible to find the file %s\n",NVERSION_HYCU,name);
   
   current_read_files[pile].address = yyin;
   yy_switch_to_buffer(yy_create_buffer(yyin,YY_BUF_SIZE));
   LP_printf(Simul->poutputs,"File %s opened\n",current_read_files[pile].name);
   free(name);

}

/* Procedure used to close the files read by the parser
* Called when the parser identifies the end of a file
*
* - closure of the file
* - decrementation of the file pile (pile--)
* - reset of the variables (name, line_nb, buffer) to the values of the last file (if it exists)
* - no more file to read if pile<0
*/

int yywrap()
{
   fclose(current_read_files[pile].address);
   yy_delete_buffer(YY_CURRENT_BUFFER);
   /*if (line_nb == 0)
      printf("\nWARNING : Empty file %s\n",current_read_files[pile].name);*/
   pile--;
   //num = 0;
    if (pile < 0)
      return 1;
   else {
      yyin = current_read_files[pile].address;

      #if defined GCC482 ||defined GCC472 ||  defined GCC471
      yy_switch_to_buffer((YY_BUFFER_STATE) current_read_files[pile].buffer);//gcc4.7
      #elif defined GCC441  || defined GCC442 ||  defined GCC443  || defined GCC444
// || defined GCC445
      yy_switch_to_buffer((yy_buffer_state*)current_read_files[pile].buffer);//gcc4.4
      #else 
// defined GCC445
      yy_switch_to_buffer(current_read_files[pile].buffer);//gcc3.x
      #endif /*test on GCC*/

      line_nb = current_read_files[pile].line_nb;
      return 0;
   }
}
