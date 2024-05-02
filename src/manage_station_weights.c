/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: manage_station_weights.c
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
////#include "global_PROSE.h"
#include "ext_PROSE.h"


// 10/05/2022 MH: function to initiate the structure holding the station weight parameters
s_weight_calc_assim *PROSE_init_calc_stat_weight()
{
  s_weight_calc_assim *pstat_weight;

  pstat_weight = new_weight_calc_assim();
  bzero((char *) pstat_weight, sizeof(s_weight_calc_assim));
  pstat_weight->weighted_stations = NO_TS;
  pstat_weight->vel_low_flow = VEL_LOW_FLOW_RIVE  ;
  pstat_weight->bact_growth_rate = BACTERIA_GROWTH_RATE_RIVE ;
  pstat_weight->bact_yield = BACTERIA_YIELD_RIVE ;
  pstat_weight->fast_bdom_low_flow = FAST_BDOM_LOW_FLOW_RIVE ;
  pstat_weight->bact_biomass = BACT_BIOMASS_RIVE;
  pstat_weight->OM_flux_threshold = OM_FLUX_THRESHOLD_RIVE;
  pstat_weight->OM_decay_rate = OM_DECAY_RATE_RIVE;
  pstat_weight->N_flows_significant = N_FLOWS_SIGNIFICANT_RIVE;

  return pstat_weight;
 
}



/* the function to create weight for observation stations based on their distance and amount of OM from OM sources  */

void Prose_create_obs_station_weight(s_carac_assim *passim, FILE *fp)
{
   
   int i,nobs,j,row;
   double sum = 0, sum_regul=0;

   //TODO: launch a function to fill weigh_calc_for_assim struct at each time step to be dynamic
   
   //input data from the already filled struct
   double velocity = passim->weight_calc_assim->vel_low_flow; //(m/s) shall be taken from calculations at upstream cell
   double yield = passim->weight_calc_assim->bact_yield; //mean of pdf for dynamic
   double miu = passim->weight_calc_assim->bact_growth_rate; // (h^-1) shall be read from bact params reading
   int N_flows = passim->weight_calc_assim->N_flows_significant; // shall be read from number of OM outflows passing the threshold of gC/s
   //int N_flows = passim->weight_calc_assim->N_flows_significant; // shall be read from number of OM outflows passing the threshold of gC/s

   double bact[] = {passim->weight_calc_assim->bact_biomass}; // gC from inflows reading
   double mod_1[] = {passim->weight_calc_assim->fast_bdom_low_flow}; //gC from inflows reading
   double pk_flows[] = {0}  ; //from the equation reading
   double decay_rate[] = {passim->weight_calc_assim->OM_decay_rate};

   /*
   double *bact = passim->weight_calc_assim->bact_biomass; // gC from inflows reading
   double *mod_1 = passim->weight_calc_assim->fast_bdom_low_flow; //gC from inflows reading
   double pk_flows[] = {0}  ; //from the equation reading
    */

   

   // depletion time and distance of each OM source
   double t_depl[N_flows];
   double dist_depl[N_flows];
   for (i = 0; i < N_flows; i++)
        {
            decay_rate[i] = bact[i] * miu / yield;
            t_depl[i] = mod_1[i] / decay_rate[i] * NMIN_HOUR_TS * NSEC_MIN_TS ;
            dist_depl[i] = t_depl[i] * velocity;

	    sum +=  dist_depl[i] * mod_1[i] ;

	    //LP_printf(Simul->poutputs, "t_dep =%f , dist_depl=%f \n" ,  t_depl[i], dist_depl[i] );
	    //LP_printf(Simul->poutputs,"sum = %f \n",sum);

        }
   
   // table to store the distance of each obs_station from each OM_source, -> 0 if the obs_station is upstream of a source
   double dist_i[passim->N_obs][N_flows];
   for (nobs = 0; nobs < passim->N_obs; nobs++)
        {
            for (j=0 ; j < N_flows; j++)
                {
		  if (passim->pobs[nobs]->pk > (pk_flows[j]*1000))
                        {
			  dist_i[nobs][j] = passim->pobs[nobs]->pk - (pk_flows[j]*1000) ;
                        }
                    else {
                            dist_i[nobs][j] = EPS_TS;
                         }
		    LP_printf(Simul->poutputs, "dist_obs_%d_nflow_%d = %f \n",nobs+1,j+1,dist_i[nobs][j]  );
                }
        }

    // arrays to store the weight of each obs_station
    for (i=0; i < passim->N_obs; i++)
     {
       passim->pobs[i]->weight_i = 0.;
       passim->pobs[i]->weight_station = 0.;
     }
   
   double weight_i[] = {0};
   double weight_i_Norm[] = {0};
   for (nobs = 0; nobs < passim->N_obs; nobs++)
        {
            for (j=0 ; j < N_flows; j++)
                {
                    if (dist_i[nobs][j] < dist_depl[j])
                        {
			  passim->pobs[nobs]->weight_i += dist_i[nobs][j] * mod_1[j];
                        }
                    else {
                           passim->pobs[nobs]->weight_i += EPS_TS;
                         }

		    LP_printf(Simul->poutputs, "Weight_obs_%d_nflow_%d (mod_1 = %f) = %f \n",nobs+1,j+1, mod_1[j] ,passim->pobs[nobs]->weight_i);                    
                }
            //sum += passim->pobs[nobs]->weight_i ;
	    //LP_printf(Simul->poutputs,"sum = %f \n",sum);
        }
   // calculating the regularized weights
    for (row = 0; row < passim->N_obs; row++)
        {
	    LP_printf(Simul->poutputs,"Wi_cumul_%d = %f \n",row+1,passim->pobs[row]->weight_i);
            //weight_i_Norm[row] = weight_i[row] / sum ; 
            //passim->pobs[row]->weight_station = weight_i_Norm[row];
	    passim->pobs[row]->weight_station = passim->pobs[row]->weight_i / sum ;

	    sum_regul += passim->pobs[row]->weight_station;

            LP_printf(Simul->poutputs,"Weight_regul @ obs_station_%d = %f \n",row+1, passim->pobs[row]->weight_station);
        }

    /*
    // normalization of the regularized values by the sum of local weights
        for (row = 0; row < passim->N_obs; row++)
        {
	    LP_printf(Simul->poutputs,"Wi_cumul_%d = %f \n",row+1,passim->pobs[row]->weight_i);
            //weight_i_Norm[row] = weight_i[row] / sum ; 
            //passim->pobs[row]->weight_station = weight_i_Norm[row];
	    passim->pobs[row]->weight_station = passim->pobs[row]->weight_station / sum_regul ;

            //LP_printf(Simul->poutputs,"Weight_norm @ obs_station_%d = %f \n",row+1, passim->pobs[row]->weight_station);
        }
    */
     
    // normalization of the regularized values by the max of local weights
    passim->weight_calc_assim->weight_max = TS_max(passim->pobs[0]->weight_station, passim->pobs[1]->weight_station );
        for (row = 0; row < passim->N_obs; row++)
        {
	    LP_printf(Simul->poutputs,"Wi_cumul_%d = %f \n",row+1,passim->pobs[row]->weight_i);
            //weight_i_Norm[row] = weight_i[row] / sum ; 
            //passim->pobs[row]->weight_station = weight_i_Norm[row];
	    passim->pobs[row]->weight_station = passim->pobs[row]->weight_station / passim->weight_calc_assim->weight_max ;

            //LP_printf(Simul->poutputs,"Weight_norm @ obs_station_%d = %f \n",row+1, passim->pobs[row]->weight_station);
        }
    

    //return weight_i_Norm[x];

}
