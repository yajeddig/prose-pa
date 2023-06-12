/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: reoxygenation.c
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
#ifdef OMP
#include <omp.h>
#endif
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

void PROSE_init_reoxy_species(s_species_ttc **p_species,int nspecies,int nele,int np,FILE *fp)
{
	char *name_rive, *name_ttc;
	int e;
	name_rive = Simul->psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[O2][0]->name;
	for(e = 0; e < nspecies; e++)
	{
		name_ttc = p_species[e]->name;
		if(strcmp(name_ttc,name_rive) == 0)
		{
			p_species[e]->plink->preoxy_ttc = TTC_init_reoxy(nele);
			p_species[e]->oxygen = YES_TS;
			break;
		}
	}
}

void Prose_cherche_debit_oxy(s_element_hyd *pele, s_BC_char_hyd **BC_char,int nworks, s_species_ttc *pspecies, int id_abs_ele, FILE *fp)
{
	int i;
	double qd = 0,qn = 0.,qtot;
	double hamont, haval,frac_debit_oxy;
	s_hydro_hyd *phyd, *phydn;
    s_element_hyd *pelen;//élément aval
	
	phyd = pele->center->hydro;
	pelen = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO];
	phydn = pelen->center->hydro;
	qtot = pele->face[X_HYD][TWO]->hydro->Q[T_HYD];
	for(i = 0; i < nworks; i++)
	{
		hamont = phyd->Z[T_HYD] - BC_char[i]->fion_t_old;
		haval = phydn->Z[T_HYD] - BC_char[i]->fion_t_old;
		if(BC_char[i]->fion_type == ZWT) { // écoulement surverse
		if((hamont > EPS_TS) && (hamont >= 1.5*haval)) // écoulement denoye
		   qd += (BC_char[i]->hydwork_param[MU]) * BC_char[i]->hydwork_param[WIDTH] * sqrt(2 * GR) * pow(hamont,1.5);
		if((hamont > EPS_TS) && haval < 1.5*haval)
		   qn += (1.5*sqrt(3.0) * BC_char[i]->hydwork_param[MU]) * BC_char[i]->hydwork_param[WIDTH] * sqrt(2 * GR) * pow(hamont-haval,0.5);
	   }
			
	}
	if(qd + qn > 0.)
	    pspecies->plink->preoxy_ttc->frac_debit_oxy[id_abs_ele] = qd/(qd + qn);	
	else
	{
		pspecies->plink->preoxy_ttc->frac_debit_oxy[id_abs_ele] = 0.;
	}
	if((qd + qn) > qtot)
	 pspecies->plink->preoxy_ttc->debit_surverse[id_abs_ele] = qtot; //debit_surverse : debit denoye
	else
	pspecies->plink->preoxy_ttc->debit_surverse[id_abs_ele] = qd + qn;
			
}


void Prose_calc_rd_rn(s_element_hyd *pele,double tempe, s_species_ttc *pspecies,int id_abs_ele, double Osat,FILE *fp)
{
    double chute;
	double rd,rn;
	double thetaT = 1.02;
	double froude,reynold,froude_jet;
	double qd;
	s_hydro_hyd *phyd, *phydn;
    s_element_hyd *pelen;//élément aval
	s_simul *psimulbio;
	
	phyd = pele->center->hydro;
	pelen = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO];
	phydn = pelen->center->hydro;
	
	//qtot = pele->face[X_HYD][TWO]->hydro->Q[T_HYD];
	qd = pspecies->plink->preoxy_ttc->debit_surverse[id_abs_ele];
	/***SW 04/12/2019 calculate Osate only for grid id_abs_ele if libseb calculated***/
	if((Simul->calc_mode[H_T] == YES_TS) && (Simul->calc_mode[SEB] == YES_TS))
	{
		psimulbio = Simul->psimul_bio[0][id_abs_ele]; 
		tempe = psimulbio->section->meteo->tempe_value;
		calc_O2_sat(tempe,psimulbio);
		Osat = psimulbio->section->compartments[WATER][0]->pspecies[O2][0]->dissolved->gas->Csat;
	}	
	chute = phyd->Z[T_HYD] - phydn->Z[T_HYD];
	if(pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->formule_rea ==  HOLLER)
	{
	rd = 1.0 + pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->holler*chute;
	rd = exp(log(rd)*pow(thetaT,tempe-20.));
	pspecies->plink->preoxy_ttc->rd[id_abs_ele] = rd;
	}
	else  // prose369 il faut definir holler
	{
		if(qd > 0.) 
		{
			froude_jet = 2.0*GR*pow(chute,3.0)*phyd->Width*phyd->Width;
			froude_jet /= (qd*qd);
			froude_jet = pow(froude_jet,0.25);
			rd=exp(pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->holler*pow(froude_jet,0.589)*pow(chute/24.,0.773));
			rd = exp(log(rd)*pow(thetaT,tempe-15.));
			pspecies->plink->preoxy_ttc->rd[id_abs_ele] = rd;
		}
		
	}
	pspecies->plink->preoxy_ttc->osat = Osat; // C'est mieux de ne pas mettre ici
	if(pspecies->plink->preoxy_ttc->frac_debit_oxy[id_abs_ele] < 1.0)
	{
		froude = pele->face[X_HYD][TWO]->hydro->Vel * pele->face[X_HYD][TWO]->hydro->Vel/(GR*pelen->face[X_HYD][ONE]->hydro->Surf/pelen->face[X_HYD][ONE]->hydro->Width);
		froude = sqrt(froude);
		froude = pow(froude,0.21);
		reynold = pelen->face[X_HYD][ONE]->hydro->Q[T_HYD]/pelen->face[X_HYD][ONE]->hydro->Width/1.143/1000000.;
        reynold = pow(reynold,0.75);
		rn = 1.0 + 1.0043 / 1000000. * froude *  reynold;
		rn = exp(log(rn)*pow(thetaT,tempe-15));
		pspecies->plink->preoxy_ttc->rn[id_abs_ele] = rn;
	 }
}
