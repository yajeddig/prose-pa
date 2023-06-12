/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: manage_link_prose.c
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
#ifdef OMP
#include <omp.h>
#endif
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
#include "global_PROSE.h"
#include "ext_PROSE.h"
/*
void PROSE_fill_var(s_species_ttc **p_species, int nspecies_ttc, int nele_ttc)
{
  //number of species for transport
  int nspecies, ns = 0,ne;
  int nele,e,nsub;
  
  char *name_ttc, *name_rive;
  nspecies = nspecies_ttc;
  nele = nele_ttc;
  
  for(ns = 0; ns < nspecies; ns ++)
  {
	  name_ttc = p_species[ns]->name;	  
	  for(e = 0; e < NSPECIES; e++)
	  {
		for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++) 
		{
			name_rive = Simul->psimul_bio[0]->section->compartments[WATER][0]->pspecies[e][nsub]->name;
			if(strcmp(name_ttc,name_rive) == 0)
			{
				for(ne = 0; ne < nele; ne++)
				{
				   p_species[ns]->plink->pvar_ttc->var[ne] = Simul->psimul_bio[ne]->section->compartments[WATER][0]->pspecies[e][nsub]->C; 
				   LP_printf(Simul->poutputs,"name = %s id_ele = %d, C = %f diff_mol = %f\n",name_rive,ne,p_species[ns]->plink->pvar_ttc->var[ne],p_species[ns]->plink->pbase_ttc->psolute_ttc->diff_mol[ne]);

				}
			}
		}			
	  }
  }
}*/

/*function used to fill the concentrations of one species for all elements*/
/*linked with C-RIVE pspecies->C at actual time step t*/
void PROSE_fill_var_one_species(s_species_ttc *pspecies, int nele, int np,FILE *fp)
{
	int ne, e, nsub,nthreads;
	char *name_ttc, *name_rive;
	//double volume;
    nthreads = Simul->psmp->nthreads;
	name_ttc = pspecies->name;
	for(e = 0; e < NSPECIES; e++)
	{
	   for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++) 
	   {
		    name_rive = Simul->psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[e][nsub]->name;
			if(strcmp(name_ttc,name_rive) == 0) // find the same species
			{
	            //#ifdef OMP
	            //omp_set_num_threads(nthreads);
	            //#pragma omp parallel for 
	            //#endif 
				for(ne = 0; ne < nele; ne++)
				{
				   pspecies->plink->pvar_ttc->var[ne] = Simul->psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[e][nsub]->C; // only one sublayer 
                    //volume = Simul->psimul_bio[ne]->section->compartments[WATER][0]->state->volume;
					//pspecies->plink->pvar_ttc->var[ne] = Simul->psimul_bio[ne]->section->compartments[WATER][0]->pspecies[e][nsub]->C*volume;
				   //LP_printf(fp,"name = %s id_ele = %d, C = %f diff_mol = %f\n",name_rive,ne,p_species[ns]->plink->pvar_ttc->var[ne],p_species[ns]->plink->pbase_ttc->psolute_ttc->diff_mol[ne]);
				}
			}
		}			
	  }	
}

void PROSE_fill_var_one_annex_species(s_species_ttc *pspecies, int nele, int nsub, int phy, int np, FILE *fp)
{
	int ne,nthreads;
	//double volume;
	nthreads = Simul->psmp->nthreads;
	//#ifndef CDA // SW 21/11/201/
	//#ifdef OMP
	//omp_set_num_threads(nthreads);
	//#pragma omp parallel for 
	//#endif
    //#endif 	
	for(ne = 0; ne < nele; ne++)
	{
		//volume = Simul->psimul_bio[ne]->section->compartments[WATER][0]->state->volume;
		pspecies->plink->pvar_ttc->var[ne] = Simul->psimul_bio[np][ne]->section->compartments[WATER][0]->pannex_var[phy][nsub]->C; //*volume;
	}
}

void PROSE_fill_var_all_annex_species(s_species_ttc ***p_phy_species, int nele, int np, FILE *fp)
{
	int phy, nsub;
	s_species_ttc *pspecies;
	
	for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
	{
		for(phy = 0; phy < 3; phy++)
		{
			pspecies = p_phy_species[nsub][phy];
			PROSE_fill_var_one_annex_species(pspecies, nele, nsub, phy, np,fp);
		}
	}
}
/*function used to fill the concentrations of all species*/
void PROSE_fill_var_all_species(s_species_ttc **p_species, int nspecies, int nele, int np,FILE *fp)
{
	int ns;
    s_species_ttc *pspecies;
	
   for(ns = 0; ns < nspecies; ns ++)
   {	
      pspecies = p_species[ns];
	  PROSE_fill_var_one_species(pspecies, nele, np,fp);
   }
}

void PROSE_update_conc_bio_one_species(s_species_ttc *pspecies, int nele,int nthreads, int np,FILE *fp)
{
	int ne, e, nsub;
	char *name_ttc, *name_rive;
	//double volume;
	name_ttc = pspecies->name;
	for(e = 0; e < NSPECIES; e++)
	{
		if(e != PHY){
	   for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++) 
	   {
		    name_rive = Simul->psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[e][nsub]->name;
			if(strcmp(name_ttc,name_rive) == 0) // find the same species
			{
	            //#ifndef CDA // SW 21/11/2018
				//#ifdef OMP
	            //omp_set_num_threads(nthreads);
	            //#pragma omp parallel for 
	            //#endif
                //#endif                
     			for(ne = 0; ne < nele; ne++)
				{
				   //volume = Simul->psimul_bio[ne]->section->compartments[WATER][0]->state->volume;
				   Simul->psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[e][nsub]->C = pspecies->plink->pvar_ttc->var[ne];//volume;  // only one sublayer 
				   Simul->psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[e][nsub]->newC = pspecies->plink->pvar_ttc->var[ne];//volume;  // only one sublayer
				   //LP_printf(fp,"name = %s id_ele = %d, C = %f diff_mol = %f\n",name_rive,ne,p_species[ns]->plink->pvar_ttc->var[ne],p_species[ns]->plink->pbase_ttc->psolute_ttc->diff_mol[ne]);
				}
			}
		}			
	  }	
	}
}
/* function rewrite in PROSE_update_phy_ctot();
void PROSE_update_conc_bio_one_annex_species(s_species_ttc *pspecies, int nele, int nsub, int phy, int np, FILE *fp)
{
	int ne,nthreads;
    //double volume;
    nthreads = Simul->psmp->nthreads;
	#ifdef OMP
	omp_set_num_threads(nthreads);
	#pragma omp parallel for 
	#endif 
	for(ne = 0; ne < nele; ne++)
   {
	  //volume = Simul->psimul_bio[ne]->section->compartments[WATER][0]->state->volume;
	  Simul->psimul_bio[np][ne]->section->compartments[WATER][0]->pannex_var[phy][nsub]->C = pspecies->plink->pvar_ttc->var[ne];//volume;  // only one sublayer 				   //LP_printf(fp,"name = %s id_ele = %d, C = %f diff_mol = %f\n",name_rive,ne,p_species[ns]->plink->pvar_ttc->var[ne],p_species[ns]->plink->pbase_ttc->psolute_ttc->diff_mol[ne]);
   }   

}

void PROSE_update_conc_bio_all_annex_species(s_species_ttc ***p_phy_species, int nele, int np, FILE *fp)
{
	int phy, nsub;
	s_species_ttc *pspecies;
	
	for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
	{
		for(phy = 0; phy < 3; phy++)
		{
			pspecies = p_phy_species[nsub][phy];
			PROSE_update_conc_bio_one_annex_species(pspecies, nele, nsub, phy, np fp);
		}
	}
}*/
/*function used to update the concentrations of all species after transport*/
void PROSE_update_conc_bio_all_species(s_species_ttc **p_species, int nspecies, int nele, int nthreads, int np,FILE *fp)
{
	int ns;
    s_species_ttc *pspecies;
	
   for(ns = 0; ns < nspecies; ns ++)
   {	
      pspecies = p_species[ns];
	  PROSE_update_conc_bio_one_species(pspecies, nele, nthreads, np,fp);
   }
}

/*function  set the porosity of water to 1*/
 
void PROSE_fill_param_base_one_species(s_species_ttc *pspecies, int nele, int np, FILE *fp)
{
	int e, nsub,ne,nthreads;
	char *name_ttc, *name_rive;
    nthreads = Simul->psmp->nthreads;	
    name_ttc = pspecies->name;
	for(e = 0; e < NSPECIES; e++)
	{
		for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++)
		{
                    if(Simul->psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[e] != NULL) // SW 22/03/2023 for TA pH DIC not defined
                    {
                                          name_rive = Simul->psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[e][nsub]->name;
		       if(strcmp(name_ttc,name_rive) == 0)
			{
	            //#ifndef CDA // SW 21/11/2018
				//#ifdef OMP
	            //omp_set_num_threads(nthreads);
	            //#pragma omp parallel for 
	            //#endif
                //#endif				
				for(ne = 0; ne < nele; ne++)
				{
					pspecies->plink->pbase_ttc->param_syst[POROSITY_TTC][ne] = 1.;
					pspecies->plink->pbase_ttc->param_syst[DISPERSIVITY_TTC][ne] = 0.; // SW 20/04/2018 a verifier
					pspecies->plink->pbase_ttc->psolute_ttc->diff_mol[ne] = Simul->psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[e][nsub]->diff_mol_ttc;
				}
			 }
                    }
		}
	}	
}

void PROSE_fill_param_base_one_annex_species(s_species_ttc *pspecies, int nele, int nsub, int np, FILE *fp)
{
	int ne,nthreads;
	nthreads = Simul->psmp->nthreads;
	//#ifndef CDA // SW 21/11/2018
	//#ifdef OMP
	//omp_set_num_threads(nthreads);
	//#pragma omp parallel for 
	//#endif
    //#endif	
	for(ne = 0; ne < nele; ne++)
	{
		pspecies->plink->pbase_ttc->param_syst[POROSITY_TTC][ne] = 1.;
		pspecies->plink->pbase_ttc->param_syst[DISPERSIVITY_TTC][ne] = 0.; // SW 20/04/2018 a verifier
		pspecies->plink->pbase_ttc->psolute_ttc->diff_mol[ne] = Simul->psimul_bio[np][ne]->section->compartments[WATER][0]->pspecies[PHY][nsub]->diff_mol_ttc;
	}
}

void PROSE_fill_param_base_all_annex_species(s_species_ttc ***p_phy_species, int nele, int np,FILE *fp)
{
	int phy,nsub;
	s_species_ttc *pspecies;
	
	for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
	{
		for(phy = 0; phy < 3; phy++)
		{
			pspecies = p_phy_species[nsub][phy];
			PROSE_fill_param_base_one_annex_species(pspecies, nele, nsub, np, fp);
		}
	}
}
void PROSE_fill_param_base_all_species(s_species_ttc **p_species, int nspecies, int nele, int np, FILE *fp)
{
	int ns;
	s_species_ttc *pspecies;
	
	for(ns = 0; ns < nspecies; ns++)
	{
		pspecies = p_species[ns];
		PROSE_fill_param_base_one_species(pspecies, nele,np,fp);
	}
}

/*
void PROSE_fill_param_base(s_species_ttc **p_species, int nspecies, int nele)
{
	int ns, ne,e,nsub;
	char *name_ttc, *name_rive;
	for(ns = 0; ns < nspecies; ns++)
	{
		name_ttc = p_species[ns]->name;
	    for(e = 0; e < NSPECIES; e++)
	    {
		   for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++) 
		   {
			name_rive = Simul->psimul_bio[0]->section->compartments[WATER][0]->pspecies[e][nsub]->name;
			if(strcmp(name_ttc,name_rive) == 0)
			{
		       for(ne = 0; ne < nele; ne++)
		       {
			    p_species[ns]->plink->pbase_ttc->param_syst[POROSITY_TTC][ne] = 1.;
		        p_species[ns]->plink->pbase_ttc->param_syst[DISPERSIVITY_TTC][ne] = 0.; // SW 20/04/2018 a verifier			
			    p_species[ns]->plink->pbase_ttc->psolute_ttc->diff_mol[ne] = Simul->psimul_bio[ne]->section->compartments[WATER][0]->pspecies[e][nsub]->diff_mol_ttc;
		       }
			}
		   }
		}		
	}
}*/


void PROSE_fill_u_one_species(s_species_ttc *pspecies, s_chyd *pchyd, double tempe, double Osat,FILE *fp)
{
	int id_abs_ele;
	int ne, r, icard,sub_icard;
	char *name_rive,*name_ttc;
	s_reach_hyd *preach;
	s_element_hyd *pele;
    s_face_hyd *pface;
	
	name_rive = Simul->psimul_bio[0][0]->section->compartments[WATER][0]->pspecies[O2][0]->name; // SW 23/10/2018 oxygene name
	name_ttc = pspecies->name;
    for(r = 0; r < pchyd->counter->nreaches; r++)
    {
	  preach = pchyd->p_reach[r];
//#ifdef OMP //SW 06/09/2018 ajout block openmp
//      int taille;
//      int nthreads;
//      nthreads=Simul->psmp->nthreads;
//      omp_set_num_threads(nthreads);

//      Simul->psmp->chunk=PC_set_chunk_size_silent(fp,preach->nele-1,nthreads);
//      taille=Simul->psmp->chunk;
//#pragma omp parallel shared(preach,nthreads,taille,pspecies) private(ne,pele,pface,icard,sub_icard,id_abs_ele)
 // {
//#pragma omp for schedule(dynamic,taille) 
//#endif		  
	  for(ne = 0; ne < preach->nele; ne++)
	  {
		pele = preach->p_ele[ne];
	    id_abs_ele = pele->id[ABS_HYD];
		//if(id_abs_ele == 150)
			//printf("id == 150\n");
		pface = pele->face[X_HYD][ONE];
		for(icard = 0; icard < NORTH_TTC; icard++)
		{
			for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++ ){

	    // SW 23/04/2018 EAST_TTC correspond to TWO_HYD, WEST_TTC correspont ONE_HYD pour v >0 aval
		if(icard == EAST_TTC) {
			
           *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][0] = pele->face[X_HYD][TWO]->hydro->Q[T_HYD];//pele->face[X_HYD][TWO]->hydro->Vel; //  5.; //SW 03/05/2018 a verifier vitesse de l'element centre ou au face ?
		   *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][1] = 0; // SW 29/05/2018 non active
           sub_icard = SUB_CARD_TTC;
		   //LP_printf(fp,"id = %d vam = %f vav = %f\n",id_abs_ele,pele->face[X_HYD][ONE]->hydro->Vel,pele->face[X_HYD][TWO]->hydro->Vel);
		   }
        else if(icard == WEST_TTC){
		  if(ne == 0 && HYD_test_sing_face(pface,UPSTREAM) == CONFLUENCE)
		  {
			  *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][sub_icard] = pface->limits[UPSTREAM]->faces[ONE][sub_icard]->hydro->Q[T_HYD];//pface->limits[UPSTREAM]->faces[ONE][sub_icard]->hydro->Vel;// 5.; //
		  }
		  else if((ne == 0) && (HYD_test_sing_face(pface,UPSTREAM) == HYDWORK)) // SW 23/10/2018 element aval du bief et type barrage pour reoxygeneration
			{
				if(strcmp(name_ttc,name_rive) == 0) // oxygene species
				{
				*pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][0] = pele->face[X_HYD][ONE]->hydro->Q[T_HYD];
				// la fraction du debit en surverse sauvegarde a 1.
				*pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][1] = 0.0;
				Prose_cherche_debit_oxy(pface->limits[UPSTREAM]->faces[ONE][0]->element[ONE], pface->BC_char,pface->nworks,pspecies,id_abs_ele,fp);
                
				sub_icard = SUB_CARD_TTC;
				Prose_calc_rd_rn(pface->limits[UPSTREAM]->faces[ONE][0]->element[ONE],tempe,pspecies,id_abs_ele,Osat,fp);
				}
				else{
				*pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][0] = pele->face[X_HYD][ONE]->hydro->Q[T_HYD];//pele->face[X_HYD][TWO]->hydro->Vel; //  5.; //SW 03/05/2018 a verifier vitesse de l'element centre ou au face ?
		        *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][1] = 0; // SW 29/05/2018 non active
                sub_icard = SUB_CARD_TTC;	
				}
			}
		  else{
           *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][0] = pele->face[X_HYD][ONE]->hydro->Q[T_HYD];//pele->face[X_HYD][ONE]->hydro->Vel; //   5.; //
		   *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][1] = 0; // SW 29/05/2018 non active
           sub_icard = SUB_CARD_TTC;		   
		//LP_printf(fp,"id_ele = %d, icard = %d, v = %f\n",id_abs_ele,icard,*pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard]);

		  }
		  }
	  }
		}
		}
// #ifdef OMP
//} /* end of parallel section */
//#endif
	}
}

void PROSE_fill_u_all_species(s_species_ttc **p_species, int nspecies, s_chyd *pchyd, double tempe, double Osat, FILE *fp)
{
	int ns;
    s_species_ttc *pspecies;
	
	for(ns = 0; ns < nspecies; ns++)
	{
		pspecies = p_species[ns];
		PROSE_fill_u_one_species(pspecies, pchyd, tempe,Osat,fp);
	}
      
}

void PROSE_fill_u_all_annex_species(s_species_ttc ***p_phy_species, s_chyd *pchyd,double tempe,double Osat, FILE *fp)
{
	int phy,nsub;
	s_species_ttc *pspecies;
	
	for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
	{
		for(phy = 0; phy < 3; phy++)
		{
			pspecies = p_phy_species[nsub][phy];
			PROSE_fill_u_one_species(pspecies, pchyd,tempe,Osat, fp);
		}
	}
}


/*void PROSE_fill_u(s_species_ttc **p_species, int nspecies_ttc, int nele_ttc, s_chyd *pchyd, FILE *fp)
{
  //number of species for transport
  int nspecies, ns = 0,ne;
  int nele,icard,icard_hyd;
  int r,id_abs_ele;

  s_reach_hyd *preach;
  s_element_hyd *pele;
  nspecies = nspecies_ttc;
  nele = nele_ttc;

  for(r = 0; r < pchyd->counter->nreaches; r++)
  {
	  preach = pchyd->p_reach[r];
	  for(ne = 0; ne < preach->nele; ne++)
	  {
		  pele = preach->p_ele[ne];
		  id_abs_ele = pele->id[ABS_HYD];
		  for(ns = 0; ns < nspecies; ns++)
		  {
			 icard_hyd = TWO;
		     for(icard = 0; icard < NB_CARD_TTC; icard++)
		     {
				 if(icard < NORTH_TTC)
				 {					                     
					*p_species[ns]->plink->pu_ttc->puface[id_abs_ele].uface[icard] = pele->face[X_HYD][icard_hyd--]->hydro->Vel; // SW 23/04/2018 WEST_TTC correspont ONE_HYD, EAST_TTC correspont TWO_HYD pour v >0
					LP_printf(fp,"id_ele = %d, icard = %d, v = %f\n",id_abs_ele,icard,*p_species[ns]->plink->pu_ttc->puface[id_abs_ele].uface[icard]);
				 }
	
		     }
	  }
	  
  }
  }
}*/

void PROSE_alloc_u_dist_one_species(s_species_ttc *pspecies, int nele, s_chyd *pchyd, FILE *fp)
{
  int ne, r, icard, id_abs_ele, sub_icard;
  s_reach_hyd *preach;
  s_element_hyd *pele;	
  //#ifndef CDA
  //#ifdef OMP //SW 06/09/2018 ajout block openmp
  //      int taille;
  //      int nthreads;
  //     nthreads=Simul->psmp->nthreads;
  //      omp_set_num_threads(nthreads);
  //
  //     Simul->psmp->chunk=PC_set_chunk_size_silent(fp,pchyd->counter->nreaches-1,nthreads);
  //      taille=Simul->psmp->chunk;
  //#pragma omp parallel shared(pchyd,nthreads,taille,pspecies) private(preach,ne,pele,icard,sub_icard,id_abs_ele)
  //  {
  //#pragma omp for schedule(dynamic,taille) 
  //#endif
  //#endif	
  for(r = 0; r < pchyd->counter->nreaches; r++)
    {
      preach = pchyd->p_reach[r];
      for(ne = 0; ne < preach->nele; ne++)
	{
	  pele = preach->p_ele[ne];
	  id_abs_ele = pele->id[ABS_HYD];
          pspecies->plink->pbase_ttc->param_syst[SIZE_TTC][id_abs_ele] = pele->length;
	  for(icard = 0; icard < NB_CARD_TTC; icard++)
	    {
	      for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++){
		pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][sub_icard] = (double *)malloc(sizeof(double));
		if(icard > WEST_TTC)
		  *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][sub_icard] = 0.;
	      }
	    }
	}
      
    }
  //#ifndef CDA 
  // #ifdef OMP
  //} /* end of parallel section */
  //#endif
  //#endif
}

void PROSE_alloc_u_dist_all_species(s_species_ttc **p_species, int nspecies, int nele, s_chyd *pchyd, FILE *fp)
{
	int ns;
	s_species_ttc *pspecies;
	
	for(ns = 0; ns < nspecies; ns++)
	{
		pspecies = p_species[ns];
		PROSE_alloc_u_dist_one_species(pspecies, nele, pchyd, fp);
	}
}

void PROSE_alloc_u_dist_all_annex_species(s_species_ttc ***p_phy_species, int nele, s_chyd *pchyd, FILE *fp)
{
	int phy,nsub;
	s_species_ttc *pspecies;
	
	for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
	{
		for(phy = 0; phy < 3; phy++)
		{
			pspecies = p_phy_species[nsub][phy];
		    PROSE_alloc_u_dist_one_species(pspecies, nele, pchyd, fp);	
		}
	}		
}
/*
void PROSE_alloc_u_dist(s_species_ttc **p_species, int nspecies_ttc, int nele_ttc, s_chyd *pchy, FILE *fp)
{
  //number of species for transport
  int nspecies, ns = 0,ne;
  int nele,icard;
  int r,id_abs_ele;
  s_reach_hyd *preach;
  s_element_hyd *pele;
  nspecies = nspecies_ttc;
  nele = nele_ttc;
  
  for(r = 0; r < pchy->counter->nreaches; r++)
  {
	  preach = pchy->p_reach[r];
	  for(ne = 0; ne < preach->nele; ne++)
	  {
		  pele = preach->p_ele[ne];
		  id_abs_ele = pele->id[ABS_HYD];
		  for(ns = 0; ns < nspecies; ns++)
		  {
			 p_species[ns]->plink->pbase_ttc->param_syst[SIZE_TTC][id_abs_ele] = pele->length;
		     for(icard = 0; icard < NB_CARD_TTC; icard++)
		     {
				 p_species[ns]->plink->pu_ttc->puface[id_abs_ele].uface[icard] = (double *)malloc(sizeof(double));
                 if(icard > WEST_TTC)
				 {
					*p_species[ns]->plink->pu_ttc->puface[id_abs_ele].uface[icard] = 0.;
				 }
			 }
		  }
	  }
  }
}*/

void PROSE_init_neighbor_one_species(s_species_ttc *pspecies, int nele, FILE *fp)
{
	int e, icard,sub_icard;
	
	pspecies->plink->pneigh_ttc = TTC_init_neigh(nele);
//#ifndef CDA	
//#ifdef OMP //SW 06/09/2018 ajout block openmp
//      int taille;
 //     int nthreads;
 //     nthreads=Simul->psmp->nthreads;
//      omp_set_num_threads(nthreads);
//
//      Simul->psmp->chunk=PC_set_chunk_size_silent(fp,nele-1,nthreads);
//      taille=Simul->psmp->chunk;
//#pragma omp parallel shared(nthreads,taille,pspecies) private(e,icard,sub_icard)
 // {
//#pragma omp for schedule(dynamic,taille) 
//#endif
//#endif
	for(e = 0; e < nele; e++)
	{
	   for(icard = 0; icard < NB_CARD_TTC; icard++)
	   {
		   for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++){
			 pspecies->plink->pneigh_ttc->pivois[e].ivois[icard][sub_icard] = NONE_TTC;
	   }
	   		pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[e].icl[icard] = NO_BOUND_TTC;

	   }
	}
//#ifndef CDA	
 //#ifdef OMP
//} /* end of parallel section */
//#endif
//#endif
}

void PROSE_fill_neighbor_all_species(s_species_ttc **p_species, int nspecies, int nele, s_chyd *pchyd, FILE *fp)
{
	int ns;
	s_species_ttc *pspecies1;
	
    for(ns = 0; ns < nspecies; ns++)
    {
     pspecies1 = p_species[ns];		
	 PROSE_init_neighbor_one_species(pspecies1, nele, fp);
	 PROSE_fill_neighbor_one_species(pspecies1, pchyd, fp);
	}
}

void PROSE_fill_neighbor_all_annex_species(s_species_ttc ***p_phy_species, int nele, s_chyd *pchyd, FILE *fp)
{
	int phy,nsub;
	s_species_ttc *pspecies;
	
	for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
	{
		for(phy = 0; phy < 3; phy++)
		{
			pspecies = p_phy_species[nsub][phy];
			PROSE_init_neighbor_one_species(pspecies, nele, fp);
	        PROSE_fill_neighbor_one_species(pspecies, pchyd, fp);
		}
	}
}
/*void PROSE_fill_neighbor_one_species(s_species_ttc *pspecies, s_chyd *pchyd, FILE *fp)
{
 
  int ne = 0;
  int r,id_abs_ele;
  
  s_reach_hyd *preach;
  s_element_hyd *pele;
   
	 for(r = 0; r < pchyd->counter->nreaches; r++)
	 {
		 preach = pchyd->p_reach[r];
		 for(ne = 0; ne < preach->nele; ne++)
		 {
			 pele = preach->p_ele[ne];
			 id_abs_ele = pele->id[ABS_HYD];
			if(id_abs_ele == 150)
			  printf("id == 150\n");
	         pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[NORTH_TTC] = NEU_FACE_TTC;
	         pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[NORTH_TTC] = 0.0;
	         pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
	         pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[SOUTH_TTC] = NEU_FACE_TTC;
	         pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[SOUTH_TTC] = 0.0;
	         pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
			 if(ne == 0 && HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM) == DISCHARGE) //  amont DIRI_CENTER_TTC
		     {
			    pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
			    pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[WEST_TTC] = DIRI_FACE_TTC;//DIRI_CENTER_TTC;
				//pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
				//pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[NORTH_TTC] = DIRI_CENTER_TTC;
				//pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
				//pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[SOUTH_TTC] = DIRI_CENTER_TTC;
				
                //break; // SW 15/05/2018 pour sortir la boucle ne nele 
			}
			if(ne == 0 && HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM) == CONFLUENCE)
			{
				pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];
				pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
				pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][1] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][1]->element[ONE]->id[ABS_HYD];
				pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
                pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->length + pele->length)/2;
                pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][1] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][1]->element[ONE]->length + pele->length)/2;
			}
			else if(ne == 0 && HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM) == DIFFLUENCE)
			{
		     pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];	      
	         pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
		     pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->length + pele->length)/2; // distance entre le centre de la maille et son voinsin.				
			}
			else if(ne == 0 && pele->face[X_HYD][ONE]->limits[UPSTREAM] != NULL && HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM) != DISCHARGE)
			{
				if(HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM) == CONF_DIFF)
					LP_printf(fp," il existe conf_diff\n");
		     pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];	      
	         pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
		     pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->length + pele->length)/2; // distance entre le centre de la maille et son voinsin.				
				
			}
	         else if(ne != 0 && pele->face[X_HYD][ONE]->element[ONE] != NULL)
	         {
		     pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->element[ONE]->id[ABS_HYD];	      
	         pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
		     pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->element[ONE]->length + pele->length)/2; // distance entre le centre de la maille et son voinsin.
	         }
			//else if(HYD_test_sing_face(pele->face[X_HYD][TWO],DOWNSTREAM) == DIFFLUENCE){
             // pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->element[TWO]->id[ABS_HYD];				
			//}
			if((ne == preach->nele -1) && (HYD_test_sing_face(pele->face[X_HYD][TWO],DOWNSTREAM) == DIFFLUENCE)){
				pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->id[ABS_HYD];
				pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
				pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][1] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][1]->element[TWO]->id[ABS_HYD];
				pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
                pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->length + pele->length)/2;
                pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][1] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][1]->element[TWO]->length + pele->length)/2;				
			}
           else if((ne == preach->nele -1) && (HYD_test_sing_face(pele->face[X_HYD][TWO],DOWNSTREAM) == CONFLUENCE)){
		     pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->id[ABS_HYD];	      
	         pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
		     pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->length + pele->length)/2;
			   
		   }	
			else if((ne != preach->nele -1) && (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]) != NULL)
			{
		     pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->id[ABS_HYD];	      
	         pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
		     pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->length + pele->length)/2;
				
			}		   
	        else if(pele->face[X_HYD][TWO]->element[TWO] != NULL)
	        {
		     pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->element[TWO]->id[ABS_HYD];	      
	         pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
		     pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->element[TWO]->length + pele->length)/2;
		     //LP_printf(fp,"id_abs_ele = %d, nvois = %d, id_vois = %d\n",id_abs_ele,p_species[ns]->plink->pneigh_ttc->nvois[id_abs_ele],p_species[ns]->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC]);
             //if(ne == 0 && HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM) == DISCHARGE)
			 //{
			//	pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[EAST_TTC] = DIRI_CENTER_TTC;
			//	pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;				 
			 //}
			}
			else // element aval
			{
			  pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
	          pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[EAST_TTC] = DIRI_FACE_TTC;	
			}

		 }
	 }	  
}*/

void PROSE_fill_neighbor_one_species(s_species_ttc *pspecies, s_chyd *pchyd, FILE *fp) // + initialize BC upstream and downstream for HT (AB 08.10.2019)
{
  /*number of species for transport*/
  int ne = 0;
  int r,id_abs_ele;
  char *temp_name;

  s_reach_hyd *preach;
  s_element_hyd *pele;
  s_face_hyd *pface;
  s_face_hyd *pfacep;
  //#ifndef CDA  
  //#ifdef OMP //SW 06/09/2018 ajout block openmp
  //        int taille;
  //        int nthreads;
  //        nthreads=Simul->psmp->nthreads;
  //       omp_set_num_threads(nthreads);
  
  //        Simul->psmp->chunk=PC_set_chunk_size_silent(fp,pchyd->counter->nreaches-1,nthreads);
  //       taille=Simul->psmp->chunk;
  //#pragma omp parallel shared(nthreads,taille,pspecies,pchyd) private(r,preach,ne,pele,pface,pfacep,id_abs_ele)
  // {
  //#pragma omp for schedule(dynamic,taille) 
  //#endif
  //#endif
  
  /*** Allows to make a verification on the specie's type ***/ 
  temp_name = (char *)malloc(109*sizeof(char));
  strcpy(temp_name, "Temp"); /*******************************/
  
  for(r = 0; r < pchyd->counter->nreaches; r++)
    {
      preach = pchyd->p_reach[r];
      pele = preach->p_ele[0]; // premier element du bief
      pface = pele->face[X_HYD][ONE];
      id_abs_ele = pele->id[ABS_HYD];
      
      if(id_abs_ele == 11)
        printf("id == 11\n");
      pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[NORTH_TTC] = NEU_FACE_TTC;
      pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[NORTH_TTC] = 0.0;
      pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
      pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[SOUTH_TTC] = NEU_FACE_TTC;
      pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[SOUTH_TTC] = 0.0;
      pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
      switch(HYD_test_sing_face(pface,UPSTREAM))
	{
	case DISCHARGE: {
	  pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
	  pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[WEST_TTC] = DIRI_FACE_TTC;//DIRI_CENTER_TTC; NEU_FACE_TTC;//
	  //if(strcmp(temp_name, pspecies->name) == 0) // AB 08.10.2019
	  //  pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[WEST_TTC] = T_0 + 20.; // To be implemented in the command file !!!!!!!!!!!!!!!!!!!!!!!
	  break;
	}
	case CONFLUENCE: {
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][1] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][1]->element[ONE]->id[ABS_HYD];
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->length + pele->length)/2;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][1] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][1]->element[ONE]->length + pele->length)/2;				 
	  break;
	}
	case CONF_DIFF: {
	  pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
	  pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[WEST_TTC] = CONF_DIFF_TTC;//DIRI_CENTER_TTC;
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][1] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][1]->element[ONE]->id[ABS_HYD];
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->length + pele->length)/2;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][1] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][1]->element[ONE]->length + pele->length)/2;				 
	  break; 
	}
	case DIFFLUENCE:{
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];	      
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->length + pele->length)/2; // distance entre le centre de la maille et son voinsin.								 
	  break;
	}
	  /*SW 24/10/2018 add HYDWORK*/
	case HYDWORK : {
	  if(pspecies->oxygen == YES_TS || strcmp(temp_name, pspecies->name) == 0) // AB 08.10.2019 // SW 04/05/2021 pourquoi barrage_TTC pour temperature ?
	    {
	      pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
	      pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[WEST_TTC] = BARRAGE_TTC;
	      pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];
	      pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	      pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->length + pele->length)/2;
	    }
	  else{
	    pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];	      
	    pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	    pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->length + pele->length)/2; // distance entre le centre de la maille et son voinsin.								 					
	  }
	  break;
	}				
	  
	default:{
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];	      
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->length + pele->length)/2; // distance entre le centre de la maille et son voinsin.								 
	}
	}
      
      if(preach->nele > 1)
	{
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->element[TWO]->id[ABS_HYD];	      
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->element[TWO]->length + pele->length)/2;			 
	  //LP_printf(fp,"id_abs_ele = %d id_neigh = %d\n",id_abs_ele,pele->face[X_HYD][TWO]->element[TWO]->id[ABS_HYD]);
	}
      
      //element au sein du bief
      for(ne = 1; ne < preach->nele -1; ne++)
	{	
	  pele = preach->p_ele[ne];
	  id_abs_ele = pele->id[ABS_HYD];
          
          //if(id_abs_ele == 2)
             //printf("debug\n");
	  // NORTH and SOUTH faces
	  pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[NORTH_TTC] = NEU_FACE_TTC;
	  pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[NORTH_TTC] = 0.0;
	  pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
	  pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[SOUTH_TTC] = NEU_FACE_TTC;
	  pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[SOUTH_TTC] = 0.0;
	  pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
	  
	  // WEST face
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->element[ONE]->id[ABS_HYD];	      
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->element[ONE]->length + pele->length)/2; // distance entre le centre de la maille et son voinsin.
	  //LP_printf(fp,"id_abs_ele = %d id_neigh = %d\n",id_abs_ele,pele->face[X_HYD][ONE]->element[ONE]->id[ABS_HYD]);
	  // EAST face			 
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->element[TWO]->id[ABS_HYD];	      
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->element[TWO]->length + pele->length)/2;
	  //LP_printf(fp,"id_abs_ele = %d id_neigh = %d\n",id_abs_ele,pele->face[X_HYD][TWO]->element[TWO]->id[ABS_HYD]);
	}
      
      pele = preach->p_ele[preach->nele-1]; // dernier element du bief
      pfacep = pele->face[X_HYD][TWO];
      id_abs_ele = pele->id[ABS_HYD];
      //if (id_abs_ele == 3774)
      //printf("DEBUB\n");
      pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[NORTH_TTC] = NEU_FACE_TTC;
      pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[NORTH_TTC] = 0.0;
      pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
      pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[SOUTH_TTC] = NEU_FACE_TTC;
      pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[SOUTH_TTC] = 0.0;
      pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
      if(preach->nele > 1){
	pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[WEST_TTC][0] = pele->face[X_HYD][ONE]->element[ONE]->id[ABS_HYD];	      
	pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[WEST_TTC][0] = (pele->face[X_HYD][ONE]->element[ONE]->length + pele->length)/2; // distance entre le centre de la maille et son voinsin.
	//LP_printf(fp,"id_abs_ele = %d id_neigh = %d\n",id_abs_ele,pele->face[X_HYD][ONE]->element[ONE]->id[ABS_HYD]);
	
      }		 
      switch(HYD_test_sing_face(pfacep,DOWNSTREAM))
	{
	case DIFFLUENCE:{
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->id[ABS_HYD];
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][1] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][1]->element[TWO]->id[ABS_HYD];
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->length + pele->length)/2;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][1] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][1]->element[TWO]->length + pele->length)/2;				
	  break; 
	}
	case CONFLUENCE:{
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->id[ABS_HYD];	      
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->length + pele->length)/2;
	  break;
	}
	case CONF_DIFF:{
	  // a verifier avec diffusion
	  //pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
	  //pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[EAST_TTC] = DIRI_FACE_TTC;//DIRI_CENTER_TTC;
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->id[ABS_HYD];
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][1] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][1]->element[TWO]->id[ABS_HYD];
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->length + pele->length)/2;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][1] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][1]->element[TWO]->length + pele->length)/2;				
	  break; 				 
	}
	case OPEN:{ // Aval ??
	  pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
	  pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[EAST_TTC] = DIRI_FACE_TTC;	
	  LP_error(fp,"open limit for id = %d\n",id_abs_ele);
	  break;				 
	}
	case WATER_LEVEL :{
	  //if(pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0] == NULL) // Aval
	  if (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->ndownstr_reaches == 0)
	    {
	      //pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
	      //pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[EAST_TTC] = DIRI_FACE_TTC;
              /* SW 23/06/2021 why we need to define NEW_FACE_TTC for the last face ?*/
	      /*if(strcmp(temp_name, pspecies->name) == 0) // AB 08.10.2019
		{
		  printf("ALLES GUT !!!^^\n");
		  pspecies->plink->pboundary_ttc->nbound[id_abs_ele]++;
		  pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[id_abs_ele].icl[EAST_TTC] = NEU_FACE_TTC;
		  pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[EAST_TTC] = 0.;
		}*/
	    }
	  else
	    {
	      pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->id[ABS_HYD];	      
	      pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	      pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->length + pele->length)/2;
	    }
	  break;
	}
	default:{
	  pspecies->plink->pneigh_ttc->pivois[id_abs_ele].ivois[EAST_TTC][0] = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->id[ABS_HYD];	      
	  pspecies->plink->pneigh_ttc->nvois[id_abs_ele]++;
	  pspecies->plink->pneigh_ttc->pdelta_L[id_abs_ele].delta_L[EAST_TTC][0] = (pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO]->length + pele->length)/2;
	}
	}		 
      
    }
  //#ifndef CDA
  //#ifdef OMP
  //  }
  //#endif
  //#endif
}	  


void PROSE_find_apport_t_all_species(s_species_ttc **p_species, int nspecies, s_chyd *pchyd, double t,int np, FILE *fp) // MH 05/12/2021 : np added to conform with PROSE_find_apport_t_one_species equation, however this functions is used no where else
{
	int ns;
	s_species_ttc *pspecies;
	
    for(ns = 0; ns < nspecies; ns++)
	{
		pspecies = p_species[ns];
		PROSE_find_apport_t_one_species(pspecies, pchyd, t,np, fp);  
	}
	
}

/*function used to find apport inflow at elements (faces), and modify the cl var pointer*/
void PROSE_find_apport_t_one_species(s_species_ttc *pspecies, s_chyd *pchyd, double t,int np, FILE *fp)  //MH: 05/12/2021 np added as an argument for b1 
{
  int r, ne, e, nsub;
  int napp, icard;
  int id_abs_ele;
  double q_inflow, q_inflow_area, h_iter, h_i;
  double val_inflow= 0;
  double toc_inflow = 0;
  //double c_conf_diff;
  //double dt;
  //double q1,q2;
  //int id1,id2;
  //double surf,volume; // SW 17/05/2018  calculer une vitesse pour apport, on prend pele->face[X_HYD][ONE]->hydro->Surf, sont face amont de l'element
  s_reach_hyd *preach;
  s_element_hyd *pele;
	
  char *name_ttc, *name_rive;

  name_ttc = pspecies->name;
  //dt = pspecies->pchronos->dt;
  //#ifdef OMP //SW 06/09/2018 ajout block openmp
  //        int taille;
  //        int nthreads;
  //       nthreads=Simul->psmp->nthreads;
  //       omp_set_num_threads(nthreads);
	
  //        Simul->psmp->chunk=PC_set_chunk_size_silent(fp,pchyd->counter->nreaches -1,nthreads);
  //       taille=Simul->psmp->chunk;
  //#pragma omp parallel shared(nthreads,taille,pspecies,name_ttc,pchyd) private(preach,ne,e,nsub,name_rive,pele,id_abs_ele,icard,napp,val_inflow,q_inflow)
  // {
  //#pragma omp for schedule(dynamic,taille) 
  //#endif
  for (r = 0; r < pchyd->counter->nreaches; r++) 
    {
      preach = pchyd->p_reach[r];
      for(ne = 0; ne < preach->nele; ne++) 
	{
	  pele = preach->p_ele[ne];
	  //LP_printf(fp,"id = %d, long = %f\n",pele->id[ABS_HYD],pele->length);
	  id_abs_ele = pele->id[ABS_HYD];
	  for(e = 0; e < NSPECIES; e++)
	    {
	      for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++)
		{
		  name_rive = Simul->psimul_bio[0][0]->section->compartments[WATER][0]->pspecies[e][nsub]->name;
		  //if(strcmp("zoo1",name_rive) == 0)
		  //printf("yes\n");
		  q_inflow_area = 0.; // SW 30/01/2020
		  if(strcmp(name_ttc,name_rive) == 0)
		    {
		      if(ne == 0 && HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM) == DISCHARGE) // apport amont DIRI_FACE_TTC
			{
				
			  if(pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[e][nsub] != NULL)
			    {
			      // MH 18/01/2022 to calculate MOD123MOP123 of upstream inflows fromm TOC*share(b1)
			      if (e == MOD || e == MOP)
				{

				// to differentiate with rivers not having macrospecies
				  	if( (pele->face[X_HYD][ONE]->pt_inflows[0]->flow_in_macrospecies[TOC] != NULL) && (Simul->passim->param_range[B1_RIVER_DA][PARAM_UP] > INFINITE_DIV_TS )  )
				  {
				    pele->face[X_HYD][ONE]->pt_inflows[0]->flow_in_macrospecies[TOC] = TS_function_t_pointed(t, pele->face[X_HYD][ONE]->pt_inflows[0]->flow_in_macrospecies[TOC],fp);
				    toc_inflow = TS_function_value_t(t, pele->face[X_HYD][ONE]->pt_inflows[0]->flow_in_macrospecies[TOC],fp);					      
				    // MOD123 or MOP123 = TOC * corresponding share
				    val_inflow = toc_inflow * PROSE_DA_mod_mop_fract(TOC, e, nsub, Simul, np, fp);
				    //val_inflow = 3;
				    //LP_printf(fp,"val_UPSTREAM = %3.2f \n",val_inflow);
				    //LP_printf(fp,"np = %i, TOC_upstream = %3.4f ,val_inflow = %3.4f \n",np, toc_inflow/83.3, val_inflow/83.3);
				    }
				else
				  {
				    pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[e][nsub] = TS_function_t_pointed(t,pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[e][nsub],fp); // SW 03/04/2020
				    val_inflow = TS_function_value_t(t,pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[e][nsub],fp);
				  }
				pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[WEST_TTC] = val_inflow ;
			        }
			      else
				{
				      
				//surf = pele->face[X_HYD][ONE]->hydro->Surf;
				//q_inflow = pele->face[X_HYD][ONE]->hydro->Q[T_HYD];
				//pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[EAST_TTC] = TS_function_value_t(t,pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[e][nsub],fp); // SW 25/04/2018 modify valcl to s_ft
				pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[e][nsub] = TS_function_t_pointed(t,pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[e][nsub],fp); // SW 03/04/2020
				pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[WEST_TTC] = TS_function_value_t(t,pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[e][nsub],fp); // SW 25/04/2018 modify valcl to s_ft
				//LP_printf(fp,"id_ele = %d, e = %d, j = %d , val = %f, q = %f\n",pele->id[ABS_HYD],e,nsub,pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[WEST_TTC],pele->face[X_HYD][ONE]->hydro->Q[T_HYD]);
				//pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[SOUTH_TTC] = TS_function_value_t(t,pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[e][nsub],fp); // SW 25/04/2018 modify valcl to s_ft
				//pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[NORTH_TTC] = TS_function_value_t(t,pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[e][nsub],fp); // SW 25/04/2018 modify valcl to s_ft

				}		       
			    }
			  else
			    {
			      //pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[EAST_TTC] = 0.0; 
			      pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[WEST_TTC] = 0.0; 
			      //pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[SOUTH_TTC] = 0.0; 
			      //pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[NORTH_TTC] = 0.0; 
			    }
			}
		      /*if(ne == 0 && HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM)== CONF_DIFF) // apport amont DIRI_FACE_TTC
			{
			      
			q1 = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->face[X_HYD][TWO]->hydro->Q[T_HYD];
			id1 = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];
			q2 = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][1]->element[ONE]->face[X_HYD][TWO]->hydro->Q[T_HYD];
			id2 = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][1]->element[ONE]->id[ABS_HYD];
			c_conf_diff = (pspecies->plink->pvar_ttc->var[id1]*q1 + pspecies->plink->pvar_ttc->var[id2]*q2)/(q1+q2);
			pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[WEST_TTC] = c_conf_diff; // SW 25/04/2018 modify valcl to s_ft
			      
			}*/				 
		      else if(pele->face[X_HYD][TWO]->element[TWO] == NULL) // element aval
			{
			  //pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[EAST_TTC] = pspecies->plink->pvar_ttc->var[id_abs_ele]; // 
			}
		      // else 
		      // {
		      for(icard = 0; icard < NELEMENT_HYD; icard++)
			{
			  pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard+2] = 0.0; // face NORTH_TTC and SOUTH_TTC
			  for(napp = 0; napp < pele->face[Y_HYD][icard]->ninflows; napp++)
			    {
			      if(pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[e][nsub] != NULL)
				{
				  pele->face[Y_HYD][icard]->pt_inflows[napp]->discharge = TS_function_t_pointed(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->discharge,fp); // SW 03/04/2020
				  q_inflow = TS_function_value_t(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->discharge,fp);
				  //if(id_abs_ele == 3392)
				  // printf("ok\n");
				  //surf = pele->face[X_HYD][ONE]->hydro->Surf;
				  q_inflow_area += q_inflow;
				  if(q_inflow > 0.)
				    {
				    // MH 06/12/2021 to calculate MOD123MOP123 of lateral inflows fromm TOC*share(b1)
				    if (e == MOD || e == MOP)
				      {
				      // to differentiate with rivers not having macrospecies 
				      if( (pele->face[Y_HYD][icard]->pt_inflows[napp]->flow_in_macrospecies[TOC] != NULL) && (Simul->passim->param_range[B1_RIVER_DA][PARAM_UP] > INFINITE_DIV_TS )  )
					{
					  pele->face[Y_HYD][icard]->pt_inflows[napp]->flow_in_macrospecies[TOC] = TS_function_t_pointed(t, pele->face[Y_HYD][icard]->pt_inflows[napp]->flow_in_macrospecies[TOC],fp);
					  toc_inflow = TS_function_value_t(t, pele->face[Y_HYD][icard]->pt_inflows[napp]->flow_in_macrospecies[TOC],fp);					      
					  //toc_inflow = 5;

					  // MOD123 or MOP123 = TOC * corresponding share
					  val_inflow = toc_inflow * PROSE_DA_mod_mop_fract(TOC, e, nsub, Simul, np, fp);
					  //val_inflow = 5 ;
					  //LP_printf(fp,"val_SIDEFLOW = %3.2f \n",val_inflow);
					  //LP_printf(fp,"np = %i, TOC_tribu = %3.4f ,val_inflow = %3.4f \n",np, toc_inflow/83.3, val_inflow/83.3);
					}
				      else
					{
					  pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[e][nsub] = TS_function_t_pointed(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[e][nsub],fp); // SW 03/04/202
					  val_inflow = TS_function_value_t(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[e][nsub],fp);
					  //LP_printf(fp,"val_inflow_2 = %3.2f \n",val_inflow);
					}
				      // val_inflow = TS_function_value_t(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[e][nsub],fp);
				    } else
				      {
				      pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[e][nsub] = TS_function_t_pointed(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[e][nsub],fp); // SW 03/04/202
				      val_inflow = TS_function_value_t(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[e][nsub],fp);
				      }
				    pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard+2] += q_inflow*val_inflow; // NEW_FACE_TTC flux impose
				    }
				  else
				    { // prelevement
				    //val_inflow = pspecies->plink->pvar_ttc->var[id_abs_ele]; //Simul->psimul_bio[id_abs_ele]->section->compartments[WATER][0]->pspecies[e][nsub]->C; 
				    //volume = pele->length*pele->center->hydro->Surf; // volume de la maille
				    pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard+2] += q_inflow; //volume; // NEW_FACE_TTC flux impose
				    //pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard+2] += 0.16*val_inflow;
				  }
				}    
			    }
				
			}
		      //}
		      //break;
		    }
		  /*** wet surf ***/ // SW 30/01/2020
		  //if(ne == 0 && r == 1)
		  //LP_printf(fp,"id_abs_ele == 0\n");
		  h_iter = pele->center->hydro->H[ITER_HYD];
		  h_i = pele->center->hydro->H[T_HYD];
		  //pspecies->plink->pparam_calc_ttc->coeff[SURF_T_TTC][id_abs_ele] = TS_function_value_t(h_i,pele->center->hydro->surf,fp); //
		  //pspecies->plink->pparam_calc_ttc->coeff[SURF_TITER_TTC][id_abs_ele] = TS_function_value_t(h_iter,pele->center->hydro->surf,fp); // pele->center->hydro->Surf;
          
		  /* SW 09/06/2021 use HL as wet surface */
		  pspecies->plink->pparam_calc_ttc->coeff[SURF_TITER_TTC][id_abs_ele] = TS_function_value_t(h_iter,pele->center->hydro->width,fp) * h_iter;
		  pspecies->plink->pparam_calc_ttc->coeff[SURF_T_TTC][id_abs_ele] = pele->center->hydro->Width * h_i;
		  
		}					
	    }			 
		
	}
	    
    }
  //#ifdef OMP
  //	   }
  //#endif	
}


void PROSE_find_apport_t_one_annex_species(s_species_ttc *pspecies, s_chyd *pchyd, double t, int nsub, int phy,FILE *fp)
{
	int r, ne;
	int napp, icard;
	int id_abs_ele;
	//double dt;
	double qapp, varapp,qapp_area, h_i,h_iter;
	//double surf,volume; // SW 17/05/2018  calculer une vitesse pour apport, on prend pele->face[X_HYD][ONE]->hydro->Surf, sont face amont de l'element
	double coef;
	//int id1,id2;
	//double q1,q2,c_conf_diff;
	s_reach_hyd *preach;
	s_element_hyd *pele;

	switch(phy){
		case PHYF : coef = PHY2PHYF; break;
		case PHYS : coef = PHY2PHYS; break;
		case PHYR : coef = PHY2PHYR; break;
		default : LP_error(fp,"error for phy coefficient\n");
	}
	
	//dt = pspecies->pchronos->dt;
//#ifdef OMP //SW 06/09/2018 ajout block openmp
//        int taille;
 //       int nthreads;
//        nthreads=Simul->psmp->nthreads;
//        omp_set_num_threads(nthreads);

//        Simul->psmp->chunk=PC_set_chunk_size_silent(fp,pchyd->counter->nreaches -1,nthreads);
//        taille=Simul->psmp->chunk;
//#pragma omp parallel shared(nthreads,taille,pspecies,pchyd,nsub) private(preach,ne,pele,id_abs_ele,icard,napp,varapp,qapp)
 // {
//#pragma omp for schedule(dynamic,taille) 
//#endif
	for (r = 0; r < pchyd->counter->nreaches; r++) 
	{
       preach = pchyd->p_reach[r];
       for(ne = 0; ne < preach->nele; ne++) 
	   {
		 pele = preach->p_ele[ne];
		 id_abs_ele = pele->id[ABS_HYD];
         if(ne == 0 && HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM)== DISCHARGE) // apport amont DIRI_FACE_TTC	
         {
			 if(pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[PHY][nsub] != NULL)
			 {	 
		       //qapp = pele->face[X_HYD][ONE]->hydro->Q[T_HYD];
                           pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[PHY][nsub] = TS_function_t_pointed(t,pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[PHY][nsub],fp);
			   pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[WEST_TTC] = coef*TS_function_value_t(t,pele->face[X_HYD][ONE]->pt_inflows[0]->app_bio[PHY][nsub],fp); 
			 }
			 else
			     pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[WEST_TTC] = 0.0; 
		 }
         /*if(ne == 0 && HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM)== CONF_DIFF) // apport amont DIRI_FACE_TTC	
         {		 
		    q1 = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->face[X_HYD][TWO]->hydro->Q[T_HYD];
			id1 = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]->id[ABS_HYD];
			q2 = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][1]->element[ONE]->face[X_HYD][TWO]->hydro->Q[T_HYD];
			id2 = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][1]->element[ONE]->id[ABS_HYD];
		    c_conf_diff = (pspecies->plink->pvar_ttc->var[id1]*q1 + pspecies->plink->pvar_ttc->var[id2]*q2)/(q1+q2);
			pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[WEST_TTC] = c_conf_diff; // SW 25/04/2018 modify valcl to s_ft
		 }*/   
		else if(pele->face[X_HYD][TWO]->element[TWO] == NULL) // element aval
		{
			//pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[EAST_TTC] = pspecies->plink->pvar_ttc->var[id_abs_ele]; // 
		}
        //else 
		//{
		   for(icard = 0; icard < NELEMENT_HYD; icard++)
		   {
			   pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard+2] = 0.0; // face NORTH_TTC and SOUTH_TTC
			   for(napp = 0; napp < pele->face[Y_HYD][icard]->ninflows; napp++)
			   {
				   if(pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[PHY][nsub] != NULL)
				   {
                                          pele->face[Y_HYD][icard]->pt_inflows[napp]->discharge = TS_function_t_pointed(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->discharge,fp);
					  qapp = TS_function_value_t(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->discharge,fp);
					  //surf = pele->face[X_HYD][ONE]->hydro->Surf;
					  if(qapp >0.)
					  {
                                          pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[PHY][nsub] = TS_function_t_pointed(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[PHY][nsub],fp);
					  varapp = TS_function_value_t(t,pele->face[Y_HYD][icard]->pt_inflows[napp]->app_bio[PHY][nsub],fp);
   
                      pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard+2] += coef*qapp*varapp; // NEW_FACE_TTC flux impose
                      }
     				  else
					  {//varapp = 0.;  // SW prelevement
				      //volume = pele->length*pele->center->hydro->Surf;
					  pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard+2] += qapp;//volume; // NEW_FACE_TTC flux impose
							  //pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard+2] += 0.16*varapp;
					  }						
						}    
				      }
					  
				   }
	  /*** wet surf ***/ // SW 30/01/2020
	  //if(ne == 0 && r == 1)
	  //LP_printf(fp,"id_abs_ele == 0\n");
	  h_iter = pele->center->hydro->H[ITER_HYD];
	  h_i = pele->center->hydro->H[T_HYD];
	  //pspecies->plink->pparam_calc_ttc->coeff[SURF_T_TTC][id_abs_ele] = TS_function_value_t(h_i,pele->center->hydro->surf,fp); //
	  //pspecies->plink->pparam_calc_ttc->coeff[SURF_TITER_TTC][id_abs_ele] = TS_function_value_t(h_iter,pele->center->hydro->surf,fp); // pele->center->hydro->Surf;
	  
          /* SW 09/06/2021 use HL as wet surface */
          pspecies->plink->pparam_calc_ttc->coeff[SURF_TITER_TTC][id_abs_ele] = TS_function_value_t(h_iter,pele->center->hydro->width,fp) * h_iter;
          pspecies->plink->pparam_calc_ttc->coeff[SURF_T_TTC][id_abs_ele] = pele->center->hydro->Width * h_i;
 
	//}				   
	//}
	}
	}
//#ifdef OMP
//  }
//#endif  
}

void PROSE_find_apport_t_all_annex_species(s_species_ttc ***p_phy_species, s_chyd *pchyd, double t, FILE *fp)
{
    int phy,nsub;
	s_species_ttc *pspecies;
	
	for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
	{
		for(phy = 0; phy < 3; phy++)
		{
			pspecies = p_phy_species[nsub][phy];
			PROSE_find_apport_t_one_annex_species(pspecies, pchyd, t, nsub, phy,fp);
		}
	}
}

/*
void PROSE_fill_cl_var_t(s_species_ttc **p_species, int nspecies_ttc, int nele_ttc, s_chyd *pchyd, double t, FILE *fp)
{
  //number of species for transport
  int nspecies, ns = 0;
  int nele;
  int id_abs_ele_amont,id_abs_ele_aval;
  int e, nsub;
  char *name_ttc, *name_rive;
  
  //s_species_ttc **p_species;
  s_reach_hyd *preach;
  //s_face_hyd *pface;
  s_element_hyd *pele_amont;
  s_element_hyd *pele_aval;
  nspecies = nspecies_ttc;
  nele = nele_ttc;
  
  name_ttc = p_species[ns]->name;
  for(ns = 0; ns < nspecies; ns++)
  {
	pele_amont = pchyd->p_reach[0]->p_ele[0]; // premier element du bief premier bief
	preach = pchyd->p_reach[pchyd->counter->nreaches - 1]; // dernier bief
    pele_aval = preach->p_ele[preach->nele-1]; //dernier element du dernier bief
	id_abs_ele_amont = pele_amont->id[ABS_HYD];
	id_abs_ele_aval = pele_aval->id[ABS_HYD];
	for(e = 0; e < NSPECIES; e++)
	{
	   for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++) 
	   {
	      name_rive = Simul->psimul_bio[0]->section->compartments[WATER][0]->pspecies[e][nsub]->name;
		  if(strcmp(name_ttc,name_rive) == 0)
		  {
			 p_species[ns]->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele_amont].valcl[WEST_TTC] = 71.48; // SW 23/04/2018 apport a faire TS_find_value_t
		     p_species[ns]->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele_aval].valcl[EAST_TTC] = 0.; //p_species[ns]->plink->pvar_ttc->var[id_abs_ele_aval];
		  }
	  }
	  
    }
  }
}*/

void PROSE_fill_surf_ttc_one_species(s_species_ttc *pspecies, s_chyd *pchyd, FILE *fp)
{
  int ne;
  int r,id_abs_ele;
  double h_iter,h_i;
  s_reach_hyd *preach;
  s_element_hyd *pele;
  //#ifdef OMP //SW 06/09/2018 ajout block openmp
  //       int taille;
  //       int nthreads;
  //       nthreads=Simul->psmp->nthreads;
  //       omp_set_num_threads(nthreads);
  
  //       Simul->psmp->chunk=PC_set_chunk_size_silent(fp,pchyd->counter->nreaches -1,nthreads);
  //       taille=Simul->psmp->chunk;
  //#pragma omp parallel shared(nthreads,taille,pspecies,pchyd) private(preach,ne,pele,id_abs_ele,h_iter)
  // {
  //#pragma omp for schedule(dynamic,taille) 
  //#endif   
  for(r = 0; r < pchyd->counter->nreaches; r++)
    {
      preach = pchyd->p_reach[r];
      
      for(ne = 0; ne < preach->nele; ne++)
	{
	  pele = preach->p_ele[ne];
	  id_abs_ele = pele->id[ABS_HYD];
	  //if(ne == 0 && r == 1)
	  //LP_printf(fp,"id_abs_ele == 0\n");
	  h_iter = pele->center->hydro->H[ITER_HYD];
	  h_i = pele->center->hydro->H[T_HYD];
	  //pspecies->plink->pparam_calc_ttc->coeff[SURF_T_TTC][id_abs_ele] = pele->center->hydro->Surf; //TS_function_value_t(h_i,pele->center->hydro->surf,fp); //
	  //pspecies->plink->pparam_calc_ttc->coeff[SURF_TITER_TTC][id_abs_ele] = pele->center->hydro->Surf; //TS_function_value_t(h_iter,pele->center->hydro->surf,fp); // pele->center->hydro->Surf; //
          
          pspecies->plink->pparam_calc_ttc->coeff[SURF_TITER_TTC][id_abs_ele] = TS_function_value_t(h_iter,pele->center->hydro->width,fp) * h_iter;
          pspecies->plink->pparam_calc_ttc->coeff[SURF_T_TTC][id_abs_ele] = pele->center->hydro->Width * h_i;

	}
    }
  //#ifdef OMP
  // }
  //#endif  
}

void PROSE_fill_surf_ttc_all_species(s_species_ttc **p_species, s_species_ttc ***p_phy_species, int nspecies, s_chyd *pchyd, FILE *fp)
{
   int e = 0;	
   int phy,nsub;
   s_species_ttc *pspecies;

	for(e = 0; e < nspecies; e++)
	{
        pspecies = p_species[e];
		PROSE_fill_surf_ttc_one_species(pspecies,pchyd,fp);
	}
	for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
	{
		for(phy = 0; phy < 3; phy++)
		{
			pspecies = p_phy_species[nsub][phy];
			PROSE_fill_surf_ttc_one_species(pspecies,pchyd,fp);
		}
	}		
	
}

void PROSE_set_iapplic_all_species(s_species_ttc **p_species, int nspecies, int nele,FILE *fp)
{
	int e;
	s_species_ttc *pspecies;
	//int iappli=3;
	for(e = 0; e < nspecies; e++){
		pspecies = p_species[e];
		pspecies->iappli_gc = e+3;//e%(IAPPLI - 5) + 3; // SW 13/06/2018 TTC_GC commence par 3
        if(pspecies->pgc != NULL)
		{
			free(pspecies->pgc->b);
			free(pspecies->pgc);
		}
		pspecies->pgc = GC_create_gc(pspecies->iappli_gc);
		pspecies->pgc->b = TTC_create_tab_double(nele);
		pspecies->pgc->appl_nb = pspecies->iappli_gc;
		//pspecies->pgc->sw_int[NS_GC]= 2;
		}
	
}

void PROSE_set_iapplic_all_annex_species(s_species_ttc ***p_phy_species, int nele,FILE *fp)
{
	int phy,nsub;
	s_species_ttc *pspecies;
	
	for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[PHY]; nsub++)
	{
		for(phy = 0; phy < 3; phy++)
		{	
	      pspecies = p_phy_species[nsub][phy];
		  pspecies->iappli_gc = 21 + nsub + phy;//nspecies(nsub + phy)%(IAPPLI - 5) + 3;
		  if(pspecies->pgc != NULL)
		  {
		    free(pspecies->pgc->b);
		    free(pspecies->pgc);
		  }
		  pspecies->pgc = GC_create_gc(pspecies->iappli_gc);
		  pspecies->pgc->b = TTC_create_tab_double(nele);		  
		}
	}
}

void Prose_cal_advflux(s_carac_ttc *pcarac_ttc,s_param_calc_ttc *pparam_calc_ttc,s_species_ttc *pspecies,int i, double *H_flux_ttc, FILE *fpout)
{
  double dist,delta_L,*var,u,coef_reoxy;
  int tcl;
  int type;
  int icard,id_neigh,sub_icard,id_neigh1,id_neigh2;
  double theta,dlip,cl_val,rhow,cw,rhowcw;
  double u_id_neigh1,u_id_neigh2;
  var=pspecies->plink->pvar_ttc->var;
  dist=pspecies->plink->pbase_ttc->param_syst[SIZE_TTC][i]; // taille de la maille
  type=pspecies->type; 
  if(type==HEAT_TTC)
		{
		rhow=pspecies->plink->pbase_ttc->pthermic[WATER_TTC]->param[RHO_TTC][i];
		cw=pspecies->plink->pbase_ttc->pthermic[WATER_TTC]->param[HEAT_CAP_TTC][i];
		rhowcw=rhow*cw;
		}
	else if(type==SOLUTE_TTC)
		{
		rhowcw=1;
		}
	for (icard=0;icard<NB_CARD_TTC;icard++) {
		for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++){ // SW 29/05/2018 pour confluence
		//u=*pspecies->plink->pu_ttc->puface[i].uface[icard];
		u=*pspecies->plink->pu_ttc->puface[i].uface[icard][sub_icard]; // SW 29/05/2018 pour confluence
		dlip=pparam_calc_ttc->pdisp_ttc->pdisp_face[i].disp_face[icard]; // dispersion + conduction
		tcl=pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[i].icl[icard];
		cl_val=pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[i].valcl[icard];

			if (icard==EAST_TTC)
			{
				if(u>0)
				{
/// WARNING A CHANGER LORSQUE AVEC LA TAILLE DES SUBFACE dist doit etre egale à la taille de la subface.
		   		//qadv= -(rhowcw*u*var[i])*dist; // SW 17/10/2018 flux sortant de la maille
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard] = -(rhowcw*u*var[i]);
				}
				else if (tcl==NO_BOUND_TTC && u<0) //face active
				{
				//id_neigh=pspecies->plink->pneigh_ttc->pivois[i].ivois[icard];
				id_neigh=pspecies->plink->pneigh_ttc->pivois[i].ivois[icard][sub_icard]; // SW 29/05/2018 pour confluence
	 			//delta_L=pspecies->plink->pneigh_ttc->pdelta_L[i].delta_L[icard]; // distance entre le centre de la maille et de son voisin
                delta_L=pspecies->plink->pneigh_ttc->pdelta_L[i].delta_L[icard][sub_icard]; // SW 29/05/2018 pour confluence
/// WARNING A CHANGER LORSQUE AVEC LA TAILLE DES SUBFACE dist doit etre egale à la taille de la subface.
				if(id_neigh != NONE_TTC) // SW 29/05/2018 pour confluence
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=(rhowcw*u*var[id_neigh]); // flux entrant de la maille // SW
				}
				else if (tcl==DIRI_FACE_TTC && u<0) //face active no use in prose-p
				{
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=(rhowcw*u*cl_val);
				}
				else //face flux impose
				{
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=0;
				}
			}
			if (icard==WEST_TTC)
			{
				if(u<0)
				{
/// WARNING A CHANGER LORSQUE AVEC LA TAILLE DES SUBFACE dist doit etre egale à la taille de la subface.
		   		//qadv= (rhowcw*u*var[i])*dist; // flux sortant de la maille
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard] = (rhowcw*u*var[i]);
				}
				else if (tcl==NO_BOUND_TTC && u>0) //face active
				{
				//if((i == 98 || i==105) && (pspecies->oxygen == YES_TS))
					//printf("yes\n");
				//id_neigh=pspecies->plink->pneigh_ttc->pivois[i].ivois[icard];
				id_neigh=pspecies->plink->pneigh_ttc->pivois[i].ivois[icard][sub_icard]; //SW 29/05/2018 pour confluence
	 			//delta_L=pspecies->plink->pneigh_ttc->pdelta_L[i].delta_L[icard]; // distance entre le centre de la maille et de son voisin
				delta_L=pspecies->plink->pneigh_ttc->pdelta_L[i].delta_L[icard][sub_icard]; //SW 29/05/2018 pour confluence
/// WARNING A CHANGER LORSQUE AVEC LA TAILLE DES SUBFACE dist doit etre egale à la taille de la subface.
				if(id_neigh != NONE_TTC) //SW 29/05/2018 pour confluence
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]= (rhowcw*u*var[id_neigh]); 
				}
				else if (tcl==DIRI_FACE_TTC && u>0) //face active
				{
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=(rhowcw*u*cl_val);
				}
				else if(tcl == CONF_DIFF_TTC)
				{
					id_neigh1 = pspecies->plink->pneigh_ttc->pivois[i].ivois[icard][0];
					id_neigh2 = pspecies->plink->pneigh_ttc->pivois[i].ivois[icard][1];
					u_id_neigh1 = *pspecies->plink->pu_ttc->puface[id_neigh1].uface[EAST_TTC][0];
					u_id_neigh2 = *pspecies->plink->pu_ttc->puface[id_neigh2].uface[EAST_TTC][0];
					pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][0]= (rhowcw*u*(var[id_neigh1]*u_id_neigh1 + var[id_neigh2]*u_id_neigh2)/(u_id_neigh1 + u_id_neigh2)); 								
				    sub_icard = SUB_CARD_TTC;
				}
				else if(tcl == BARRAGE_TTC)
				{
					if(pspecies->oxygen == YES_TS)
					{
						id_neigh=pspecies->plink->pneigh_ttc->pivois[i].ivois[icard][sub_icard];
						coef_reoxy = TTC_reoxygeneration_rhs(pspecies->plink->preoxy_ttc->rd[i], pspecies->plink->preoxy_ttc->rn[i], pspecies->plink->preoxy_ttc->frac_debit_oxy[i], fpout);
						//u -= pspecies->plink->preoxy_ttc->debit_surverse[i]*coef_reoxy; //implicite theta = 1.
						//pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard] = (rhowcw*u*var[id_neigh]*(1 - coef_reoxy));
						//pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard] -= pspecies->plink->preoxy_ttc->debit_surverse[i]*coef_reoxy*rhowcw*var[id_neigh];
						pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard] = (u - pspecies->plink->preoxy_ttc->debit_surverse[i]*coef_reoxy)*rhowcw*var[id_neigh];
						pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard] += pspecies->plink->preoxy_ttc->debit_surverse[i]*coef_reoxy*pspecies->plink->preoxy_ttc->osat*rhowcw;
						//pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard] += u*coef_reoxy*pspecies->plink->preoxy_ttc->osat*rhowcw;
						sub_icard = SUB_CARD_TTC;
					}
                                        if(strcmp("Temp", pspecies->name) == 0) // SW 05/05/2021 for temperature, it is the same than WEST_TTC when u > 0.
                                        {
                                            id_neigh=pspecies->plink->pneigh_ttc->pivois[i].ivois[icard][sub_icard];
                                            if(id_neigh != NONE_TTC)
				                pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]= (rhowcw*u*var[id_neigh]);
                                        }
				}
				else //face flux impose
				{
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=0;
				}
			}
			if (icard==NORTH_TTC)
			{
				if(u>0)
				{
/// WARNING A CHANGER LORSQUE AVEC LA TAILLE DES SUBFACE dist doit etre egale à la taille de la subface.
		   		pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=(rhowcw*u*var[i]);
				}
				else if (tcl==NO_BOUND_TTC && u<0) //face active
				{
				//id_neigh=pspecies->plink->pneigh_ttc->pivois[i].ivois[icard];
				id_neigh=pspecies->plink->pneigh_ttc->pivois[i].ivois[icard][sub_icard]; //SW 29/05/2018 pour confluence
	 			//delta_L=pspecies->plink->pneigh_ttc->pdelta_L[i].delta_L[icard]; // distance entre le centre de la maille et de son voisin
				delta_L=pspecies->plink->pneigh_ttc->pdelta_L[i].delta_L[icard][sub_icard]; //SW 29/05/2018 pour confluence
/// WARNING A CHANGER LORSQUE AVEC LA TAILLE DES SUBFACE dist doit etre egale à la taille de la subface.
				if(id_neigh != NONE_TTC) //SW 29/05/2018 pour confluence
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=(rhowcw*u*var[id_neigh]);
				}
				else if (tcl==DIRI_FACE_TTC && u<0) //face active
				{
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=(rhowcw*u*cl_val);
				sub_icard = SUB_CARD_TTC; // SW 29/05/2018 no confluence pour DIRI_FACE_TTC
				}
				else //face flux impose
				{
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=0;
				}
			}

			if (icard==SOUTH_TTC)
			{
				if(u<0)
				{
/// WARNING A CHANGER LORSQUE AVEC LA TAILLE DES SUBFACE dist doit etre egale à la taille de la subface.
		   		pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=(rhowcw*u*var[i]);
				}
				else if (tcl==NO_BOUND_TTC && u>0) //face active
				{
				//id_neigh=pspecies->plink->pneigh_ttc->pivois[i].ivois[icard];
				id_neigh=pspecies->plink->pneigh_ttc->pivois[i].ivois[icard][sub_icard]; //SW 29/05/2018 pour confluence
	 			//delta_L=pspecies->plink->pneigh_ttc->pdelta_L[i].delta_L[icard]; // distance entre le centre de la maille et de son voisin
				delta_L=pspecies->plink->pneigh_ttc->pdelta_L[i].delta_L[icard][sub_icard]; //SW 29/05/2018 pour confluence
/// WARNING A CHANGER LORSQUE AVEC LA TAILLE DES SUBFACE dist doit etre egale à la taille de la subface.
				if(id_neigh != NONE_TTC) //SW 29/05/2018 pour confluence
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=(rhowcw*u*var[id_neigh]);
				}
				else if (tcl==DIRI_FACE_TTC && u>0) //face active
				{
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=(rhowcw*u*cl_val);
				sub_icard = SUB_CARD_TTC; // SW 29/05/2018 no confluence pour DIRI_FACE_TTC
				}
				else //face flux impose
				{
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][sub_icard]=0;
				sub_icard = SUB_CARD_TTC; // SW 29/05/2018 no confluence pour DIRI_FACE_TTC
				}
			}
		if(tcl == NEU_FACE_TTC)	
		{
                    if(type == SOLUTE_TTC)
                    { 
			cl_val=pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[i].valcl[icard];
			if(cl_val > 0)
			{
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][0] = cl_val*rhowcw;				
			}
			else if(cl_val < 0) // flux sortant
			{
				pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][0] = cl_val*var[i]*rhowcw;
			}
			sub_icard = SUB_CARD_TTC;
                    }
                    else if(type == HEAT_TTC) // SW 05/05/2021 add here for heat flux
                    {
                        if(icard==NORTH_TTC && Simul->calc_mode[SEB] == YES_TS) // SW 05/05/2021 need to minus atmospheric heat flux which is calculated in PROSE_calc_energy_init_and_flux_from_seb()
                            cl_val = pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[i].valcl[icard] - H_flux_ttc[i]; // apport ou rejet unit is K m2/s
                        
                        pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[i].face[icard][0] = cl_val*rhowcw*dist; // unit is J/s
                        sub_icard = SUB_CARD_TTC;
                    }
		}
		}//fin sub_icard
	}
	
}


void Prose_cal_tot_advflux(s_carac_ttc *pcarac_ttc,s_param_calc_ttc *pparam_calc_ttc,s_species_ttc *pspecies, double *H_flux_ttc, FILE *fpout)
{
    int i, nele;
	nele=pcarac_ttc->count[NELE_TTC];
    for(i = 0; i < nele; i++)
	{
      Prose_cal_advflux(pcarac_ttc,pparam_calc_ttc,pspecies,i, H_flux_ttc, fpout);
    }      
}

void PROSE_calc_mb_one_species(s_species_ttc *pspecies, s_total_mb *pmb, int i, int nele, int np, double dt, FILE *fp)
{
	int ne, e, nsub,icard,sub_icard;
	char *name_ttc, *name_rive;
	//double volume;
    //nthreads = Simul->psmp->nthreads;
	name_ttc = pspecies->name;
	double mass_entrant,mass_entrant1;
	double mass_sortant,mass_sortant1;
	double mass_lateral_inflow, mass_reoxy,coef_reoxy;
	double flux_mass;

	s_lp_pk_hyd *ppk;
	s_element_hyd *pele;
	s_reach_hyd *preach;
	s_chyd *pchyd;	
	int nr,ns,id_neigh;
	int find_upstream_grid, find_downstream_grid;
	int limits_nr, limits_ne;
	find_upstream_grid = NO_TS;
	find_downstream_grid = NO_TS;

	for(e = 0; e < NSPECIES; e++)
	{
	   for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++) 
	   {
		    name_rive = Simul->psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[e][nsub]->name;
			if(strcmp(name_ttc,name_rive) == 0) // find the same species
			{
	            //#ifdef OMP
	            //omp_set_num_threads(nthreads);
	            //#pragma omp parallel for 
	            //#endif 
				mass_entrant = 0.;
				mass_sortant = 0.;
				mass_lateral_inflow = 0;
				mass_reoxy = 0.;
				/*** SW 06/12/2019 add mass_balance output for the defined domaine by user ***/
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
						for(ne = 0; ne < preach->nele; ne++) 
						{
							pele = preach->p_ele[ne];					
							if ((pele->center->pk >= ppk->pk_up) && (pele->center->pk <= ppk->pk_down)) 
							{
								ns = pele->id[ABS_HYD]; // id_abs_ele 
								for(icard=NORTH_TTC;icard<NB_CARD_TTC;icard++) 
								{
									for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++)
									{
										mass_lateral_inflow +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[icard][sub_icard];
										//if(flux_mass > 0.)	
											//mass_lateral_inflow += flux_mass;
										//else if(flux_mass < 0.)
											//mass_sortant += flux_mass;
									}
								}
							}
							if((find_upstream_grid == NO_TS) || (find_downstream_grid == NO_TS)) // SW 19/03/2020 upstream element inside a reach
							{
								if((ne != 0) && (pele->center->pk >= ppk->pk_up) && (pele->face[X_HYD][ONE]->element[ONE]->center->pk < ppk->pk_up))
								{
									mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][0]; // upstream inflow
									//find_upstream_grid = YES_TS;
									mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][1]; // upstream inflow
                                                                         //LP_printf(fp,"id_ele = %d, e = %d, nsub = %d fin = %f\n",ns,e,nsub,mass_entrant);
									/* we check if the upstream point is located on a branch*/
									if(preach->limits[DOWNSTREAM]->nupstr_reaches > 1) // we have at least two branches here
									{
										
										for(limits_nr = 0; limits_nr < preach->limits[DOWNSTREAM]->nupstr_reaches; limits_nr++)
										{
											if(preach->id != preach->limits[DOWNSTREAM]->p_reach[UPSTREAM][limits_nr]->id)
											{
												limits_ne = preach->limits[DOWNSTREAM]->p_reach[UPSTREAM][limits_nr]->p_ele[preach->limits[DOWNSTREAM]->p_reach[UPSTREAM][limits_nr]->nele - 1]->id[ABS_HYD];
												mass_entrant -= pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[limits_ne].face[EAST_TTC][0]; // EAST is negatif
											}
										}
									}
								}
								if((ne != preach->nele-1) && (pele->center->pk <= ppk->pk_down) && (pele->face[X_HYD][TWO]->element[TWO]->center->pk > ppk->pk_down))
								{
									mass_sortant += pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[EAST_TTC][0];
									//find_downstream_grid == YES_TS;
                                    if(preach->limits[UPSTREAM]->ndownstr_reaches > 1)
									{
										for(limits_nr = 0; limits_nr < preach->limits[UPSTREAM]->ndownstr_reaches; limits_nr++)
										{
											if(preach->id != preach->limits[UPSTREAM]->p_reach[DOWNSTREAM][limits_nr]->id)
											{
												limits_ne = preach->limits[UPSTREAM]->p_reach[DOWNSTREAM][limits_nr]->p_ele[0]->id[ABS_HYD];
												mass_sortant -= pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[limits_ne].face[WEST_TTC][0]; // WEST is positif
											}
										}										
									}										
								}
								if((ne == 0) && (pele->center->pk >= ppk->pk_up))
								{
									if(pele->face[X_HYD][ONE]->limits[UPSTREAM]->type == DISCHARGE) // first gird of the model
									{
										mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][0]; // upstream inflow
                                                                                //LP_printf(fp,"id_ele = %d, e = %d, nsub = %d fin = %f\n",ns,e,nsub,mass_entrant);
										//find_upstream_grid = YES_TS;										
									}
									else if(pele->face[X_HYD][ONE]->limits[UPSTREAM]->type == CONFLUENCE)
									{
										
										if(pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][ONE]->element[ONE]->center->pk > ppk->pk_up)
										{
											if(strcmp(pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][ONE]->element[ONE]->reach->river,ppk->river) != 0) // another river simulated, like Marne
											{
												id_neigh = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][ONE]->element[ONE]->id[ABS_HYD];
												mass_lateral_inflow += pspecies->plink->pvar_ttc->var[id_neigh] * pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][ONE]->element[ONE]->face[X_HYD][TWO]->hydro->Q[T_HYD];
											}
										}										
										if(pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][TWO]->element[ONE]->center->pk > ppk->pk_up)
										{
											if(strcmp(pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][TWO]->element[ONE]->reach->river,ppk->river) != 0) // another river simulated, like Marne
											{
												id_neigh = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][TWO]->element[ONE]->id[ABS_HYD];
												mass_lateral_inflow += pspecies->plink->pvar_ttc->var[id_neigh] * pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][TWO]->element[ONE]->face[X_HYD][TWO]->hydro->Q[T_HYD];
											}
										}
                                        if((pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][ONE]->element[ONE]->center->pk < ppk->pk_up) && (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][TWO]->element[ONE]->center->pk < ppk->pk_up))
										{
										    mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][0]; 
										    mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][1]; 											
										}											
									}
									else if((pele->face[X_HYD][ONE]->limits[UPSTREAM]->type == HYDWORK) && (e == O2)) // reoxygeneration
									{
				                            id_neigh=pspecies->plink->pneigh_ttc->pivois[ns].ivois[WEST_TTC][0];
						                    coef_reoxy = TTC_reoxygeneration_rhs(pspecies->plink->preoxy_ttc->rd[ns], pspecies->plink->preoxy_ttc->rn[ns], pspecies->plink->preoxy_ttc->frac_debit_oxy[ns], fp);
						                     
												
										mass_reoxy +=  pspecies->plink->preoxy_ttc->debit_surverse[ns]*coef_reoxy*pspecies->plink->preoxy_ttc->osat - 
										pspecies->plink->preoxy_ttc->debit_surverse[ns]*coef_reoxy*pspecies->plink->pvar_ttc->var[id_neigh];
									}									
									else if(pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][ONE]->element[ONE]->center->pk < ppk->pk_up)
									{
									    if(preach->limits[DOWNSTREAM]->nupstr_reaches > 1) // we have at least two branches here
									    {
										     //int limits_nr, limits_ne;
										    for(limits_nr = 0; limits_nr < preach->limits[DOWNSTREAM]->nupstr_reaches; limits_nr++)
										    {
											   if(preach->id != preach->limits[DOWNSTREAM]->p_reach[UPSTREAM][limits_nr]->id)
											   {
												    limits_ne = preach->limits[DOWNSTREAM]->p_reach[UPSTREAM][limits_nr]->p_ele[preach->limits[DOWNSTREAM]->p_reach[UPSTREAM][limits_nr]->nele - 1]->id[ABS_HYD];
												    mass_entrant -= pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[limits_ne].face[EAST_TTC][0]; // EAST is negatif
											   }
										    } 
									     }
										 //else
										 //{
											 
								        mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][0]; 
										mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][1]; 											 
										 //}
									}
								}
								if((find_downstream_grid == NO_TS) && (ne == preach->nele-1) && (pele->center->pk <= ppk->pk_down))
								{
                                    if(pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][ONE]->element[TWO]->center->pk > ppk->pk_down)
									{
										if(preach->limits[UPSTREAM]->ndownstr_reaches > 1)
										{
											for(limits_nr = 0; limits_nr < preach->limits[UPSTREAM]->ndownstr_reaches; limits_nr++)
											{
												if(preach->id != preach->limits[UPSTREAM]->p_reach[DOWNSTREAM][limits_nr]->id)
												{
													limits_ne = preach->limits[UPSTREAM]->p_reach[DOWNSTREAM][limits_nr]->p_ele[0]->id[ABS_HYD];
													mass_sortant -= pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[limits_ne].face[WEST_TTC][0]; // WEST is positif
												}
											}										
										}									
									    //else
											mass_sortant += pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[EAST_TTC][0];	
									}										
								}
								}
						}
						nr++;
					}
					// sum of all grid in defined domaine
					pmb->mbspecies[WATER][FIN_RIVE][e][nsub] += mass_entrant * dt; // SW 04/06/2019 pour nstep != 1
					pmb->mbspecies[WATER][FIN_LATERAL][e][nsub] +=  mass_lateral_inflow * dt;
					pmb->mbspecies[WATER][FOUT][e][nsub] += mass_sortant * dt;
					pmb->mbspecies[WATER][REOX][e][nsub] += mass_reoxy * dt;
				}
                /*** SW 06/12/2019 add mass_balance output for the defined domaine by user END***/
				else {				
				for(ne = 0; ne < nele; ne++) // for the entire domaine
				{			//mass_entrant1 = 0.;
				            //mass_sortant1 = 0.;
	                for (icard=0;icard<NB_CARD_TTC;icard++) {
		               for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++){
                         flux_mass =  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ne].face[icard][sub_icard];
                         if(flux_mass > 0.)	
						 {
                            mass_entrant += flux_mass;
							//mass_entrant1 += flux_mass;
						 }
                         else if(flux_mass < 0.)
						 {
                            mass_sortant += flux_mass;
							//mass_sortant1 += flux_mass;
						 }
					   }
					}
					//if(e == 0 && nsub == 0)
                      //fprintf(Simul->fphy,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ne,pspecies->plink->pvar_ttc->var_old[ne],mass_entrant1,mass_sortant1,pspecies->plink->pvar_ttc->var[ne],
				  //Simul->psimul_bio[np][ne]->section->compartments[WATER][0]->state->volume,Simul->pchyd->p_reach[ne/40]->p_ele[ne%40]->face[X_HYD][ONE]->hydro->Q[T_HYD],
				  //Simul->pchyd->p_reach[ne/40]->p_ele[ne%40]->face[X_HYD][TWO]->hydro->Q[T_HYD],Simul->pchyd->p_reach[ne/40]->p_ele[ne%40]->center->hydro->Surf);				
				}
				//pmb->mbspecies[WATER][FIN_RIVE][e][nsub] = mass_entrant * dt;
				//pmb->mbspecies[WATER][FOUT][e][nsub] = mass_sortant * dt;
				pmb->mbspecies[WATER][FIN_RIVE][e][nsub] += mass_entrant * dt; // SW 04/06/2019 pour nstep != 1
				pmb->mbspecies[WATER][FOUT][e][nsub] += mass_sortant * dt;
				//fclose(Simul->fphy);
			}
			}
		}			
	  }	
}


/* void PROSE_calc_mb_one_species(s_species_ttc *pspecies, s_total_mb *pmb, int i, int nele, int np, double dt, FILE *fp)
{
	int ne, e, nsub,icard,sub_icard;
	char *name_ttc, *name_rive;
	//double volume;
    //nthreads = Simul->psmp->nthreads;
	name_ttc = pspecies->name;
	double mass_entrant,mass_entrant1;
	double mass_sortant,mass_sortant1;
	double flux_mass;

	s_lp_pk_hyd *ppk;
	s_element_hyd *pele;
	s_reach_hyd *preach;
	s_chyd *pchyd;	
	int nr,ns;

	for(e = 0; e < NSPECIES; e++)
	{
	   for(nsub = 0; nsub < Simul->counter_bio->nsubspecies[e]; nsub++) 
	   {
		    name_rive = Simul->psimul_bio[np][0]->section->compartments[WATER][0]->pspecies[e][nsub]->name;
			if(strcmp(name_ttc,name_rive) == 0) // find the same species
			{
	            //#ifdef OMP
	            //omp_set_num_threads(nthreads);
	            //#pragma omp parallel for 
	            //#endif 
				mass_entrant = 0.;
				mass_sortant = 0.; */
				/*** SW 06/12/2019 add mass_balance output for the defined domaine by user ***/
				/*if(Simul->lp_pk[i] != NULL) // output pk defined
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
						for(ne = 0; ne < preach->nele; ne++) 
						{
							pele = preach->p_ele[ne];					
							if ((pele->center->pk >= ppk->pk_up) && (pele->center->pk <= ppk->pk_down)) 
							{
								ns = pele->id[ABS_HYD]; // id_abs_ele 
								for(icard=0;icard<NB_CARD_TTC;icard++) 
								{
									for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++)
									{
										flux_mass =  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[icard][sub_icard];
										if(flux_mass > 0.)	
											mass_entrant += flux_mass;
										else if(flux_mass < 0.)
											mass_sortant += flux_mass;
									}
								}
							}
						}
						nr++;
					}
					// sum of all grid in defined domaine
					pmb->mbspecies[WATER][FIN_RIVE][e][nsub] += mass_entrant * dt; // SW 04/06/2019 pour nstep != 1
					pmb->mbspecies[WATER][FOUT][e][nsub] += mass_sortant * dt;
				}
                /*** SW 06/12/2019 add mass_balance output for the defined domaine by user END***/
				/*else {				
				for(ne = 0; ne < nele; ne++) // for the entire domaine
				{			//mass_entrant1 = 0.;
				            //mass_sortant1 = 0.;
	                for (icard=0;icard<NB_CARD_TTC;icard++) {
		               for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++){
                         flux_mass =  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ne].face[icard][sub_icard];
                         if(flux_mass > 0.)	
						 {
                            mass_entrant += flux_mass;
							//mass_entrant1 += flux_mass;
						 }
                         else if(flux_mass < 0.)
						 {
                            mass_sortant += flux_mass;
							//mass_sortant1 += flux_mass;
						 }
					   }
					}
					//if(e == 0 && nsub == 0)
                      //fprintf(Simul->fphy,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",ne,pspecies->plink->pvar_ttc->var_old[ne],mass_entrant1,mass_sortant1,pspecies->plink->pvar_ttc->var[ne],
				  //Simul->psimul_bio[np][ne]->section->compartments[WATER][0]->state->volume,Simul->pchyd->p_reach[ne/40]->p_ele[ne%40]->face[X_HYD][ONE]->hydro->Q[T_HYD],
				  //Simul->pchyd->p_reach[ne/40]->p_ele[ne%40]->face[X_HYD][TWO]->hydro->Q[T_HYD],Simul->pchyd->p_reach[ne/40]->p_ele[ne%40]->center->hydro->Surf);				
				}
				//pmb->mbspecies[WATER][FIN_RIVE][e][nsub] = mass_entrant * dt;
				//pmb->mbspecies[WATER][FOUT][e][nsub] = mass_sortant * dt;
				pmb->mbspecies[WATER][FIN_RIVE][e][nsub] += mass_entrant * dt; // SW 04/06/2019 pour nstep != 1
				pmb->mbspecies[WATER][FOUT][e][nsub] += mass_sortant * dt;
				//fclose(Simul->fphy);
			}
			}
		}			
	  }	
} */

void PROSE_calc_mb_all_species(s_carac_ttc *pcarac_ttc, int nspecies, double t, double dt, int np, double *H_flux_ttc, FILE *fp)
{
	int ns,nele,i;
	s_species_ttc *pspecies;
	s_param_calc_ttc *pparam_calc_ttc;
	nele = pcarac_ttc->count[NELE_TTC];
    s_total_mb *pmb;	
    for(ns = 0; ns < nspecies; ns++)
	{
		pspecies = pcarac_ttc->p_species[ns];
		pparam_calc_ttc = pspecies->plink->pparam_calc_ttc;
		Prose_cal_tot_advflux(pcarac_ttc,pparam_calc_ttc,pspecies, H_flux_ttc, fp);
		/*** SW 06/12/2019 add for loop for mutiple mb output***/
		for (i = 0; i < Simul->counter_bio->nmb; i++) 
		{
			pmb = Simul->total_mb[np][i];
			if((t >= pmb->t[BEGINNING]) && (t <= pmb->t[END]))
			{ 		
				PROSE_calc_mb_one_species(pspecies, pmb, i, nele, np, dt, fp);
			}
		}
	}	
}

void PROSE_init_temp_specie(s_carac_ttc *pcarac_heat_ttc, int nele)
{
  pcarac_heat_ttc->p_species = (s_species_ttc **)malloc(sizeof(s_species_ttc*)); // Initialization of pcarac_ttc is done through the yacc input.
  pcarac_heat_ttc->p_species[0] = TTC_init_species(nele); // Only one specie allocated as one considers the heat transport equation only.
  pcarac_heat_ttc->p_species[0]->name = (char *)malloc(109*sizeof(char));
  strcpy(pcarac_heat_ttc->p_species[0]->name, "Temp");
  printf("Name: %s -> variable will be treated!\n", pcarac_heat_ttc->p_species[0]->name);
  pcarac_heat_ttc->p_species[0]->media_type = WATER_TTC;
  pcarac_heat_ttc->p_species[0]->calc_process[DIFFUSION_TTC] = NO_TS;
  pcarac_heat_ttc->p_species[0]->calc_process[ADVECTION_TTC] = YES_TS;
  pcarac_heat_ttc->p_species[0]->theta = Simul->pcarac_heat_ttc->theta;
}

void PROSE_init_heat_transport(s_species_ttc *pspecies, s_chyd *pchyd, int nele, double t, FILE *fp)
{
  int e, icard, sub_icard, id_abs_ele;
  //s_reach_hyd *preach;
  //s_element_hyd *pele;
  //s_face_hyd *pface;
  double var;
  s_link_ttc *plink;
  //s_ft *init_T;
  //int init_from_file;

  plink = pspecies->plink; 
  //init_T = Simul->init_T;
  //init_from_file = Simul->init_from_file;
  
  Simul->init_T = TS_browse_ft(Simul->init_T,BEGINNING); // SW 05/12/2019 for init_temp
  
  for(e = 0; e < nele; e++)
    {
      /*** init_bounds ***/
      for(icard = 0; icard < NB_CARD_TTC; icard++)
	{
	  for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++)
	    pspecies->plink->pneigh_ttc->pivois[e].ivois[icard][sub_icard] = NONE_TTC;
	  
	  pspecies->plink->pboundary_ttc->ptype_bound_ttc->picl[e].icl[icard] = NO_BOUND_TTC;
	}
      
      /*** init_temp ***/
      if (Simul->init_from_file == YES_TS)
	{
	  //var = TS_function_value_t(e,init_T,fp);
	  //plink->pvar_ttc->var[e] = T_0 + var;
	  if(Simul->init_T != NULL)
	  {
		  id_abs_ele = (int)Simul->init_T->t;
		  plink->pvar_ttc->var[id_abs_ele] = T_0 + Simul->init_T->ft;
		  Simul->init_T = Simul->init_T->next;
	  }
	  else
		  LP_error(fp,"number of values error. Check number of lines in init_T file.\n");
	}
      else
		  plink->pvar_ttc->var[e] = T_0 + 20.;
      //LP_printf(fp,"Initial temp of cell [%d]: %lf\n",e,plink->pvar_ttc->var[e]);

      /*** init_paramsystem ***/
      pspecies->plink->pbase_ttc->param_syst[POROSITY_TTC][e] = 1.;
      pspecies->plink->pbase_ttc->param_syst[DISPERSIVITY_TTC][e] = 0.; // SW 20/04/2018 a verifier

      /*** init_base_parameters ***/ 
      plink->pbase_ttc->pthermic[SOLID_TTC]->param[LAMBDA_TTC][e]=LAMBDAW_TTC;
      plink->pbase_ttc->pthermic[SOLID_TTC]->param[RHO_TTC][e]=RHOW_TTC;
      plink->pbase_ttc->pthermic[SOLID_TTC]->param[HEAT_CAP_TTC][e]=HEAT_CAPW_TTC;
      //Creating param_therm_water!!!!!!!!
      TTC_default_pthermic(plink->pbase_ttc,nele);
    }
	
  /*** build_neighbors + build_bounds ***/
  PROSE_fill_neighbor_one_species(pspecies, pchyd, fp);

  /*** build_distances + init_velocity ***/
  PROSE_alloc_u_dist_one_species(pspecies, nele, pchyd, fp);
  /*** build wet surface SW 19/05/2021 ***/
  PROSE_fill_surf_ttc_one_species(pspecies,pchyd,fp);
}

void PROSE_update_heat_transport(s_species_ttc *pspecies, s_chyd *pchyd, double *H_flux_ttc, int nele, double t, FILE *fp)
{
  int r, ne, icard, sub_icard;
  int id_abs_ele;
  int napp;
  double  qapp, varapp, dist, qapp_area, interpol_t;
  s_reach_hyd *preach;
  s_element_hyd *pele;
  s_face_hyd *pface;
  double h_i, h_iter;
  //if(t > 1464.5*86400)
   //printf("debug t = %f\n",t);
  for(r = 0; r < pchyd->counter->nreaches; r++)
    {
      preach = pchyd->p_reach[r];
      for(ne = 0; ne < preach->nele; ne++)
	{
	  pele = preach->p_ele[ne];
	  id_abs_ele = pele->id[ABS_HYD];
	  pface = pele->face[X_HYD][ONE];
	  qapp_area = 0.;
	  for(icard = 0; icard < NB_CARD_TTC; icard++)
	    {
	      switch (icard)
		{
		  /*** flow rate ***/
		  // SW 23/04/2018 EAST_TTC correspond to TWO_HYD, WEST_TTC correspont ONE_HYD pour v >0 aval
		case EAST_TTC :
		  *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][0] = pele->face[X_HYD][TWO]->hydro->Q[T_HYD];//pele->face[X_HYD][TWO]->hydro->Vel; //  5.; //SW 03/05/2018 a verifier vitesse de l'element centre ou au face ?
		  *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][1] = 0; // SW 29/05/2018 non active
		  //LP_printf(fp,"id = %d vam = %f vav = %f\n",id_abs_ele,pele->face[X_HYD][ONE]->hydro->Vel,pele->face[X_HYD][TWO]->hydro->Vel);
		  //printf("%lf %d %d\n",*pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][0],id_abs_ele,icard);
		  break;
		case WEST_TTC :
		  if(ne == 0 && HYD_test_sing_face(pface,UPSTREAM) == CONFLUENCE)
		    {
		      for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++ )
			*pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][sub_icard] = pface->limits[UPSTREAM]->faces[ONE][sub_icard]->hydro->Q[T_HYD];//pface->limits[UPSTREAM]->faces[ONE][sub_icard]->hydro->Vel;// 5.; //
		    }
		  else if((ne == 0) && (HYD_test_sing_face(pface,UPSTREAM) == HYDWORK)) // SW 23/10/2018 element aval du bief et type barrage pour reoxygeneration
		    {
		      *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][0] = pele->face[X_HYD][ONE]->hydro->Q[T_HYD];//pele->face[X_HYD][TWO]->hydro->Vel; //  5.; //SW 03/05/2018 a verifier vitesse de l'element centre ou au face ?
		      *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][1] = 0; // SW 29/05/2018 non active
		    }
		  else
		    {
		      *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][0] = pele->face[X_HYD][ONE]->hydro->Q[T_HYD];//pele->face[X_HYD][ONE]->hydro->Vel; //   5.; //
		      *pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard][1] = 0; // SW 29/05/2018 non active
		      //LP_printf(fp,"id_ele = %d, icard = %d, v = %f\n",id_abs_ele,icard,*pspecies->plink->pu_ttc->puface[id_abs_ele].uface[icard]);		
		    }
		  
		  /*** boundary condition at entrances ***/		  
		  if(ne == 0 && HYD_test_sing_face(pele->face[X_HYD][ONE],UPSTREAM) == DISCHARGE) // apport amont DIRI_FACE_TTC
		    {
		      if(pele->face[X_HYD][ONE]->pt_inflows[0]->temperature != NULL)
			{
			  pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard] = T_0 + TS_function_value_t(t,pele->face[X_HYD][ONE]->pt_inflows[0]->temperature,fp);
			}
		      else
			{
			  // Coller un message d'erreur... avec temp. par défaut
			  pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard] = T_0 + 20.;
			}
		    }
		  break;
		default : // face NORTH_TTC and SOUTH_TTC
		  pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard] = 0.0;
		  for(napp = 0; napp < pele->face[Y_HYD][icard-2]->ninflows; napp++) // on va rester sur la face SUD par raison de commodité
		    { 
		      if(pele->face[Y_HYD][icard-2]->pt_inflows[napp]->temperature != NULL)
			{
			  //printf("apport thermique\n");
			  dist = pspecies->plink->pbase_ttc->param_syst[SIZE_TTC][id_abs_ele];
			  qapp = TS_function_value_t(t,pele->face[Y_HYD][icard-2]->pt_inflows[napp]->discharge,fp);
			  if(qapp > 0.)
			    {
			      pele->face[Y_HYD][icard-2]->pt_inflows[napp]->temperature = TS_function_t_pointed(t,pele->face[Y_HYD][icard-2]->pt_inflows[napp]->temperature,fp); // SW 15/03/2021
                              interpol_t = TS_function_value_t(t,pele->face[Y_HYD][icard-2]->pt_inflows[napp]->temperature,fp);
                              //if((int)interpol_t == CODE_TS) // SW 10/03/2021 CODE_TS = -9999, provided by user in inflows
                              //{
                                 //LP_printf(Simul->poutputs, "default water temperature %f is used for inflow %s \n", Simul->default_t_inflows,
                                   //pele->face[Y_HYD][icard-2]->pt_inflows[napp]->name);
                                // interpol_t = Simul->default_t_inflows;
                              //}
                              if(interpol_t < 0.)
                                  LP_error(fp, "inflow %s temperature < 0., interpol_t = %f\n",pele->face[Y_HYD][icard-2]->pt_inflows[napp]->name, interpol_t);
                              varapp = T_0 + interpol_t;
			      pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard] += varapp*qapp / dist; // NEW_FACE_TTC flux impose
			      qapp_area += qapp;
			    }
			  //else if(qapp < 0.)
			    //{   /* SW 2704/2021 add here prelevement. However if there are some influents and prelevements in one element, it donesn't work anymore*/
                            //    pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard] += qapp / dist;
                            //    qapp_area += qapp;
                            //    LP_printf(fp, "prelevement %s in 1st t = %f qapp = %f\n", pele->face[Y_HYD][icard-2]->pt_inflows[napp]->name, t, qapp);
			   // }
			}
                        else
                         {
                            /*** SW 06/10/2020 if no temperature is provided, we can use a general temperature defined by user in biology section or calculate a
                                temperature using average, delay etc. defined bu user ***/
                          dist = pspecies->plink->pbase_ttc->param_syst[SIZE_TTC][id_abs_ele];
			  qapp = TS_function_value_t(t,pele->face[Y_HYD][icard-2]->pt_inflows[napp]->discharge,fp);
                          // no prelevement
                          if(qapp > 0.)
                          {
                            if(Simul->calc_mode[RIVE] == YES_TS) //SW 26/01/2021
                            {
                                s_fmeteo *ptempe;
                                ptempe = Simul->psimul_bio[0][0]->section->meteo->temperature;
                      
                                if(ptempe->calc == YES_RIVE) // SW here YES_RIVE == 1 in C-RIVE
                                   varapp = T_0;
                              //varapp = T_0 + double_period(Simul->chronos->d_d,ptempe->mean,ptempe->amplitude,ptempe->delay,1); // need to add d_d in chronos
                                else if(ptempe->function_meteo != NULL) // general temperature provided
                                {
                                   ptempe->function_meteo = TS_function_t_pointed(t,ptempe->function_meteo,Simul->poutputs);
	                           varapp = T_0 + TS_function_value_t(t,ptempe->function_meteo,Simul->poutputs);
                                }
                                else
                                  //LP_warning(Simul->poutputs, "no temperature provided for inflows %s. 20 �C used ?!!!\n",pele->face[Y_HYD][icard-2]->pt_inflows[napp]->name);
                                  varapp = T_0 + 20; // SW 10/03/2021 a voir
			    }
                            else
                            {
                                LP_warning(Simul->poutputs, "no temperature provided for inflows %s. 20 �C used ?!!!\n",pele->face[Y_HYD][icard-2]->pt_inflows[napp]->name);  
                                varapp = T_0 + 20;
                                
                            }
                            pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard] += varapp*qapp / dist; // DIRI_FACE_TTC flux impose
			    qapp_area += qapp;                            
                          }
                          
                          /* SW 2704/2021 for prelevement, pele->face[Y_HYD][icard-2]->pt_inflows[napp]->temperature == NULL in input.y : 2188*/
                          else if( qapp < 0.)
                          { 
                             pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard] += qapp * pspecies->plink->pvar_ttc->var[id_abs_ele] / dist; //explicit
                             qapp_area += qapp;
                             //LP_printf(fp, "prelevement %s in 2nd t = %f qapp = %f\n", pele->face[Y_HYD][icard-2]->pt_inflows[napp]->name, t, qapp);
                             
                          }   

                         }
		    }
		  
		  /*** heat exchange budget ***/
		  // In the main.
		  if (icard == NORTH_TTC && Simul->calc_mode[SEB] == YES_TS)
		    {
		      pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard] += H_flux_ttc[id_abs_ele];
                      //if(H_flux_ttc[id_abs_ele] < 0.)
		          //LP_printf(fp,"%lf %d \n",pspecies->plink->pboundary_ttc->pval_bc_ttc->pvalcl[id_abs_ele].valcl[icard],id_abs_ele);
		    }
		  break;
		}    
	    }
	  
	  /*** wet surf ***/
	  //if(ne == 0 && r == 1)
	  //LP_printf(fp,"id_abs_ele == 0\n");
	  h_iter = pele->center->hydro->H[ITER_HYD];
	  h_i = pele->center->hydro->H[T_HYD];
	  
          //pspecies->plink->pparam_calc_ttc->coeff[SURF_T_TTC][id_abs_ele] = TS_function_value_t(h_i,pele->center->hydro->surf,fp);
          //pspecies->plink->pparam_calc_ttc->coeff[SURF_TITER_TTC][id_abs_ele] = TS_function_value_t(h_iter,pele->center->hydro->surf,fp);
          
          /* SW 09/06/2021 using HL as wet surface, rectangular section */
          pspecies->plink->pparam_calc_ttc->coeff[SURF_TITER_TTC][id_abs_ele] = TS_function_value_t(h_iter,pele->center->hydro->width,fp) * h_iter;
          pspecies->plink->pparam_calc_ttc->coeff[SURF_T_TTC][id_abs_ele] = pele->center->hydro->Width * h_i;
	}
    }
}
