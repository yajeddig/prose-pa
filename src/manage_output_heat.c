/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: manage_output_heat.c
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

void PROSE_create_files_eb_heat(s_mbheat_mb **pheatmb, FILE *flog) // SW 04/05/2021
{
  /* Loop indexes */
  int nmb, imb, p;
  s_mbheat_mb *pmb;
  
  nmb = pheatmb[0]->nmb;

  for(imb = 0; imb < nmb; imb++)
  {
      pmb = pheatmb[imb];
      char filename[STRING_LENGTH_MB];
      char pkamont[10];
      char pkaval[10];
      
      sprintf(filename,"%s/energy_balance_heat/",getenv("RESULT"));
      strcat(filename,"mass_balance_heat_");
      strcat(filename,pmb->name);
      strcat(filename,"_");
      strcat(filename,pmb->pk_mb_heat->river);
      sprintf(pkamont, "%4.2f", pmb->pk_mb_heat->pk_up / 1000.);
      sprintf(pkaval, "%4.2f", pmb->pk_mb_heat->pk_down / 1000.);
      strcat(filename,"_");
      strcat(filename,pkamont);
      strcat(filename,"_");
      strcat(filename,pkaval);
      strcat(filename,"\0");
      printf("%s\n",filename);
      if ((pmb->fpmb = fopen(filename,"w")) == NULL) 
          LP_error(flog, "problem when opening the file %s\n",filename);

      fprintf(pmb->fpmb,"Date;t;");

      for(p = FIN_LATERAL_MB; p < NHEAT_TERMS_MB; p++)
          fprintf(pmb->fpmb,"%s;", MB_name_process(p));
      
      fprintf(pmb->fpmb,"\n");
  }
}


char *PROSE_name_T(int param)
{
  char *name;

   switch(param) {
   case TW_IO : {name = strdup("TW");break;}
   case TA_IO : {name = strdup("TA");break;}
   default : {name = strdup("DEFAULT, unknown parameters in PROSE_name_T(int param)");break;}
   }

   return name;
}



void PROSE_print_pl_header_tmp(FILE *fout,s_output_hyd *pout_hyd)
{
  int e, nsub;
  
  for (e = 0; e < NTEMP_IO; e++){
    if (pout_hyd->pout->tempvar[e] == YES_TS){
            fprintf(fout,"%s ",PROSE_name_T(e));
	}
  }
  
  fprintf(fout,"\n");
}


void PROSE_long_profile_heat(s_carac_seb **p_surf_heat_flux, double t, s_output_hyd *result, s_inout_set_io *pinout, s_chyd *pchyd, FILE *fp)
{
  s_lp_pk_hyd *ppk;
  char name[MAXCHAR_PROSE];
  //int pk0_entier,pkf_entier,ligne;
  int k,pkf_atteint,nb;
  //int premier_passage;
  int lx,lx0;
  //double t_entier;
  double pk,pk0,pkf,pklim;
  s_element_hyd *pele;
  s_reach_hyd *preach;
  s_file_io *fic_water;
  //s_file_io *fic_vase;
  int i;
 
  for (i = 0; i < result->npk; i++)
    {
      nb = 0; // SW 25/01/2018 il faut intialiser, sinon la ligne 573 nb++
      //premier_passage = NO_TS;
      ppk = result->lp_pk[i];
      //LP_printf(fp,"reach_nb = %d\n",ppk->reach_nb[0]);
      if (ppk->reach_nb[0] == -1)
	HYD_calculate_output_domain(ppk,LONG_PROFILE,pchyd,fp);
      
      /* Initialisations */
      pk0 = ppk->pk_up;
      pkf = ppk->pk_down;
      preach = NULL;
      
      if (ppk->reach_nb[0] >= 0)
	{
	  preach = pchyd->p_reach[ppk->reach_nb[0]];
	  pklim = preach->limits[DOWNSTREAM]->pk;
	}
    
      pele = ppk->pele;
      //    ligne = nb = 0;
      pkf_atteint = NO_TS;
      k = 1;
      pk = pele->center->pk;
      
      /* Ouverture des fichiers :
       * Plusieurs cas :
       * - Soit un fichier unique qui contient tout 
       * (resultat->fic->adresse != NULL) ;
       * - Soit triés en plusieurs fichiers ...
       * - Fichiers au formats GNUPLOT et fichiers au format OPENDX 
       */
      //t_entier = (int)(t*result->time_unit);
      //pk0_entier = (int)(pk0/1000.0);
      //pkf_entier = (int)(pkf/1000.0);
      
      /* Fichiers au format GNUPLOT, triés par PK et par date */
      /* Ouverture à chaque nouvelle date.
       * Contiendra, à chaque date et pour chaque PK, les informations
       * de positons curviligne (PK) et les valeurs moyennes des variables
       * demandées. */
      //ppk->fic = (s_file_io *)calloc(1,sizeof(s_file_io));
      fic_water = (s_file_io *)calloc(1,sizeof(s_file_io));
      //fic_vase = (s_file_io *)calloc(1,sizeof(s_file_io));
      
      
      sprintf(name,"%s/longitudinal_profiles/lp_temp_t%4.2f_pk%4.2f_pk%4.2f_%s_%d.txt",pinout->name_outputs,t*result->pout->time_unit,pk0/1000.0,pkf/1000.,ppk->river,nb);
      if ((fic_water->address = fopen(name,"w")) == 0) 
	LP_error(fp,"file %s\n",name);	
      fic_water->name = strdup(name);	
      
      fprintf(fic_water->address,"#File generated by %s%4.2f\n",CODE_NAME,NVERSION_PROSE);
      fprintf(fic_water->address,"# Temperature Longitudinal profile (mean transversal values) between PK %4.2f and %4.2f at time %f in %s\n",pk0/1000.0,pkf/1000.0,t*result->pout->time_unit,name_layer(0));
      fprintf(fic_water->address,"#PK POINT(XC YC) "); // SW 04/06/2021 add X Y coordinates 19/10/2022 add POINT for SIG
      
      //HYD_print_variables(ppk->fic->address,result);
      PROSE_print_pl_header_tmp(fic_water->address,result);//NF 12/10/2020
      
      lx = lx0 = 0;
      
      while ((preach != NULL) && (pk < pkf)) {
	
	if ((preach->limits[DOWNSTREAM]->pk < pklim) || (pkf_atteint == YES_TS)) {
	  if (nb == 0)
	    lx0 = lx;
	  nb++;
	  /* Les données de chaque branche sont réparties dans différents 
	   * fichiers. On ouvre donc un nouveau fichier.
	   */
	  fclose(fic_water->address);
	  sprintf(name,"%s/longitudinal_profiles/lp_temp_t%4.2f_pk%4.2f_pk%4.2f_%s_%d.txt",pinout->name_outputs,t*result->pout->time_unit,pk0/1000.,pkf/1000.,ppk->river,nb);
	  if ((fic_water->address = fopen(name,"w")) == 0) 
	    LP_error(fp,"file %s\n",name);	
	  fic_water->name = strdup(name);
	}	  
	
	while (pele != NULL) {
	  /* Pour éviter les retours en arrière */
	  if (pk <= preach->limits[DOWNSTREAM]->pk) {
	    /* Format GNUPLOT */
	    fprintf(fic_water->address,"%f;",pele->center->pk / 1000.); // SW 19/10/2022  add ; for importing points into SIG

            /* SW 04/06/2021 add pk coordinates*/
            fprintf(fic_water->address,"POINT(%f %f);",pele->center->description->xc, pele->center->description->yc);

	    PROSE_print_average_result_heat(pele, p_surf_heat_flux, pele->center->pk, Simul->pchyd, ppk->river, result, fic_water->address); // SW 26/01/2021 add river name
	    lx++;
	  }
	  pk += pele->length;
	  if (pk > pkf) {
	    pkf_atteint = YES_TS;
	    break;
	  }
	  pele = pele->face[X_HYD][TWO]->element[TWO];
	}
	
	if (preach->limits[DOWNSTREAM]->pk >= pklim)
	  pklim = preach->limits[DOWNSTREAM]->pk;
	if ((ppk->reach_nb[k] >= 0) &&(k < pchyd->counter->nreaches)){
	  preach = pchyd->p_reach[ppk->reach_nb[k++]];
	  pele = preach->p_ele[0];
	  pk = preach->limits[UPSTREAM]->pk;
	}
	else
	  preach = NULL;
      }    
      if (nb == 0)
	lx0 = lx;
      
      fclose(fic_water->address);
      
      //if (result->graphics == GNUPLOT)
      //HYD_print_GNUPLOT_lp(result,pinout,ppk,nb,t,fp);
    }
}

void  PROSE_print_average_result_heat(s_element_hyd *pele, s_carac_seb **p_surf_heat_flux, double pk, s_chyd *pchyd, char *river,s_output_hyd *result, FILE *pfic)
{
  int e, yes_no;
  int id_abs, id_abs_n,n;
  int phyfsr;
  double dx_bio;
  double val_e, val_n;
  s_element_hyd *pelen;
  double val;
  
  id_abs = pele->id[ABS_HYD];
  
  /* SW 26/01/2021 check upstream or downstream elements */
  n = HYD_find_river(river,TRANSV_PROFILE,pchyd,Simul->poutputs);
  if((fabs(pk - pchyd->upstr_reach[n]->limits[UPSTREAM]->pk) > EPS_TS) && (fabs(pk - pchyd->downstr_reach[0]->limits[DOWNSTREAM]->pk) > EPS_TS)) { // pk is the upstream pk
    if(pk > pele->center->pk)
    {
      pelen = pele->face[X_HYD][TWO]->element[TWO]; //element aval
      if(pelen != NULL)
	dx_bio = (pk - pele->center->pk)/(pelen->center->pk - pele->center->pk);
      else
	{
	  pelen = pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[TWO][0]->element[TWO];
          dx_bio = (pk - pele->center->pk)/(pelen->center->pk - pele->center->pk);		  
	}
      id_abs_n = pelen->id[ABS_HYD];
    }
    else
    {
      pelen = pele->face[X_HYD][ONE]->element[ONE]; //element amont
      if(pelen != NULL) // SW 18/05/2019 face libre
	dx_bio = (pele->center->pk - pk)/(pele->center->pk - pelen->center->pk);
      else // SW 18/05/2019 face limits
	{
	  pelen = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[ONE][0]->element[ONE]; // element amont
	  dx_bio = (pele->center->pk - pk)/(pele->center->pk - pelen->center->pk);
	}   
      id_abs_n = pelen->id[ABS_HYD];
    }
  
  
  dx_bio = dx_bio > 1.0 ? 1.0 : dx_bio;
  dx_bio = dx_bio < 0.0 ? 0.0 : dx_bio;
  }
  else
    {
       id_abs_n = id_abs; /* SW 26/01/2021 */
       dx_bio = 0.;
    }

  // SW 28/05/2018 pour l'instant one sublayer
  for(e = 0; e < NTEMP_IO; e++){
    if(result->pout->tempvar[e] == YES_TS){
      switch(e) {
      case TW_IO :
	val_e = Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[id_abs] - T_0;//NF 12/10/2020 les sorties ne peuvent être qu'en °C
	val_n = Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[id_abs_n] - T_0;//NF 12/10/2020 Les sorties ne peuvent être qu'en °C
	break;
      case TA_IO :
	val_e = p_surf_heat_flux[id_abs]->pH_inputs->meteo[T_AIR] - T_0;//NF 12/10/2020 Les sorties ne peuvent être qu'en °C
	val_n = p_surf_heat_flux[id_abs_n]->pH_inputs->meteo[T_AIR] - T_0;//NF 12/10/2020 Les sorties ne peuvent être qu'en °C
	break;
      }
      val = (1.0 - dx_bio) * val_e + dx_bio * val_n;
      if(isnan(val))
	fprintf(pfic,"%e;",val_e / result->pout->tempvar_unit[e]);
      else
	fprintf(pfic,"%e;",val / result->pout->tempvar_unit[e]);		 
      yes_no = YES_TS;
      // }
    }
  }
  fprintf(pfic,"\n");
  fflush(pfic);
}

/* Creates the files in which the transv profile (time series at a given location) temperature outputs will be printed */
void PROSE_transv_profile_format_temp(s_chyd *pchyd,s_output_hyd ***p_outputs,  FILE *fp) 
{
  /* Loop indexes */
  int layer,j,i;
  int npk, nts;
  s_ts_pk_hyd *ppk;

    
  nts = pchyd->counter->nts;
  Simul->transv_profil_temp = (FILE ***)calloc(nts,sizeof(FILE **));      	
  for (j = 0; j < nts; j++) {
    npk = p_outputs[TRANSV_PROFILE][j]->npk;  
    Simul->transv_profil_temp[j] = (FILE **)calloc(npk,sizeof(FILE *));			
    for (i = 0; i < npk; i++) {
      ppk = p_outputs[TRANSV_PROFILE][j]->ts_pk[i];
      char filename[MAXCHAR_PROSE];
      sprintf(filename,"%s/time_series/ts_temp_pk%4.2f_%s_branch%d.txt",getenv("RESULT"),ppk->pk / 1000.,ppk->river,ppk->branch_nb);	      
      if ((Simul->transv_profil_temp[j][i] = fopen(filename,"w")) == NULL) 
	LP_error(fp,"problem when opening the file %s\n",filename);
      fprintf(Simul->transv_profil_temp[j][i],"#File generated by %s%4.2f\n",CODE_NAME,NVERSION_PROSE);
      fprintf(Simul->transv_profil_temp[j][i], "#Date day "); /* SW 07/01/2021 print date */
      //fprintf(Simul->transv_profil_temp[j][i], "#time ");
      PROSE_print_pl_header_tmp(Simul->transv_profil_temp[j][i],p_outputs[TRANSV_PROFILE][j]);
    }
  }
  
}


//NF 12/10/2020 printing Temp profiles
void PROSE_print_outputs_formats_temp(s_output_hyd ***p_outputs,s_inout_set_io *pinout,s_chyd *pchyd,FILE *fp)
{ 
  if (pinout->calc[TRANSV_PROFILE] == YES_TS) {
    PROSE_transv_profile_format_temp(pchyd,p_outputs,fp);
  }
  //NF12/10/2020 C'est ici qu'il faudra creer les energy balance
  
  //SW 04/05/2021 add here energy balance
  if(Simul->calc_mode[EB_HEAT] == YES_TS)
       PROSE_create_files_eb_heat(Simul->mb_heat, Simul->poutputs);
}

//NF 12/10/2020
void PROSE_transv_profile_temp(double t,s_output_hyd *result,s_carac_seb **p_surf_heat_flux, int i,  FILE *fp)
{
  s_element_hyd *pele;
  s_ts_pk_hyd *ppk;
  int j;
  //int id_abs;

  for (j = 0; j < result->npk; j++) {
    ppk = result->ts_pk[j];
    pele = ppk->pele;
    if (pele != NULL) {
      /* SW 07/01/2021 print date to ouput file */
      TS_print_date(Simul->chronos->pd[CUR_CHR],Simul->date_format,Simul->transv_profil_temp[i][j], Simul->poutputs);
      fprintf(Simul->transv_profil_temp[i][j], " ");
      TS_print_time(Simul->chronos->pd[CUR_CHR],Simul->transv_profil_temp[i][j], Simul->poutputs); 
      fprintf(Simul->transv_profil_temp[i][j], ";");

      fprintf(Simul->transv_profil_temp[i][j],"%f;",t * result->pout->time_unit);
      PROSE_print_average_result_heat(pele,p_surf_heat_flux, ppk->pk, Simul->pchyd, ppk->river, result,Simul->transv_profil_temp[i][j]);
    }
  }
}


//NF 12/10/2020
void PROSE_print_outputs_heat(double t,s_carac_seb **p_surf_heat_flux,s_output_hyd ***p_outputs,s_inout_set_io *pinout,s_chyd *pchyd,s_chronos_CHR *chronos,FILE *fp)
{
  double t0;
  s_output_hyd *result;
  s_mbheat_mb *pmb;
  int i;

  
  if(Simul->calc_mode[H_T] == YES_TS){
  
    if (pinout->calc[TRANSV_PROFILE] == YES_TS) {
      
      for (i = 0; i < pchyd->counter->nts; i++) {
	result = p_outputs[TRANSV_PROFILE][i];
	t0 = result->pout->t_out[CUR_IO];
	
	if ((t >= result->pout->t_out[BEGINNING]) && (t <= t0) && 
	    (t + chronos->dt > t0) && (t <= result->pout->t_out[END]))
	  {  
	    PROSE_transv_profile_temp(t,result,p_surf_heat_flux,i,fp);
	
	  }
      }
    }
  
    if ((pinout->calc[LONG_PROFILE] == YES_TS))
      {
	for (i = 0; i < pchyd->counter->nlp; i++)
	  {
	    result = p_outputs[LONG_PROFILE][i];
	    t0 = result->pout->t_out[CUR_IO];
	    if ((t >= result->pout->t_out[BEGINNING]) && (t <= t0) && (t + chronos->dt > t0) && (t <= result->pout->t_out[END]))
	      {		
		PROSE_long_profile_heat(p_surf_heat_flux, t, result, pinout, pchyd, fp);
	      }
	  }
      }
    
   //SW 04/05/2021 add here energy balance print
  if (Simul->calc_mode[EB_HEAT] == YES_TS) {
    
    for (i = 0; i < Simul->mb_heat[0]->nmb; i++) {
    pmb = Simul->mb_heat[i];
    if((t >= pmb->t[BEGINNING]) && (t <= pmb->t0) && 
       (t + chronos->dt > pmb->t0)&& (t <= pmb->t[END]))
    {		
	Prose_print_energy_balance_heat(t, pmb, fp);
	pmb->t0 += pmb->ndt;
    }

    }
  }
  }
}

/* SW 05/05/2021 printing of energy balance */
void Prose_print_energy_balance_heat(double t, s_mbheat_mb *pmb, FILE *fp)
{
    int p;

    /* SW 06/07/2023 print date to ouput file */
    TS_print_date(Simul->chronos->pd[CUR_CHR],Simul->date_format,pmb->fpmb, Simul->poutputs);
    fprintf(pmb->fpmb, " ");
    TS_print_time(Simul->chronos->pd[CUR_CHR],pmb->fpmb, Simul->poutputs); 
    fprintf(pmb->fpmb, ";");
    
    fprintf(pmb->fpmb, "%f;",t / pmb->time_unit);
    pmb->mbheat[ERROR_MB] = pmb->mbheat[EEND_MB] - pmb->mbheat[EINIT_MB];

    for(p = FIN_LATERAL_MB; p < EINIT_MB; p++)
    {
        pmb->mbheat[ERROR_MB] -= pmb->mbheat[p]; // unit is J

        /* in Prose-p, standard unit is �E (mmol) for energy, here unit of energy budget is J (in Joule) */
        /* so firstly, Joule is convertied to �E (using J_UE * 0.36), then convertied in output unit asked by user (for example in cal) in print step */

        pmb->mbheat[p] *= J_UE * 0.36;
        fprintf(pmb->fpmb, "%.6f;",pmb->mbheat[p] * pmb->unit_mb_heat); // SW 06/07/2023 replace %e\t by %.6f;
        pmb->mbheat[p] = 0.;
    }
    if(fabs(pmb->mbheat[ERROR_MB]) > 1000)
        LP_warning(Simul->poutputs, "energy balance for %s is not correct, error = %f J. Please check the two pks are located on main river branch.\n", pmb->name, pmb->mbheat[ERROR_MB]);

    for(p = EINIT_MB; p < NHEAT_TERMS_MB; p++)
    {
        /* in Prose-p, standard unit is �E (mmol) for energy, here unit of energy budget is J (in Joule) */
        /* so firstly, Joule is convertied to �E (using J_UE * 0.36), then convertied in output unit asked by user (for example in cal) in print step */

        pmb->mbheat[p] *= J_UE * 0.36;
        fprintf(pmb->fpmb, "%.6f;",pmb->mbheat[p] * pmb->unit_mb_heat); // SW 06/07/2023 replace %e\t by %.6f;
        pmb->mbheat[p] = 0.;
        
    }
    fprintf(pmb->fpmb, "\n");
}
/* Prinitng of the temp_init file, which can be used to initialize upcoming simulations */
void PROSE_print_final_temp(s_output_hyd ***p_outputs, s_inout_set_io *pinout, s_chyd *pchyd, FILE *fp) 
{
  /* Name of the hyd_init file */
  char temp_init_path[MAXCHAR_PROSE];
  /* temp_init file */
  FILE *temp_init;
  /* Element, face, reach_nb */
  int e,f,r;
  
  if (p_outputs[FINAL_STATE][0]->fic == NULL)
    sprintf(temp_init_path,"%s/hyd_init_T_%d",pinout->name_outputs,pchyd->counter->niter);
  else
    sprintf(temp_init_path,"%s/%s_T",pinout->name_outputs,p_outputs[FINAL_STATE][0]->fic->name);
  
  if ((temp_init = fopen(temp_init_path,"w")) == 0)
    LP_error(fp,"Impossible to open Tfinal state file\n");
  
  for (r = 0; r < pchyd->counter->nreaches; r++)
    for (e = 0; e < pchyd->p_reach[r]->nele; e++)
      LP_printf(temp_init,"%d %f\n",pchyd->p_reach[r]->p_ele[e]->id[ABS_HYD],Simul->pcarac_heat_ttc->p_species[0]->plink->pvar_ttc->var[ pchyd->p_reach[r]->p_ele[e]->id[ABS_HYD] ] - T_0);
  
  fclose(temp_init);
}


/*

//NF 11/10/2020 printing Temp profiles taken from libhyd, not the good one, should check the print bio longitudinal profile. If in PROSE, then add directely the counter_bio and the test on Simul->calc_mod[H_T] for printing TW_IO, and calc_MOD[SEB] for TA_IO;

void PROSE_create_GNUPLOT_heat_lp(int i,s_output_hyd ***p_outputs,s_inout_set_io *pinout,FILE *fp)
{
  int j,k;
  s_lp_pk_hyd *ppk;
  char name_gpfile[MAXCHAR_PROSE];
  

  

  for (k = 0; k < NTEMP_IO; k++) {
    
    if (p_outputs[LONG_PROFILE][i]->pout->tempvar[k] == YES_TS) {
      
      for (j = 0; j < p_outputs[LONG_PROFILE][i]->npk; j++) {
	
	ppk = p_outputs[LONG_PROFILE][i]->lp_pk[j];
	
	ppk->gnu_fic[k] = new_file_io();
	
	sprintf(name_gpfile,"%s/gnuplot/lp_pk%4.2f_pk%4.2f_%s_%s.gp",
		pinout->name_outputs,ppk->pk_up / 1000.,ppk->pk_down / 1000.,ppk->river,HYD_name_hyd_var(k));
	ppk->gnu_fic[k]->name = strdup(name_gpfile);

	if ((ppk->gnu_fic[k]->address = fopen(ppk->gnu_fic[k]->name,"w")) == 0)
	  LP_error(fp,"Impossible to open longitudinal profile output file %s\n",ppk->gnu_fic[k]->name);
	
	fprintf(ppk->gnu_fic[k]->address,"#File generated by %s%4.2f\n",CODE_NAME,NVERSION_PROSE);
	fprintf(ppk->gnu_fic[k]->address,"set term x11\n");
	fprintf(ppk->gnu_fic[k]->address,"set autoscale\n");
	fprintf(ppk->gnu_fic[k]->address,"set key below\n");
	fprintf(ppk->gnu_fic[k]->address,"set grid\n");
	fprintf(ppk->gnu_fic[k]->address,"set xlabel \"pk\"\n");
	fprintf(ppk->gnu_fic[k]->address,"set xtics %f,%f,%f\n",ppk->pk_up / 1000.,(ppk->pk_down - ppk->pk_up) / 10000.,ppk->pk_down / 1000.);
	fprintf(ppk->gnu_fic[k]->address,"set xrange [%d:%d]\n\n",(int)(ppk->pk_up / 1000.),(int)(ppk->pk_down / 1000. + 1.));
	
	fclose(ppk->gnu_fic[k]->address);
      }
    }
  }
}


//NF 11/10/2020 printing Temp profiles
void PROSE_heat_long_profile_format(s_chyd *pchyd,s_output_hyd ***p_outputs,s_inout_set_io *pinout,FILE *fp)
{
  int i;

  for (i = 0; i < pchyd->counter->nlp; i++) {
    if (p_outputs[LONG_PROFILE][i]->graphics == GNUPLOT) {

      HYD_create_GNUPLOT_heat_lp(i,p_outputs,pinout,fp);
    }
  }
}

*/


/* SW 04/05/2021 calculation of energy at t-1 (EINIT), and flux SW, LW, LH, SH at t, so called after calculation of libseb. Unit is J */
void PROSE_calc_flux_from_seb(double t, double dt, s_species_ttc *pspecies_temp, s_mbheat_mb **mb_heat, s_chyd *pchyd, s_carac_seb **p_shf, FILE *fp)
{
   int imb, nmb, nr, ne, ns;
   double vol, cp_water, rho_water;
   s_mbheat_mb *pmb;
   s_element_hyd *pele;
   s_reach_hyd *preach;

   for(imb = 0; imb < mb_heat[0]->nmb; imb++)
   {
       pmb = mb_heat[imb];
       if((t >= pmb->t[BEGINNING]) && (t <= pmb->t[END]))
       {
           if(pmb->pk_mb_heat != NULL) // output pk defined
	   {
               nr = 0;
               if(pmb->pk_mb_heat->reach_nb[0] == -1)
	           HYD_calculate_output_domain(pmb->pk_mb_heat,MASS_BALANCE,pchyd,fp);
               preach = NULL;
               
               while((pmb->pk_mb_heat->reach_nb[nr] >= 0) &&(nr < pchyd->counter->nreaches))
               {
                   preach = pchyd->p_reach[pmb->pk_mb_heat->reach_nb[nr]];
                   for(ne = 0; ne < preach->nele; ne++)
                   {
                       pele = preach->p_ele[ne];
                       if ((pele->center->pk >= pmb->pk_mb_heat->pk_up) && (pele->center->pk <= pmb->pk_mb_heat->pk_down))
                       {
                           ns = pele->id[ABS_HYD];
                            
                           pmb->mbheat[SW_MB] += p_shf[ns]->pH_budget[FLUX_SW] * pele->center->hydro->Width * pele->length * dt;
                           pmb->mbheat[LW_MB] += p_shf[ns]->pH_budget[FLUX_LW] * pele->center->hydro->Width * pele->length * dt;
                           pmb->mbheat[LH_MB] += p_shf[ns]->pH_budget[FLUX_LH] * pele->center->hydro->Width * pele->length * dt;
                           pmb->mbheat[SH_MB] += p_shf[ns]->pH_budget[FLUX_SH] * pele->center->hydro->Width * pele->length * dt;
                       }

                   }
                   nr++;
               } //end while
           }
       }
   } //end for
}

/* SW 04/05/2021 calculation of energy at t-1 (EINIT), and flux SW, LW, LH, SH at t, so called after calculation of libseb. Unit is J */
void PROSE_calc_energy_init(double t, double dt, s_species_ttc *pspecies_temp, s_mbheat_mb **mb_heat, s_chyd *pchyd, FILE *fp)
{
   int imb, nmb, nr, ne, ns;
   double vol, cp_water, rho_water;
   s_mbheat_mb *pmb;
   s_element_hyd *pele;
   s_reach_hyd *preach;

   for(imb = 0; imb < mb_heat[0]->nmb; imb++)
   {
       pmb = mb_heat[imb];
       if((t >= pmb->t[BEGINNING]) && (t <= pmb->t[END]))
       {
           if(pmb->pk_mb_heat != NULL) // output pk defined
	   {
               nr = 0;
               if(pmb->pk_mb_heat->reach_nb[0] == -1)
	           HYD_calculate_output_domain(pmb->pk_mb_heat,MASS_BALANCE,pchyd,fp);
               preach = NULL;
               
               while((pmb->pk_mb_heat->reach_nb[nr] >= 0) &&(nr < pchyd->counter->nreaches))
               {
                   preach = pchyd->p_reach[pmb->pk_mb_heat->reach_nb[nr]];
                   for(ne = 0; ne < preach->nele; ne++)
                   {
                       pele = preach->p_ele[ne];
                       if ((pele->center->pk >= pmb->pk_mb_heat->pk_up) && (pele->center->pk <= pmb->pk_mb_heat->pk_down))
                       {
                           ns = pele->id[ABS_HYD];
                           cp_water = pspecies_temp->plink->pbase_ttc->pthermic[SOLID_TTC]->param[HEAT_CAP_TTC][ns]; //4185 J/kg/K
                           rho_water = pspecies_temp->plink->pbase_ttc->pthermic[SOLID_TTC]->param[RHO_TTC][ns]; //1000 kg/m3
                           //vol = pele->length * pele->center->hydro->H[ITER_HYD] * TS_function_value_t(pele->center->hydro->H[ITER_HYD],pele->center->hydro->width,fp); //m3 use volume at t - dt, since var[ns] is before transport here
                           vol = pele->length * pspecies_temp->plink->pparam_calc_ttc->coeff[SURF_T_TTC][ns]; // SW 19/05/2021 SURF_T_TTC has not been updated
                           /* here unit is J */
                           //LP_printf(Simul->poutputs, "ns1 = %d debug\n", ns);
                           pmb->mbheat[EINIT_MB] += cp_water * rho_water * vol * pspecies_temp->plink->pvar_ttc->var[ns]; 
                           
                       }

                   }
                   nr++;
               } //end while
           }
       }
   } //end for
}

/* SW 05/05/2021 calculation of one energy balance (i) for temeprature */
void PROSE_calc_one_eb_temperature_and_energy_end(s_species_ttc *pspecies, s_mbheat_mb *pmb, int i, double dt, FILE *fp)
{
	int ne, e, nsub,icard,sub_icard;
	double mass_entrant = 0.;
	double mass_sortant = 0.;
	double mass_lateral_inflow = 0.;
        double rhowcw, vol, Surf;

	s_lp_pk_hyd *ppk;
	s_element_hyd *pele;
	s_reach_hyd *preach;
	s_chyd *pchyd;	
	int nr,ns,id_neigh;
	int find_upstream_grid, find_downstream_grid;
	int limits_nr, limits_ne;
	find_upstream_grid = NO_TS;
	find_downstream_grid = NO_TS;

 
	/*** SW 06/12/2019 add mass_balance output for the defined domaine by user ***/
	if(pmb->pk_mb_heat != NULL) // output pk defined
	{
	    ppk = pmb->pk_mb_heat;
	    nr = 0;	
	    pchyd = Simul->pchyd;
	    if(ppk->reach_nb[0] == -1)
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

                        /* calculation of end energy after transport and seb, unit is J */
                        rhowcw = pspecies->plink->pbase_ttc->pthermic[WATER_TTC]->param[RHO_TTC][ns] * pspecies->plink->pbase_ttc->pthermic[WATER_TTC]->param[HEAT_CAP_TTC][ns];
                        //vol = pele->length * pele->center->hydro->Surf; //m3 
                        /* SW 19/05/2021 */
                        Surf = pspecies->plink->pparam_calc_ttc->coeff[SURF_T_TTC][ns];
                        vol = pele->length * Surf;
                        pmb->mbheat[EEND_MB] += rhowcw * vol * pspecies->plink->pvar_ttc->var[ns];
                        //LP_printf(Simul->poutputs, "ns2 = %d debug\n", ns); 
                        /* end calculation of end energy after transport and seb */

		        for(icard=NORTH_TTC;icard<NB_CARD_TTC;icard++) 
		        {
			    for(sub_icard = 0; sub_icard < SUB_CARD_TTC; sub_icard++)
			        mass_lateral_inflow +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[icard][sub_icard];
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

                            /* we check if the upstream point is located on a branch */
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
                                        rhowcw = pspecies->plink->pbase_ttc->pthermic[WATER_TTC]->param[RHO_TTC][id_neigh] * pspecies->plink->pbase_ttc->pthermic[WATER_TTC]->param[HEAT_CAP_TTC][id_neigh];
                                        mass_lateral_inflow += rhowcw * pspecies->plink->pvar_ttc->var[id_neigh] * pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][ONE]->element[ONE]->face[X_HYD][TWO]->hydro->Q[T_HYD];
				    }
				}											
				if(pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][TWO]->element[ONE]->center->pk > ppk->pk_up)
				{
				    if(strcmp(pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][TWO]->element[ONE]->reach->river,ppk->river) != 0) // another river simulated, like Marne
				    {
				        id_neigh = pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][TWO]->element[ONE]->id[ABS_HYD];
                                        rhowcw = pspecies->plink->pbase_ttc->pthermic[WATER_TTC]->param[RHO_TTC][id_neigh] * pspecies->plink->pbase_ttc->pthermic[WATER_TTC]->param[HEAT_CAP_TTC][id_neigh];
				        mass_lateral_inflow += rhowcw * pspecies->plink->pvar_ttc->var[id_neigh] * pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][TWO]->element[ONE]->face[X_HYD][TWO]->hydro->Q[T_HYD];
				    }
				}
                                if((pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][ONE]->element[ONE]->center->pk < ppk->pk_up) && (pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][TWO]->element[ONE]->center->pk < ppk->pk_up))
			        {
		                    mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][0]; 
				    mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][1]; 											
				}											
			    }									
			    else if(pele->face[X_HYD][ONE]->limits[UPSTREAM]->faces[X_HYD][ONE]->element[ONE]->center->pk < ppk->pk_up)
			    {
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
				mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][0]; 
				mass_entrant +=  pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[WEST_TTC][1]; 											 
										 
			     }
			}
			if((find_downstream_grid == NO_TS) && (ne == preach->nele-1) && (pele->center->pk <= ppk->pk_down))
			{
                            /* SW 23/06/2021 just one main branch is possible */
                            if(pele->center->pk == ppk->pk_down)
                                mass_sortant += pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[EAST_TTC][0];
                            /* SW 23/06/2021 there is a problem to check downstream point */
                            /*if(pele->face[X_HYD][TWO]->limits[DOWNSTREAM]->faces[X_HYD][ONE]->element[TWO]->center->pk > ppk->pk_down)
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
									 
				mass_sortant += pspecies->plink->pflux_ttc->padv_ttc->padv_face_ttc[ns].face[EAST_TTC][0];	
			    }*/										
			}
		    } // end first if(find_upstream_grid == NO_TS)
		} // end for ne
		nr++;
	    }
	    // sum of all grid in defined domaine
	    pmb->mbheat[FIN_MB] += mass_entrant * dt;  /* here unit is J */
	    pmb->mbheat[FIN_LATERAL_MB] +=  mass_lateral_inflow * dt; /* here unit is J */
	    pmb->mbheat[FOUT_MB] += mass_sortant * dt; /* here unit is J */
	}
	else
            LP_error(fp, "no extent pk defined for energy balance %s\n", pmb->name);
  	
}

/* SW 05/05/2021 calculation of all energy balances, EEND, FIN_LATERAL, FIN, FOUT for temeprature */
void PROSE_calc_all_ebs_temperature_and_energy_end(s_carac_ttc *pcarac_ttc, double t, double dt, double *H_flux_ttc, FILE *fp)
{
	int ns, i, nmb;
	s_species_ttc *pspecies;
	s_param_calc_ttc *pparam_calc_ttc;
        s_mbheat_mb *pmb;	
	
        nmb = Simul->mb_heat[0]->nmb;
        pspecies = pcarac_ttc->p_species[0];
	pparam_calc_ttc = pspecies->plink->pparam_calc_ttc;
	Prose_cal_tot_advflux(pcarac_ttc,pparam_calc_ttc,pspecies,H_flux_ttc, fp);
	
        /***for loop for mutiple mb output***/
	for (i = 0; i < nmb; i++) 
	{
	    pmb = Simul->mb_heat[i];
	    if((t >= pmb->t[BEGINNING]) && (t <= pmb->t[END]))
	        PROSE_calc_one_eb_temperature_and_energy_end(pspecies, pmb, i, dt, fp);
	}	
}

