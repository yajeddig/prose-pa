/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: manage_meteo.c
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
#include <strings.h>
//#include <libprint.h>
#include <libprint.h>
#include <omp.h>
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

long gregorian_calendar_to_jd(int y, int m, int d) // date 01/01/1850 gives "0"
{
  y+=8000;
  if(m<3) { y--; m+=12; }
  return (y*365) + (y/4) - (y/100) + (y/400) - 1200820 + (m*153+3)/5 - 92 + d - 1 - 2396759 ;
}

/* Allows to convert time data read in the main command file into */
/* the time extrema of the simulation. These time extrema will be */
/* essential to the extraction of safran's meteo data.            */

double *time_extrema(s_chronos_CHR *chronos) // in decimal days, must be initialized with "variable[NEXTREMA_TS]".
{
  double *julian_day;
  long jd_ref, jd_beg, jd_end;
  double jd_ref_f, jd_beg_f, jd_end_f; 
  double b, e, bd, ed;
  int d, m, y, h;

  julian_day = (double*)malloc(NEXTREMA_TS*sizeof(double));

  b = chronos->t[BEGINNING_CHR];
  e = chronos->t[END_CHR];

  d = chronos->day_d;
  m = chronos->month + 1;
  y = chronos->year[BEGINNING];
  h = chronos->hour_h;

  jd_ref = gregorian_calendar_to_jd(y, m, d); // Exact julian day
  jd_ref_f = jd_ref + h/24.; // Julian day with decimals

  printf("date_of_reference: %lf\n",jd_ref_f);

  
  bd = b/(24.*3600.); // Supplementary days with decimals from ref. day
  jd_beg_f = jd_ref_f + bd; // Julian beginning day with decimals

  printf("date_ini: %lf\n",jd_beg_f);

  ed = e/(24.*3600.); // Supplementary days with decimals from ref. day
  jd_end_f = jd_ref_f + ed; // Julian ending day with decimals
  
  printf("date_end: %lf\n",jd_end_f);

  julian_day[BEGINNING_TS] = jd_beg_f;
  julian_day[END_TS] = jd_end_f;

  return julian_day;
}

long safran_cell_count(long n_reach, s_rts **p_rts)
{
  long h, i, j;
  long a = 100000, b;
  long counter = n_reach;

  for (i=0; i<n_reach; i++)
    if (p_rts[i]->id_meteo < a)
	a = p_rts[i]->id_meteo;

  for (h=0; h<n_reach; h++)
    {
      j = 0;
      b = 100000;
      
      for (i=0; i<n_reach; i++)
	{
	  if (p_rts[i]->id_meteo > a && p_rts[i]->id_meteo < b)
	    b = p_rts[i]->id_meteo;
	  
	  if (p_rts[i]->id_meteo == a)
	    j = j+1;
	}
      
      if (a == b)
	break;
      
      counter = counter - (j - 1);
      a = b;
    }
  
  return counter;
}

/* Reads and stores formatted data files, from "fp" to "p_met".
   First line: Jour / Heure / Identifiants des mailles Safran (celles que l'on avait dans le fichier des correspondances précédemment)
   Following lines: Value of the Jour / Value of the Heure / Value of the stored meteo parameter for each safran's cell */

/* SW 04/02/2021 quel est le sens d'ajouter t_ini et t_fin qui ne sont pas utilis� dans la fonction ? */
void PROSE_Patm_reading_in_inputy(int check, FILE *fp, long n_reach, long n_safran, int met_param, double t_ini, double t_fin, s_rts **p_rts, s_met **p_met)
{
  char str[1000] = "", strvoid1[255] = "", strvoid2[255] = "";
 
  long i, j, k, l;
  long icheck = 0;
  int endoffile;
  double Patm;
  double dd, hh;
  
  s_ft **a, **b;
  long *V_idmet;
  int n_safran_seine = N_SEINE_SAFRAN; // Useful when n_safran is smaller than the entire count of met cells over the Seine between Ablon and Poses !!!
  int id_meteo, position[n_safran_seine];
  double buffer;
  
  
  /*** Particular case: no time chains for the mean atm pressure P_atm over the safran's cells. ***/
  
  l = 0;
  j = 0;
  
  endoffile=fscanf(fp,"%lf\n", &Patm);
  
  if (check == 0)
    {
      while(endoffile!=EOF)
	{
	  for (i=0; i<n_reach; i++)
	    if(p_rts[i]->id_meteo == l)
	      {
		p_met[j]->id_meteo = l;
		p_met[j]->pmeteo[met_param] = TS_init_ft();
		p_met[j]->pmeteo[met_param]->ft = Patm; // ft does not vary with the time: mean absolute value (P = f(z))
		//printf("%ld %lf\n",p_met[j]->id_meteo,p_met[j]->pmeteo[met_param]->ft);
		j++;
		break;
	      } 
	  endoffile=fscanf(fp,"%lf\n", &Patm);
	  l++;
	}
    }
  else
    {
      while(endoffile!=EOF)
	{
	  for (i=0; i<n_safran; i++)
	    if(p_met[i]->id_meteo == l)
	      {
		p_met[i]->pmeteo[met_param] = TS_init_ft();
		p_met[i]->pmeteo[met_param]->ft = Patm; // ft does not vary with the time: mean absolute value (P = f(z))
		//printf("%ld %lf\n",p_met[i]->id_meteo,p_met[i]->pmeteo[met_param]->ft);
	      } 
	  endoffile=fscanf(fp,"%lf\n", &Patm);
	  l++;
	}
    }
}

void PROSE_safransreading_in_inputy(int check, FILE *fp, long n_reach, long n_safran, int met_param, double t_ini, double t_fin, s_rts **p_rts, s_met **p_met)
{
  char str[1000] = "", strvoid1[255] = "", strvoid2[255] = "";
 
  long i, j, k, l;
  long icheck = 0;
  int endoffile;
  double Patm;
  double dd, hh;
  
  s_ft **a, **b;
  long *V_idmet;
  int n_safran_seine = N_SEINE_SAFRAN; // Useful when n_safran is smaller than the entire count of met cells over the Seine between Ablon and Poses !!!
  int id_meteo, position[n_safran_seine];
  double buffer;

  
  /*** Time chains for the other variables: T air, H sw, H lw, Humidity, Wind velocity. ***/
  
  a = (s_ft **)malloc(sizeof(s_ft*)*n_safran);
  b = (s_ft **)malloc(sizeof(s_ft*)*n_safran);

  
  /*************************************/
  /*** Bloc de lecture de l'en-tête. ***/
  /*************************************/
  
  if (check == 0) // I. If done for the first time ! (->!!! Not to fill in id_meteo twice or with different indices !!!<-)
    {
      k = 0;
      fscanf(fp, "%s %s", strvoid1, strvoid2);
      for (i=0; i<n_safran_seine-1; i++)
	{
	  fscanf(fp, " %d", &id_meteo);
	  if (k<n_safran)
	    {
	      for (j=0; j<n_reach; j++)
		if (p_rts[j]->id_meteo == id_meteo)
		  {
		    p_met[k]->id_meteo = id_meteo;
		    position[k] = i;
		    k++;
		    break;
		  }
	    }
	}
      fscanf(fp, " %d\n", &id_meteo);
      if (k<n_safran)
	for (j=0; j<n_reach; j++)
	  if (p_rts[j]->id_meteo == id_meteo)
	    {
	      p_met[k]->id_meteo = id_meteo;
	      position[k] = i;
	      break;
	    } 
      
      /*for (i=0; i<n_safran; i++)
	printf("%lu %lu\n", i, p_met[i]->id_meteo);*/
    }
  else // II. For the next readings..!
    {
      V_idmet = (long*)malloc(n_safran*sizeof(long));
      fscanf(fp, "%s %s", strvoid1, strvoid2);
      
      k = 0;
      for (i=0; i<n_safran_seine-1; i++)
	{
	  fscanf(fp, " %d", &id_meteo);
	  if (k<n_safran)
	    {
	      for (j=0; j<n_safran; j++)
		if (p_met[j]->id_meteo == id_meteo)
		  {
		    V_idmet[k] = id_meteo;
		    position[k] = i;
		    k++;
		    break;
		  }
	    }
	}
      fscanf(fp, " %d\n", &id_meteo);
      if (k<n_safran)
	for (j=0; j<n_safran; j++)
	  if (p_met[j]->id_meteo == id_meteo)
	    {
	      V_idmet[k] = id_meteo;
	      position[k] = i;
	      break;
	    }
      
      /*for (i=0; i<n_safran; i++)
	printf("V_idmet à la lecture (indice i en premier lieu) : %lu %lu\n", i, V_idmet[i]);*/
    }
  
  
  /*****************************************************************************/
  /*** Bloc d'atteinte de la bonne date. Revient à atteindre la bonne ligne. ***/
  /*****************************************************************************/
  
  for (i=0; i<100000000; i++)
    {
      fscanf(fp, "%lf %lf", &dd, &hh);
      dd = dd + hh/24.;
      if (dd >= t_ini - 1./24.)
	{
	  //printf("Starts on %lf (JD) with",dd);
	  break;
	}
      
      fscanf(fp, "%*[^\n]\n", NULL);
    }
  
  if (check == 0) // I. If done for the first time ! (->!!! Not to fill in id_meteo twice or with different indices !!!<-)
    { 
      /******************************************************************/
      /*** I.1 Bloc de lecture de la première ligne de données météo. ***/
      /******************************************************************/
      
      for (i=0; i<n_safran_seine-1; i++)
	{
	  fscanf(fp, " %lf", &buffer);
	  for (j=0; j<n_safran; j++)
	    if (position[j] == i)
	      {
		p_met[j]->pmeteo[met_param] = TS_init_ft();
		p_met[j]->pmeteo[met_param]->t = dd;
		p_met[j]->pmeteo[met_param]->ft = buffer;
		/*if (j==n_safran-1)
		  printf(" %lf\n", p_met[j]->pmeteo[met_param]->ft);
		else
		printf(" %lf", p_met[j]->pmeteo[met_param]->ft);*/
		
		a[j] = p_met[j]->pmeteo[met_param];
		break;
	      }
	}	  
      //Pour la dernière maille safran sur le Seine, retour à la ligne.
      fscanf(fp, " %lf\n", &buffer);
      for (j=0; j<n_safran; j++)
	if (position[j] == i)
	  {
	    p_met[j]->pmeteo[met_param] = TS_init_ft();
	    p_met[j]->pmeteo[met_param]->t = dd;
	    p_met[j]->pmeteo[met_param]->ft = buffer;
	    //printf(" %lf\n", p_met[j]->pmeteo[met_param]->ft);
	    
	    a[j] = Simul->p_met[j]->pmeteo[met_param];
	    break;
	  }
      
      
      /****************************************************/
      /*** I.2 Bloc de lecture des autres données met'. ***/
      /****************************************************/
      
      for (l=1; l<100000000; l++)
	{
	  fscanf(fp, "%lf %lf", &dd, &hh);
	  dd = dd + hh/24;
	  //printf("%lf", dd);
	  
	  for (i=0; i<n_safran_seine-1; i++)
	    {
	      fscanf(fp, " %lf", &buffer);
	      for (j=0; j<n_safran; j++)
		if (position[j] == i)
		  {
		    b[j] = TS_init_ft();
		    b[j]->t = dd;
		    b[j]->ft = buffer;
		    /*if (j==n_safran-1)
		      printf(" %lf\n", b[j]->ft);
		    else
		    printf(" %lf", b[j]->ft);*/
		    
		    a[j]->next = b[j];
		    b[j]->prev = a[j];
		    a[j] = b[j];
		    break;
		  }
	    }
	  //Pour la dernière maille safran sur la Seine, retour à la ligne.
	  fscanf(fp, " %lf\n", &buffer);
	  for (j=0; j<n_safran; j++)
	    if (position[j] == i)
	      {
		b[j] = TS_init_ft();
		b[j]->t = dd;
		b[j]->ft = buffer;
		//printf(" %lf\n", b[j]->ft);
		
		a[j]->next = b[j];
		b[j]->prev = a[j];
		a[j] = b[j];
		break;
	      }
	  
	  if(dd >= t_fin)
	    {
	      //printf("Ends on %lf (JD)\n",dd);
	      break;
	    }
	}
    }
  else // II. For the next readings..! Cases where id_meteo have already been stored.
    {
      /*******************************************************************/
      /*** II.1 Bloc de lecture de la première ligne de données météo. ***/
      /*******************************************************************/
      
      for (i=0; i<n_safran_seine-1; i++)
	{
	  fscanf(fp, " %lf", &buffer);
	  //printf(" %lf", buffer);
	  for (k=0; k<n_safran; k++)
	    if (position[k] == i)
	      {
		for (j=0; j<n_safran; j++)
		  {
		    if (V_idmet[k] == Simul->p_met[j]->id_meteo)
		      {
			p_met[j]->pmeteo[met_param] = TS_init_ft();
			p_met[j]->pmeteo[met_param]->t = dd;
			p_met[j]->pmeteo[met_param]->ft = buffer;
			/*if (k==n_safran-1)
			  printf(" %lf\n", p_met[j]->pmeteo[met_param]->ft);
			else
			printf(" %lf", p_met[j]->pmeteo[met_param]->ft);*/
			
			a[k] = Simul->p_met[j]->pmeteo[met_param];
			break;
		      }
		  }
		break;
	      }
	}
      //Pour la dernière maille safran sur le Seine, retour à la ligne.
      fscanf(fp, " %lf\n", &buffer);
      //printf(" %lf\n", buffer);
      for (k=0; k<n_safran; k++)
	if (position[k] == i)
	  {
	    for (j=0; j<n_safran; j++)
	      if (V_idmet[k] == Simul->p_met[j]->id_meteo)
		{
		  p_met[j]->pmeteo[met_param] = TS_init_ft();
		  p_met[j]->pmeteo[met_param]->t = dd;
		  p_met[j]->pmeteo[met_param]->ft = buffer;
		  //printf(" %lf\n", Simul->p_met[j]->pmeteo[met_param]->ft);
		  
		  a[k] = Simul->p_met[j]->pmeteo[met_param];
		  break;
		    }
	    break;
	  }
      
      
      /*****************************************************/
      /*** II.2 Bloc de lecture des autres données met'. ***/
      /*****************************************************/
      
      for (l=1; l<100000000; l++)
	{
	  fscanf(fp, "%lf %lf", &dd, &hh);
	  //printf("%lf %lf", dd, hh);
	  dd = dd + hh/24;
	  //printf("%lf", dd);
	  
	  for (j=0; j<n_safran_seine-1; j++)
	    {
	      fscanf(fp, " %lf", &buffer);
	      //printf(" %lf", buffer);
	      for (k=0; k<n_safran; k++)
		if (position[k] == j)
		  {
		    b[k] = TS_init_ft();
		    b[k]->t = dd;
		    b[k]->ft = buffer;
		    /*if (k==n_safran-1)
		      printf(" %lf\n", b[k]->ft);
		    else
		    printf(" %lf", b[k]->ft);*/
		    
		    a[k]->next = b[k];
		    b[k]->prev = a[k];
		    a[k] = b[k];
		    break;
		  }
	    }
	  //Pour la dernière maille safran, retour à la ligne.
	  fscanf(fp, " %lf\n", &buffer);
	  //printf(" %lf\n", buffer);
	  for (k=0; k<n_safran; k++)
	    if (position[k] == j)
	      {
		b[k] = TS_init_ft();
		b[k]->t = dd;
		b[k]->ft = buffer;
		printf(" %lf\n", b[k]->ft);
		
		a[k]->next = b[k];
		b[k]->prev = a[k];
		a[k] = b[k];
		break;
	      }
	  
	  if(dd >= t_fin)
	    {
	      //printf("%lf\n",dd);
	      break;
	    }
	}
      
      free(V_idmet);
    }
}

int *PROSE_link_icell_to_imet_SEB_old_AB(int nele, int nreach, s_rts **p_rts, s_reach_hyd **p_reach)
{
  int *id_meteo;
  int i, j, k, l;
  int nele_reach;
  int nthreads, tid, chunk;
  
  id_meteo = (int *)calloc(nele, sizeof(int));
  
  chunk = CHUNKSIZE;
  /*#pragma omp parallel shared(nele,nreach,p_rts,p_reach,nthreads,chunk) private(i,j,k,l,tid)
    {
    tid = omp_get_thread_num();
    if (tid == 0)
      {
      nthreads = omp_get_num_threads();
      //printf("Number of threads = %d\n", nthreads);
      }
      //printf("Thread %d starting...\n",tid); 
      
    #pragma omp for schedule(dynamic,chunk)*/
  for (i=0; i<nele; i++)
    {
      for (j=0; j<nreach; j++)
	{
	  nele_reach = p_reach[j]->nele;
	  for (k=0; k<nele_reach; k++)
	    {
	      if (i == p_reach[j]->p_ele[k]->id[ABS_HYD])
		{ 
                  //if(i == 1025)
                    // printf("debug\n");
		  for (l=0; l<nreach; l++)
		    {
		      if (strcasecmp(p_reach[j]->limits[UPSTREAM]->name, p_rts[l]->pid_reach->amont)==0 && strcasecmp(p_reach[j]->limits[DOWNSTREAM]->name, p_rts[l]->pid_reach->aval)==0 && p_reach[j]->branch_nb==p_rts[l]->pid_reach->voie)
			{
			  id_meteo[i] = p_rts[l]->id_meteo;
			  //printf("icell(abs) %d imeteo %d tid used: %d\n", i, id_meteo[i],tid);
			  break;
			}
		      else
			if (l==nreach-1)
			  printf("%s %s %d %d\n",p_reach[j]->limits[UPSTREAM]->name,p_reach[j]->limits[DOWNSTREAM]->name,p_reach[j]->branch_nb,i);
		    }
		  break;
		}
	    }
	  if (id_meteo[i] != 0)
	    break;
	}
      // }
    }
  
  for (j=0; j<nreach; j++)
    {
      free(p_rts[j]->pid_reach->amont);
      free(p_rts[j]->pid_reach->aval);
      free(p_rts[j]->pid_reach);
      free(p_rts[j]);
    }
  free(p_rts);
  
  return id_meteo;
}

int *PROSE_link_icell_to_imet_SEB(int nele, int n_meteocell, int nreach, s_rts **p_rts, s_reach_hyd **p_reach, s_met **p_met)
{
  int *id_meteo;
  int i, j, k, l, nmet;
  int nele_reach;
  int nthreads, tid, chunk;
  
  id_meteo = (int *)calloc(nele, sizeof(int));
  
  chunk = CHUNKSIZE;
  /*#pragma omp parallel shared(nele,nreach,p_rts,p_reach,nthreads,chunk) private(i,j,k,l,tid)
    {
    tid = omp_get_thread_num();
    if (tid == 0)
      {
      nthreads = omp_get_num_threads();
      //printf("Number of threads = %d\n", nthreads);
      }
      //printf("Thread %d starting...\n",tid); 
      
    #pragma omp for schedule(dynamic,chunk)*/
  for (i=0; i<nele; i++)
    {
      for (j=0; j<nreach; j++)
	{
	  nele_reach = p_reach[j]->nele;
	  for (k=0; k<nele_reach; k++)
	    {
	      if (i == p_reach[j]->p_ele[k]->id[ABS_HYD])
		{ 
                  //if(i == 1025)
                    // printf("debug\n");
		  for (l=0; l<nreach; l++)
		    {
		      if (strcasecmp(p_reach[j]->limits[UPSTREAM]->name, p_rts[l]->pid_reach->amont)==0 && strcasecmp(p_reach[j]->limits[DOWNSTREAM]->name, p_rts[l]->pid_reach->aval)==0 && p_reach[j]->branch_nb==p_rts[l]->pid_reach->voie)
			{
			  //id_meteo[i] = p_rts[l]->id_meteo;

                          /* SW 14/06/2021 on stocke directement nmet, qui sera utilis� dans PROSE_manage_metinput_SEB */
                          for (nmet=0; nmet<n_meteocell; nmet++)
	                  { 
	                      if (p_rts[l]->id_meteo == p_met[nmet]->id_meteo)
	                      {
	                          id_meteo[i] = nmet;
                                  break;
	                      }
	                  }
                          
			  //printf("icell(abs) %d imeteo %d tid used: %d\n", i, id_meteo[i],tid);
			  break;
			}
		      else
			if (l==nreach-1)
			  printf("%s %s %d %d\n",p_reach[j]->limits[UPSTREAM]->name,p_reach[j]->limits[DOWNSTREAM]->name,p_reach[j]->branch_nb,i);
		    }
		  break;
		}
	    }
	  if (id_meteo[i] != 0)
	    break;
	}
      // }
    }
  
  for (j=0; j<nreach; j++)
    {
      free(p_rts[j]->pid_reach->amont);
      free(p_rts[j]->pid_reach->aval);
      free(p_rts[j]->pid_reach);
      free(p_rts[j]);
    }
  free(p_rts);
  
  return id_meteo;
}

//double *interpolate_multidim_kro(double t, s_ft **p_, int dim) // to be optimised: should be read only once and stored in a big table for each meteo input
void PROSE_interpolate_multidim_kro(double t, s_ft **p_, double *meteo, int dim) // SW 14/06/2021 add meteo as an argument to avoid malloc/free in a time loop
{
  int k;
  //double *rslt;
  //rslt = (double *)malloc(sizeof(double)*dim); // SW 27/01/2021 malloc in time t loop, so a memery leak, attention !!!! Need to be freed later, see void PROSE_free_inputs_meteo_all(int n_cell, s_carac_seb **p_shf, FILE *fp)
  
  for (k=0; k<dim; k++)
    {
      p_[k] = TS_function_t_pointed(t,p_[k],Simul->poutputs);
      //rslt[k] = TS_function_value_t(t,p_[k],Simul->poutputs);
      meteo[k] = TS_function_value_t(t,p_[k],Simul->poutputs); // SW 14/06/2021 add meteo as an argument to avoid malloc/free in a time loop
      //LP_printf(Simul->poutputs, "t = %f,val %ld = %f", t, k, meteo[k]);
    }
    //LP_printf(Simul->poutputs, "\n");
  //return rslt;
}

/* SW 28/01/2021 */
void PROSE_free_inputs_meteo_all(int n_cell, s_carac_seb **p_shf)
{
  int num_cell;
  
  for(num_cell = 0; num_cell < n_cell; num_cell++)
  {
    free(p_shf[num_cell]->pH_inputs->meteo);
    p_shf[num_cell]->pH_inputs->meteo = NULL;
  }
}

/* SW 14/06/2021 */
void PROSE_alloc_inputs_meteo_all(int n_cell, s_carac_seb **p_shf, int dim)
{
  int num_cell;
  double *meteo;

  for(num_cell = 0; num_cell < n_cell; num_cell++)
  {
    p_shf[num_cell]->pH_inputs->meteo = (double *)malloc(sizeof(double)*dim);
    bzero((char *)p_shf[num_cell]->pH_inputs->meteo, sizeof(double));
  }
}


void PROSE_manage_metinput_SEB_old_AB(long n_cell, long n_meteocell, double t, s_carac_seb **p_shf, s_met **p_met, int *id_meteo)
{
  long i, k, v;
  int nthreads, chunk;
  
  nthreads = Simul->psmp->nthreads;
  chunk = PC_set_chunk_size_silent(Simul->poutputs,n_cell,nthreads);
  omp_set_num_threads(nthreads);
  
#pragma omp parallel for schedule(dynamic,chunk) shared(n_cell,n_meteocell) private(k)
  for (i=0; i<n_cell; i++)
    {
      //printf("[%02lu] Task stated with thread %d\n", i, omp_get_thread_num());
      for (k=0; k<n_meteocell; k++)
	{ 
	  //if (k == n_meteocell-1 && id_meteo[i] != p_met[k]->id_meteo)
	    //printf("id_met: %d i_cell : %ld\n",id_meteo[i], i);
	  if (id_meteo[i] == p_met[k]->id_meteo)
	    {
	      //if (k == 0 && i == 4587)
	      //printf("debug\n");
	      //p_shf[i]->pH_inputs->meteo = interpolate_multidim_kro(t, p_met[k]->pmeteo, NMETEO);
              PROSE_interpolate_multidim_kro(t, p_met[k]->pmeteo, p_shf[i]->pH_inputs->meteo, NMETEO); // SW 14/06/2021 add meteo as an argument to avoid malloc/free in time loop
	      break;
	    }
	}
    }

}
  
void PROSE_manage_metinput_SEB(long n_cell, long n_meteocell, double t, s_carac_seb **p_shf, s_met **p_met, int *id_meteo)
{
  int i, l, v;
  int nthreads, chunk;
  
  nthreads = Simul->psmp->nthreads;
  chunk = PC_set_chunk_size_silent(Simul->poutputs,n_cell,nthreads);
  omp_set_num_threads(nthreads);
  
  #pragma omp parallel for schedule(dynamic,chunk) shared(n_cell) private(l)
  for (i=0; i<n_cell; i++)
  {
      l = id_meteo[i];
      //LP_printf(Simul->poutputs, "debug, i = %d l = %d\n", i, l);
      //LP_printf(Simul->poutputs, "debug, leaving %p %p\n",p_met[l]->pmeteo,p_shf[i]->pH_inputs->meteo);
      //PROSE_interpolate_multidim_kro(t, p_met[l]->pmeteo, p_shf[i]->pH_inputs->meteo, NMETEO); // SW 14/06/2021 add meteo as an argument to avoid malloc/free in time loop
      PROSE_interpolate_multidim_kro(t, p_met[l]->pmeteo, p_shf[i]->pH_inputs->meteo, T_SOIL); //SW 13/09/2022 we added T_SOIL in numerate NMETEO, not needed in prose

  }
}

double adapt_surfheatflux_to_bc_for_ttc(double rho, double cp, double d, double surf, double Htot)
{
  //printf("%lf %lf %lf %lf\n",surf,rho,cp,d);
  return surf*Htot/(rho*cp*d);
}

void PROSE_update_meteo_for_HT(long n_cell, long n_meteocell, double t_abs, s_chyd *pchyd, s_carac_seb **p_shf, s_met **p_met, int *id_meteo, double *H_flux_ttc)
{
  int nthreads, chunk;
  int r, ne, i, icard, sub_icard;
  int id_abs_ele;
  s_reach_hyd *preach;
  s_element_hyd *pele;
  s_face_hyd *pface;
  double d, H_tot;
  int num_met;
  FILE *outfile;
  
  /* Calcul du flux brut */
  PROSE_manage_metinput_SEB(n_cell, n_meteocell, t_abs, p_shf, p_met, id_meteo);
  
  nthreads = Simul->psmp->nthreads;
  chunk = PC_set_chunk_size_silent(Simul->poutputs,n_cell,nthreads);
  omp_set_num_threads(nthreads);
  
#pragma omp parallel for schedule(dynamic,chunk) shared(pchyd) private(ne,preach,pele,i,d,H_tot)
  for(r = 0; r < pchyd->counter->nreaches; r++)
    {
      preach = pchyd->p_reach[r];
      for(ne = 0; ne < preach->nele; ne++)
	{
	  pele = preach->p_ele[ne];
	  i = pele->id[ABS_HYD];
	  
	  //d = pele->center->hydro->H[T_HYD]; // mean depth at the cell's centre
           d = pele->center->hydro->Surf / pele->center->hydro->Width; //SW 06/05/2021 
	   H_tot = SEB_H_flux_tot(p_shf[i], i, preach->limits[0]->name); // debug SW 25/01/2021 I think the mean depth is calculated like this. pele->center->hydro->H[T_HYD] is water depth but not mean depth.
	   H_flux_ttc[i] = adapt_surfheatflux_to_bc_for_ttc(Simul->pcarac_heat_ttc->p_species[0]->plink->pbase_ttc->pthermic[SOLID_TTC]->param[RHO_TTC][i], Simul->pcarac_heat_ttc->p_species[0]->plink->pbase_ttc->pthermic[SOLID_TTC]->param[HEAT_CAP_TTC][i], d, pele->center->hydro->Surf, H_tot);
	  /*if (fabs(H_flux_ttc[i]) > 0.1)
	    {
	      outfile = fopen("/home/user/Bureau/WORKING_REP/OUTPUT_PROSE/Out_meteo_test/Temp_check.txt", "a");  // written in the context of a personal test !!! careful !
	      if ( outfile == NULL ) { // error checking with fopen 
		printf("Unable to open file."); 
		exit(1);
	      }
	      //num_met = fprintf(outfile, "ind: %d Temp eau: %lf Temp air: %lf H_tot: %lf, H_atmoflux_ttc: %lf\n", i, p_shf[i]->pH_inputs->Tw, p_shf[i]->pH_inputs->meteo[T_AIR], H_tot, H_flux_ttc[i]);
	      num_met = fprintf(outfile, "ind: %d Temp eau: %lf Temp air: %lf H_tot: %lf, surf: %lf, d: %lf\n", i, p_shf[i]->pH_inputs->Tw, p_shf[i]->pH_inputs->meteo[T_AIR], pele->center->hydro->Surf, d);
	      if( num_met < 0 ) {  // fprintf returns number of values written
		printf("fprintf = %i: fprint error\n",num_met);
		exit(1);
	      }
	      fclose(outfile);
	      }*/
	  //free(p_shf[i]->pH_inputs->meteo);
	}
    }
  //}
}
