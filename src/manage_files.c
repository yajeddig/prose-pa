/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: manage_files.c
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
//#include "global_PROSE.h"
#include "ext_PROSE.h"


/*void read_bin(int nj,int din,int df,int type,int njtot,FILE *fplog)
{

  int i,j,it;
  float val;
  s_ft **pft1,**pft2;
  FILE *fp;
  char name[STRING_LENGTH];

  pft1=((s_ft **) malloc(Simul->counter->n[MTO]*sizeof(s_ft *)));
  pft2=((s_ft **) malloc(Simul->counter->n[MTO]*sizeof(s_ft *)));

  for(i=0;i<Simul->counter->n[MTO];i++)
     pft1[i]=p_mto[i]->pkro[type];

  sprintf(name,"%s.%d%d.dir",Simul->settings->pref[type],din,df);
  fp=LP_openr_file(name,Simul->poutputs,ERR_LP);
  fclose(fp);
  fp=fopen(name,"rb");
  if(fp==NULL)
    print_error("Impossible to open file %s\n",name);


  for(j=0;j<nj;j++){
	//	{
    //printf("%dto%d: j %d lu",din,df,i);
    for (it=0;it<Simul->counter->n[MTO];it++)
      {
	fread(&val,sizeof(float),1,fp);
	pft2[it]=create_ft();
	pft2[it]->t=j+njtot;
	pft2[it]->ft=val;

	if(pft1[it]!=NULL)
	  pft1[it]=chain_ft(pft1[it],pft2[it]);
	else
	  pft1[it]=pft2[it];
	fprintf(fplog," %d  = %f",it,val);

      }
//     printf("\n");
  }
  for (it=0;it<Simul->counter->n[MTO];it++)
    p_mto[it]->pkro[type]=pft1[it];

  fclose(fp);
  }*/

/*int lire_one_line(FILE *fp,float *v,int *jj)
{
  int endoffile;
  int i,j;
  int dd,mm,yyyy;

  endoffile=fscanf(fp,"%2d/%2d/%4d\n",&dd,&mm,&yyyy);
  *jj=TS_calculate_jour_julien_j(yyyy,INITIAL_YEAR_JULIAN_DAY_TS,mm-1,dd,Simul->poutputs);//NF 15/5/2015 initialement --mm
//  LP_printf(Simul->poutputs,"%2d/%2d/%4d,%d,",dd,mm,yyyy,*jj);
  for(i=0;i<Simul->counter->nrecpluie;i++){
    endoffile=fscanf(fp," \t%f",&v[i]);
//    LP_printf(Simul->poutputs,"v[%d]=%f,",i,v[i]);
  }
//  LP_printf(Simul->poutputs,"\n");
  endoffile=fscanf(fp,"\n");
  /* for(j=0;j<Simul->counter->nrecpluie;j++)
    printf("v[%d] = %f,",j,v[j]);
    printf("\nfin lire_one_line\n");*/
/*
  return endoffile;

}*/

char *file_path(char *name)
{
  FILE *fp;
  char *new_name;
  int i;

  for (i = 0; i <folder_nb; i++) {
    new_name = (char *)calloc(strlen(name_out_folder[i])+strlen(name)+2,sizeof(char));
    sprintf(new_name,"%s/%s",name_out_folder[i],name);
    fp = fopen(new_name,"r");
    if (fp != NULL)//{
      // fclose(fp);
      break;
    //}
    free(new_name);
  }

  if(i>=folder_nb)
    LP_error(Simul->poutputs,"impossible de localiser %s\n",new_name);

  return new_name;

}

