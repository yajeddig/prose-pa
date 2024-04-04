/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: calc_enkf.c
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
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <libprint.h>
#ifdef OMP
#include <omp.h>
#endif
#include "time_series.h"
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
//#include "global_PROSE.h"
#include "ext_PROSE.h"

// calculate the error covariance matrix of observation, number of observation site (m) is known
// diagonal matrix, with error N(0, var), mXm, R
s_matrix_la *Prose_calc_obs_error_cov_matrix_enkf(s_carac_assim *passim, FILE *fp)
{

    s_matrix_la *obs_err_cov;
    int nr = 0,nobs;
    double sigma;
    
    obs_err_cov = new_matrix_la();
    obs_err_cov->row_size = obs_err_cov->col_size = passim->num_t_obs;

    obs_err_cov->data = (MATRIX_TYPE **) malloc(obs_err_cov->row_size * sizeof(MATRIX_TYPE *));
    for(nobs = 0; nobs < passim->N_obs; nobs++)
    {
        if(passim->pobs[nobs]->answer_obs == YES_TS)
        {
            obs_err_cov->data[nr] = (MATRIX_TYPE *) malloc(obs_err_cov->col_size * sizeof(MATRIX_TYPE));
            memset(obs_err_cov->data[nr], 0, obs_err_cov->col_size * sizeof(MATRIX_TYPE));
        
            sigma = passim->error_obs_sigma * passim->pobs[nobs]->Obs; // error_obs is proportional to the value of observation
            //LP_printf(fp,"nr = %d s = %f pk = %f obs = %f sigma = %f\n",nr,passim->error_obs_sigma,passim->pobs[nobs]->pk,passim->pobs[nobs]->Obs, sigma);
            obs_err_cov->data[nr][nr] = pow(sigma,2);
            nr++;
        }
    }    

    return obs_err_cov;
}

s_matrix_la *Prose_calc_obs_error_cov_matrix_direct_enkf(s_carac_assim *passim, FILE *fp)
{

    s_matrix_la *obs_err_cov;
    int nr = 0,nobs, nc;
    double sigma;
    
    s_matrix_la *obs_epsilon, *obs_epsilon_transpose;


    //obs_err_cov = new_matrix_la();
    //obs_err_cov->row_size = obs_err_cov->col_size = passim->num_t_obs;

    obs_epsilon = new_matrix_la();
    obs_epsilon->row_size = passim->num_t_obs;
    obs_epsilon->col_size = passim->N_particules;

    //obs_err_cov->data = (MATRIX_TYPE **) malloc(obs_err_cov->row_size * sizeof(MATRIX_TYPE *));
    obs_epsilon->data = (MATRIX_TYPE **) malloc(obs_epsilon->row_size * sizeof(MATRIX_TYPE *));

    for(nobs = 0; nobs < passim->N_obs; nobs++)
    {
        if(passim->pobs[nobs]->answer_obs == YES_TS)
        {
            //obs_err_cov->data[nr] = (MATRIX_TYPE *) malloc(obs_err_cov->col_size * sizeof(MATRIX_TYPE));
            //memset(obs_err_cov->data[nr], 0, obs_err_cov->col_size * sizeof(MATRIX_TYPE));

            obs_epsilon->data[nr] = (MATRIX_TYPE *) malloc(obs_epsilon->col_size * sizeof(MATRIX_TYPE));

            //sigma = passim->error_obs_sigma * passim->pobs[nobs]->Obs; // error_obs is proportional to the value of observation
            //LP_printf(fp,"nr = %d s = %f pk = %f obs = %f sigma = %f\n",nr,passim->error_obs_sigma,passim->pobs[nobs]->pk,passim->pobs[nobs]->Obs, sigma);
            //obs_err_cov->data[nr][nr] = pow(sigma,2);
            
            for(nc = 0; nc < obs_epsilon->col_size; nc++)
                obs_epsilon->data[nr][nc] = passim->pobs[nobs]->epsilon[nc];            
            nr++;
        }
    }    

    if(nr == 1)
        LP_printf(fp, "only one observation\n");

    obs_epsilon_transpose = LA_transpose_matrix(obs_epsilon, fp);
    obs_err_cov = LA_calc_mutiply_matrix(obs_epsilon, obs_epsilon_transpose, fp);
    LA_calc_mutiply_constant(obs_err_cov, 1./(passim->N_particules - 1), fp); //observation covariance matrix

    LA_free_matrix(obs_epsilon, fp);
    LA_free_matrix(obs_epsilon_transpose, fp);

    return obs_err_cov;
}


//calculte prediction covariance matrix (Y-Y*)(Y-Y*)^T/(n-1), observations are already identified
s_matrix_la *Prose_calc_predict_matrix_subtract_obs_enkf(s_carac_assim *passim, s_simul ***psimul_bio, int nparticules, FILE *fp)
{
    s_matrix_la *predict;
    s_carac_obs_assim *pobs;
    int nr, nc,nobs, check_nobs = 0;

    predict = new_matrix_la();
    predict->row_size = passim->num_t_obs;
    predict->col_size = nparticules;

    predict->data = (MATRIX_TYPE **) malloc(predict->row_size * sizeof(MATRIX_TYPE *));
    
    for(nr  = 0; nr < predict->row_size; nr++)
        predict->data[nr] = (MATRIX_TYPE *) malloc(predict->col_size * sizeof(MATRIX_TYPE));
    
    for(nobs = 0; nobs < passim->N_obs; nobs++)
    {
        pobs = passim->pobs[nobs];
        if(pobs->answer_obs == YES_TS)
        {      
            for(nc = 0; nc < predict->col_size; nc++)
                predict->data[check_nobs][nc] = pobs->difference[nc] * -1; // predict Y - Y*
            check_nobs++;
            
        }

    }
     
    if(check_nobs != predict->row_size)
        LP_error(fp,"In FILE %s FUNCTION %s LINE %d, number of observation is not correct, check pobs->answer_obs, check_nobs = %d row_size = %d\n", __FILE__, __FUNCTION__, __LINE__, check_nobs, predict->row_size);

    return predict;
}


//calculte prediction covariance matrix (Y-E(Y))(Y-E(Y))^T/(n-1), observations are already identified
s_matrix_la *Prose_calc_predict_matrix_subtract_mean_enkf(s_carac_assim *passim, s_simul ***psimul_bio, int nparticules, FILE *fp)
{
    s_matrix_la *predict,*predict_tmp;
    
    s_carac_obs_assim *pobs;
    int nr, nc,nobs, check_nobs = 0;
    double nuY, sigmaY;
    //MATRIX_TYPE *mean_by_row;

    predict = new_matrix_la();
    predict->row_size = passim->num_t_obs;
    predict->col_size = nparticules;

    predict->data = (MATRIX_TYPE **) malloc(predict->row_size * sizeof(MATRIX_TYPE *));
    
    for(nr  = 0; nr < predict->row_size; nr++)
        predict->data[nr] = (MATRIX_TYPE *) malloc(predict->col_size * sizeof(MATRIX_TYPE));
    
    for(nobs = 0; nobs < passim->N_obs; nobs++)
    {
        pobs = passim->pobs[nobs];
        if(pobs->answer_obs == YES_TS)
        {      
            for(nc = 0; nc < predict->col_size; nc++)
            {
                //predict->data[check_nobs][nc] = psimul_bio[nc][pobs->id_ele_obs]->section->compartments[WATER][0]->pspecies[pobs->var][pobs->num]->C; // predict Y
                // 19/01/2022 test model error
                //sigmaY = psimul_bio[nc][pobs->id_ele_obs]->section->compartments[WATER][0]->pspecies[pobs->var][pobs->num]->C * 0.1; 
                //nuY = Prose_generate_random_norm_param(sigmaY);
                nuY = 0.;
                predict->data[check_nobs][nc] = psimul_bio[nc][pobs->id_ele_obs]->section->compartments[WATER][0]->pspecies[pobs->var][pobs->num]->C + nuY;
            }
            check_nobs++;
            
        }

    }
     
    if(check_nobs != predict->row_size)
        LP_error(fp,"In FILE %s FUNCTION %s LINE %d, number of observation is not correct, check pobs->answer_obs, check_nobs = %d row_size = %d\n", __FILE__, __FUNCTION__, __LINE__, check_nobs, predict->row_size);

    //mean_by_row = Prose_calc_mean_by_row_matrix(predict, fp);
    predict_tmp = LA_subtract_mean_by_row_matrix(predict, fp); // [Y - E(Y)]
    
    //free(mean_by_row);
    LA_free_matrix(predict, fp);

    return predict_tmp;
}

s_matrix_la *Prose_calc_cross_cov_enkf(s_matrix_la *predict_subtract_mean_trans, s_carac_assim *passim, s_simul ***psimul_bio, int nparticules, FILE *fp)
{
    s_matrix_la *cross_cov_tmp,*cross_cov_subtract_mean, *cross_cov, *cross_cov2;
    
    s_carac_obs_assim *pobs;
    int nparam;
    int nr, nc,nobs, check_nobs = 0;
    MATRIX_TYPE *mean_by_row;

    cross_cov_tmp = new_matrix_la();
    cross_cov_tmp->row_size = passim->num_Of_assimilated_param;
    //LP_printf(fp, "num of da = %d\n", passim->num_Of_assimilated_param);
    cross_cov_tmp->col_size = nparticules;
    
    // filling parameter values
    //for(nparam = 0; nparam < NPARAMDA; nparam++)
    //{
    //    if(passim->param_yesOrno[nr] == YES_TS)
    //        PROSE_extraction_parameters(psimul_bio, nparticules, nparam, fp);
    //}
    
    
    cross_cov_tmp->data = (MATRIX_TYPE **) malloc(cross_cov_tmp->row_size * sizeof(MATRIX_TYPE *));
    nparam = 0;

    for(nr  = 0; nr < NPARAMDA; nr++)
    {
        // SW 25/01/2022 check if parameter is assimilated
        if(passim->param_yesOrno[nr] == YES_TS)
        {
            cross_cov_tmp->data[nparam] = (MATRIX_TYPE *) malloc(cross_cov_tmp->col_size * sizeof(MATRIX_TYPE));
            for(nc = 0; nc < cross_cov_tmp->col_size; nc++)
            {
                cross_cov_tmp->data[nparam][nc] = passim->param[nr][nc];
                
                /* SW 13/01/2022 test with gaussian values varG */
                //cross_cov_tmp->data[nparam][nc] = passim->paramGaussian[nr][nc];
                //if(nr == MU_BACT_DA)
                    //LP_printf(fp, "nparam = %s np = %d val = %.12f\n",PROSE_name_param(nr),nc,passim->param[nr][nc]);
            }
            nparam++;
        }
    }
    
    //mean_by_row = Prose_calc_mean_by_row_matrix(cross_cov_tmp, fp);
    cross_cov_subtract_mean = LA_subtract_mean_by_row_matrix(cross_cov_tmp, fp); // [X - E(X)]
    //LP_printf(fp,"[X - E(X)]:\n");
    //Prose_print_matrix_type(cross_cov_subtract_mean,  MATRIX_A, fp);
    //LP_printf(fp,"[Y - E(Y)]^T:\n");
    //Prose_print_matrix_type(predict_subtract_mean_trans,  MATRIX_A, fp);
    cross_cov = LA_calc_mutiply_matrix(cross_cov_subtract_mean, predict_subtract_mean_trans, fp);
    //LP_printf(fp,"cross_cov:\n");
    
    LA_calc_mutiply_constant(cross_cov, 1./(nparticules - 1), fp); //prediction covariance matrix
    
    //cross_cov2 = Prose_cacl_cov_matrix_unknown_true(cross_cov_subtract_mean, predict_subtract_mean_trans, nparticules, fp);

    //LP_printf(fp,"cross_cov:\n");
    //Prose_print_matrix_type(cross_cov,  MATRIX_A, fp);

    //LP_printf(fp,"cross_cov2:\n");
    //Prose_print_matrix_type(cross_cov2,  MATRIX_A, fp);

    //free(mean_by_row);
    LA_free_matrix(cross_cov_tmp, fp);
    LA_free_matrix(cross_cov_subtract_mean, fp);

    return cross_cov;
}

// calculation of covariance matrix for an ensemble (N), sum([Y - E(Y)][Y - E(Y)]^T)/(N-1) or sum([X - E(X)][Y - E(Y)]^T)/(N-1)
//left matrix dimension is P X N or m X N; right matrix dimension is N X m
s_matrix_la *Prose_cacl_cov_matrix_unknown_true(s_matrix_la *left_matrix, s_matrix_la *right_matrix, int nparticules, FILE *fp)
{

    int nr, nc;
    int num_threads, taille;
    s_matrix_la *temp_l; // P X 1 for 1 member
    s_matrix_la *temp_r; // 1 X m for 1 member
    s_matrix_la *temp_product;
    int nrl;
    
    s_matrix_la *sum;

    temp_l = new_matrix_la();
    temp_l->row_size = left_matrix->row_size;
    temp_l->col_size = 1;

    temp_l->data = (MATRIX_TYPE **) malloc(temp_l->row_size * sizeof(MATRIX_TYPE *));
    for(nr = 0; nr < temp_l->row_size; nr++)
        temp_l->data[nr] = (MATRIX_TYPE *) malloc(temp_l->col_size * sizeof(MATRIX_TYPE));

    temp_r = new_matrix_la();
    temp_r->row_size = 1;
    temp_r->col_size = right_matrix->col_size;

    temp_r->data = (MATRIX_TYPE **) malloc(temp_r->row_size * sizeof(MATRIX_TYPE *));
    for(nr = 0; nr < temp_r->row_size; nr++)
        temp_r->data[nr] = (MATRIX_TYPE *) malloc(temp_r->col_size * sizeof(MATRIX_TYPE));

    sum = new_matrix_la();
    sum->row_size = left_matrix->row_size;
    sum->col_size = right_matrix->col_size;

    sum->data = (MATRIX_TYPE **) malloc(sum->row_size * sizeof(MATRIX_TYPE *));
    for(nr = 0; nr < temp_l->row_size; nr++)
    {
        sum->data[nr] = (MATRIX_TYPE *) malloc(sum->col_size * sizeof(MATRIX_TYPE));
        memset(sum->data[nr], 0, sum->col_size * sizeof(MATRIX_TYPE));
    }



/*#ifdef OMP
          num_threads = Simul->num_threads_par;
          taille = PC_set_chunk_size_silent(fp,left_matrix->col_size,num_threads);
          omp_set_num_threads(num_threads);
          LP_printf(fp,"num threads = %d taille = %d\n",num_threads,taille);
#pragma omp parallel for schedule(dynamic,taille) shared(sum) private(temp_l, temp_r)
#endif*/
    for(nc = 0; nc < nparticules; nc++)
    {
        for(nrl = 0; nrl < temp_l->row_size; nrl++)
            temp_l->data[nrl][0] = left_matrix->data[nrl][nc];
        
        for(nrl = 0; nrl < temp_r->col_size; nrl++)
            temp_r->data[0][nrl] = right_matrix->data[nc][nrl];

        temp_product = LA_calc_mutiply_matrix(temp_l, temp_r, fp);
        //Prose_print_matrix_type(temp_product,  MATRIX_A, fp);
        sum = LA_addition_matrix_replace(sum, temp_product, fp);
        //LP_printf(fp,"nc = %d sum:\n",nc);
        //Prose_print_matrix_type(sum,  MATRIX_A, fp);
        LA_free_matrix(temp_product, fp);
    }

    LA_calc_mutiply_constant(sum, 1./(nparticules-1), fp);
    LA_free_matrix(temp_l, fp);
    LA_free_matrix(temp_r, fp);
    
   return sum;

}

s_matrix_la *Prose_calc_gain_enkf(s_carac_assim *passim, s_simul ***psimul_bio, int nparticules, FILE *fp)
{

    s_matrix_la *predict_subtract_mean, *predict_subtract_mean_trans, *predict_cov, *cross_cov,*obs_err_cov;
    s_matrix_la *sum_right, *gain_kalman;
    
    
    //calculate prediction covariance matrix sum([Y - E(Y)][Y - E(Y)]^T)/(N-1), with N is the size of ensemble
    //LP_printf(fp,"debug3 enkf\n");

    predict_subtract_mean = Prose_calc_predict_matrix_subtract_mean_enkf(passim, psimul_bio, nparticules, fp); // m X N [Y - E(Y)]

    //predict_subtract_mean = Prose_calc_predict_matrix_subtract_obs_enkf(passim, psimul_bio, nparticules, fp); // m X N [Y - Y*] test

    //LP_printf(fp,"debug4 enkf\n");

    predict_subtract_mean_trans = LA_transpose_matrix(predict_subtract_mean, fp); // [Y-E(Y)]^T N X m ou [Y - Y*]^T

    //LP_printf(fp,"debug5 enkf\n");

    predict_cov = LA_calc_mutiply_matrix(predict_subtract_mean, predict_subtract_mean_trans, fp);

    //LP_printf(fp,"predict cov :\n");
    //Prose_print_matrix_type(predict_cov,  MATRIX_A, fp);

    LA_calc_mutiply_constant(predict_cov, 1./(nparticules-1), fp); //prediction covariance matrix size of m X m with m number of observation

    //Prose_print_matrix_type(predict_cov,  MATRIX_A, fp);
    //predict_cov2 = Prose_cacl_cov_matrix_unknown_true(predict_subtract_mean, predict_subtract_mean_trans, nparticules, fp);
    //LP_printf(fp,"predict cov2 :\n");
    //Prose_print_matrix_type(predict_cov2,  MATRIX_A, fp);
    //LP_printf(fp,"debug7 enkf\n");

    //calculate cross covariation matrix [X - E(X)][Y - E(Y)]^T/(N-1), with N is the size of ensemble X = paramGaussain or param
    cross_cov = Prose_calc_cross_cov_enkf(predict_subtract_mean_trans, passim, psimul_bio, nparticules, fp); // size of P X N, with P number of parameters
 
    //LP_printf(fp,"cross_cov :\n");
    //Prose_print_matrix_type(cross_cov,  MATRIX_A, fp);

    //calculate observation error covariance matrix, no perturbation of observation
    //obs_err_cov = Prose_calc_obs_error_cov_matrix_enkf(passim, fp); // size of m X m with m number of observation

    obs_err_cov = Prose_calc_obs_error_cov_matrix_direct_enkf(passim, fp); // size of m X m with m number of observation

    //LP_printf(fp,"observation error cov matrix \n");
    //Prose_print_matrix_type(obs_err_cov,  MATRIX_A, fp);

    //calculate inverse of predict_cov + cross_cov
    sum_right = LA_addition_matrix(predict_cov, obs_err_cov, fp);

    //LP_printf(fp,"debug10 enkf\n");
    //Prose_print_matrix_type(sum_right,  MATRIX_A, fp);

    LA_LU_decomposition_L(sum_right, fp);

    //LP_printf(fp,"Matrix L :\n");
    //Prose_print_matrix_type(sum_right,  MATRIX_L, fp);
    //LP_printf(fp,"Matrix U :\n");
    //Prose_print_matrix_type(sum_right,  MATRIX_U, fp);

    LA_inverse_matrix(sum_right, fp); // size of m X m

    //LP_printf(fp,"debug11 enkf\n");

    //calculate kalman gain
    gain_kalman = LA_calc_mutiply_matrix(cross_cov, sum_right->mat_inv, fp); // size of P X m


    //LP_printf(fp,"debug12 enkf\n");
    LA_free_matrix(predict_subtract_mean, fp);
    LA_free_matrix(predict_subtract_mean_trans, fp);
    LA_free_matrix(predict_cov, fp);
    LA_free_matrix(obs_err_cov, fp);
    LA_free_matrix(cross_cov, fp);
    LA_free_matrix_all(sum_right, fp);

    return gain_kalman;
}

// differences are already calculated
void Prose_update_parameter_enkf(s_matrix_la *gain_kalman, s_carac_assim *passim, s_simul ***psimul_bio, int nele, FILE *fp)
{

    s_matrix_la *difference;
    s_matrix_la *result_upt;
    s_carac_obs_assim *pobs;
    int np, nobs, nr = 0;
    int nparam;
    double val, valG, valP;

    difference = new_matrix_la();
    difference->row_size = passim->num_t_obs;
    difference->col_size = passim->N_particules;
    
    difference->data = (MATRIX_TYPE **) malloc(difference->row_size * sizeof(MATRIX_TYPE *));
    for(nobs = 0; nobs < passim->N_obs; nobs++)
    {
        pobs = passim->pobs[nobs];
        if(pobs->answer_obs == YES_TS)
        {
            difference->data[nr] = (MATRIX_TYPE *) malloc(difference->col_size * sizeof(MATRIX_TYPE));
            // for loop to assign data
            for(np = 0; np < passim->N_particules; np++)
                difference->data[nr][np] = pobs->difference[np]; //(Y* - HY)
            nr++;
            
        }
    }

    result_upt = LA_calc_mutiply_matrix(gain_kalman, difference, fp); // size of P X N 
    nr = 0;

    for(nparam = 0; nparam < NPARAMDA; nparam++)
    {
        // SW 25/01/2022 check if parameter is assimilated
        if(passim->param_yesOrno[nparam] == YES_TS)
        {
            for(np = 0; np < passim->N_particules; np++)
            {
                    
            
                //if(nparam == MU_BACT_DA)
                  //  LP_printf(fp,"np = %d nparam = %s val = %.12f result_upt->data[nr][np] = %.12f difference = %f, passim->param[nparam][np] = %.12f\n", np, PROSE_name_param(nparam), val, result_upt->data[nr][np], difference->data[2][np], passim->param[nparam][np]);             
                
                val = passim->param[nparam][np] + result_upt->data[nr][np]; // X_t+1 = X_t + K(Y* - HY)
	        if(val > passim->param_range[nparam][PARAM_UP]) // out of maximum
		    val = passim->param_range[nparam][PARAM_UP];
	        else if(val < passim->param_range[nparam][PARAM_DOWN]) // out of min
		    val = passim->param_range[nparam][PARAM_DOWN];          

                passim->param[nparam][np] = val;
                // SW 13/01/2022
                //passim->paramGaussian[nparam][np] = valG;
                
	        //PROSE_assign_parameters_val(psimul_bio, nele, np, nparam, val, fp);
            }
            nr++;
        }
    }

    LA_free_matrix(difference, fp);
    LA_free_matrix(result_upt, fp);

}

void Prose_perturbation_and_assign_param_enkf(s_simul ***psimul_bio, int nele, s_carac_assim *passim, int nparam, FILE *fp)
{
    int np;
    double val;
    double variance, a = 0.97, sigma; // SW 27/01/2022
    double h = sqrt(1 - pow(a,2));
    double mean_param, mean_kernel;

    // SW 27/01/2022 read paper Moradkhani200_AWR_EnKF for understanding perturbation a and h
    //mean_kernel = a*param_i + (1 - a)*mean_param; sigma = h*sqrt(variance)
    // param_forcast_i is sampled from N(mean_kernel, sigma)


    // generate forcast parameter for each ensemble member dans assign it to C-RIVE
    for(np = 0; np < passim->N_particules; np++)
    {
        //mean_kernel = a*passim->param[nparam][np] + (1 - a)*mean_param;
        //sigma = h*sqrt(variance);
        //val = mean_kernel + Prose_generate_random_norm_param(sigma);
        
        sigma = (passim->param_range[nparam][PARAM_UP] - passim->param_range[nparam][PARAM_DOWN]) * passim->s_percent[nparam]; // SW 01/06/2023 add nparam for s_percent
        val = passim->param[nparam][np] + Prose_generate_random_norm_param(sigma);

        if(val > passim->param_range[nparam][PARAM_UP]) // out of maximum
            val = passim->param_range[nparam][PARAM_UP];
	else if(val < passim->param_range[nparam][PARAM_DOWN]) // out of min
            val = passim->param_range[nparam][PARAM_DOWN];
 
        passim->param[nparam][np] = val;
        PROSE_assign_parameters_val(psimul_bio, nele, np, nparam, val, fp);
        //if(nparam == MU_BACT_DA)
        //    LP_printf(fp,"in perturbation np = %d nparam = %s val = %.12f passim->param[nparam][np] = %.12f\n", np, PROSE_name_param(nparam), val, passim->param[nparam][np]);
   
    }  
}

void Prose_processus_enkf(s_simul ***psimul_bio, int nele, s_carac_assim *passim, int nparticules, double t,double t_da, FILE *fp)
{
    int answer_obs, nparam, np;
    s_matrix_la *gain_kalman;

    //CHR_begin_timer();
    // caclculate difference obs simul
    if(fabs(t_da - t) < EPS_TS) // SW 04/04/2024 data assimilation time step
        answer_obs = Prose_calc_difference_obs_simul(passim, psimul_bio, nparticules, t, Simul->poutputs);
    
    if(answer_obs == YES_TS)
    {
        //LP_printf(fp,"debug1 enkf\n");
        gain_kalman = Prose_calc_gain_enkf(passim, psimul_bio, nparticules, fp);
        //LP_printf(fp,"debug2 enkf\n");
        Prose_update_parameter_enkf(gain_kalman, passim, psimul_bio, nele, fp);

        // print analyzed parameter values
        //LP_printf(fp, "print parameter at t= %f\n",t);
        for(nparam = 0; nparam < NPARAMDA; nparam++)
	{
            // SW 25/01/2022 check if parameter is assimilated
            if(passim->param_yesOrno[nparam] == YES_TS)
            {
                PROSE_print_parameters(psimul_bio, nparticules, nparam, t, Simul->poutputs);
                Prose_perturbation_and_assign_param_enkf(psimul_bio, nele, passim, nparam, fp);
            }
	}
        // print outputs bio
        //LP_printf(fp, "print concentration at t= %f\n",t);
        for(np = 0; np < nparticules; np++)
        {
	    
            // print concentrations
	    PROSE_print_outputs_bio(t,psimul_bio[np],Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos,np,Simul->poutputs);

            //print mass_balance
	    if(Simul->calc_mode[MB_BIO] == YES_TS)
	    
		    PROSE_print_mb_bio_domaine(psimul_bio[np], t, Simul->chronos->dt, nele, np, Simul->poutputs);				  
	}

        PROSE_reset_no_answer_obs(passim,Simul->poutputs);
        LA_free_matrix(gain_kalman, fp);
        //Simul->clock->time_spent[SOLVE_RIVE] += CHR_end_timer();
    }
    else // no observation at time t
    {
        // print analyzed parameter values
        for(nparam = 0; nparam < NPARAMDA; nparam++)
	{
            // SW 25/01/2022 check if parameter is assimilated
            if(passim->param_yesOrno[nparam] == YES_TS)
                PROSE_print_parameters(psimul_bio, nparticules, nparam, t, Simul->poutputs);
	}

        // print outputs bio
        for(np = 0; np < nparticules; np++)
        {
	    
            // print concentrations
	    PROSE_print_outputs_bio(t,psimul_bio[np],Simul->outputs,Simul->pinout,Simul->pchyd,Simul->chronos,np,Simul->poutputs);

            //print mass_balance
	    if(Simul->calc_mode[MB_BIO] == YES_TS)
		    PROSE_print_mb_bio_domaine(psimul_bio[np], t, Simul->chronos->dt, nele, np, Simul->poutputs);				  
	}
        
    }
}
