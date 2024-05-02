/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: functions_simul_PROSE.h
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

/*in manage_simulation */

s_simul_PROSE *PROSE_init_simulation();
void PROSE_init_output_simulation(int , FILE *);

/*manage_link_prose */
//void PROSE_fill_var(s_species_ttc **, int , int );
//void PROSE_fill_param_base(s_species_ttc ** , int , int );
void PROSE_alloc_u_dist(s_species_ttc ** , int , int , s_chyd *, FILE *);
void PROSE_fill_u(s_species_ttc **, int , int , s_chyd *, FILE *);
void PROSE_fill_neighbor(s_species_ttc **, int , int , s_chyd *, FILE *);
//void PROSE_fill_cl_var_t(s_species_ttc **, int , int , s_chyd *, double , FILE *);
/*manage_print.c*/
void PROSE_print_mats(s_species_ttc *,FILE *);
void print_LHS_RHS(s_gc *);
void PROSE_fill_sparse(s_gc *,double *,	void *);
void PROSE_allocate_spmatrix(int, int, int,FILE *);
void PROSE_allocate_spmatrix_phy(int, int, FILE *);

/*in manage_link_prose.c*/
void PROSE_fill_var_one_species(s_species_ttc *, int , int ,FILE *);
void PROSE_fill_var_all_species(s_species_ttc **, int , int , int, FILE *);
void PROSE_fill_param_base_one_species(s_species_ttc *, int , int,FILE *);
void PROSE_fill_param_base_all_species(s_species_ttc **, int , int , int, FILE *);
void PROSE_fill_u_one_species(s_species_ttc *, s_chyd *, double, double,FILE *);
void PROSE_fill_u_all_species(s_species_ttc **, int , s_chyd *, double, double,FILE *);
void PROSE_alloc_u_dist_one_species(s_species_ttc *, int , s_chyd *, FILE *);
void PROSE_alloc_u_dist_all_species(s_species_ttc **, int , int , s_chyd *, FILE *);
void PROSE_init_neighbor_one_species(s_species_ttc *, int , FILE *);
void PROSE_fill_neighbor_all_species(s_species_ttc **, int , int, s_chyd *, FILE *);
void PROSE_fill_neighbor_one_species(s_species_ttc *, s_chyd *, FILE *);
void PROSE_find_apport_t_all_species(s_species_ttc **, int , s_chyd *, double,int , FILE *);
void PROSE_find_apport_t_one_species(s_species_ttc *, s_chyd *, double, int , FILE *);
void PROSE_update_conc_bio_one_species(s_species_ttc *, int ,int, int,FILE *);
void PROSE_update_conc_bio_all_species(s_species_ttc **, int , int , int,int,FILE *);
void PROSE_fill_var_one_annex_species(s_species_ttc *, int , int , int , int, FILE *);
void PROSE_fill_var_all_annex_species(s_species_ttc ***, int , int, FILE *);
//void PROSE_update_conc_bio_one_annex_species(s_species_ttc *, int , int , int , FILE *);
//void PROSE_update_conc_bio_all_annex_species(s_species_ttc ***, int , FILE *);
void PROSE_fill_param_base_one_annex_species(s_species_ttc *, int , int , int, FILE *);
void PROSE_fill_param_base_all_annex_species(s_species_ttc ***, int, int, FILE *);
void PROSE_fill_u_all_annex_species(s_species_ttc ***, s_chyd *, double,double, FILE *);
void PROSE_alloc_u_dist_all_annex_species(s_species_ttc ***, int , s_chyd *, FILE *fp);
void PROSE_fill_neighbor_all_annex_species(s_species_ttc ***, int , s_chyd *, FILE *);
void PROSE_find_apport_t_one_annex_species(s_species_ttc *, s_chyd *, double , int , int ,FILE *);
void PROSE_find_apport_t_all_annex_species(s_species_ttc ***, s_chyd *, double , FILE *);
void PROSE_fill_surf_ttc_one_species(s_species_ttc *, s_chyd *, FILE *);
void PROSE_fill_surf_ttc_all_species(s_species_ttc **, s_species_ttc ***, int , s_chyd *, FILE *);
void PROSE_set_iapplic_all_species(s_species_ttc **, int ,int, FILE *);
void PROSE_set_iapplic_all_annex_species(s_species_ttc ***, int, FILE *);

void Prose_cal_advflux(s_carac_ttc *,s_param_calc_ttc *,s_species_ttc *,int , double *, FILE *);
void Prose_cal_tot_advflux(s_carac_ttc *,s_param_calc_ttc *,s_species_ttc *,double *, FILE *);
void PROSE_calc_mb_one_species(s_species_ttc *, s_total_mb *, int ,int , int , double , FILE *);
void PROSE_calc_mb_all_species(s_carac_ttc *, int , double , double , int , double *, FILE *);

void PROSE_init_temp_specie(s_carac_ttc *, int);
void PROSE_init_heat_transport(s_species_ttc *, s_chyd *, int, double, FILE *);
void PROSE_update_heat_transport(s_species_ttc *, s_chyd *, double *, int, double, FILE *);

/*in manage_meteo.c*/
double *time_extrema(s_chronos_CHR *);
long safran_cell_count(long, s_rts **);
void PROSE_Patm_reading_in_inputy(int, FILE *, long, long, int, double, double, s_rts **, s_met **);
void PROSE_safransreading_in_inputy(int, FILE *, long, long, int, double, double, s_rts **, s_met **);
int *PROSE_link_icell_to_imet_SEB(int, int, int, s_rts **, s_reach_hyd **, s_met **);
void PROSE_interpolate_multidim_kro(double, s_ft **, double *, int);
void PROSE_manage_metinput_SEB(long, long, double, s_carac_seb **, s_met **, int *);
void PROSE_update_meteo_for_HT(long, long, double, s_chyd *, s_carac_seb **, s_met **, int *, double *);
void PROSE_free_inputs_meteo_all(int , s_carac_seb **);
void PROSE_alloc_inputs_meteo_all(int , s_carac_seb **, int );

/*in manage_file.c*/
char *file_path(char *);

/*in solve_prose.c*/
void PROSE_ttc_all_species(s_carac_ttc *, s_chyd *, double ,double , int, double,double,FILE *);
void PROSE_ttc_one_species(s_carac_ttc *,s_param_calc_ttc * ,s_species_ttc * ,double , double , FILE *);
void PROSE_steady_solve(s_carac_ttc *,s_param_calc_ttc *,s_species_ttc *, int,FILE *);
void PROSE_solve_gc(s_gc *, int );
void PROSE_steady_solve_sparse(s_carac_ttc *, s_param_calc_ttc *, s_species_ttc *, double *, void *,FILE *);
void PROSE_ttc_one_species_sparse(s_carac_ttc *,s_param_calc_ttc *,s_species_ttc *,double ,double ,double *, void *,int,FILE *);
void PROSE_one_section(s_simul *, double , double, int, FILE *);
void PROSE_all_sections(s_simul **, double , double , int , int, FILE *);

/*in manage_bio_prose.c*/
void PROSE_init_param_reaction_state_comp(s_simul **, int , FILE *);
void PROSE_init_O2_conc_sat(s_simul **, double, int, FILE *);
void PROSE_link_hyd_rive_hydro(s_simul **, double , s_chyd *, int,FILE *);
void PROSE_link_hyd_rive_geom(s_simul **, s_chyd *, int, FILE *);
void PROSE_calc_interface_water_sed(s_simul **,s_chyd *, int,FILE *);
void PROSE_update_phy_ctot(s_simul **, int , int, FILE *);
void PROSE_calc_volume_water_all_sections(s_simul **, int , double , int, FILE *);
void PROSE_calc_hyd(double ,s_simul *);

/*calc_outputs.c*/
void PROSE_create_files_mb_bio(int, int, FILE *);
void PROSE_transv_profile_format_bio(s_chyd *,s_output_hyd ***,int, FILE *); // SW 26/04/2018
void PROSE_calc_mb_minit_all_sections(s_simul **, double , double , int , int, FILE *);
void PROSE_calc_mb_sections(double , s_simul **, int, int, FILE *);
void PROSE_calc_mb_mend_all_sections(s_simul **, double , double , int , int, FILE *);
void PROSE_print_variables(FILE *,s_output_hyd *);
void PROSE_print_outputs_formats_bio(s_output_hyd ***,s_inout_set_io *,s_chyd *,int , FILE *);
void PROSE_print_outputs_bio(double ,s_simul **,s_output_hyd ***,s_inout_set_io *,s_chyd *,s_chronos_CHR *,int, FILE *);
void PROSE_transv_profile_bio(double ,s_output_hyd *, s_simul **, int , int , int, FILE *);
void PROSE_print_average_result_bio(s_element_hyd *, s_simul **, s_chyd *, int , double , char *, s_output_hyd *,FILE *);
void PROSE_print_conc_init(int , int, s_inout_set_io *,FILE *);
void PROSE_print_conc_volume_final(int , s_inout_set_io *,FILE *);
void PROSE_print_mb_bio_domaine(s_simul **, double , double , int , int , FILE *);
void PROSE_create_files_weights(int , FILE *);
void PROSE_create_files_parameters(int , FILE *);
void PROSE_print_parameters(s_simul ***, int ,int , double , FILE *);
void PROSE_print_weights(double , int , FILE *);
s_total_mb * Prose_copy_total_mb(s_total_mb *);
void PROSE_print_resampling_size(double , int , FILE *);
void PROSE_print_extract_parameters(s_simul ***, int ,int , double , FILE *);
void PROSE_long_profile(s_simul **, double , s_output_hyd *, int , s_inout_set_io *, s_chyd *, FILE *);
//void PROSE_print_outputs_long(s_simul **,double t,s_output_hyd ***,s_inout_set_io *,s_chyd *,s_chronos_CHR *,FILE *);
void PROSE_create_output_extents_bio(int ,s_lp_pk_hyd *);
void PROSE_print_final_temp(s_output_hyd ***,s_inout_set_io *,s_chyd *,FILE *);
void PROSE_update_outputs_hyd_bio_heat_cur_io(double ,s_output_hyd ***,s_inout_set_io *,s_chyd *,s_chronos_CHR *, int , FILE *);
void PROSE_print_tube_mesh_hyd(double , s_output_tube_type  *, s_chyd *, s_def_tub ****, FILE *);
void PROSE_print_sediment_variables(FILE *,s_output_hyd *);



//calc_base_assimilation.c
void Prose_init_rand_param(s_simul ***, s_carac_assim *, int , FILE *);
int Prose_calc_normalized_weights_and_sample_size(s_carac_assim *, FILE *);
void Prose_fill_re_sample_elim_dupli(s_carac_assim *, FILE *);
void Prose_reinitialization_weights(s_carac_assim *, FILE *);
void Prose_resample_particules(s_carac_assim *, FILE *); 
void Prose_duplications_particules(int , int , FILE *);
void Prose_duplication_conc_one_specie(s_species *, s_species *, FILE *);
void Prose_duplication_conc_one_annnexe_specie(s_annex_var *, s_annex_var *, FILE *);
void Prose_duplication_param_one_specie(s_species *, s_species *, FILE *);
void Prose_duplication_mass_balance_one_specie(s_species *, s_species *, FILE *);
void Prose_duplication_mass_balance_one_annexe_specie(s_annex_var *, s_annex_var *, FILE *);
void Prose_duplication_total_mass_balance( int , int ,int , int , int , FILE *);
void Prose_duplication_total_mass_balance_annexe(int , int ,int , int , int , FILE *);
s_matrix_la* Prose_cov_matrix(s_carac_assim *, FILE *);
s_matrix_la* Prose_inv_cov_matrix(s_matrix_la *, int , FILE *);
double Prose_det_cov_matrix(s_matrix_la *, FILE *);
s_matrix_la* Prose_calc_matrix_produt(s_carac_assim *, s_matrix_la *, int , int , FILE *);
void Prose_calc_weight_all_particules(s_carac_assim *, FILE *);
void Prose_perturbation_parameters(s_simul ***, s_carac_assim *, int , int , FILE *);
void Prose_weights_to_weights_prev(s_carac_assim *, int , FILE *);
void Prose_determine_trophic_state(s_simul ***, s_carac_assim *, FILE *);
void PROSE_extraction_parameters(s_simul ***, int ,int , FILE *);
void Prose_perturb_phy_param(s_simul ***, s_carac_assim *, int , FILE *);
void Prose_perturb_bact_param(s_simul ***, s_carac_assim *, int , FILE *);
void PROSE_assign_epsilon_parameters(s_simul ***, int ,int , int , double , FILE *);
void PROSE_assign_parameters_val(s_simul ***, int ,int , int , double , FILE *);
double PROSE_density_normal(double , double , double , FILE *); // SW 20/05/2022
double PROSE_sum_perturbation_densties_productions(s_carac_assim *, FILE *); // SW 20/05/2022
void PROSE_new_weights_after_perturbation(s_carac_assim *, FILE *); // SW 20/05/2022

//manage_output_heat.c
void PROSE_transv_profile_format_temp(s_chyd *pchyd,s_output_hyd ***p_outputs,  FILE *fp) ;//NF 12/10/2020
void PROSE_print_outputs_heat(double t, s_carac_seb **, s_output_hyd ***, s_inout_set_io *, s_chyd *, s_chronos_CHR *, FILE *);//AB 12/10/2020
void PROSE_long_profile_heat(s_carac_seb **, double, s_output_hyd *, s_inout_set_io *, s_chyd *, FILE *);//AB 12/10/2020
void PROSE_print_average_result_heat(s_element_hyd *, s_carac_seb **, double, s_chyd *, char *, s_output_hyd *, FILE *);//AB 12/10/2020
void PROSE_print_outputs_formats_temp(s_output_hyd ***p_outputs,s_inout_set_io *pinout,s_chyd *pchyd,FILE *fp);//NF 12/10/2020
void PROSE_calc_energy_init_and_flux_from_seb(double, double, s_species_ttc *, s_mbheat_mb **, s_chyd *, s_carac_seb **, FILE *);
void PROSE_create_files_eb_heat(s_mbheat_mb **, FILE *);
void PROSE_calc_one_eb_temperature_and_energy_end(s_species_ttc *, s_mbheat_mb *, int , double , FILE *);
void PROSE_calc_all_ebs_temperature_and_energy_end(s_carac_ttc *, double , double , double *, FILE *fp);
void Prose_print_energy_balance_heat(double , s_mbheat_mb *, FILE *);
void PROSE_calc_energy_init(double , double , s_species_ttc *, s_mbheat_mb **, s_chyd *, FILE *);
void PROSE_calc_flux_from_seb(double , double , s_species_ttc *, s_mbheat_mb **, s_chyd *, s_carac_seb **, FILE *);

//manage_assimilation
void Prose_create_obs_points(s_carac_assim *, s_carac_obs_assim *);
s_carac_obs_assim *Prose_chain_ts_obs(s_carac_obs_assim *,s_carac_obs_assim *);
s_carac_obs_assim *Prose_browse_ts_obs(s_carac_obs_assim *,int );
void Prose_find_obs_id_ele(s_carac_obs_assim *,s_chyd *,FILE *);
int Prose_calc_difference_obs_simul(s_carac_assim *passim, s_simul ***, int , double , FILE *);
void Prose_init_assimilation(s_carac_assim *, int , FILE *);
void PROSE_reset_no_answer_obs(s_carac_assim *, FILE *);
//itos.c*/
char *PROSE_name_param(int );
#ifdef CHR_CHECK
void PROSE_calculate_check_time(int , FILE *);
#endif

//reoxygenation
void PROSE_init_reoxy_species(s_species_ttc **,int ,int ,int ,FILE *);
void Prose_cherche_debit_oxy(s_element_hyd *, s_BC_char_hyd **,int , s_species_ttc *, int , FILE *);
void Prose_calc_rd_rn(s_element_hyd *,double , s_species_ttc *,int , double ,FILE *);


// in manage_stochastic_param_PROSE.c
void PROSE_toc_to_mod_mop_fract(int itype, s_simul_PROSE *Simul, s_inflow_hyd *pinflow, FILE *fp);  //MH 13/09/2021
double PROSE_DA_mod_mop_fract(int itype, int e, int nsub, s_simul_PROSE *Simul, int np, FILE *fp); // MH 05/12/2021
void PROSE_toc_to_mod_mop_fract_var_b1(int itype, s_simul_PROSE *Simul, s_inflow_hyd *pinflow, FILE *fp);  //MH 23/05/2022


// in manage_station_weights.c
s_weight_calc_assim *PROSE_init_calc_stat_weight();
void Prose_create_obs_station_weight(s_carac_assim *, FILE *);

// in randnorm.c
double Prose_normsinv(double );
double Prose_generate_random_norm_param(double );
double Prose_generate_random_param(double , double );

/*in calc_enkf */
s_matrix_la *Prose_calc_obs_error_cov_matrix_enkf(s_carac_assim *, FILE *);
s_matrix_la *Prose_calc_predict_matrix_subtract_mean_enkf(s_carac_assim *, s_simul ***, int, FILE *);
s_matrix_la *Prose_calc_cross_cov_enkf(s_matrix_la *, s_carac_assim *, s_simul ***, int, FILE *);
s_matrix_la *Prose_calc_gain_enkf(s_carac_assim *, s_simul ***, int , FILE *);
void Prose_update_parameter_enkf(s_matrix_la *, s_carac_assim *, s_simul ***, int, FILE *);
void Prose_processus_enkf(s_simul ***, int , s_carac_assim *, int, double ,  FILE *);
s_matrix_la *Prose_cacl_cov_matrix_unknown_true(s_matrix_la *, s_matrix_la *, int , FILE *);
s_matrix_la *Prose_calc_predict_matrix_subtract_obs_enkf(s_carac_assim *, s_simul ***, int , FILE *);
s_matrix_la *Prose_calc_obs_error_cov_matrix_direct_enkf(s_carac_assim *, FILE *);
void Prose_perturbation_and_assign_param_enkf(s_simul ***, int , s_carac_assim *, int , FILE *);


