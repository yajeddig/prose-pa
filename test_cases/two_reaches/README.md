A simple case study to test the ProSe-PA program. 
The case study consists of two reaches and each reach has a length of 30 km. 

The Cmds_files directory:
    1. Here are the command files, the initialization data (init_XXXX), the vitural observation data (o2_obs_XXXX) \
      of dissolved oxygen, the ranges of parameters (param_range), and the maping of reach to SAFRAN meteo data (forcing_to_reach_corresponding.txt).

    2. The command file <steady_case_study.COMM>, runs a steady state simultion. The hydraulics data (Q and Z),\
      the water temperature, the concentrations are calculated. These data are used to initilization the transient\
      state simulation.

    3. The command file <transient_case_study_da.COMM>, runs a transient state simulation with the data assimilation of DO,\
      the calculation of water temperature by libseb using SAFRAN meteo data. The calculated water temperature is used in biogeochemical simulation.
    
    4. !!***!! Before runing the case study, please check if all input folders path are correct in the command files
   
    5. To run the test, ./prose-pa0.7x command_filename log_filename

The bathymetry, Reaches, and Singularities directories: 
    1. Here are the cross sections, reaches, singularities data

The inflows directory: 
    1. Here are the inflows data, UPSTREAM_INFLOW or EFFLUENT (see user guide hydraulics section)

The inflows directory: 
    1. Here are the inflows data, UPSTREAM_INFLOW or EFFLUENT (see user guide hydraulics section)

The param_bio and layers directory:
    1. param_bio: definition of biogeochemical species
    2. layers: properties of simulated layers, and initialization of concentrations (water and sediment)

The DATA_METEO directory:
    1. Here are the hourly SAFRAN meteorological data

