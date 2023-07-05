<br />
<div align="left">
    <h2 align="left">ProSe-PA</h2>
  <p align="left">
    The ProSe-PA software simulates the hydro-biogeochemical functioning of rivers, particularly heavily urbanised rivers, and streams. The sofware can
    operate in two modes: direct calculation or data assimilation.
    In direct calculation mode, based on a semi-implicit Eulerian numerical scheme, the software simulates the functioning of the water column in contact with a benthic compartment made up of unconsolidated sediments and periphyton (librive library). It can be used to simulate the anthropisation of environments, through the explicit representation of developments such as navigation dams, sluice gates and river navigation, as well as discharges into the environment, such as those from wastewater treatment plants or combined sewer overflows. The software explicitly simulates the growth of micro-organisms in the water column and in the benthic compartment, enabling the carbon, oxygen and nutrient (nitrogen, phosphrus, silica) cycles associated with these biological processes to be quantified (librive library). Water temperature is also simulated by the software (libseb library), as are particulate and dissolved exchanges between the water column and the benthic compartment. The software can simulate 1D, pseudo-2D hydraulics of river and streams (discharge, water height) using the libhyd library. The advection-dispersion  process is simulated using libttc library.
     In data assimilation mode, ProSe-PA includes two filters for assimilating high frequency dissolved oxygen data. These two filters are a particle filter and the ensemble Kalman filter.
    <br />
    <br />
    Software developed at the Centre for geosciences and geoengineering, Mines Paris/ARMINES, PSL University, Fontainebleau, France.
    <br />
    <br />
    <a href="https://gitlab.com/prose-pa/prose-pa-/blob/main/notice/notice_CRIVE_aout2012.pdf"><strong>Explore the user guide »</strong></a>
    <br />    <br />
    <strong>Contributors</strong>
    <br />
     Shuaitao WANG, Lauriane VILMIN, Aur�lien BORDET, Masihullah HASANYAR, Thomas ROMARY, Nicolas FLIPO
    <br />
    <br />
    <strong>Project manager</strong>
    <br />
    Nicolas FLIPO
    <br />
    <br />
    <strong>Contact</strong>
    <br />
    Nicolas FLIPO <a href="mailto:nicolas.flipo@minesparis.psl.eu">nicolas.flipo@minesparis.psl.eu</a>
    <br />
    <br />
    <strong>Citation</strong>
    <br />
    Wang, S., Flipo, N., Romary, T.. (2019). Oxygen data assimilation for estimating micro-organism communities' parameters in river systems. Water Research, 165, 115021. doi:10.1016/j.watres.2019.115021
    Flipo, N., Even, S., Poulin, M., Tusseau-Vuillemin, M-H., Ameziane, T., Dauta, A. (2004). Biogeochemical modelling at the river scale: plankton and periphyton dynamics (Grand Morin case study, France). Ecol. Model. 176(3-4), 333-347. doi:10.1016/j.ecolmodel.2004.01.012. 
     Even S, Poulin M, Garnier J, Billen G, Servais P, Chesterikoff A (1998). River ecosystem modelling: application of the PROSE model to the Seine River (France). Hydrobiologia, 373, pp. 27-45. doi: 10.1023/A:1017045522336
    Vilmin, L, Aissa, N., Garnier, J., Billen, G., Mouchel, J-M., Poulin, M., Flipo, N. (2015). Impact of hydro-sedimentary processes on the dynamics of soluble reactive phosphorus in the Seine River. Biogeochemistry, 122, 229-251. doi:10.1007/s10533-014-0038-3
    Wang, S., Flipo, N., Romary, T. (2023). Which filter for data assimilation in water quality models? Focus on oxygen reaeration and heterotrophic bacteria activity. Journal of Hydrology, 620, 129423. doi:10.1016/j.jhydrol.2023.129423
    </p>
</div>

ProSe-PA couples multiple C-ANSI libraries (also under EPL v2.0 licence), available on GitLab *via* the following groups URLs :
* https://gitlab.com/ghydro
* https://gitlab.com/gtransp
* https://gitlab.com/gutil

## Copyright

[![License](https://img.shields.io/badge/License-EPL_2.0-blue.svg)](https://opensource.org/licenses/EPL-2.0)

&copy; 2023 Contributors to the ProSe-PA software.

*All rights reserved*. This software and the accompanying materials are made available under the terms of the Eclipse Public License (EPL) v2.0 which accompanies this distribution, and is available at http://www.eclipse.org/legal/epl-v20.html.

## How to install ?

**Prerequisite** : For a successfull installation, first, make sure the following packages are installed on your system : **git**, **make**, **flex**, **bison**, **gfortran** and **gcc**.


**If you want to install ProSe-PA from the main branch of ProSe-PA and all librairies, use:**

| Command                        | Installation type                                                                                                                                                                | 
|:--------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------| 
| `./make_prose-pa.sh all`                                  | Re-installs and compiles each ProSe-PA library at `$LIB_HYDROSYSTEM_PATH` location. Compiles prose-pa0.x at the `.` location.              |
| `./make_prose-pa.sh PATH`                                 | Installs and compiles each ProSe-PA library at `$PATH` location. Compiles prose-pa0.x at `.` location.                                     |
| `./make_prose-pa.sh PATH update`                          | Compiles prose-pa0.x at `.` location, assuming that all libraries needed for compiling are located at `$PATH`, and recompiles them.      |
| `./make_prose-pa.sh`                                      | Re-installs and compiles each library at `./.` location. Compiles prose-pa0.x at `.` location. equivalent to   `./make_prose-pa.sh ./`                                        |


**If you want to install CaWaQS from a branch :**

First, please make sure that :
1. the name of the branch is the same in each library. If the branch doesn't exist in the librairy, then the main branch is used as default.
2. the branch does exist in the prose-pa project. To do so, use the command  `git checkout <branch_name> `.

Then, choose your compilation option :

| Command                                                 | Installation type                                                                                                                                                                 | 
|:--------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------| 
| `./make_prose-pa_from_branches.sh -b <branch_name>`       | Re-installs and compiles each library at `$LIB_HYDROSYSTEM_PATH` location, using `<branch_name>` as a branch name if found, `main` otherwise. Compiles prose-pa0.x at `.` location. |
| `./make_prose-pa_from_branches.sh -b <branch_name> PATH`  | Installs each library and recompiles them at `$PATH` location using `<branch_name>` as a branch name if found, `main` otherwise. Compiles prose-pa0.x at `.` location.              |
| `./make_prose-pa_from_branches.sh all`                    | Calls `make_prose-pa.sh all`. Re-installs and compiles each library at `$LIB_HYDROSYSTEM_PATH` location. Compiles prose-pa0.x at `.` location.                                        |
| `./make_prose-pa_from_branches.sh PATH`                   | Calls `make_prose-pa.sh PATH` : Installs each library and recompiles them at `$PATH` location. Compiles prose-pa0.x at `.` location.                                                  |
| `./make_prose-pa_from_branches.sh`                        | Calls `make_prose-pa.sh` : Re-installs and compiles each library at `./.` location. Compiles prose-pa0.x at `.` location.                                                             |


