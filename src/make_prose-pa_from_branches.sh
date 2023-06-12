#-------------------------------------------------------------------------------
# 
# SOFTWARE NAME: ProSe-PA
# FILE NAME: make_pprose_from_branches.sh
# BRANCH NAME: main
# 
# CONTRIBUTORS: Shuaitao WANG, Lauriane VILMIN, Aurélien BORDET, Masihullah HASANYAR, Thomas ROMARY, Nicolas FLIPO
#
# PROJECT MANAGER: Nicolas FLIPO
# 
# SOFTWARE BRIEF DESCRIPTION: The ProSe-PA is a software for simulating the hydro-biogeochemical functioning of rivers, particularly heavily urbanised rivers, and streams. The sofware can
# operate in two modes: direct calculation or data assimilation.
#
# In direct calculation mode, based on a semi-implicit Eulerian numerical scheme, the software simulates the functioning of the water column in contact with a benthic compartment made up of unconsolidated
# sediments and periphyton (librive library). It can be used to simulate the anthropisation of environments, through the explicit representation of developments such as navigation dams, sluice
# gates and river navigation, as well as discharges into the environment, such as those from wastewater treatment plants or combined sewer overflows.
# The software explicitly simulates the growth of micro-organisms in the water column and in the benthic compartment, enabling the carbon, oxygen and nutrient (nitrogen, phosphrus, silica) cycles
# associated with these biological processes to be quantified (librive library). Water temperature is also simulated by the software (libseb library), as are particulate and dissolved exchanges
# between the water column and the benthic compartment. The software can simulate 1D, pseudo-2D hydraulics of river and streams (discharge, water height) using the libhyd library. The advection-dispersion 
# process is simulated using libttc library.
# 
# In data assimilation mode, ProSe-PA includes two filters for assimilating high frequency dissolved oxygen data. These two filters are a particle filter and the ensemble Kalman filter.
#
# ANSI C software developed at the Geosciences and geoengineering Department, joint research center of Mines Paris-PSL and ARMINES, Fontainebleau, France. The code is based on the coupling 
# of 12 libraries developed also at the Geosciences and geoengineering Department, mostly in ANSI C: libprint, libts, libpc, libchronos, libio, libhyd, libtube, libttc, librive, libseb, libmb, scripts.
#
# CITATION: 
# Wang, S., Flipo, N., Romary, T.. (2019). Oxygen data assimilation for estimating micro-organism communities' parameters in river systems. Water Research, 165, 115021. doi:10.1016/j.watres.2019.115021
# Flipo, N., Even, S., Poulin, M., Tusseau-Vuillemin, M-H., Ameziane, T., Dauta, A. (2004). Biogeochemical modelling at the river scale: plankton and periphyton dynamics (Grand Morin case study, France).
#    Ecol. Model. 176(3-4), 333-347. doi:10.1016/j.ecolmodel.2004.01.012. 
# Even S, Poulin M, Garnier J, Billen G, Servais P, Chesterikoff A (1998). River ecosystem modelling: application of the PROSE model to the Seine River (France). Hydrobiologia, 373, pp. 27-45.
#    doi: 10.1023/A:1017045522336
# Vilmin, L, Aissa, N., Garnier, J., Billen, G., Mouchel, J-M., Poulin, M., Flipo, N. (2015). Impact of hydro-sedimentary processes on the dynamics of soluble reactive phosphorus in the Seine River.
#    Biogeochemistry, 122, 229-251. doi:10.1007/s10533-014-0038-3
# Wang, S., Flipo, N., Romary, T. (2023). Which filter for data assimilation in water quality models? Focus on oxygen reaeration and heterotrophic bacteria activity. Journal of Hydrology, 620, 129423. 
#    doi:10.1016/j.jhydrol.2023.129423
# 
# COPYRIGHT: (c) 2023 Contributors to the ProSe-PA software. 
# CONTACT: Nicolas FLIPO <nicolas.flipo@minesparis.psl.eu>
#          
# 
# All rights reserved. This software and the accompanying materials
# are made available under the terms of the Eclipse Public License v2.0
# which accompanies this distribution, and is available at
# http://www.eclipse.org/legal/epl-v20.html
# 
#------------------------------------------------------------------------------*/

#NF 31/7/2021 updated for branching compilation
chmod +x *.sh
DIR_exe=`pwd`
CODE="all"
BNAME=""

if [ -z "$LIB_HYDROSYSTEM_PATH" ]
then
    LIB_HYDROSYSTEM_PATH=$HOME/Programmes/LIBS/
fi

PATH_INST="$LIB_HYDROSYSTEM_PATH"



if [ $# -eq 0 ] ; then
    ./make_pprose.sh
    exit
elif [ $# -eq 1 ] ; then
    ./make_pprose.sh $1
    exit
else
    if [ $1 != "-b" ] ; then

	echo "unknown option" $1
	exit
    elif [ $# -eq 3 ] ; then
	PATH_INST="$3"
    fi
    BNAME="$2"
fi

if [ ! -d "$PATH_INST" ] ; then 
    mkdir -p $PATH_INST
fi

if [ ! -f "$PATH_INST/delete_links.sh" ] ; then
    echo "No links yet.... Creating them...."
else
    echo "First deleting old links"
    $PATH_INST/delete_links.sh ./
fi

echo "installing scripts in" $PATH_INST
cd $PATH_INST
if [ -d "$PATH_INST/scripts" ] ; then
    rm -rf scripts
fi
svn checkout http://svn.geosciences.fontainebleau.ensmp.fr/repos/scripts

DIR_SHELLS="$PATH_INST/scripts/trunk/install_lib/"
cd $DIR_SHELLS/
chmod 755 *
echo "./create_links.sh" $PATH_INST
./create_links.sh $PATH_INST

cd $DIR_exe
cp Makefile Makefile_tmp
LIST_DEP=`awk -F= -f $PATH_INST/print_dependencies.awk Makefile`
echo "Dependencies of prosepa :" $LIST_DEP

cd $PATH_INST
echo "./clean_install_prosepa2.sh gcc"  $PATH_INST -b $BNAME
./clean_install_prosepa2.sh gcc $PATH_INST -b $BNAME


for i in $LIST_DEP
do
ACR=`$PATH_INST/acronyme.sh $i`
INCL=`$PATH_INST/acronyme.sh $i INCL_`
echo "ACR" $ACR
echo "INCL" $INCL

cd $DIR_exe
VERSION=`$PATH_INST/get_version.sh $i $ACR $PATH_INST $BNAME`
echo "Formatting Makefile_tmp ACR=" $ACR "VERSION=" $VERSION 
awk -F= -f $PATH_INST/write_dep_version.awk -v ACR=$ACR -v VERSION=$VERSION Makefile_tmp
mv awk.out Makefile_tmp

INCL_PATH=`$PATH_INST/get_incl.sh $i $PATH_INST $BNAME`
awk -F= -f $PATH_INST/write_incl.awk -v ACR=$INCL PATH=$INCL_PATH Makefile_tmp
mv awk.out Makefile_tmp
done #for i in $LIST_DEP

awk -F= -f $PATH_INST/modify_lib_path.awk -v path=$PATH_INST Makefile_tmp
mv awk.out Makefile_tmp

awk -F= -f ./bname.awk -v bname=$BNAME Makefile_tmp
mv awk.out Makefile_tmp

cd $DIR_exe
make -f Makefile_tmp clean
make -f Makefile_tmp all

cd $PATH_INST
delete_links.sh $PATH_INST
