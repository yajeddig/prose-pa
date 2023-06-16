#-------------------------------------------------------------------------------
# 
# SOFTWARE NAME: ProSe-PA
# FILE NAME: Makefile
# BRANCH NAME: main
# 
# CONTRIBUTORS: Shuaitao WANG, Lauriane VILMIN, Aur�lien BORDET, Masihullah HASANYAR, Thomas ROMARY, Nicolas FLIPO
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

PROGRAMME = $(NAME)$(VERSION)
PATH_INST=${LIB_HYDROSYSTEM_PATH}


NAME = prose-pa
SHELL = /bin/bash

VERSION=0.76

SOURCES.c = manage_simulation.c manage_link_prose.c manage_print.c solve_prose.c manage_bio_prose.c\
              calc_outputs.c calc_base_assimilation.c manage_assimilation.c itos.c reoxygenation.c\
		manage_meteo.c manage_files.c manage_output_heat.c manage_matrix.c matrix_operation.c calc_enkf.c \
              manage_stochastic_param_PROSE.c manage_station_weights.c randnum.c

MAIN.c =   main_PROSE.c

SOURCES.l = lexical.l

SOURCES.y = input.y

SOURCES.h =   

SOURCES.f = 

CC=gcc
#CC=g++
#CC=/usr/local/gcc-4.7.2/bin/g++ #for numerix which has gcc3.3 by default
#CC=${MPI_DIR}/bin/mpicxx
#CC=${MPI_DIR}/bin/mpicc
#CC=mpicc
#CC=mpicxx
FC=gfortran

# VERSION of GCC is needed to compile 
GCC_VER=$(shell $(CC) -v 2>&1|tail -1|awk '{print $$3}'|awk -F "." '{print $$1$$2$$3}')


COMPILO=$(shell echo $(CC) | awk -F "/" '{print $$NF}')

YACC= bison -y 
LEX= flex
YFLAGS= -dv 
LFLAGS= -v

# Configure dependencies
# libprint and libts
# Versions
V_PC=0.05
V_LP=1.24
V_TS=1.74
V_GC=0.10
V_IO=0.12
V_CHR=0.13
#V_SPA=0.02
V_HYD=0.34
V_TTC=0.09
V_RIVE=0.27
V_SEB=0.06
V_TUB=0.06
V_MB=0.01


# Paths for Include
INCL_GC=-I$(PATH_INST)/libgc/src
INCL_PC=-I$(PATH_INST)/libpc/src
INCL_TS=-I$(PATH_INST)/libts/src
INCL_LP=-I$(PATH_INST)/libprint/src
INCL_IO=-I$(PATH_INST)/libio/src/
INCL_CHR=-I$(PATH_INST)/libchronos/src/
#INCL_SPA=-I$(PATH_INST)/libspa/src/
INCL_HYD=-I$(PATH_INST)/libhyd/src/
INCL_TTC=-I$(PATH_INST)/libttc/src
INCL_RIVE=-I$(PATH_INST)/c-rive/src
INCL_SEB=-I$(PATH_INST)/libseb/src
INCL_TUB=-I$(PATH_INST)/libtube/src
INCL_PROSE=-I$(HOME)/Programmes/prose/src
INCL_SPASR=-I$(PATH_INST)/libgc/src/sparse_11_7_2011
INCL_GSL=-I/usr/$(PREFIX_GSL)/include
INCL_MB=-I$(PATH_INST)/libmb/src
# Attention, l'ordre est important!

LIBS= -L$(PATH_INST) $(PATH_INST)/libgc/src/sparse_11_7_2011/sparse.a $(LIBGSL)  -lrive$(V_RIVE)_$(COMPILO) -lttc$(V_TTC)_$(COMPILO) -lseb$(V_SEB)_$(COMPILO) -lhyd$(V_HYD)_$(COMPILO) -lmb$(V_MB)_$(COMPILO) -ltube$(V_TUB)_$(COMPILO) -lchronos$(V_CHR)_$(COMPILO) -lio$(V_IO)_$(COMPILO)  -lgc$(V_GC)_$(COMPILO) -lts$(V_TS)_$(COMPILO)  -lpc$(V_PC)_$(COMPILO)   -lprint$(V_LP)_$(COMPILO) -lc -lm

#LIBS= -L$(PATH_INST) $(PATH_INST)/libgc/src/sparse_11_7_2011/sparse.a  -lrive$(V_RIVE)_$(COMPILO) -lttc$(V_TTC)_$(COMPILO) -lhyd$(V_HYD)_$(COMPILO) -lchronos$(V_CHR)_$(COMPILO) -lio$(V_IO)_$(COMPILO)  -lgc$(V_GC)_$(COMPILO) -lts$(V_TS)_$(COMPILO)  -lpc$(V_PC)_$(COMPILO)   -lprint$(V_LP)_$(COMPILO) -lc -lm

OPTCOMP = -g
#OPTW= -Wall
OPTW= -Wno-write-strings
#OPTD = -DOMP -fopenmp -DCHR_CHECK
#OPTD = -DTS_DEBUG
#OPTD = -DSO_OLD

# This option indicates the version of GCC to the pre-compiler to properly compile the lexical.l
CFLAGSD = -DGCC$(GCC_VER) #Mettre en commentaire pour certaine version de lex
OPTOMP= -fopenmp
#OPTD = -DOMP $(OPTOMP) -DCOUPLED_RIVE
OPTD = -DOMP $(OPTOMP) -DCOUPLED_RIVE -DCDA
CFLAGS_GLOBAUX= $(OPTCOMP) $(OPTW) 
CFLAGS= $(CFLAGS_GLOBAUX) $(OPTD)  $(CFLAGSD) 


#INCLUDES = -I/usr/local/gcc-4.7.2/include -I/usr/local/gcc-4.7.2/lib/gcc/x86_64-unknown-linux-gnu/4.7.2/include/#for system using gcc older than 4.x
#INCLUDES = -I/usr/local/gcc-4.7.2/include -I/usr/local/gcc-4.7.2/lib/gcc/x86_64-unknown-linux-gnu/4.7.2/include/#for system using gcc older than 4.x #for numerix with petsc

#CPP : Preprocessor, gere les inclure.
#C'est ici qu'il faut indiquer les chemins d'acces aux .h
CPPFLAGS=  $(INCL_GC) $(INCL_LP) $(INCL_TS) $(INCL_PC) $(INCL_IO) $(INCL_CHR) $(INCL_HYD) $(INCL_TTC) $(INCL_RIVE) $(INCL_SEB) $(INCL_TUB) $(INCL_MB) $(INCL_SPASR) $(INCL_PROSE)
#CPPFLAGS=  $(INCL_GC) $(INCL_LP) $(INCL_TS) $(INCL_PC) $(INCL_IO) $(INCL_CHR) $(INCL_HYD) $(INCL_TTC) $(INCL_RIVE) $(INCL_SPASR) $(INCL_PROSE)


# Variables d��duites des pr��c��dentes:
OBJ1=$(SOURCES.y:.y=.o) $(SOURCES.l:.l=.o)
OBJ2=$(SOURCES.c:.c=.o)
MAIN=$(MAIN.c:.c=.o)
MAINFLU=$(MAINFLU.c:.c=.o)
OBJ4=$(LIBFLU.c:.c=.o)
OBJ =  $(OBJ1) $(OBJ2) $(MAIN)
OBJ_LIB = $(OBJ1) $(OBJ2) $(OBJ4)

OBJFLU =  $(OBJ1) $(OBJ2) $(MAINFLU)

MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
MKFILE_DIR := $(dir $(MKFILE_PATH))


#$(PROGRAMME) :  $(OBJ) 
#	$(LINK.c) -o $(PROGRAMME)  $(OBJ) $(LIBS)

#$(PROGRAMME) :  $(OBJ) chkopts  
#	@echo 'Linking the .o'
#	-$(CLINKER) $(CFLAGS) -o $(PROGRAMME) $(OBJ) ${PETSC_KSP_LIB} $(LIBS) 
#	-$(CLINKER) $(LIBS) -o $(PROGRAMME) $(OBJ) ${PETSC_KSP_LIB}

$(PROGRAMME) :  $(OBJ)  
	#$(CC) $(CFLAGS) $(CPPFLAGS)  -o $(PROGRAMME)  $(OBJ)  $(LIBS) 
	$(FC) $(CFLAGS) -o $(PROGRAMME) $(OBJ)  $(LIBS) 

# Pour nettoyer quand on change de machine:
clean:
	@echo 'Cleaning old .o'
	rm -f $(OBJ) $(LINTFILES) $(PROGRAMME) core y.* apoub*

default :$(PROGRAMME)

all : $(PROGRAMME)

# For compiling the program and remove the source files
allpub: $(PROGRAMME)
	rm *.c *.o *.y *.l *.h Makefile Makefile_tmp
	sudo rm /usr/local/bin/prose 
	sudo ln -s $(MKFILE_DIR)/$(NAME)$(VERSION) /usr/local/bin/prose
	rm -rf $(MKFILE_DIR)/../../branches

llib :   $(OBJ_LIB) chkopts  
	ar rv libHyCu$(VERSION).a $(OBJ_LIB) 
#	ar -q $(HYCU_DIR)/libFakeHyCu.a $(HYCU_DIR)/libFakeHyCu.o

llibpub :   $(OBJ_LIB) chkopts  
	ar rv libHyCu$(VERSION).a $(OBJ_LIB) 
	rm *.c *.o *.y *.l *.h Makefile

# Pas s��r que ��a marche toujours pour Flumy depuis qu'il n'y a plus PetsC

flutest :   $(OBJFLU)  
	$(CC) $(CFLAGS) $(CPPFLAGS)  -o $(PROGRAMME)-test  $(OBJFLU)  $(LIBS)
	@echo 'To check HyCu try > make test'

test : flutest
	./$(PROGRAMME)-test

# Suppression du $(RM) par d��faut :
.y.c :
	$(YACC.y) $< -o $*.c
#.y.o :
#	$(YACC.y) $<
#	$(COMPILE.c) -o $@ y.tab.c
%.o : %.c
	@echo 'Creating $(@)...'
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $(basename $@).c -o $@ 


print : 
	@echo "**************************************************************"
	@echo "**************************************************************"
	@echo "Settings for the Compilation of version $(VERSION) of $(PROGRAMME) in $(MKFILE_DIR)"
	@echo "With MAIN FILE $(MAIN.c)"
	@echo "$(CC) version : $(GCC_VER)"
	@echo "Compiler extension : $(COMPILO)"
	@echo "INCLUDES sat to $(INCLUDES)"
	@echo "compiling options $(CFLAGS)"
	@echo "with inclusion of  $(CPPFLAGS)"
	@echo "and static libraries $(LIBS)"
	@echo "**************************************************************"
	@echo "**************************************************************"
