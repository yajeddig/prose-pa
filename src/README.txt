Installation automatique de ProSe-PA

> make_pprose.sh all
Reinstalle et compile Toutes les libraries dans $LIB_HYDROSYSTEM_PATH. Compile prose-pa0.32 dans ./

> make_pprose.sh PATH (à regarder pour libgc où le chemin d'accès au headers n'est pas bon)
Installe toutes les librairies et les recompile dans PATH Compile prose-pa0.32 dans PATH

> make_pprose.sh PATH update
Compile prose-pa0.32 dans ./ en cherchant les librairies dans PATH et en les y recompilant

###### installation out school ###############
> make_pprose_outschool.sh PATH reinstall login
Reinstalle les librairies dans $PATH, et compile prosepa dans ./
> make_pprose_outschool.sh PATH reinstall login -b <branch_name>
Reinstalle les librairies dans $PATH, et compile prosepa dans ./ a partir des branches <branch_name>


##### installation from branches ######

If you want to install ProSe-P from LIBS branches then be sure that the name of the branch is the same in all libs and use

> make_pprose_from_branches.sh all
calls  make_pprose.sh all
Reinstalle et compile Toutes les libraries dans $LIB_HYDROSYSTEM_PATH. Compile prosepX.XX dans ./

> make_pprose_from_branches.sh PATH
calls  make_pprose.sh PATH
Installe toutes les librairies et les recompile dans PATH Compile prosepX.XX dans ./

> make_pprose_from_branches.sh
calls  make_pprose.sh 
Reinstalle et compile Toutes les libraries dans ./. Compile prosepX.XX dans ./


> make_pprose_from_branches.sh -b <branch_name>
Reinstalle et compile Toutes les libraries dans $LIB_HYDROSYSTEM_PATH en utilisant le nom de branche s'il existe et trunk sinon. Compile prosepX.XX dans ./


> make_pprose_from_branches.sh -b <branch_name> PATH
Installe toutes les librairies et les recompile dans PATH en utilisant le nom de branche s'il existe et trunk sinon. Compile prosepX.XX dans ./
