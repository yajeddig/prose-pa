# Guide de compilation — ProSe-PA v0.78

## Contexte

ProSe-PA est un modèle hydro-biogéochimique développé par Mines Paris (Geosciences), couplant hydraulique 1D, transport de solutés et biogéochimie (modèle RIVE). Le code C-ANSI dépend de 12 bibliothèques externes hébergées sur GitLab, dont les branches `main` évoluent indépendamment **sans gestion de versions compatible** (pas de tags, pas de lock file).

> **Problème principal** : `libttc` (transport/transfert) a subi un refactoring structurel majeur entre v0.13 et v0.33 (suppression de `plink`, éclatement de `s_species_ttc` en sous-structures `s_settings_ttc`, `s_carac_ttc`, etc.). ProSe-PA v0.78 est écrit pour l'API v0.13. Cloner `main` HEAD provoque **~324 erreurs de compilation**.

## Prérequis système

```bash
# Ubuntu 22.04+ / Debian 12+
sudo apt-get install gcc gfortran flex bison make
```

## Étape 1 — Cloner les sources

```bash
WORK_DIR=$HOME/prose-pa-build
mkdir -p $WORK_DIR/LIBS && cd $WORK_DIR

# ProSe-PA
git clone https://gitlab.com/prose-pa/prose-pa.git src

# 12 bibliothèques dépendantes
cd $WORK_DIR/LIBS
git clone https://gitlab.com/gutil/libprint.git
git clone https://gitlab.com/gutil/libts.git
git clone https://gitlab.com/gutil/libpc.git
git clone https://gitlab.com/gutil/libio.git
git clone https://gitlab.com/gutil/libgc.git
git clone https://gitlab.com/gutil/libchronos.git
git clone https://gitlab.com/ghydro/libhyd.git
git clone https://gitlab.com/gtransp/libttc.git
git clone https://gitlab.com/gtransp/c-rive.git
git clone https://gitlab.com/gtransp/libseb.git
git clone https://gitlab.com/gutil/libtube.git
git clone https://gitlab.com/gutil/libmb.git
```

## Étape 2 — Rollback libttc à v0.13

**C'est l'étape critique.** Sans elle, la compilation échoue.

```bash
cd $WORK_DIR/LIBS/libttc
git checkout 36b8f006
```

> Ce commit (2022-09-01) correspond au dépôt initial de libttc v0.13, seule version compatible avec l'API utilisée par ProSe-PA v0.78.

## Étape 3 — Compiler les bibliothèques

L'ordre de compilation respecte les dépendances entre bibliothèques. Pour chaque lib, adapter le `Makefile` :

- `PATH_INST` → chemin absolu vers `$WORK_DIR/LIBS`
- Remplacer les chemins `branches/C302/src` ou `branches/main/src` par `src`

```bash
export PATH_INST=$WORK_DIR/LIBS

# Ordre de compilation (dépendances croissantes)
LIBS_ORDER="libprint libts libpc libio libgc libchronos libttc libhyd c-rive libseb libtube libmb"

for lib in $LIBS_ORDER; do
    echo "=== Compiling $lib ==="
    cd $PATH_INST/$lib/src

    # Adapter les chemins dans le Makefile
    sed -i "s|PATH_INST=.*|PATH_INST=$PATH_INST|" Makefile
    sed -i 's|/branches/C302/src|/src|g' Makefile
    sed -i 's|/branches/main/src|/src|g' Makefile
    sed -i 's|/trunk/src|/src|g' Makefile

    make clean && make lib
done
```

### Cas particulier : libgc (sparse)

Le sous-répertoire `sparse_11_7_2011` nécessite une compilation séparée :

```bash
cd $PATH_INST/libgc/src/sparse_11_7_2011
gcc -c -O2 *.c
ar rcs sparse.a *.o
```

## Étape 4 — Compiler ProSe-PA

```bash
cd $WORK_DIR/src
```

Adapter le `Makefile` avec les versions effectivement compilées :

| Variable | Valeur |
|----------|--------|
| `PATH_INST` | `$WORK_DIR/LIBS` |
| `V_TTC` | `0.13` |
| `V_HYD` | valeur dans `libhyd/src/Makefile` (ex: `0.44`) |
| `V_GC` | idem (ex: `0.13`) |
| `V_IO` | idem (ex: `0.21`) |
| `V_RIVE` | idem (ex: `0.37`) |

> **Astuce** : les versions sont lisibles via `grep "^V_" $PATH_INST/<lib>/src/Makefile`

```bash
# Adapter le Makefile de prose-pa
sed -i "s|PATH_INST=.*|PATH_INST=$PATH_INST|" Makefile
sed -i "s|INCL_PROSE=.*|INCL_PROSE=-I$WORK_DIR/src|" Makefile
# Ajuster chaque V_XXX selon les versions compilées

make clean && make all
```

Résultat attendu : binaire `prose-pa0.78` (~2.7 Mo).

```bash
./prose-pa0.78
# > ProSe-PA0.78 CommandFileName DebugFileName
```

## Résumé des versions compatibles

```
prose-pa    v0.78   main HEAD       ✅
libttc      v0.13   commit 36b8f006 ⚠️  ROLLBACK OBLIGATOIRE
libhyd      v0.44   main HEAD       ✅
libgc       v0.13   main HEAD       ✅
libio       v0.21   main HEAD       ✅
c-rive      v0.37   main HEAD       ✅
libprint    v1.24   main HEAD       ✅
libts       v1.74   main HEAD       ✅
libpc       v0.05   main HEAD       ✅
libchronos  v0.13   main HEAD       ✅
libseb      v0.09   main HEAD       ✅
libtube     v0.08   main HEAD       ✅
libmb       v0.01   main HEAD       ✅
```

## Avertissements

- Les warnings `new_function` et `new_simulation redefined` (entre libts/c-rive et libhyd) sont bénins — macros identiques définies dans deux headers.
- La compilation utilise `gfortran` pour le link final (dépendance Fortran dans libgc/sparse).
- Le flag `-DGCC1330` est défini automatiquement si `gcc --version` ≥ 13.x.
- Le flag `-DCOUPLED_RIVE` active le couplage biogéochimique. Sans lui, seul le transport thermique fonctionne.
