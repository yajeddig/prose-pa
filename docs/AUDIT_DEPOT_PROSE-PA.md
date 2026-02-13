# Audit qualité du dépôt ProSe-PA v0.78

> **Date** : Février 2025  
> **Périmètre** : dépôt GitLab `prose-pa/prose-pa` + 12 bibliothèques dépendantes (`gutil/*`, `ghydro/*`, `gtransp/*`)  
> **Méthode** : analyse statique du code source, des Makefiles, de l'historique git et de la structure des repos

---

## 1. Structure du projet

```
prose-pa/
├── src/                  ← TOUT le code source (flat, aucun sous-dossier)
│   ├── 20 fichiers .c
│   ├── 13 fichiers .h
│   ├── input.y + lexical.l   (parser Bison/Flex)
│   ├── 22 fichiers .o         ⚠️ objets compilés versionnés
│   ├── prose-pa0.78           ⚠️ binaire versionné
│   ├── Makefile_tmp           ⚠️ fichier temporaire versionné
│   ├── make_prose-pa.sh~      ⚠️ backup éditeur versionné
│   └── README.txt             ⚠️ obsolète (réfère à v0.32)
├── LIBS/                 ← 12 dépendances clonées manuellement (pas de submodule)
├── post-processing/      ← vide (3× .gitkeep uniquement)
├── README.md             ← description générale + instructions d'installation
└── LICENSE               ← Eclipse Public License v2.0
```

### Constats

- **Arborescence plate** — 33 fichiers source dans un seul répertoire `src/`, sans séparation logique (`core/`, `io/`, `da/`, `include/`). Pour 21 000 lignes de code, c'est la limite de maintenabilité.
- **Artefacts de build versionnés** — 22 fichiers `.o`, le binaire `prose-pa0.78`, `Makefile_tmp` et `make_prose-pa.sh~` sont présents dans le repo. Cause : **aucun `.gitignore`** dans le projet.
- **`post-processing/` vide** — 3 fichiers `.gitkeep`, zéro script de post-traitement. Le README.md mentionne un lien vers une notice PDF (`notice_CRIVE_aout2012.pdf`) dont l'URL est cassée (typo dans le chemin GitLab : `prose-pa-/` au lieu de `prose-pa/`).
- **Aucun dossier** `tests/`, `examples/`, `doc/`, `benchmarks/`.

---

## 2. Documentation

| Élément | État | Commentaire |
|---|---|---|
| README.md (racine) | ✅ Existe | Bon contenu : description, citations, install, licence |
| README.txt (src/) | ⚠️ Obsolète | Réfère à `prose-pa0.32`, noms de scripts anciens (`make_pprose.sh`) |
| Guide utilisateur | ❌ Lien cassé | Pointe vers un PDF absent sur GitLab |
| Documentation technique / API | ❌ Absente | Pas de Doxygen, pas de `Doxyfile` |
| Documentation d'architecture | ❌ Absente | Aucun diagramme, aucune description du couplage entre libs |
| CHANGELOG | ❌ Absent | |
| CONTRIBUTING | ❌ Absent | |
| Format du fichier de commande | ❌ Non documenté | 461 keywords dans le lexer, aucune doc associée |
| Commentaires en-tête | ✅ 32/33 fichiers | Header EPL standardisé avec contributeurs et citations |
| Commentaires inline | ⚠️ 12% dead code | 1 584 lignes commentées `//...` = code mort conservé, pas de la documentation utile |

**Verdict** : le README est correct pour un premier contact, mais il n'existe **aucune documentation technique exploitable** — ni pour un utilisateur (format des fichiers d'entrée, paramètres attendus), ni pour un développeur (architecture interne, API des bibliothèques, conventions de couplage). La seule notice utilisateur référencée date de 2012 et son URL est cassée.

---

## 3. Système de branches et versioning

### État des branches par dépôt

| Repo | Branches | Tags |
|---|---|---|
| `prose-pa/prose-pa` | `main`, `SW` | ∅ |
| `gtransp/libttc` | `main`, `C302`, `C345`, `C350` | `c-345-033` |
| `ghydro/libhyd` | `main`, `C302`, `C345`, `C350`, `19inflows_BCK` | `c-345-042` |
| `gutil/libgc` | `main`, `C302`, `C345`, `C350` | `c-345-013` |
| `gutil/libio` | `main`, `C302`, `C345`, `C350` | `c-345-018` |
| `gtransp/c-rive` | `main`, `SW` | ∅ |
| `gutil/libprint` | `main`, `C302`, `C345`, `C350` | `c-345-124` |
| `gtransp/libseb` | `main`, `C302`, `C345`, `C350` | `c-345-009` |
| `gutil/libts` | `main`, `C302`, `C345`, `C350` | `c-345-174` |
| `gutil/libpc` | `main`, `C302`, `C345`, `C350` | `c-345-005` |
| `gutil/libchronos` | `main`, `C302`, `C345`, `C350` | `c-345-013` |
| `gutil/libtube` | Privé (404) | — |
| `gutil/libmb` | Privé (404) | — |

### Problèmes critiques

- **Aucun tag sur prose-pa** — impossible de relier un commit à un numéro de version. Le `0.78` est hardcodé dans `Makefile` et `PROSE.h`, pas dans un tag git.

- **Nomenclature de branches incohérente** — les bibliothèques `gutil/ghydro/gtransp` utilisent `C302`, `C345`, `C350` (identifiants de projets CaWaQS), tandis que prose-pa et c-rive utilisent `SW` (initiales d'un contributeur). Aucune convention partagée.

- **Pas de correspondance branches prose ↔ branches libs** — le script `make_prose-pa_from_branches.sh` suppose le même nom de branche dans tous les repos. En pratique, `SW` n'existe que dans prose-pa et c-rive, pas dans libttc (qui a `C302`). Résultat : le script clone `main` de libttc (v0.33), qui est **incompatible** avec prose-pa `main` (attend v0.13). C'est la cause des 324 erreurs de compilation constatées (voir `BUILD_GUIDE.md`).

- **Pas de git submodules** — les 12 bibliothèques sont clonées indépendamment par les scripts shell. Aucun mécanisme de pinning (commit SHA, tag, submodule). La reproductibilité dépend entièrement de l'état de `main` au moment du `git clone`.

- **Tags sur les libs nommés `c-345-0XX`** → liés au projet CaWaQS v3.45, pas à ProSe-PA. Inutilisables pour figer une version compatible ProSe-PA.

- **2 repos privés** — `libtube` et `libmb` renvoient 404 pour les utilisateurs non authentifiés, ce qui empêche la compilation sans accès préalable.

---

## 4. Reproductibilité du build

| Critère | État |
|---|---|
| Verrouillage des versions de dépendances | ❌ Aucun (ni submodule, ni lock file, ni tag) |
| Version compatible documentée | ❌ Nulle part (avant notre `BUILD_GUIDE.md`) |
| Build reproductible | ❌ `git clone` + `make` donne un résultat différent selon la date |
| Système de build | ⚠️ `Makefile` brut avec chemins hardcodés |
| CI/CD | ❌ Aucun — pas de `.gitlab-ci.yml` |
| Chemin hardcodé dans Makefile | `INCL_PROSE=-I$(HOME)/Programmes/prose/src` — spécifique à la machine du mainteneur |

Le script de build `make_prose-pa.sh` repose sur une variable d'environnement `$LIB_HYDROSYSTEM_PATH` (par défaut `$HOME/Programmes/LIBS/`) et un ensemble de scripts `awk` utilitaires clonés depuis `gutil/scripts`. L'ensemble forme une chaîne fragile qui n'a jamais été testée hors de l'environnement du labo d'origine.

---

## 5. Qualité du code source

### Métriques

| Métrique | Valeur | Appréciation |
|---|---|---|
| LOC total (src/) | 21 103 | Modéré |
| Plus gros fichier | `input.y` — 5 172 lignes | ⚠️ Monolithe : parser + init + allocation + I/O |
| `main_PROSE.c` | 939 lignes | ⚠️ Fonction `main()` de ~900 lignes — devrait être découpée |
| `manage_link_prose.c` | 2 357 lignes | Couplage prose↔libttc, très sensible aux changements d'API |
| Dead code commenté | 12% (1 584 lignes) | ⚠️ Fragments `//.../` laissés partout, bruit significatif |
| TODO/FIXME/HACK | 1 seul | Suspect : soit le code est considéré "fini", soit personne ne les pose |
| Warnings compilation | ~30 | Redéfinitions macros, sprintf overflow, format args — non traités |
| Conventions de nommage | Mixtes | `PROSE_xxx()`, `Prose_xxx()`, `calc_xxx()` — incohérent |
| Gestion mémoire | `malloc` sans `free` visible dans `main()` | Fuites probables (acceptable pour un run unique) |
| Tests unitaires | 0 | Aucun test automatisé |
| Cas-tests de référence | 0 | Aucun jeu de données d'exemple fourni |

### Observations qualitatives

- **`input.y` (5 172 lignes)** est un fichier Bison monolithique qui cumule le parsing du fichier de commande, l'allocation mémoire des structures, l'initialisation des paramètres physiques et l'ouverture des fichiers de sortie. Toute modification du format d'entrée impacte potentiellement l'ensemble.

- **`main_PROSE.c` (939 lignes)** contient une seule fonction `main()` qui gère l'initialisation, les deux modes (steady/transient), la boucle temporelle, l'assimilation de données et les sorties. Le contrôle de flux repose sur une cascade de `if(Simul->calc_mode[X] == YES_TS)` imbriqués sur 3-4 niveaux.

- **12% de dead code** sous forme de lignes commentées (`//for(np=0...`, `//PROSE_fill_var...`, `//fpno3 = fopen...`). Ce code mort rend la lecture difficile et masque la logique réelle.

- **Warnings ignorés** — les macros `new_function()` et `new_simulation()` sont redéfinies identiquement dans `libts`, `c-rive` et `libhyd`. Les `sprintf` vers des buffers trop petits (`char s[4]` pour un entier) sont des bugs latents.

---

## 6. Synthèse — Notation par axe

| Axe | Note /5 | Justification |
|---|---|---|
| **Arborescence projet** | 2/5 | Flat, pas de séparation logique, artefacts de build versionnés |
| **Documentation** | 1.5/5 | README correct, mais zéro doc technique/utilisateur accessible |
| **Branching & versioning** | 1/5 | Pas de tags prose-pa, pas de submodules, incompatibilités silencieuses entre repos |
| **Reproductibilité build** | 1/5 | Chemins hardcodés, pas de lock, pas de CI, `main` = cible mouvante |
| **Qualité code** | 2.5/5 | Headers corrects, structure modulaire raisonnable, mais dead code, main monolithique |
| **Tests & validation** | 0/5 | Aucun test, aucun exemple, aucun benchmark |

---

## 7. Recommandations

### Immédiat (avant toute contribution)

1. **Ajouter un `.gitignore`** — exclure `*.o`, `prose-pa*`, `*.tmp`, `*~`, `*.orig`, `input.c`, `lexical.c`
2. **Supprimer les artefacts** du repo (`git rm --cached *.o prose-pa0.78 Makefile_tmp`)
3. **Figer les versions des bibliothèques** via un fichier `VERSIONS.lock` (a minima) ou `git submodule` :
   ```
   # VERSIONS.lock — commits compatibles avec prose-pa 0.78
   libttc      36b8f006  # v0.13 — OBLIGATOIRE (main=v0.33 incompatible)
   libhyd      HEAD      # v0.44
   libgc       HEAD      # v0.13
   libio       HEAD      # v0.21
   c-rive      HEAD      # v0.37
   libprint    HEAD      # v1.24
   libts       HEAD      # v1.74
   libpc       HEAD      # v0.05
   libchronos  HEAD      # v0.13
   libseb      HEAD      # v0.09
   libtube     HEAD      # v0.08
   libmb       HEAD      # v0.01
   ```

### Court terme

4. **Ajouter un `.gitlab-ci.yml` minimal** — un job `compile` qui clone les libs aux bonnes versions et lance `make`. Si ça casse, on le sait immédiatement.
5. **Migrer vers CMake** — remplace les 13 Makefiles avec chemins hardcodés par un système portable et introspectable. Gestion native des dépendances entre libs.
6. **Corriger le chemin hardcodé** dans `Makefile` : `INCL_PROSE=-I$(HOME)/Programmes/prose/src` → relatif ou paramétrable.

### Moyen terme

7. **Inclure un cas-test de référence** — un jeu de données réduit (ex. tronçon de Seine, quelques pas de temps) avec les résultats attendus, permettant de valider toute recompilation.
8. **Refactorer `main_PROSE.c`** — extraire les phases en fonctions : `init_phase()`, `steady_solve()`, `transient_loop()`, `da_step()`.
9. **Refactorer `input.y`** — séparer le parsing grammatical de l'initialisation des structures et de l'allocation mémoire.
10. **Nettoyer le dead code** — supprimer les 1 584 lignes commentées, se fier à git pour l'historique.
11. **Documenter le format du fichier de commande** — les 461 keywords du lexer ne sont documentés nulle part. Un utilisateur ne peut pas construire un cas d'étude sans exemple fonctionnel ou assistance directe.

---

## 8. Contexte

Ce dépôt est un **outil de recherche mono-équipe** (Geosciences, Mines Paris) publié en open source sous EPL 2.0. Il n'a pas été conçu pour être compilé par des utilisateurs extérieurs sans assistance directe du mainteneur (Nicolas Flipo). L'absence de réponse du mainteneur depuis 2 ans rend le projet effectivement inutilisable sans rétro-ingénierie, comme en témoigne le travail nécessaire pour simplement obtenir une compilation fonctionnelle (identification du commit compatible libttc v0.13 par analyse de l'historique git — voir `BUILD_GUIDE.md`).
