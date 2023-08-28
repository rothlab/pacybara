#!/bin/bash
# Copyright (C) 2021  Jochen Weile, Roth Lab
#
# This file is part of BarseqPro.
#
# BarseqPro is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BarseqPro is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with BarseqPro.  If not, see <https://www.gnu.org/licenses/>.

CONDA=1
DEPENDENCIES=1
PREFIX="${HOME}/bin/"
#fail on error, even within pipes; disallow use of unset variables, enable history tracking
set -euo pipefail +H

#helper function to print usage information
usage () {
  cat << EOF

install.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

Installs Pacybara
Usage: install.sh [-c|--skipCondaEnv] [-d|--skipDependencies] [-s|--skipAll]
                  [-p|--PREFIX <DIRECTORY>] 

-c|--skipCondaEnv        : skip conda environment setup
-d|--skipDependencies    : skip installing R package dependencies
-s|--skipAll             : skips all conda and dependencies setup. 
                           Only installs script files. This is eqivalent
                           to install.sh -c -d .
-p|--PREFIX  : Sets the installation directory. Defaults to $PREFIX

EOF
 exit $1
}

#Parse Arguments
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      usage 0
      shift
      ;;
    -c|--skipCondaEnv)
      CONDA=0
      shift
      ;;
    -d|--skipDependencies)
      DEPENDENCIES=0
      shift
      ;;
    -s|--skipAll)
      CONDA=0
      DEPENDENCIES=0
      shift
      ;;
    -p|--prefix)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        PREFIX=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    --) # end of options indicates that only positional arguments follow
      shift
      PARAMS="$PARAMS $@"
      eval set -- ""
      ;;
    -*|--*=) # unsupported flags
      echo "ERROR: Unsupported flag $1" >&2
      usage 1
      ;;
    *) # positional parameter
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
#reset command arguments as only positional parameters
eval set -- "$PARAMS"

checkPath() {
  VIABLE=""
  for LOCATION in $(echo $PATH|tr : "\n"|sort|uniq); do
    if [[ -w ${LOCATION} ]]; then
      if [[ -z $VIABLE ]]; then
        VIABLE=$LOCATION
      else
        VIABLE="${VIABLE} ; $LOCATION"
      fi
    fi
  done
  echo "$VIABLE"
}

#helper function to check if user wants to proceed with installation.
checkProceed() {
  sleep 1
  printf "\nAre you sure you want to proceed anyway? [y/n]: "
  read -r ANSWER
  if ! [[ $ANSWER == "y" || $ANSWER == "yes" ]]; then
    exit 1
  fi
}

#Check that PREFIX is writable and listed in PATH
if ! [[ -w $PREFIX ]]; then
  printf "\033[0;31mERROR: You don't have permissions to write to your chosen installation directory (${PREFIX}).

Please use the --PREFIX option to choose a more appropriate target directory. \033[0m\n">&2

  exit 1

elif ! [[ $PATH == *${PREFIX%/}* ]]; then 

  printf "\033[1;33mWARNING: Your chosen installation directory (${PREFIX}) is not listed in your \$PATH variable, which means that your command shell would be unable to find it there.

Please use the --PREFIX option to choose a more appropriate target directory. For example, here are some writable locations listed in your \$PATH that could work: $(checkPath)

Alternatively, you could add $PREFIX to your \$PATH definition in your .bashrc file.\033[0m\n"
  checkProceed
fi

#Check for the presence of clusterutil
if [[ -z $(command -v "waitForJobs.sh") ]] ; then
  printf "\033[1;33mWARNING: clusterutil does not appear to be installed!
You will not be able to use HPC multiplexing without it.\033[0m\n"
  checkProceed
fi
#check whether conda, mamba or micromamba are installed
if [[ -n $(command -v conda) ]]; then
  CONDAMGR=conda
elif [[ -n $(command -v mamba) ]]; then
  CONDAMGR=mamba
elif [[ -n $(command -v micromamba) ]]; then
  CONDAMGR=micromamba
elif [[ -n "$CONDAENV" ]]; then
  printf "\033[0;31mERROR: No conda installation was found!\033[0m\n">&2
  exit 1
fi
#check that we're in the conda base environment
if ! [[ -z $CONDA_DEFAULT_ENV || $CONDA_DEFAULT_ENV == "base" ]]; then
  printf "\033[0;31mERROR: You must run this installer from your conda 'base' environment!\033[0m\n">&2
  exit 1
fi

#install conda environment
if [[ $CONDA == 1 ]]; then
  #Check for presence of conda
  if [[ -z $(command -v "conda") ]] ; then
    printf "\033[0;31mERROR: conda is not installed! Unable to proceed!
Please either install anaconda/miniconda or use the --skipCondaEnv argument.\033[0m\n">&2
    exit 1
  fi
  if ${CONDAMGR} env list|grep -q pacybara; then
    printf "\033[1;33mWARNING: Conda environment 'pacybara' already exists. Skipping conda setup.\033[0m\n">&2
    sleep 1
  else
    ${CONDAMGR} env create -f pacybara_env.yml
  fi
fi

if [[ $DEPENDENCIES == 1 ]]; then
  #if the pacybara conda environment was used, then we need to install the packages 
  # into that environment's R installation
  if [[ $CONDA == 1 ]]; then
    #so we need to activate the environment first
    if ${CONDAMGR} env list|grep -q pacybara; then
      source "$CONDA_PREFIX/etc/profile.d/${CONDAMGR}.sh"
      ${CONDAMGR} activate pacybara
    else 
      printf "\033[0;31mERROR: The pacybara conda environment could not be found!\033[0m\n">&2
      exit 1
    fi
  #otherwise we install it into the default R installation.
  elif [[ -z $(command -v Rscript) ]] ; then
    printf "\033[0;31mERROR: Conda setup was skipped, but no R installation was found in the conda environment!
If you really want to skip the conda setup, please manually install all dependencies first.\033[0m\n">&2
    exit 1
  fi
  #install required R packages
  Rscript -e '
    cran.needed <- c("remotes","argparser","hash","bitops","pbmcapply")
    cran.missing <- setdiff(cran.needed,rownames(installed.packages()))
    if (length(cran.missing) > 0) {
      install.packages(cran.missing,repos="https://cloud.r-project.org/")
    }
    library(remotes)
    github.needed <- data.frame(user=c("jweile","jweile","VariantEffect"),package=c("yogitools","yogiseq","hgvsParseR"))
    github.missing <- which(!(github.needed$package %in% rownames(installed.packages())))
    if (length(github.missing) > 0) {
      toInstall <- apply(github.needed[github.missing,],1,paste,collapse="/")
      invisible(lapply(toInstall,install_github))
    }
    total <- length(cran.missing)+length(github.missing)
    if (total == 0) {
      cat("All required R packages are already present!\n")
    } else {
      cat("Installed ",total," R packages!\n")
    }
  '
fi

cp -v src/*sh src/*R src/*.py ${PREFIX}

echo "Installation complete!"
