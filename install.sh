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
    -p|--PREFIX)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        PREFIX=1
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

#Check for the presence of clusterutil
if [[ -z $(command -v "waitForJobs.sh") ]] ; then
  echo "Warning: clusterutil does not appear to be installed!">&2
  echo "You will not be able to use HPC multiplexing without!">&2
  sleep 1
fi
#check that we're in the conda base environment
if ! [[ -z $CONDA_DEFAULT_ENV || $CONDA_DEFAULT_ENV == "base" ]]; then
  echo "ERROR: You must run this installer from your conda 'base' environment!">&2
  exit 1
fi

#install conda environment
if [[ $CONDA == 1 ]]; then
  #Check for presence of conda
  if [[ -z $(command -v "conda") ]] ; then
    echo "ERROR: conda is not installed! Unable to proceed!">&2
    echo "Please either install anaconda/miniconda or use the --skipCondaEnv argument.">&2
    exit 1
  fi
  if conda env list|grep -q pacybara; then
    echo "Warning: Conda environment 'pacybara' already exists. Skipping conda setup.">&2
    sleep 1
  else
    conda env create -f pacybara_env.yml
  fi
fi

if [[ $DEPENDENCIES == 1 ]]; then
  #check if R is installed
  if [[ -z $(command -v Rscript) ]] ; then
    #if it's not found, try to activiate the new conda environment
    if conda env list|grep -q pacybara; then
      source "$CONDA_PREFIX/etc/profile.d/conda.sh"
      conda activate pacybara
    else 
      echo "ERROR: Neither an R installation, nor the 'pacybara' conda environment were found!">&2
      exit 1
    fi
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
  '
fi

cp -v src/*sh src/*R src/*.py ${PREFIX}

echo "Installation complete!"
