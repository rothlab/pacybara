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
if [[ -z $(command -v "conda") ]] ; then
  echo "Warning: conda is not installed! Proceeding without.">&2
  sleep 1
fi

if [[ $CONDA == 1 ]]; then
  conda env create -f myfile.yaml
fi

if [[ $DEPENDENCIES == 1 ]]; then
  if conda env list|grep -q pacybara; then
    conda activate pacybara
  else 
    echo "Warning: Pacybara conda environment not found. Proceeding without.">&2
  fi
  Rscript -e '
    install.packages(c("remotes","argparser","hash","bitops","pbmcapply"),repos="https://cloud.r-project.org/")
    library(remotes)
    install_github("jweile/yogitools")
    install_github("jweile/yogiseq")
    install_github("VariantEffect/hgvsParseR")
  '
fi

cp -v src/*sh src/*R src/*.py ${PREFIX}

echo "Installation complete!"
