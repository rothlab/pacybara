#!/bin/bash
# Copyright (C) 2021, 2022  Jochen Weile, The Roth Lab
#
# This file is part of Pacybara.
#
# Pacybara is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Pacybara is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Pacybara.  If not, see <https://www.gnu.org/licenses/>.

#This script depends on submitjob.sh and waitForJobs.sh !

BASEDIR=~/projects/pacybara/workspace/pacybara
CONDAENV=pacybara
CONDAARG="--conda $CONDAENV"
BLACKLIST=""

#helper function to print usage information
usage () {
  cat << EOF

submitPacybara.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

A convenient job submission script of pacybara
Usage: submitPacybara.sh [-d|--basedir <DIR>]
  [-b|--blacklist <BLACKLIST>] [-c|--conda <ENV>] <PARAMS>

<PARAMS>       : A barseq parameter sheet file
-d|--basedir   : The base directory in which the output directory will be created. Default: $BASEDIR
-b|--blacklist : An optional comma-separated blacklist of nodes to avoid.
-c|--conda     : Conda environment to activate. Default: $CONDAENV

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
    -b|--blacklist)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        BLACKLIST=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -c|--conda)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        CONDAARG="--conda $2"
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -d|--basedir)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        BASEDIR=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    --) # end of options indicates that the main command follows
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

#parameter file
PARAMETERS="$(realpath "$1")"
if ! [[ -r "$PARAMETERS" ]]; then
  echo "$PARAMETERS not found or unreadable!">&2
  exit 1
fi

if [[ -z $BLACKLIST ]]; then
  BLARG=""
else
  BLARG="--blacklist $BLACKLIST"
fi

#helper function to extract relevant sections from a parameter file
extractParamSection() {
  INFILE="$1"
  SECTION="$2"
  TMPDIR="$3"
  # mkdir -p tmp
  case "$SECTION" in
    ARGUMENTS) TMPFILE=$(mktemp -p "$TMPDIR");;
    *SEQUENCE*) TMPFILE=$(mktemp -p "$TMPDIR" --suffix=.fasta);;
    SAMPLE) TMPFILE=$(mktemp -p "$TMPDIR" --suffix=.tsv);;
    *) die "Unrecognized section selected!";;
  esac
  RANGE=($(grep -n "$SECTION" "$INFILE"|cut -f1 -d:))
  sed -n "$((${RANGE[0]}+1)),$((${RANGE[1]}-1))p;$((${RANGE[1]}))q" "$INFILE">"$TMPFILE"
  echo "$TMPFILE"
}

#load environment variables from parameter sheet
TMPDIR=$(mktemp -d -p ./)
source $(extractParamSection "$PARAMETERS" 'ARGUMENTS' "$TMPDIR")
rm -r "$TMPDIR"

echo "Submitting job..."

mkdir -p "$WORKSPACE"
cd "$WORKSPACE"
submitjob.sh -n pacybara -c 12 -m 24G -t 7-00:00:00 \
	-l pacybara.log -e pacybara.log \
	-q guru $CONDAARG $BLARG -- \
	pacybara.sh $BLARG "$PARAMETERS"

echo "Logfile: ${WORKSPACE}/pacybara.log"

echo "Done!"
