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

INFILE=$1
# REFSEQ="LPL_cds.fa"
REFSEQ=$2
OUTFILE=$(echo "$INFILE"|sed -r "s/\\.txt.gz$/_varcalls.txt/")

mkdir -p chunks
zcat "$INFILE"|split -l 1000 - chunks/chunk
for CHUNK in $(ls chunks/*); do
    ID=$(basename $CHUNK)
    submitjob.sh -t 00:30:00 -n "$ID" -- ./pbVarcall_worker.sh "$CHUNK" "$REFSEQ"
done

waitForJobs.sh

cat chunks/*_varcall.txt>$OUTFILE&&rm chunks/*

submitjob.sh -t 00:30:00 -n translate -- pbVarcall_translate.R $OUTFILE parameters.json

waitForJobs.sh
