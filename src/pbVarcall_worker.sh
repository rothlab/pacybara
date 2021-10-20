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
REFSEQ=${2:-"LDLR_cds.fa"}

DEFAULTIFS=$IFS

callVars() {
    ALNFILE=$1
    BARCODE=$2
    OUTFILE=$3
    LLEN=60
    # VARFILE=$(echo "$ALNFILE"|sed -r "s/\\.aln$/.vc/")
    #read all lines from the file into an array
    mapfile -t LINES<"$ALNFILE"
    #find the line number with FASTA headers (">")
    HEADERS=($(grep -n ">" $1|cut -f 1 -d :))
    printf "${BARCODE}\t">>$OUTFILE
    #the number of sequence alignment lines
    NLINES=$((${HEADERS[1]}-2))
    INSCOUNT=0
    MCOUNT=0
    #i iterates over the lines in the alignment
    for (( i = 0; i < $NLINES; i++ )); do
        ALTLINE=${LINES[$(($i+1))]}
        REFLINE=${LINES[$(($i+${HEADERS[1]}))]}
        # LLEN=${#REFLINE}
        # echo "$(($i+1)) $(($i+${HEADERS[1]})) - $LLEN"
        #j iterates over bases in line
        for (( j = 0; j < $LLEN; j++ )); do
            #if the read differes from the reference at this base
            if [[ ${ALTLINE:$j:1} != ${REFLINE:$j:1} ]]; then
                #prepend a semicolon if this isn't the first mutation in here
                if [[ $MCOUNT > 0 ]]; then
                    printf ";">>$OUTFILE
                fi
                #if the reference base is a "-", it's an insertion
                if [[ ${REFLINE:$j:1} == "-" ]]; then
                    ((INSCOUNT++))
                    POS=$(($LLEN*$i + $j+1 - $INSCOUNT))
                    printf "$(($POS-1))_${POS}ins${ALTLINE:$j:1}">>$OUTFILE
                #if the read base is a "-", it's a deletion
                elif [[ ${ALTLINE:$j:1} == "-" ]]; then
                    POS=$(($LLEN*$i + $j+1 - $INSCOUNT))
                    printf "${POS}del">>$OUTFILE
                #otherwise it's a substitution
                else
                    POS=$(($LLEN*$i + $j+1 - $INSCOUNT))
                    printf "${POS}${REFLINE:$j:1}>${ALTLINE:$j:1}">>$OUTFILE
                fi
                ((MCOUNT++))
            fi
        done
    done
    if [[ $MCOUNT == 0 ]]; then
        printf "=">>$OUTFILE
    fi
    echo "">>$OUTFILE
}

processChunk() {
    CHUNKFILE=$1
    REFSEQ=${2:-"LDLR_cds.fa"}
    mkdir -p alignments

    # OUTFILE=$(echo "$CHUNKFILE"|sed -r "s/\\.txt$/_varcall.txt/")
    OUTFILE=${CHUNKFILE}_varcall.txt
    while IFS= read -r LINE; do
        IFS=$'\t' read -ra SEQS <<<"$LINE"
        IFS=$DEFAULTIFS
        BC=${SEQS[0]}
        # echo "Processing barcode: $BC"
        ALNOUT=alignments/${BC}.aln
        echo "${SEQS[1]}"|needle -filter -bsequence $REFSEQ -outfile "$ALNOUT" -auto -aformat3 fasta &&\
        callVars "$ALNOUT" "$BC" "$OUTFILE" && rm "$ALNOUT"
    done < "$CHUNKFILE"
    rm "$CHUNKFILE"
}

processChunk $INFILE $REFSEQ
