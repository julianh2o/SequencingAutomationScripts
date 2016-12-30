#!/bin/bash

DB=tsa_nt
PREVIOUS=previous.txt
RESULTS_FOLDER=results
OUTPUT_FOLDER=./fastaupdates
MAFFT_OUTPUT_FOLDER=./mafftoutput
MAFFT_RTF_OUTPUT_FOLDER=./mafftoutput_rtf

QUERY_LIST=queries.fasta
TAXA_QUERIES=`cat taxa.txt | sed 's/\/\/.*//g'`

echo "$TAXA_QUERIES" | while read taxa
do
    TAXA_IDENT=`echo -n $taxa | awk -F'###' '{print $1}'`
    QUERY=`echo -n $taxa | awk -F'###' '{print $2}'`
    E=`echo -n $taxa | awk -F'###' '{print $3}'`

    extractfastaseq.py $QUERY_LIST -l | while read fasta
    do
        HEAD=`echo -n $fasta | awk -F'###' '{print $1}'`
        SEQ=`echo -n $fasta | awk -F'###' '{print $2}'`
        FASTA_IDENT=`echo -n $HEAD | awk '{print $1}' | tr -d '>'`

        QUERY_NAME="$TAXA_IDENT"_"$FASTA_IDENT"

        echo "Blasting: $QUERY_NAME"
        mkdir -p "$RESULTS_FOLDER/$QUERY_NAME"
        cd "$RESULTS_FOLDER/$QUERY_NAME"

        echo -e "$HEAD\n$SEQ" > QUERY.fasta
        touch $PREVIOUS
        diffblast.sh $DB "$QUERY" QUERY.fasta $PREVIOUS $OUTPUT_FOLDER $MAFFT_OUTPUT_FOLDER $MAFFT_RTF_OUTPUT_FOLDER

        cd ../..

    done

done
