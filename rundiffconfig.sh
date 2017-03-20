#!/bin/bash

#Look at rundiffsample.config for an example of how to configure this
#Run with rundiffconfig.sh myconfig.config

CONFIG_FILE=$1
source $1

echo "PREVIOUS_FILE: $PREVIOUS_FILE"
echo "ECUTOFF:       $ECUTOFF"
echo "RESULT_FOLDER: $RESULT_FOLDER"
echo "RESULT_SEQ:    $RESULT_SEQ"
echo "RESULT_MAFFT:  $RESULT_MAFFT"
echo "RESULT_RTF:    $RESULT_RTF"
echo "BLAST_CMD:     $BLAST_CMD"

#echo "QUERY:         $QUERY"
#echo "DATABASE:      $DATABASE"
#echo "INCLUDE:       $INCLUDE"
#echo "EXCLUDE:       $EXCLUDE"

SEQUENCE_OUTPUT=$RESULT_FOLDER/$RESULT_SEQ
MAFFT_OUTPUT=$RESULT_FOLDER/$RESULT_MAFFT
RTF_OUTPUT=$RESULT_FOLDER/$RESULT_RTF

echo

if [ ! -f INITIAL_BLAST.blast ]; then
    echo "Executing: " $BLAST_CMD
    $BLAST_CMD > INITIAL_BLAST.blast
fi

echo "Executing: " diffblast.sh INITIAL_BLAST.blast $PREVIOUS_FILE $SEQUENCE_OUTPUT $MAFFT_OUTPUT $RTF_OUTPUT $ECUTOFF
diffblast.sh INITIAL_BLAST.blast $PREVIOUS_FILE $SEQUENCE_OUTPUT $MAFFT_OUTPUT $RTF_OUTPUT $ECUTOFF

rm INITIAL_BLAST.blast
