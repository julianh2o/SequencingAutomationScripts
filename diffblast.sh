#!/bin/bash

DB=$1
ENTREZ=$2
INPUT_FASTA=$3
PREVIOUS=$4
OUTPUT_FOLDER=$5
MAFFT_OUTPUT_FOLDER=$6
MAFFT_RTF_OUTPUT_FOLDER=$7
ECUTOFF=$8 #1e-100

if [ ! -f INITIAL_BLAST.blast ]; then
    tblastn -remote -db $DB -entrez_query "$ENTREZ" -query $INPUT_FASTA -outfmt 5 > INITIAL_BLAST.blast
fi

viewblast.py list INITIAL_BLAST.blast -f '{Hit_accession}' | sed 's/^[ \t]*//;s/[ \t]*$//' | sort | uniq > CURRENT_ACCESSIONS.txt
cat $PREVIOUS | sed 's/\/\/.*//g' | sed 's/^[ \t]*//;s/[ \t]*$//' | sort | uniq > PREVIOUS_SORTED.txt
comm -13 PREVIOUS_SORTED.txt CURRENT_ACCESSIONS.txt | tr -d '\t' > NEW_ACCESSIONS.txt

NEW_COUNT=`wc -l NEW_ACCESSIONS.txt | awk '{print $1}'`
PREVIOUS_ACCESSIONS=`wc -l PREVIOUS_SORTED.txt | awk '{print $1}'`
echo "Found $NEW_COUNT new accessions ($PREVIOUS_ACCESSIONS previous)"

#rm -rf $OUTPUT_FOLDER
#rm -rf $MAFFT_OUTPUT_FOLDER
#rm -rf $MAFFT_RTF_OUTPUT_FOLDER

mkdir -p $OUTPUT_FOLDER
mkdir -p $MAFFT_OUTPUT_FOLDER
mkdir -p $MAFFT_RTF_OUTPUT_FOLDER

cat NEW_ACCESSIONS.txt | while read ACC; do
    [ -z "$ACC" ] && continue #skip blank lines

    FASTA_NAME=$ACC
    FASTA_PATH=$OUTPUT_FOLDER/$FASTA_NAME.txt
    MAFFT_FASTA_PATH=$MAFFT_OUTPUT_FOLDER/$FASTA_NAME.fasta
    MAFFT_RTF_PATH=$MAFFT_RTF_OUTPUT_FOLDER/$FASTA_NAME.rtf

    INDEX=`viewblast.py list INITIAL_BLAST.blast -f '{Hit_accession} {Hit_num}' | grep $ACC | awk '{print $2}'`
    E=`viewblast.py info $INDEX INITIAL_BLAST.blast -f '{e}'`
    TOOWEAK=`awk 'BEGIN { print ('$E'>'$ECUTOFF')? "1" : "0" }'`

    if [[ "$TOOWEAK" = "1" ]]; then
        echo "Skipping weak result: $ACC E=$E"
        continue #skip weak results
    fi

    if [ ! -f $FASTA_PATH ]; then
        echo "Downloading into $FASTA_NAME (e: $E)";
        HEADER=`viewblast.py info $INDEX INITIAL_BLAST.blast -f '{Hit_accession} {Hit_def} {Hit_hsps/Hsp/Hsp_bit-score} {e} {Hit_hsps/Hsp/Hsp_identity} {Hit_hsps/Hsp/Hsp_align-len}'`
        FRAME=`viewblast.py info $INDEX INITIAL_BLAST.blast -f '{Hit_hsps/Hsp/Hsp_hit-frame}'`
        COMMENT=`viewblast.py info $INDEX INITIAL_BLAST.blast -f '{Hit_def} bit_s:{Hit_hsps/Hsp/Hsp_bit-score} s:{Hit_hsps/Hsp/Hsp_score} e:{Hit_hsps/Hsp/Hsp_evalue} ident:{Hit_hsps/Hsp/Hsp_identity}'`
        NOHEADER=`fetchaccession.py $ACC | fastafromtraces.sh | tail -n -1`
        echo -e ">$HEADER\n$NOHEADER" | translate.py -r $FRAME > $FASTA_PATH

        MAFFTFASTA=`cat $INPUT_FASTA $FASTA_PATH`
        echo "$MAFFTFASTA" > TMP_MAFFT_INPUT.fasta
        mafft TMP_MAFFT_INPUT.fasta > $MAFFT_FASTA_PATH 2> /dev/null
        ANALYSIS=`alignmentanalysis.py $MAFFT_FASTA_PATH`
        cat $MAFFT_FASTA_PATH | formatmafft.py --header "$ANALYSIS" -o $MAFFT_RTF_PATH

        echo "$ACC //$COMMENT" >> $PREVIOUS
    fi
done
