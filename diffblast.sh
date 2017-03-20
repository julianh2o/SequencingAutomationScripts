#!/bin/bash

QUERY_FASTA=$1
BLAST_FILE=$2
PREVIOUS=$3
OUTPUT_FOLDER=$4
MAFFT_OUTPUT_FOLDER=$5
MAFFT_RTF_OUTPUT_FOLDER=$6
ECUTOFF=$7 #1e-100

viewblast.py list $BLAST_FILE -f '{Hit_accession}' | sed 's/^[ \t]*//;s/[ \t]*$//' | sort | uniq > CURRENT_ACCESSIONS.txt
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

DATE=`date +%Y-%m-%d`

echo "//[$DATE] Running diffblast: fasta:$QUERY_FASTA previous:$PREVIOUS ecut:$ECUTOFF output:$OUTPUT_FOLDER,$MAFFT_OUTPUT_FOLDER,$MAFFT_RTF_OUTPUT_FOLDER" >> $PREVIOUS

viewblast.py list $BLAST_FILE -f '{Hit_accession}|{Hit_num}|{e}|{Hit_def}|{Hit_hsps/Hsp/Hsp_bit-score}|{Hit_hsps/Hsp/Hsp_identity}|{Hit_hsps/Hsp/Hsp_align-len}|{Hit_hsps/Hsp/Hsp_hit-frame}|{cover}' | while read LINE; do
    ACC=`echo $LINE | awk -F'|' '{print $1}'`
    INDEX=`echo $LINE | awk -F'|' '{print $2}'`
    E=`echo $LINE | awk -F'|' '{print $3}'`
    DESCRIPTION=`echo $LINE | awk -F'|' '{print $4}'`
    HSP_BIT_SCORE=`echo $LINE | awk -F'|' '{print $5}'`
    HSP_IDENTITY=`echo $LINE | awk -F'|' '{print $6}'`
    HSP_ALIGN_LEN=`echo $LINE | awk -F'|' '{print $7}'`
    HSP_HIT_FRAME=`echo $LINE | awk -F'|' '{print $8}'`
    COVER=`echo $LINE | awk -F'|' '{print $9}'`
    ISNEW=`cat NEW_ACCESSIONS.txt | grep $ACC`

#echo "ACC: $ACC"
#echo "index: $INDEX"
#echo "e: $E"
#echo "desc: $DESCRIPTION"
#echo "bit score: $HSP_BIT_SCORE"
#echo "ident: $HSP_IDENTITY"
#echo "align len: $HSP_ALIGN_LEN"
#echo "hit frame: $HSP_HIT_FRAME"
#echo "cover: $COVER"

    [ -z $ISNEW ] && continue

    FASTA_NAME=$ACC
    FASTA_PATH=$OUTPUT_FOLDER/$FASTA_NAME.fasta
    MAFFT_FASTA_PATH=$MAFFT_OUTPUT_FOLDER/$FASTA_NAME.fasta
    MAFFT_RTF_PATH=$MAFFT_RTF_OUTPUT_FOLDER/$FASTA_NAME.rtf

    TOOWEAK=`awk 'BEGIN { print ('$E'>'$ECUTOFF')? "1" : "0" }'`

    if [[ "$TOOWEAK" = "1" ]]; then
        echo "Skipping weak result: $ACC E=$E"
        continue #skip weak results
    fi

    if [ ! -f $FASTA_PATH ]; then
        echo "Downloading $FASTA_NAME (e: $E)";
        HEADER="$ACC $DESCRIPTION $HSP_HIT_SCORE $E $HSP_IDENTITY $HSP_ALIGN_LEN"
        NOHEADER=`fetchaccession.py $ACC | fastafromtraces.sh | tail -n +1`

        echo -e ">$HEADER\n$NOHEADER" | translate.py -r $HSP_HIT_FRAME > $FASTA_PATH
        if [ $? -ne 0 ]; then
            echo "Translation failed for $ACC, skipping result!"
            continue;
        fi

        MAFFTFASTA=`cat $QUERY_FASTA $FASTA_PATH`
        echo "$MAFFTFASTA" > TMP_MAFFT_INPUT.fasta
        mafft TMP_MAFFT_INPUT.fasta > $MAFFT_FASTA_PATH 2> /dev/null

        ANALYSIS=`alignmentanalysis.py $MAFFT_FASTA_PATH`
        cat $MAFFT_FASTA_PATH | formatmafft.py --header "$ANALYSIS" -o $MAFFT_RTF_PATH

        IDENT=`echo -n "$ANALYSIS" | tail -n 1 | awk -F'\t' '{printf("%d%%",100*($4/$3))}'`
        COMMENT="[$DATE] $DESCRIPTION e:$E cover:$COVER ident:$IDENT"
        echo "$ACC //$COMMENT" >> $PREVIOUS
    fi
done

#clean up
rm -f CURRENT_ACCESSIONS.txt NEW_ACCESSIONS.txt PREVIOUS_SORTED.txt TMP_MAFFT_INPUT.fasta
