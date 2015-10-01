#!/bin/bash

LAB_TRINITY=~/blastdemo/databases/FASTA_2015/Lab_Trinity.fasta
NEOFL_TRINITY=~/blastdemo/databases/FASTA_2015/Neofl_Trinity.fasta

ACCESSION=`head -n 1 accessions.txt`

rm -rf $ACCESSION || true
mkdir $ACCESSION
cd $ACCESSION

#Download sequence by accession number
fetchaccession.py $ACCESSION > 1_HSP_SEQUENCE.fasta

#Run the initial blast
tblastn -db $LAB_TRINITY -query 1_HSP_SEQUENCE.fasta -outfmt 5 > 2_TBLASTN_RESULT.blast

#Get the first result from the blast
RESULT_NUM=1
RESULT_EVALUE=`viewblast.py info $RESULT_NUM 2_TBLASTN_RESULT.blast -f '{e}'`
RESULT_IDENTIFIER=`viewblast.py info $RESULT_NUM 2_TBLASTN_RESULT.blast -f '{Hit_def}' | awk -F' ' '{print $1}'`
RESULT_CONTIG=`cat 2_TBLASTN_RESULT.blast | viewblast.py info $RESULT_NUM -f '{Hit_def}'`

#Get the full contig from the original fasta file and translate it
cat $LAB_TRINITY | extractfastacontig.py $RESULT_CONTIG | translate.py > 3_TRANSLATED_RESULT.fasta

flybase.py blastp org/translation 3_TRANSLATED_RESULT.fasta > 4_FLYBASE_RESULT.blast

FLYBASE_EVALUE1=`viewblast.py info 1 4_FLYBASE_RESULT.blast -f '{e}'`
FLYBASE_ACCESSION1=`viewblast.py info 1 4_FLYBASE_RESULT.blast -f '{Hit_accession}'`
FLYBASE_NCBI1=`viewblast.py info 1 4_FLYBASE_RESULT.blast -f '{Hit_def}' | perl -nle 'print "$1" if (/GB_protein:(.*?),/)'`
FLYBASE_HEADER1=`fetchaccession.py $FLYBASE_NCBI1 | head -n 1`

FLYBASE_EVALUE2=`viewblast.py info 1 4_FLYBASE_RESULT.blast -f '{e}'`
FLYBASE_ACCESSION2=`viewblast.py info 2 4_FLYBASE_RESULT.blast -f '{Hit_accession}'`
FLYBASE_NCBI2=`viewblast.py info 2 4_FLYBASE_RESULT.blast -f '{Hit_def}' | perl -nle 'print "$1" if (/GB_protein:(.*?),/)'`
FLYBASE_HEADER2=`fetchaccession.py $FLYBASE_NCBI2 | head -n 1`

echo -e "Accession\tResult Num\tResult evalue\tresult identifier\tflybase1 evalue\tflybase1 accession\tflybase1 ncbi accession\tflybase1 ncbi header\tflybase2 evalue\tflybase2 accession\tflybase2 ncbi accession\tflybase2 ncbi header"
echo -e "$ACCESSION\t$RESULT_NUM\t$RESULT_EVALUE\t$RESULT_IDENTIFIER\t$FLYBASE_EVALUE1\t$FLYBASE_ACCESSION1\t$FLYBASE_NCBI1\t$FLYBASE_HEADER1\t$FLYBASE_EVALUE2\t$FLYBASE_ACCESSION2\t$FLYBASE_NCBI2\t$FLYBASE_HEADER2"


