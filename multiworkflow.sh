#!/bin/bash

LAB_TRINITY=~/blastdemo/databases/FASTA_2015/Lab_Trinity.fasta
NEOFL_TRINITY=~/blastdemo/databases/FASTA_2015/Neofl_Trinity.fasta

echo -e "Accession\tResult Num\tResult evalue\tresult identifier\tflybase1 evalue\tflybase1 accession\tflybase1 ncbi accession\tflybase1 ncbi header\tflybase2 evalue\tflybase2 accession\tflybase2 ncbi accession\tflybase2 ncbi header"

rm -rf WorkDir || true
mkdir WorkDir
cd WorkDir

for ACCESSION in `cat ../accessions.txt`;
do
	mkdir $ACCESSION
	cd $ACCESSION

	#Download sequence by accession number
	fetchaccession.py $ACCESSION > 1_HSP_SEQUENCE.fasta

	#Run the initial blast
	>&2 echo "Running initial blast: $ACCESSION"
	tblastn -db $LAB_TRINITY -query 1_HSP_SEQUENCE.fasta -outfmt 5 > 2_TBLASTN_RESULT.blast

	#Get the first result from the blast
	#viewblast.py list 2_TBLASTN_RESULT.blast -f '{Hit_id}'
	viewblast.py list 2_TBLASTN_RESULT.blast -f '{Hit_id}' | head -n 20 | nl | while read line; do
		RESULT_NUM=`echo $line | awk -F' ' '{print $1}'`
		RESULT_EVALUE=`viewblast.py info $RESULT_NUM 2_TBLASTN_RESULT.blast -f '{e}'`
		RESULT_IDENTIFIER=`viewblast.py info $RESULT_NUM 2_TBLASTN_RESULT.blast -f '{Hit_id}'`
		RESULT_CONTIG=`cat 2_TBLASTN_RESULT.blast | viewblast.py info $RESULT_NUM -f '{Hit_def}'`

		mkdir "RESULT_$RESULT_NUM"
		cd "RESULT_$RESULT_NUM"

		echo "$RESULT_CONTIG" > a_RESULT_CONTIG.txt
		cat $LAB_TRINITY | extractfastacontig.py $RESULT_CONTIG > b_ORIGINAL_SEQ.txt

		#Get the full contig from the original fasta file and translate it
		cat $LAB_TRINITY | extractfastacontig.py $RESULT_CONTIG | translate.py > 3_TRANSLATED_RESULT.fasta

		>&2 echo "Running flybase blast: $RESULT_NUM $RESULT_IDENTIFIER"
		flybase.py blastp org/translation 3_TRANSLATED_RESULT.fasta > 4_FLYBASE_RESULT.blast 2> /dev/null

		FLYBASE_EVALUE1=`viewblast.py info 1 4_FLYBASE_RESULT.blast -f '{e}'`
		FLYBASE_ACCESSION1=`viewblast.py info 1 4_FLYBASE_RESULT.blast -f '{Hit_accession}'`
		FLYBASE_NCBI1=`viewblast.py info 1 4_FLYBASE_RESULT.blast -f '{Hit_def}' | perl -nle 'print "$1" if (/GB_protein:(.*?)[,; ]/)'`

		FLYBASE_HEADER1=""
		if [[ ! -z  $FLYBASE_NCBI1  ]]
		then
			FLYBASE_HEADER1=`fetchaccession.py $FLYBASE_NCBI1 | head -n 1`
		fi

		FLYBASE_EVALUE2=`viewblast.py info 2 4_FLYBASE_RESULT.blast -f '{e}'`
		FLYBASE_ACCESSION2=`viewblast.py info 2 4_FLYBASE_RESULT.blast -f '{Hit_accession}'`
		FLYBASE_NCBI2=`viewblast.py info 2 4_FLYBASE_RESULT.blast -f '{Hit_def}' | perl -nle 'print "$1" if (/GB_protein:(.*?)[,; ]/)'`

		FLYBASE_HEADER2=""
		if [[ ! -z  $FLYBASE_NCBI2  ]]
		then
			FLYBASE_HEADER2=`fetchaccession.py $FLYBASE_NCBI2 | head -n 1`
		fi

		echo -e "$ACCESSION\t$RESULT_NUM\t$RESULT_EVALUE\t$RESULT_IDENTIFIER\t$FLYBASE_EVALUE1\t$FLYBASE_ACCESSION1\t$FLYBASE_NCBI1\t$FLYBASE_HEADER1\t$FLYBASE_EVALUE2\t$FLYBASE_ACCESSION2\t$FLYBASE_NCBI2\t$FLYBASE_HEADER2"

		cd ..
	done

	cd ..
done
