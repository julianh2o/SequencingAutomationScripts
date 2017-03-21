#!/bin/bash

DATABASE=/home/lenzlab/home-data-lenzlab/Alaska/Collection_Timeline/Assemblies/Cap3_Combos/Wk0-7/wk0-7_cap.fasta
QUERIES=/home/lenzlab/opt-data-lenzlab/Illumina_2011/Calfi_GoM_Assembly/Longest/c-finmarchicus_96k_ref.fasta

mkdir -p blasts
Q=1
extractfastaseq.py $QUERIES -l | while read line
do
	echo "$line" | sed "s/###/\n/" > query.fasta
	HEADER=`head -n 1 query.fasta`

	BLAST_FILE="blasts/$Q.blast"
	if [[ ! -f $BLAST_FILE ]]; then
		tblastx -db $DATABASE -query query.fasta -outfmt 5 > $BLAST_FILE
	fi

	HITS=`viewblast.py list $BLAST_FILE -e '1e-100' -f '{Hit_def}'`
	if [[ -z  $HITS  ]]; then
		E=`viewblast.py info 1 $BLAST_FILE -f '{e}'`
		echo "[$Q] NO HITS [$E] $HEADER"
	else
		echo "[$Q]" $HEADER '####' `xargs "$HITS"`
	fi

	((Q++))
done
