#!/usr/bin/bash

#seq 10 prints every number between 1 and 10
#extractfastaseq.py mafftout.txt -c counts the entries in your file
#tail -n +2 chops off the "1" element, we dont need to remap 1 to 1
seq `extractfastaseq.py aligned.fasta -c` | tail -n +2 |  while read line
do
    echo "remapping $line"
    #this line performs the remapping using -d $line to set the destination to the appropriate index
    #it outputs to a file named region_2.tsv
    remapregion.py regions.tsv aligned.fasta -d $line > region_$line.tsv
    formatmafft.py aligned.fasta -o region_$line.rtf -m region_$line.tsv -r $line
done
