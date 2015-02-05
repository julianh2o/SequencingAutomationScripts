http://www.ncbi.nlm.nih.gov/protein/263970?report=fasta

#Find all stored sequences that look like our opsin (in calanus finmarchigus)
./blast.py --entrez="txid6837[orgn]" tblastn tsa_nt opsin.fasta > out.blast

#Show the first 5
./viewblast.py out.blast list -n 5

#Get the nucleotide sequence of one of the records that look like our opsin
./viewblast.py out.blast contig 1

#Translate the first result to an amino acid sequence using the best frame
./viewblast.py out.blast contig 1 | ./translate.py

#Extract the reading frame
./viewblast.py out.blast contig 1 | ./translate.py -f

#Blast the experiment's sequence back against the database to verify that we got the right one
./viewblast.py out.blast contig 1 | ./translate.py -f | ./blast.py tblastn nr | ./viewblast.py list -n 5

