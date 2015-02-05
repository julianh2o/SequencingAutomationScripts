Setup
-----
First of all, these scripts require python3. Download it if you don't already have access to it: https://www.python.org/downloads/

Furthermore, these scripts require the requests module, your best bet is to install them via pip: 
python3 -m pip install requests

or in some cases just pip install requests

Examples
--------
Note: all of these commands are formatted to run on a unix/linux machine.. for a windows friendly command line example see below...

http://www.ncbi.nlm.nih.gov/protein/263970?report=fasta

#Find all stored sequences that look like our opsin (in Calanus finmarchicus)
./blast.py --entrez="txid6837[orgn]" tblastn tsa_nt opsin.fasta > out.blast

#Show the first 5
./viewblast.py list -n 5 out.blast

#Get the nucleotide sequence of one of the records that look like our opsin
./viewblast.py contig 1 out.blast

#Translate the first result to an amino acid sequence using the best frame
./viewblast.py contig 1 out.blast | ./translate.py

#Extract the reading frame
./viewblast.py contig 1 out.blast | ./translate.py -f

#Blast the experiment's sequence back against the database to verify that we got the right one
./viewblast.py contig 1 out.blast | ./translate.py -f | ./blast.py tblastn nr | ./viewblast.py list -n 5


#All together now..
./blast.py --entrez="txid6837[orgn]" tblastn tsa_nt opsin.fasta | ./viewblast.py contig 1 | ./translate.py -f | ./blast.py tblastn nr | ./viewblast.py list -n 5


#On a windows machine
blast.py --entrez=txid6837[orgn] tblastn tsa_nt opsin.fasta|viewblast.py contig 1|translate.py -f|blast.py tblastn nr|viewblast.py list -n 5

