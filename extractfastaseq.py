#!/usr/bin/env python
import sys
import argparse
from argparse import RawTextHelpFormatter

help="""
Extracts sequences from a fasta file

Outputs a single fasta entry with the header that contains a string
    extractfastaseq.py input.fasta TR3733

Outputs the first fasta sequence in the file
    extractfastaseq.py input.fasta -n 1

Outputs all entries EXCEPT the first fasta sequence in the file
    extractfastaseq.py input.fasta -r 1

Outputs the number of fasta sequences in the given file (may take time to run on large files)
    extractfastaseq.py input.fasta -c

Outputs one fasta entry per line with "###" delimiting the header from the sequence
eg: >TR3733|c0_g1_i5len=6587path=[###MDDSRVGSPNGSLDGGVI..
    extractfastaseq.py input.fasta -l

Iterates over each fasta sequence in your file
    extractfastaseq.py input.fasta -l | while read line
    do
        HEAD=`echo -n $line | awk -F'###' '{print $1}'
        SEQ=`echo -n $line | awk -F'###' '{print $2}'

        #do something with the fasta here
    done

Iterate through the numbers between 2 and the number of fasta sequences
(useful for remapping regions)
    seq `extractfastaseq.py input.fasta -c` | while read n
    do
        echo $n
    done
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('fasta', type=argparse.FileType('r'), nargs="?", help='the FASTA database')
parser.add_argument('accession', nargs="?", default=sys.stdin, help='The accession number to find')
parser.add_argument('-n', dest="extractByIndex", type=int, default=None,help='Return the sequence by index number')
parser.add_argument('-r', dest="removeByIndex", type=int, default=None,help='Returns a fasta file with one of the sequences removed')
parser.add_argument('-c', dest="count", action="store_const", const=True, default=False,help='Return the number of sequences in the fasta file')
parser.add_argument('-l', dest="list", action="store_const", const=True, default=False,help='List the fasta output one per line: [header]###[sequence]')

def nth(iterable, n, default=None):
    "Returns the nth item or a default value"
    return next(islice(iterable, n, None), default)

def fasta(f):
    i = 0;
    accum = None;
    for line in f:
        line = line.strip()
        i+=1;
        if (line.startswith(">")):
            if (accum is not None):
                yield accum;
            accum = (line,"");
            continue;
        if (accum is not None):
            accum = (accum[0],accum[1]+line);

    if (accum is not None): yield accum;


def main():
    args = parser.parse_args()
    if (args.list):
        for head,seq in fasta(args.fasta):
            print("%s###%s" % (head,seq));
        sys.exit(0)
    if (args.count):
        i=1
        for head,seq in fasta(args.fasta):
            i+=1
        print(i-1)
        sys.exit(0)
    if (args.removeByIndex):
        i=1
        for head,seq in fasta(args.fasta):
            if (i != args.removeByIndex):
                print("%s\n%s" % (head,seq));
            i+=1
        sys.exit(0)
    if (args.extractByIndex):
        i=1
        for head,seq in fasta(args.fasta):
            if (i == args.extractByIndex):
                print("%s\n%s" % (head,seq));
                sys.exit(0)
            i+=1
    if (args.accession):
        accession = args.accession
        if not isinstance(accession,str):
            accession = accession.read()
        accession = accession.strip()
        for head,seq in fasta(args.fasta):
            if (accession in head):
                print("%s\n%s" % (head,seq));
                sys.exit(0)

if __name__ == "__main__":
    main()
