#!/usr/bin/env python
import sys
import argparse
from argparse import RawTextHelpFormatter

help="""
Extracts a sequence from a fasta file by comparing to the header
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('fasta', type=argparse.FileType('r'), nargs="?", help='the FASTA database')
parser.add_argument('accession', nargs="?", default=sys.stdin, help='The accession number to find')
parser.add_argument('-n', dest="extractByIndex", type=int, default=None,help='Return the sequence by index number')
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
    if (args.count):
        i=1
        for head,seq in fasta(args.fasta):
            i+=1
        print(i)
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
