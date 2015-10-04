#!/usr/bin/env python
import sys
import argparse
from argparse import RawTextHelpFormatter

help="""
Extracts a sequence from a fasta file by comparing to the header
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('fasta', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the FASTA database')

args = parser.parse_args()

i = 0;
accum = None;
header = sys.stdin.read().strip();
for line in args.fasta:
    i+=1;
    if (line.startswith(">")):
        if (accum is not None):
            break;
        if (header in line):
            accum = "";
    if (accum is not None):
        accum += line;

if (accum is not None): print(accum,end="");
