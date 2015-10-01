#!/usr/local/bin/python3
import sys
import argparse

parser = argparse.ArgumentParser(description='Extracts a sequence from a fasta file by comparing to the header')
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
