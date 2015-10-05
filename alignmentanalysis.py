#!/usr/bin/env python

import time
import sys
import os.path
import argparse
import colorama
import re
from argparse import RawTextHelpFormatter
from extractfastaseq import fasta
from formatmafft import isSimilar
from viewblast import formatTable

help="""
Reports useful statistics on a mafft alignment
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('fasta', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the MAFFT aligned FASTA sequence(s)')
parser.add_argument('-p', "--percent", action="store_const", const=True, default=False,help='Displays the values as a percent of the length of the sequence')

def asPercent(a,b):
    return "{0:.2f}%".format(a/b * 100) 

def main():
    args = parser.parse_args()

    stats = [];
    reference = None;
    for head,seq in fasta(args.fasta):
        if reference is not None:
            matching = 0;
            similar = 0;
            different = 0;
            for i,c in enumerate(seq):
                if c == reference[i]:
                    matching += 1;
                elif isSimilar(c,reference[i]):
                    similar += 1;
                else:
                    different += 1;

        if (head.__len__() > 20): head = head[:20]
        seqlen = seq.__len__()
        if (args.percent):
            sstats = {
                "head":head,
                "alignedlength":seqlen,
                "sequencelength":asPercent(seqlen - seq.count("-"),seqlen)
            }
            if reference is not None:
                sstats.update({
                    "matching":asPercent(matching,seqlen),
                    "similar":asPercent(similar,seqlen),
                    "different":asPercent(different,seqlen),
                });
        else:
            sstats = {
                "head":head,
                "alignedlength":seqlen,
                "sequencelength":seqlen - seq.count("-"),
            }
            if reference is not None:
                sstats.update({
                    "matching":matching,
                    "similar":similar,
                    "different":different,
                });

        stats.append(sstats);
        if reference is None:
            reference = seq

    tabular = []
    fields = ["head","alignedlength","sequencelength","matching","similar","different"]
    tabular.append(fields)
    for stat in stats:
        row = []
        for f in fields:
            a = ""
            if f in stat:
                a = stat[f]
            row.append(a)
        tabular.append(row)
    formatTable(tabular)
        

if __name__ == "__main__":
    main()
