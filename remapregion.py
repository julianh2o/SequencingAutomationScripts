#!/usr/bin/env python

import time
import sys
import os.path
import argparse
import colorama
import re
from formatmafft import parseRegionInfo
from formatmafft import parseFasta
from argparse import RawTextHelpFormatter

help="""
Converts a region definition file from one reference to another

Remaps the region to the second sequence using the first as a reference
    remapregion.py regions.tsv aligned.fasta > newreference.tsv

Remmaps the region to the third sequence using the first as a reference
    remapregion.py regions.tsv aligned.fasta -d 3 > newreference.tsv

Remaps the region to the first sequence's frame using the third sequence as a reference
    remapregion.py regions.tsv aligned.fasta -d 1 -s 3 > newreference.tsv

Remaps from the first sequence to the second
Then generates newreference.fasta that excludes the first sequence
Finally formats the alignment using the new regions (the same amino acids will be highlighted, but the original reference is not used)
    remapregion.py regions.tsv aligned.fasta > newregions.tsv
    extractfastaseq.py aligned.fasta -l | tail -n +2 | sed 's/###/\\n/' > newreference.fasta
    formatmafft.py newreference.fasta -o output.rtf -m newregions.txt

"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('regions', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='The region file to translate')
parser.add_argument('fasta', type=argparse.FileType('r'), help='The aligned fasta file to use')
parser.add_argument('-s', '--source', type=int,default=1,help='[default: 1] The index of the old reference sequence (first sequence is 1)')
parser.add_argument('-d', '--dest', type=int,default=2,help='[default: r] The index of the new reference sequence (first sequence is 1)')

class Bunch(dict):
    __getattr__, __setattr__ = dict.get, dict.__setitem__

def findIndex(index,sequence):
    offset = 0;
    for i,c in enumerate(sequence):
        if (c == "-"): offset += 1;
        if i-offset >= index:
            return i
    return None

def translateLocation(location,source,dest):
    index = findIndex(location,source);
    offset = dest[:index].count("-")
    if (index == offset): return None
    newloc = index-offset
#print(location,index,offset,newloc);
    return newloc

def main():
    args = parser.parse_args()

    fasta = parseFasta(args.fasta);
    regions = parseRegionInfo(args.regions);

    source = fasta[args.source-1][1]
    dest = fasta[args.dest-1][1]

    newregions = [];
    for region in regions:
        start,end,feature_type,description,color_number = region;
        tstart = translateLocation(start,source,dest);
        tend = translateLocation(end,source,dest);
        if (tstart == None):
            tstart = 0;
        if (tend == None):
            rangeinfo = "-";
        else:
            rangeinfo = "%s-%s" % (tstart, tend)
        if (color_number == None): color_number = "";
        print("%s\t%s\t%s\t%s" % (rangeinfo,feature_type,description,color_number));


if __name__ == "__main__":
    main()

