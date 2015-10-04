#!/usr/bin/env python

import time
import sys
import os.path
import argparse
import colorama
import re
from formatmafft import parseRegionInfo
from formatmafft import parseFasta

help="""
Converts a region definition file from one reference to another
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('regions', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='The region file to translate')
parser.add_argument('fasta', type=argparse.FileType('r'), help='The aligned fasta file to use')
parser.add_argument('-s', '--source', type=int,default=1,help='The index of the old reference sequence (first sequence is 1)')
parser.add_argument('-d', '--dest', type=int,default=2,help='The index of the new reference sequence (first sequence is 1)')

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

