#!/usr/bin/env python

import requests
import urllib
import re
import time
import sys
import os.path
import argparse
import xml.etree.ElementTree as ET
from argparse import RawTextHelpFormatter

help="""
This script lets you extract results and information from blast results stored as XML

The format parameter can be any strings. The following sequences can be used to insert record information:
    {Hit_num} {Hit_id} {Hit_def} {Hit_accession} {Hit_len} {Hit_hsps/Hsp/Hsp_num} {Hit_hsps/Hsp/Hsp_bit-score} {Hit_hsps/Hsp/Hsp_score}
    {Hit_hsps/Hsp/Hsp_evalue} {Hit_hsps/Hsp/Hsp_query-from} {Hit_hsps/Hsp/Hsp_hit-from} {Hit_hsps/Hsp/Hsp_hit-to}
    {Hit_hsps/Hsp/Hsp_query-frame} {Hit_hsps/Hsp/Hsp_hit-frame} {Hit_hsps/Hsp/Hsp_identity} {Hit_hsps/Hsp/Hsp_positive}
    {Hit_hsps/Hsp/Hsp_gaps} {Hit_hsps/Hsp/Hsp_align-len} {Hit_hsps/Hsp/Hsp_qseq} {Hit_hsps/Hsp/Hsp_hseq} {Hit_hsps/Hsp/Hsp_midline}

of particular use are
    {Hit_num} - replaced with the result number starting at 1
    {Hit_id} - replaecd with the hit ID. This is often the accession number or contig information but is not necessarily consistent
    {Hit_def} - replaecd with most of the remaining fasta header, again, this seems to be dependent on the database
    {Hit_accession} - replaced with the accession number as known by the database, if blasting against NCBI this is a GB accession
    {e} - replaced with the evalue of the first "sub hit"
    {cover} - replaced with the query coverage of all "sub hits"

Show the accession number and the evalue of the top 5 results separated by commas
    viewblast.py list out.blast -n 5 -f '{Hit_accession},{e}'

Show the hit number and the accession number of the top 5 results separated by space
    viewblast.py list out.blast -n 5 -f '{Hit_num} {Hit_accession} {e}'

The tab character has special meaning here. Instead of inserting a tab, this creates a equalized column view for display
    viewblast.py list out.blast -n 5 -f '{Hit_num}\t{Hit_accession}\t{e}'

Returns the "{Hit_accession}" of the first result
    viewblast.py info 1 out.blast

Returns only the fasta header "Hit_def" field of the first result
    viewblast.py info 1 out.blast -f '{Hit_def}'

Returns a list of all result accession numbers from the blast
    viewblast.py list out.blast -f '{Hit_accession}'

This only works with an NCBI blast, it uses the result Hit_accession number to then look up the contig on GenBank
    viewblast.py contig 1 out.blast

Only show results who's E-value is less than or equal to the provided evalue
    viewblast.py list out.blast -e '1e-50'

Output the list of HSPs from the 2nd hit
    viewblast.py hsp 2 result.blast

Format the list of HSPs from the 2nd hit
    viewblast.py hsp 2 result.blast -f '{Hsp_num} {Hsp_score} {Hsp_evalue} {Hsp_hit-frame}'

Get help on the various options
    viewblast.py info -h
    viewblast.py contig -h
    viewblast.py hsp -h
    viewblast.py list -h

You can of course, pipe in a blast result directly. This command gets the first 5 results of the blast
tblastn -db tsa_nt -query query.txt -remote | viewblast.py list -n 5
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
command_subparser = parser.add_subparsers(help="View mode")

info_parser = command_subparser.add_parser("info",help="View a single result");
info_parser.add_argument("index", type=int);
info_parser.add_argument("-f","--format", help="Format of the response", default="{Hit_accession}");
info_parser.add_argument('blast', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the BLAST results file in XML')

hsp_parser = command_subparser.add_parser("hsp",help="List Hsp information from hit");
hsp_parser.add_argument("index", type=int);
hsp_parser.add_argument("-f","--format", help="Format of the response", default="[{Hsp_num}] score:{Hsp_score} e:{Hsp_evalue} frame:{Hsp_hit-frame}");
hsp_parser.add_argument('blast', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the BLAST results file in XML')

contig_parser = command_subparser.add_parser("contig",help="View a result's contig");
contig_parser.add_argument("index", type=int);
contig_parser.add_argument('blast', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the BLAST results file in XML')

list_parser = command_subparser.add_parser("list",help="List the results");
list_parser.add_argument("-f","--format", help="Format of the response", default="{Hit_id}\t{Hit_def}\t{Hit_hsps/Hsp/Hsp_evalue}\t{Hit_accession}");
list_parser.add_argument("-n","--max", type=int, help="Number of records to show");
list_parser.add_argument("-e","--ecutoff", type=str, help="Cut off at at a given E value");
list_parser.add_argument('blast', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the BLAST results file in XML')


def fetchContig(accession):
    r = requests.get("http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&dopt=fasta&val=%s" % accession);
    return r.text.strip();

def formatTable(output):
    lengths = [0] * len(output[0]);
    for row in output:
        i = 0;
        for d in row:
            d = str(d);
            if (len(d) > lengths[i]): lengths[i] = len(d);
            i+=1;
    
    for row in output:
        i = 0;
        for d in row:
            d = str(d);
            print(d.ljust(lengths[i]+3),end='');
            i += 1;
        print();

def doFormat(formatString,root,othermapping):
    for k,v in othermapping.items():
        formatString = formatString.replace("{"+k+"}",v);
    def replacementFunction(match):
        return root.find(match.group(1)).text;
    formatString = re.sub("{(.*?)}",replacementFunction,formatString);
    return formatString;

def addCoverage(coverage,seq):
    inserted = False;
    for i,cov in enumerate(coverage):
        startsBefore = seq[0] < cov[0];

        if (startsBefore):
            inserted = True;
            coverage.insert(i,seq);
            break;
    if (not inserted): coverage.append(seq);

    i = -1;
    while(i < len(coverage) - 2):
        i += 1;

        a = coverage[i];
        b = coverage[i+1];

        if (a[1] < b[0]): continue #no coalescing required, disparate
        if (a[1] > b[1]):
            coverage.pop(i+1);
            i -= 1;
            continue;

        coverage[i] = (a[0],b[1]);
        coverage.pop(i+1);
        i -= 1;

def calculateCoverage(hsps):
    coverage = [];
    for hsp in hsps:
        seq = (int(hsp.find("Hsp_query-from").text),int(hsp.find("Hsp_query-to").text))
        addCoverage(coverage,seq);
    
    coverageCount = 0;
    for seq in coverage:
        diff = seq[1] - seq[0] + 1;
        coverageCount += diff;
    
    return coverageCount;

def calculateMultiframe(ecutoff,hsps):
    frameOfHit = None;
    for hsp in hsps:
        e = float(hsp.find("Hsp_evalue").text);
        hitFrame = hsp.find("Hsp_hit-frame").text;
        if (e <= float(ecutoff)):
            if (frameOfHit is not None):
                if (frameOfHit != hitFrame):
                    return True;
            frameOfHit = hitFrame;
    return False;

def formatAlignment(args,root,alignment):
    queryLength = int(root.find("BlastOutput_iterations/Iteration/Iteration_query-len").text);
    isMultiframe = calculateMultiframe(args.ecutoff,alignment.findall("Hit_hsps/Hsp")) if hasattr(args,"ecutoff") else False
    multiframe = "MULTIFRAME" if isMultiframe else "";
    hspCount = len(alignment.find("Hit_hsps").getchildren());
    coverCount = calculateCoverage(alignment.findall("Hit_hsps/Hsp"));
    coverage = "{:.0f}".format(100*float(coverCount) / float(queryLength)) + "%";
    hitid = alignment.find("Hit_id").text;
    hitdef = alignment.find("Hit_def").text;
    accession = alignment.find("Hit_accession").text;
    evalue = float(alignment.find("Hit_hsps/Hsp/Hsp_evalue").text);
    fstring = args.format.replace("\\t","\t");
    otherParams = {
        "e":str(evalue),
        "cover":coverage,
        "hspCount":str(hspCount),
        "multiframe":multiframe
    }
    fstring = doFormat(fstring,alignment,otherParams);
    return fstring

def hsp(args):
    results = args.blast.read();

    if (len(results.strip()) == 0): exit(1);
    rootNode = ET.fromstring(results);
    alignments = rootNode.findall("BlastOutput_iterations/Iteration/Iteration_hits/Hit")
    if (len(alignments) == 0):
        sys.exit(0)
    i = 1;
    for alignment in alignments:
        if (i == args.index):
            hsps = alignment.findall("Hit_hsps/Hsp");
            for hsp in hsps:
                print(doFormat(args.format.replace("\\t","\t"),hsp,{}));
            exit(0);
        i += 1
hsp_parser.set_defaults(func=hsp);

            
def info(args):
    results = args.blast.read();

    if (len(results.strip()) == 0): exit(1);
    rootNode = ET.fromstring(results);
    alignments = rootNode.findall("BlastOutput_iterations/Iteration/Iteration_hits/Hit")
    if (len(alignments) == 0):
        sys.exit(0)
    i = 1;
    for alignment in alignments:
        if (i == args.index):
            print(formatAlignment(args,rootNode,alignment));
            exit(0);
        i += 1
info_parser.set_defaults(func=info);

def contig(args):
    results = args.blast.read();

    if (len(results.strip()) == 0): exit(1);
    rootNode = ET.fromstring(results);
    alignments = rootNode.findall("BlastOutput_iterations/Iteration/Iteration_hits/Hit")
    if (len(alignments) == 0):
        sys.exit(0)
    i = 1;
    for alignment in alignments:
        if (i == args.index):
            accession = alignment.find("Hit_accession").text;
            contig = fetchContig(accession);
            print(contig);
            exit(0);
        i += 1
contig_parser.set_defaults(func=contig);

def list(args):
    results = args.blast.read();

    if (len(results.strip()) == 0): exit(1);
    rootNode = ET.fromstring(results);
    alignments = rootNode.findall("BlastOutput_iterations/Iteration/Iteration_hits/Hit")
    if (len(alignments) == 0):
        sys.exit(0)

    output = [];
    i = 0;
    for alignment in alignments:
        evalue = float(alignment.find("Hit_hsps/Hsp/Hsp_evalue").text);

        if (args.ecutoff):
            if (not float(evalue) <= float(args.ecutoff)): continue;

        fstring = formatAlignment(args,rootNode,alignment);
        output.append(fstring.split("\t"));

        if (i == args.max): break;
        i+=1;
    if (len(output)): formatTable(output);
list_parser.set_defaults(func=list);


def main():
    args = parser.parse_args();
    args.func(args);

if __name__ == "__main__":
    main();
