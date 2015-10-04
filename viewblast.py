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
View the XML blast results
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
command_subparser = parser.add_subparsers(help="View mode")

info_parser = command_subparser.add_parser("info",help="View a single result");
info_parser.add_argument("index", type=int);
info_parser.add_argument("-f","--format", help="Format of the response", default="{Hit_accession}");
info_parser.add_argument('blast', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the BLAST results file in XML')

contig_parser = command_subparser.add_parser("contig",help="View a result's contig");
contig_parser.add_argument("index", type=int);
contig_parser.add_argument('blast', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the BLAST results file in XML')

list_parser = command_subparser.add_parser("list",help="List the results");
list_parser.add_argument("-f","--format", help="Format of the response", default="{Hit_id}\t{Hit_def}\t{Hit_hsps/Hsp/Hsp_evalue}\t{Hit_accession}");
list_parser.add_argument("-n","--max", type=int, help="Number of records to show");
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

def doFormat(formatString,alignment):
    formatString = formatString.replace("{e}","{Hit_hsps/Hsp/Hsp_evalue}");
    def replacementFunction(match):
        return alignment.find(match.group(1)).text;
    formatString = re.sub("{(.*?)}",replacementFunction,formatString);
    return formatString;

            
def info(args):
    results = args.blast.read();

    rootNode = ET.fromstring(results);
    alignments = rootNode.findall("BlastOutput_iterations/Iteration/Iteration_hits/Hit")
    i = 1;
    for alignment in alignments:
        hitid = alignment.find("Hit_id").text;
        hitdef = alignment.find("Hit_def").text;
        accession = alignment.find("Hit_accession").text;
        evalue = float(alignment.find("Hit_hsps/Hsp/Hsp_evalue").text);
        if (i == args.index):
            fstring = args.format;
            fstring = doFormat(fstring,alignment);
            print(fstring);
            exit(0);
        i += 1
info_parser.set_defaults(func=info);

def contig(args):
    results = args.blast.read();

    rootNode = ET.fromstring(results);
    alignments = rootNode.findall("BlastOutput_iterations/Iteration/Iteration_hits/Hit")
    i = 1;
    for alignment in alignments:
        hitid = alignment.find("Hit_id").text;
        hitdef = alignment.find("Hit_def").text;
        accession = alignment.find("Hit_accession").text;
        evalue = float(alignment.find("Hit_hsps/Hsp/Hsp_evalue").text);
        if (i == args.index):
            contig = fetchContig(accession);
            print(contig);
            exit(0);
        i += 1
contig_parser.set_defaults(func=contig);

def list(args):
    results = args.blast.read();
    rootNode = ET.fromstring(results);
    alignments = rootNode.findall("BlastOutput_iterations/Iteration/Iteration_hits/Hit")

    output = [];
    i = 0;
    for alignment in alignments:
        hitid = alignment.find("Hit_id").text;
        hitdef = alignment.find("Hit_def").text;
        accession = alignment.find("Hit_accession").text;
        evalue = float(alignment.find("Hit_hsps/Hsp/Hsp_evalue").text);
        fstring = args.format;
        fstring = doFormat(fstring,alignment);
        output.append(fstring.split("\t"));
        if (i == args.max): break;
        i+=1;
    formatTable(output);
list_parser.set_defaults(func=list);


def main():
    args = parser.parse_args();
    args.func(args);

main();
