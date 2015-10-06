#!/usr/bin/env python

import requests
import sys
import time
import argparse
from argparse import RawTextHelpFormatter

help="""
Fetch the fastsa entry from NCBI for the given accession number

Loads the fasta file from NCBI for a given accession number
    fetchaccession.py AAN12345

The accession number can also be provided by piping
    echo "AAN12345" | fetchaccession.py
"""

#http://stackoverflow.com/questions/19966519/argparse-optional-stdin-argument
class StreamType(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        try:
            setattr(namespace, self.dest, values.readline().strip())
        except AttributeError:
            setattr(namespace, self.dest, values)

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument("accession", nargs='?', default=sys.stdin, action=StreamType)

args = parser.parse_args()

def fetchContig(accession):
    r = None;
    while (r == None or "Resource temporarily unavailable" in r.text):
        if (r is not None): time.sleep(3)
        r = requests.get("http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&dopt=fasta&val=%s" % accession);
    return r.text.strip();

print(fetchContig(args.accession));

