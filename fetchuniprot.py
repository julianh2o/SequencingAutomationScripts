#!/usr/bin/python

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
Fetch region information from swissprot for coloring alignments

List all region information that swissprot has on a specific accession number
    fetchuniprot.py P35500

Filter out only the transmembrane regions
    fetchuniprot.py P35500 -m transmem

Filter out only the repeats
    fetchuniprot.py P35500 -m repeat

Accession number can be piped to this command
    echo "P35500" | fetchuniprot.py -m repeat
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
parser.add_argument("-m","--match",default=None, help='Match only a specific type of section (string match)');


args = parser.parse_args();

swissprotUrl = "http://www.uniprot.org/uniprot/%s.xml" % args.accession
r = requests.get(swissprotUrl);

ns = {'ns': 'http://uniprot.org/uniprot'};
rootNode = ET.fromstring(r.text);
sections = rootNode.findall("ns:entry/ns:feature",ns)
for section in sections:
    feature_type = section.attrib["type"]
    description = section.attrib["description"]
    start = section.find("ns:location/ns:begin",ns);
    end = section.find("ns:location/ns:end",ns);
    if (start is not None): start = start.attrib["position"];
    if (end is not None): end = end.attrib["position"];
    if (not args.match or args.match in feature_type):
        print("%s-%s\t%s\t%s\t" % (start,end,feature_type,description));

