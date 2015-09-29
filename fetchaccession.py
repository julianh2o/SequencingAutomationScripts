#!/usr/local/bin/python3
import requests
import sys
import time

def fetchContig(accession):
    r = None;
    while (r == None or "Resource temporarily unavailable" in r.text):
        if (r is not None): time.sleep(3)
        r = requests.get("http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&dopt=fasta&val=%s" % accession);
    return r.text.strip();

print(fetchContig(sys.argv[1]));

