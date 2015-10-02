#!/usr/bin/env python

import requests
import urllib
import re
import time
import sys
import os.path
import argparse

parser = argparse.ArgumentParser(description='Performs a blast using NCBI\'s web tool')
parser.add_argument('program', help='program: blastn blastp blastx tblastn tblastx')
parser.add_argument('database', help='database: nr refseq_rna refseq_genomic chromosome est gss htgs pat pdb alu dbsts Whole_Genome_Shotgun_contigs tsa_nt rRNA_typestrains/prokaryotic_16S_ribosomal_RNA')
parser.add_argument('fasta', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the FASTA sequence(s)')

def extractBlastInfo(responseText):
    info = {};
    for m in re.finditer("<!--.*?QBlastInfoBegin(.*?)QBlastInfoEnd.*?-->",responseText,re.DOTALL):
        content = m.group(1)
        for mm in re.finditer("\W+(.*?)\W*?=\W*?(.*)",content):
            info[mm.group(1)] = mm.group(2);
    return info;

def startBlast(program,database,query):
    params = {
        "program" : program,
        "database" : database,
        "sequence" : urllib.parse.quote(query),
        "org":"dmel",
        "tax" : "drosophila",
        "outputformat":"XML"
    }

    r = requests.post("http://flybase.org/blast/index.html",params=params,allow_redirects=False);
    url = r.headers["location"];
    return url;

def waitForCompletion(jobUrl):
    i = 0;
    while(True):
        r = requests.get(jobUrl);
        if ("text/xml" in r.headers["content-type"]):
            return r.text;

        i += 1
        sleeptime = 5*i;
        if (sleeptime > 30): sleeptime = 30;
        time.sleep(sleeptime);

def main():
    args = parser.parse_args()

    fasta = args.fasta.read()
    sys.stderr.write("Program: %s\nDatabase: %s\nQuery: \n%s" % (args.program,args.database,fasta));
    jobUrl = startBlast(args.program,args.database,fasta);
    sys.stderr.write("JobURL: %s\n" % jobUrl);
    time.sleep(5);
    results = waitForCompletion(jobUrl);

    with open("_latest.blast","w") as f:
        f.write(results);
    print(results.strip());
    sys.stderr.write("View in browser: %s\n" % jobUrl);

main();
