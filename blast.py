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
parser.add_argument('--rid', help='retrieve results for an existing blast by rid')
parser.add_argument('--entrez', help='an entrez_query to filter results')
parser.add_argument('fasta', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the FASTA sequence(s)')

def extractBlastInfo(responseText):
    info = {};
    for m in re.finditer("<!--.*?QBlastInfoBegin(.*?)QBlastInfoEnd.*?-->",responseText,re.DOTALL):
        content = m.group(1)
        for mm in re.finditer("\W+(.*?)\W*?=\W*?(.*)",content):
            info[mm.group(1)] = mm.group(2);
    return info;

def startBlast(program,database,entrez,query):
    params = {
        "CMD" : "Put",
        "PROGRAM" : program,
        "DATABASE" : database,
        "QUERY" : urllib.parse.quote(query)
    }

    if (entrez):
        params["ENTREZ_QUERY"] = entrez;

    r = requests.post("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi",params=params);
    info = extractBlastInfo(r.text);
    return info['RID'], int(info['RTOE']);

def waitForCompletion(rid):
    i = 0;
    while(True):
        r = requests.get("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=%s" % rid);
        info = extractBlastInfo(r.text);
        status = info["Status"]
        if (status == "READY"): return True;
        if (status == "FAILED"): return False;
        if (status == "UNKNOWN"): return False;

        i += 1
        sleeptime = 5*i;
        if (sleeptime > 30): sleeptime = 30;
        time.sleep(sleeptime);

def getResults(rid):
    r = requests.get("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID=%s" % rid);
    return r.text;

def main():
    args = parser.parse_args()

    if (not args.rid):
        fasta = args.fasta.read()
        sys.stderr.write("Program: %s\nDatabase: %s\nEntrez: %s\nQuery: \n%s" % (args.program,args.database,args.entrez,fasta));
        rid,rtoe = startBlast(args.program,args.database,args.entrez,fasta);
        rid = rid.strip()
        sys.stderr.write("RID: %s\n" % rid);
        sys.stderr.write("Estimated time: %d seconds\n" % rtoe);
        time.sleep(rtoe);
        waitForCompletion(rid);
    else:
        rid = args.rid;
    results = getResults(rid)
    with open("_latest.blast","w") as f:
        f.write(results);
    print(results.strip());
    sys.stderr.write("View in browser: http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=HTML&RID=%s\n" % rid);

main();
