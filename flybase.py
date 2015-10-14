#!/usr/bin/env python

import requests
import urllib
import re
import time
import sys
import os.path
import argparse
from argparse import RawTextHelpFormatter

help="""
Performs a blast using Flybase\'s web tool

Note: this program uses "dmel" and "drosphilia" as org and tax respectively

Perform a blast and get XML result
    flybase.py blasn org/translation query.fasta > result.blast

Perform a blast using piped data
    extractfastaseq.py input.fasta -n 1 | flybase.py blasn org/translation > result.blast

Perform a blast and immediately return the first result info
    flybase.py blasn org/translation query.txt | viewblast.py into 1
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('program', help='program: blastn blastp blastx tblastn tblastx')
parser.add_argument('database', help='database')
parser.add_argument('fasta', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the FASTA sequence(s)')

def extractBlastInfo(responseText):
    info = {};
    for m in re.finditer("<!--.*?QBlastInfoBegin(.*?)QBlastInfoEnd.*?-->",responseText,re.DOTALL):
        content = m.group(1)
        for mm in re.finditer("\W+(.*?)\W*?=\W*?(.*)",content):
            info[mm.group(1)] = mm.group(2);
    return info;

def pretty_print_POST(req):
    """
    At this point it is completely built and ready
    to be fired; it is "prepared".

    However pay attention at the formatting used in 
    this function because it is programmed to be pretty 
    printed and may differ from the actual request.
    """
    print('{}\n{}\n{}\n\n{}'.format(
        '-----------START-----------',
        req.method + ' ' + req.url,
        '\n'.join('{}: {}'.format(k, v) for k, v in req.headers.items()),
        req.body.decode("ascii"),
    ))


def startBlast(program,database,query):
    params = {
        "program" : ('',program),
        "database" : ('',database),
        "seqfile" : ('sequence.txt',query),
#"sequence" : ('',''),
#"advoptionsopen" : ('','false'),
#"matrix" : ('','BLOSUM62'),
#"otheropts" : ('',''),
#"geneticcodes" : ('','1'),
#"wordsize" : ('',''),
#"_lowcomp" : ('',''),
#"numalignments" : ('','25'),
#"numdescriptions" : ('','25'),
#"_legacyengine" : ('',''),
#"_megablast" : ('',''),
#"expect" : ('','10'),
        "org": ('',"dmel"),
        "tax" : ('',"drosophila"),
        "outputformat":('',"XML")
    }

    req = requests.Request("POST","http://flybase.org/blast/index.html",files=params,headers={
        'User-Agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/45.0.2454.101 Safari/537.36'
    });
    prepared = req.prepare()
    s = requests.Session()
    r = s.send(prepared,allow_redirects=False);
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
