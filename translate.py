#!/usr/bin/env python

from __future__ import print_function
import re
import sys
import fileinput
import argparse
import util
from argparse import RawTextHelpFormatter
import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

help="""
Translate DNA/RNA sequence into a protein sequence

Translates the sequence and displays only the translation with the longest unbroken sequence
Note: a warning is generated for any unbroken sequences that are longer than 50% of the longest
    translate.py sequence.fasta

Translates the sequence and displays all possible translations
    translate.py sequence.fasta -a

Traslates the sequence and then extracts the open reading frame
    translate.py sequence.fasta -f

Translates the sequence but wraps to 80 characters instead of 60
    translate.py sequence.fasta --wrap 80

Translates the sequence and does not wrap the output
    translate.py sequence.fasta --wrap 0

Translates the sequence and echos a warning to stderr for any sequences that are longer than 100 amino acids
    translate.py sequence.fasta -w 100

Translates the sequence and echos a warning to stderr for any sequences that are longer than 25% of the longest sequence
    translate.py sequence.fasta -w 25%

Translates the sequence provided on stdin
    cat sequence.fasta | translate.py
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument("-a", dest="showAll", action="store_const", const=True, default=False,help='show all possible translation')
parser.add_argument("-r", "--readingframe", type=int, default=None,help='Choose a specific reading frame only')
parser.add_argument('file', metavar="FILE", type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='input sequence')
parser.add_argument('-f', dest="readingFrameOnly", action="store_const", const=True, default=False,help='show only the largest reading frame')
parser.add_argument('-w', "--warn",dest="warnLength", default="50%",help='warn to stderr if there are any reading frames other than the longest longer than this number (or percent of the longest)')
parser.add_argument('-g', "--noregen",dest="regenerateFromThreePrime", action="store_const", const=False, default=True,help='Use the same codon groupings rather than regenerating the reading frame from the 3\' end (default is compatible with expasy and blast, with flag is compatible with emboss)')
parser.add_argument('-o', "--orderByFrame",dest="orderByFrame", action="store_const", const=True, default=False,help='Order by frame index instead of match (used with -a)')
parser.add_argument("--wrap",dest="wrap", type=int, default="60",help='Wrap the fasta to this width')

mapping = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
    "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

def trimTo3(sequence):
    extraCodons = len(sequence) % 3
    if (extraCodons == 0): return sequence;
    seq = sequence[:-extraCodons];
    return seq;

#independentReverseFrame means that we recreate the frame from the 3' end without trimming the overhanging characters
#use this to be more compatible with ncbi's blast tool
def genPermutation(sequence,frameIndex,reverse,independentReverseFrame):
    rframeIndex = (len(sequence) - frameIndex) % 3;
    if (reverse):
        sequence = complement(sequence);
        sequence = sequence[::-1]
        if (not independentReverseFrame): frameIndex = rframeIndex
        sequence = sequence[frameIndex:];
    else:
        sequence = sequence[frameIndex:];

    return sequence;

def permute(sequence,independentReverseFrame):
    permutations = [];
    a = genPermutation(sequence,0,False,independentReverseFrame);
    b = genPermutation(sequence,1,False,independentReverseFrame);
    c = genPermutation(sequence,2,False,independentReverseFrame);
    rca = genPermutation(sequence,0,True,independentReverseFrame);
    rcb = genPermutation(sequence,1,True,independentReverseFrame);
    rcc = genPermutation(sequence,2,True,independentReverseFrame);

    permutations.append(("1",a));
    permutations.append(("2",b));
    permutations.append(("3",c));
    permutations.append(("-1",rca));
    permutations.append(("-2",rcb));
    permutations.append(("-3",rcc));
    return permutations;

def complement(sequence):
    compl = "";
    for c in sequence:
        if (c == "A"): compl += "U";
        if (c == "U"): compl += "A";
        if (c == "C"): compl += "G";
        if (c == "G"): compl += "C";
    return compl;

def subdivide(seq):
    i = 0;
    out = "";
    for c in seq:
        if (i%3 == 0): out += " ";
        out += c;
        i+=1
    return out

def translate(sequence):
    translated = "";
    i=0;
    while(i<len(sequence)):
        if (len(sequence) < i+3 or "N" in sequence[i:i+3]):
            translated += "X";
        else:
            translated += mapping[sequence[i:i+3]];
        i+=3;
    return translated;

def detectDNA(sequence):
    u = 0;
    t = 0;
    for c in sequence:
        if (c is 'U'): u += 1;
        if (c is 'T'): t += 1;
    if (u > 0 and t > 0):
        return None;
    if (u > 0):
        return False;
    if (t > 0):
        return True;

    return None

def dnaToRna(sequence):
    return sequence.replace("T","U");

def countWithoutStop(translated):
    count = 0;
    maxCount = 0;
    for c in translated:
        if (c == '*'): count = 0;
        count += 1;
        if (count > maxCount): maxCount = count;
    return maxCount;

def extractLongestReadingFrame(translated):
    frames = []
    split = translated.split("*");
    for index in range(split.__len__()):
        tok = split[index]
        if (index == 0):
            frames.append((tok.__len__(),index,tok));
        elif ("M" in tok):
            first = tok.index("M")
            frame = tok[first:];
            frames.append((frame.__len__(),index,frame));

    def getKey(item):
        return item[0];
    frames.sort(key=getKey,reverse=True);

    for size,index,frame in frames:
        return frame;

def main():
    args = parser.parse_args()
    sequence = args.file.read()
    header,_,remain = sequence.partition("\n");
    if (header.startswith(">")):
        sequence = remain;
    else:
        header = None;
            
    sequence = sequence.replace("\n","").upper();
    isDna = detectDNA(sequence);
    if (isDna is None):
        eprint("ERROR: Could not detect DNA/RNA");
        sys.exit(1);

    if (isDna):
        sequence = dnaToRna(sequence);

#print("IN SEQ         ",sequence)
#print("IN SEQ         ",complement(sequence));
    allTranslated = [];
    for name,seq in permute(sequence,args.regenerateFromThreePrime):
#print("translating",name.ljust(2),subdivide(seq));
        translated = translate(seq);
        unbrokenLength = countWithoutStop(translated);
        if args.readingframe:
            if int(name) == args.readingframe:
                allTranslated.append((unbrokenLength,name,translated));
        else:
            allTranslated.append((unbrokenLength,name,translated));

    def getKey(item):
        return item[0];
    if (not args.orderByFrame): allTranslated.sort(key=getKey,reverse=True);
    
    for trans in allTranslated[1:]:
        unbrokenLength, name, translated = trans;
        longest = allTranslated[0][0]
        warnSize = args.warnLength.strip();
        if args.warnLength.endswith("%"):
            warnSize = longest * float(warnSize.strip("%"))/100
        warnSize = int(warnSize)
        if (unbrokenLength >= warnSize):
            name,_,_ = header.partition(" ");
            name = name.strip(">");
            sys.stderr.write("Translate: long frame detected in %s length:%d (longest: %d)\n" % (name,unbrokenLength,allTranslated[0][0]));

    if (args.readingFrameOnly):
        if (header): print(header,"[longest reading frame only]");
        print(util.wrap(extractLongestReadingFrame(allTranslated[0][2]),args.wrap));
    elif (not args.showAll):
        if (header): print(header);
        print(util.wrap(allTranslated[0][2],args.wrap));
    else:
        #showall
        for trans in allTranslated:
            unbrokenLength, name, translated = trans;
            print("Frame %s (max length %d)" % (name,unbrokenLength));
            if (header): print(header)
            print(translated);
            print();

main();
