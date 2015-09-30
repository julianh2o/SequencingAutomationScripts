#!/usr/local/bin/python3

import re
import sys
import fileinput
import argparse

parser = argparse.ArgumentParser(description='Translate DNA/RNA sequence into a protein sequence')
parser.add_argument("-a", dest="showAll", action="store_const", const=True, default=False,help='show all possible translation')
parser.add_argument('file', metavar="FILE", type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='input sequence')
parser.add_argument('-f', dest="readingFrameOnly", action="store_const", const=True, default=False,help='show only the largest reading frame')
parser.add_argument('-w', "--warn",dest="warnLength", type=int, default=100,help='warn to stderr if there are any reading frames other than the longest longer than this number')

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

def genPermutation(sequence,frameIndex,reverse):
    rframeIndex = (len(sequence) - frameIndex) % 3;
    if (reverse):
        sequence = complement(sequence);
        sequence = sequence[::-1]
        sequence = sequence[rframeIndex:];
    else:
        sequence = sequence[frameIndex:];

    return sequence;

def permute(sequence):
    permutations = [];
    a = genPermutation(sequence,0,False);
    b = genPermutation(sequence,1,False);
    c = genPermutation(sequence,2,False);
    rca = genPermutation(sequence,0,True);
    rcb = genPermutation(sequence,1,True);
    rcc = genPermutation(sequence,2,True);

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

def translate(sequence):
    translated = "";
    i=0;
    while(i<len(sequence)):
        if (len(sequence) < i+3):
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
    longest = None;
    for m in re.finditer("M.*?\*",translated):
        if (not longest or len(longest) < len(m.group(0))):
            longest = m.group(0);
    return longest;

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
        print("ERROR: Could not detect DNA/RNA");

    if (isDna):
        sequence = dnaToRna(sequence);

    allTranslated = [];
    for name,seq in permute(sequence):
        translated = translate(seq);
        unbrokenLength = countWithoutStop(translated);
        allTranslated.append((unbrokenLength,name,translated));

    def getKey(item):
        return item[0];
    allTranslated.sort(key=getKey,reverse=True);
    
    for trans in allTranslated[1:]:
        unbrokenLength, name, translated = trans;
        if (unbrokenLength >= args.warnLength):
            name,_,_ = header.partition(" ");
            name = name.strip(">");
            sys.stderr.write("Translate: long frame detected in %s length:%d (longest: %d)\n" % (name,unbrokenLength,allTranslated[0][0]));

    if (args.readingFrameOnly):
        if (header): print(header,"[longest reading frame only]");
        print(extractLongestReadingFrame(allTranslated[0][2]));
    elif (not args.showAll):
        if (header): print(header);
        print(allTranslated[0][2]);
    else:
        #showall
        for trans in allTranslated:
            unbrokenLength, name, translated = trans;
            print("Frame %s (max length %d)" % (name,unbrokenLength));
            if (header): print(header)
            print(translated);
            print();

main();
