#!/usr/bin/env python

import time
import sys
import os.path
import argparse
import colorama
import re
from argparse import RawTextHelpFormatter

help="""
Formats a MAFFT fasta
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('fasta', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the MAFFT aligned FASTA sequence(s)')
parser.add_argument('-o', "--output", type=argparse.FileType('w'), nargs="?", default=sys.stdout, help='Output file')
parser.add_argument('-f', "--fastaout", action="store_const", const=True, default=False, help="Output in a fasta format")

parser.add_argument("--wrap",dest="wrap", type=int, default="60",help='Wrap the fasta to this width')
parser.add_argument('-C', "--colors", nargs="?", default=None, help='Provide custom colors in the format "255,0,255 0,255,255"')

parser.add_argument('-s', dest="similar", action="store_const", const=True, default=False,help='Show similar proteins')
parser.add_argument('-a', dest="coloraminoacids", action="store_const", const=True, default=False,help='Color amino acids according to their group')

parser.add_argument('-m', "--motifs", type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='Provide a regions file for coloration')
parser.add_argument('-c', "--capitalization", type=argparse.FileType('r'), nargs="?", help='Provides a region file for capitalization')
parser.add_argument('-b', "--bold", type=argparse.FileType('r'), nargs="?", help='Provides a region file for bolding')
parser.add_argument('-u', "--underline", type=argparse.FileType('r'), nargs="?", help='Provides a region file for underline')
parser.add_argument('-i', "--italics", type=argparse.FileType('r'), nargs="?", help='Provides a region file for italics')

#('green',         0,  255,   0 ),
#('darkblue',     0,    0, 128 ),
#('teal',          0,  128, 128 ),
#('violet',      128,    0, 128 ),
#('darkred',    128,    0,   0 ),
#('darkyellow', 128,  128,   0 ),
#('darkgrey',   128,  128, 128 ),
#('grey',        192,  192, 192 )

def parseRegionInfo(regionLines):
    regions = []
    for line in regionLines:
        line = line.strip()
        if (line.__len__() == 0): continue

        split = line.split("\t")
        if (split[0].strip() == "-"):
            start = None;
            end = None;
        else:
            start,_,end = split[0].partition("-")
            start = int(start);
            end = int(end);
        region_type = split[1]
        region_description = split[2]
        color_number = None
        if (len(split) >= 4): color_number = int(split[3]);
        regions.append((start,end,region_type,region_description,color_number))

    def getKey(item):
        if item[0] is None: return 0;
        return item[0]
    regions.sort(key=getKey)
    return regions

def getRegion(regions,index):
    for r in regions:
        if (r[0] == None): continue;
        if (r[0]-1 > index):
            return None
        if r[0]-1 <= index and r[1]-1 >= index:
            return r
    return None

def appendRegionStyle(style,regions,index):
    for rstyle,rlist in regions:
        region = getRegion(rlist,index);
        if (region == None): continue;
        if (not isinstance(rstyle,list)):
            rstyle = [rstyle];
        cindex = rlist.index(region)
        if (region[4] is not None): cindex = region[4];
        cindex = cindex % len(rstyle);
        style = style.add(rstyle[cindex]);
    return style;

def parseFasta(fasta):
    header = None
    sequence = ""
    sequences = []
    for line in fasta:
        line = line.strip()
        if (line.startswith(">")):
            if (header is not None):
                sequences.append((header,sequence))
            header = line
            sequence = ""
        else:
            sequence += line
    if (header is not None):
        sequences.append((header,sequence))
    return sequences


similarAminos="FLIMVATP SYCWQNG HRK DE".lower()
def isSimilar(a,b):
    for group in similarAminos.split(" "):
        if (a.lower() in group and b.lower() in group): return True
    return False;

def appendMatchStyle(style,c,orig,ss,isReference):
    if (c == "-"): return style
    if (ss.similar is not None and (isReference or (isSimilar(c,orig) and c != orig))):
        style = style.add(ss.similar)
    elif (c != orig):
        style = style.add(ss.nonmatch)
    return style

def appendAminoAcidColoring(style,c,ss):
    if (c == "-"): return style
    for i,group in enumerate(similarAminos.split(" ")):
        if (c.lower() in group):
            style = style.add(ss.amino[i])
    return style

def doAlignmentOutput(sequences,regions,ss,write,args):
    longestHeader = None
    for header,sequence in sequences:
        if (longestHeader is None or header.__len__() > longestHeader):
            longestHeader = header.__len__()
    if (longestHeader > 20): longestHeader = 20

    i = 0
    sequenceLength = sequences[0][1].__len__()
    while(i<sequenceLength):
        isReference = True;
        for head,seq in sequences:
            head = head[:min(head.__len__(),longestHeader)]
            write(head.ljust(longestHeader+3))
            for l in range(args.wrap):
                if (i+l >= sequenceLength): break
                c = seq[i+l]
                style = Style(ss.default)
                offset = sequences[0][1][:i+l].count("-")
                style = appendRegionStyle(style,regions,i+l-offset)
                style = appendMatchStyle(style,c,sequences[0][1][i+l],ss,isReference)
                if (args.coloraminoacids):
                    style = appendAminoAcidColoring(style,c,ss);
                write(c,style)
            write("\n")
            isReference = False;
        write("\n")
        i+= args.wrap

######################
def doFastaOutput(sequences,regions,ss,write,args):
    sequenceLength = sequences[0][1].__len__()
    isReference = True;
    for head,seq in sequences:
        write(head+"\n")
        offset = 0;
        i = 0
        while(i<sequenceLength):
            l = 0;
            while(l < args.wrap):
                if (i >= sequenceLength): break
                c = seq[i]
                style = Style(ss.default)
                if (ss.capitalization): style["capitalize"] = False;
                offset = sequences[0][1][:i].count("-")
                style = appendRegionStyle(style,regions,i-offset)
                style = appendMatchStyle(style,c,sequences[0][1][i],ss,isReference)
                if (args.coloraminoacids):
                    style = appendAminoAcidColoring(style,c,ss);
                i+=1
                if (c == "-"): continue;
                write(c,style)
                l+=1;
            write("\n")
        isReference = False;
######################

class Style(dict):
    def __init__(self,dic={}):
        self["bold"] = None
        self["italics"] = None
        self["underline"] = None
        self["capitalize"] = None
        self["foreground"] = None; #black
        self["background"] = None; #white

        self.update(dic)
    def add(self,addstyle):
        filtered = {k: v for k, v in addstyle.items() if v != None}
        self.update(filtered)
        return self;



def rtfStyle(a,style):
    a = a.replace("\n","\\\n")
    style["foreground"] = style["foreground"] or 0;
    style["background"] = style["background"] or 1;
    if (style["capitalize"] is not None):
        if (style["capitalize"]):
            a = a.upper()
        else:
            a = a.lower()

    out = ""
    if (style["italics"]): out += r"\i"
    if (style["bold"]): out += r"\b"
    if (style["underline"] and a != "-"): out += r"\ulw"
    out += r"\cf%d" % style["foreground"]
    out += "\cb%d \highlight%d" % (style["background"],style["background"])
    out += a
    if (style["underline"] and a != "-"): out += r"\ul0"
    if (style["bold"]): out += r"\b0"
    if (style["italics"]): out += r"\i0"
    return out

class Bunch(dict):
    __getattr__, __setattr__ = dict.get, dict.__setitem__

def parseColorString(colorString):
    if ("\n" in colorString):
        tokens = colorString.split("\n");
    else:
        tokens = colorString.split(" ");
    colors = []
    for tok in tokens:
        tok,_,_ = tok.partition("#");
        tok = tok.strip()
        if (tok.__len__() == 0): continue;
        colors.append([int(numeric_string) for numeric_string in tok.split(",")])
    return colors

def generateStyle(rtf,args):
    style = Bunch()
    if (rtf):
        style.default = Style();
        style.similar = None;
        style.nonmatch = Style({"bold":True})
        if args.capitalization: style.default["capitalize"]=False;
        if args.similar:
            style.default["foreground"]=2
            style.similar = Style({"foreground":0})
            style.nonmatch = Style({"bold":True,"foreground":0})

        colorString = "255,255,0 0,255,255 70,228,70 255,0,255 255,0,0 100,100,255"
        if (args.colors is not None):
            if (os.path.exists(args.colors)):
                with open(args.colors,"r") as f:
                    colorString = f.read()
            else:
                colorString = args.colors;
        
        reservedColors = parseColorString("255,255,255 150,150,150")
        aminoColors = parseColorString("0,255,0 255,255,0 255,0,0 0,0,255")
        colorList = parseColorString(colorString)
        style.region = [];
        style.amino = [];
        i = 1;
        colorTable = ""
        for c in reservedColors:
            colorTable += r"\red%d\green%d\blue%d;" % (c[0],c[1],c[2])
            i+=1
        for c in aminoColors:
            colorTable += r"\red%d\green%d\blue%d;" % (c[0],c[1],c[2])
            style.amino.append(Style({"foreground":i}));
            i+=1
        for c in colorList:
            colorTable += r"\red%d\green%d\blue%d;" % (c[0],c[1],c[2])
            style.region.append(Style({"background":i}));
            i+=1
        style.colorTable = colorTable;
    else:
        style.region = ["yellow", "green", "cyan", "magenta", "red"]
        style.nonmatch = "blue"
    return style

def main():
    global out
    args = parser.parse_args()

    sequences = parseFasta(args.fasta)

    rtfOutput = args.output and args.output.name.endswith(".rtf")
    ss = generateStyle(rtfOutput,args)

    regions = []
    if args.motifs:
        regions.append((ss.region,(parseRegionInfo(args.motifs.read().split("\n")))));

    if args.capitalization:
        regions.append((Style({"capitalize":True}),(parseRegionInfo(args.capitalization.read().split("\n")))));

    if args.bold:
        regions.append((Style({"bold":True}),(parseRegionInfo(args.bold.read().split("\n")))));

    if args.underline:
        regions.append((Style({"underline":True}),(parseRegionInfo(args.underline.read().split("\n")))));

    if args.italics:
        regions.append((Style({"italics":True}),(parseRegionInfo(args.italics.read().split("\n")))));

    if (rtfOutput):
        def write(a,style=Style()):
            global out
            out += rtfStyle(a,style)

        out = r"{\rtf1\ansi\ansicpg1252\cocoartf1348\cocoasubrtf170\n{\fonttbl\f0\fswiss\fcharset0 Courier New;}"
        out += r"{\colortbl;%s}" % ss.colorTable
        out += r"\margl1440\margr1440\vieww12600\viewh7800\viewkind0"
        out += r"\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural"
        out += r"\f0\fs24"

        if (args.fastaout):
            doFastaOutput(sequences,regions,ss,write,args)
        else:
            doAlignmentOutput(sequences,regions,ss,write,args)

        out += "}"

        args.output.write(out)

    else:
        coloramaFgColorMapping = {
            "black":colorama.Fore.BLACK,
            "red":colorama.Fore.RED,
            "green":colorama.Fore.GREEN,
            "yellow":colorama.Fore.YELLOW,
            "blue":colorama.Fore.BLUE,
            "magenta":colorama.Fore.MAGENTA,
            "cyan":colorama.Fore.CYAN,
            "white":colorama.Fore.WHITE
        }

        coloramaBgColorMapping = {
            "black":colorama.Back.BLACK,
            "red":colorama.Back.RED,
            "green":colorama.Back.GREEN,
            "yellow":colorama.Back.YELLOW,
            "blue":colorama.Back.BLUE,
            "magenta":colorama.Back.MAGENTA,
            "cyan":colorama.Back.CYAN,
            "white":colorama.Back.WHITE
        }

        def write(a,fg="black",bg="white"):
            print(coloramaFgColorMapping[fg]+coloramaBgColorMapping[bg]+a+colorama.Fore.RESET+colorama.Back.RESET)

        if (args.fastaout):
            doFastaOutput(sequences,regions,ss,write,args)
        else:
            doAlignmentOutput(sequences,regions,ss,write,args)

if __name__ == "__main__":
    main()
