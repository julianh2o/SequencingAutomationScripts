#!/usr/bin/env python

import time
import sys
import os.path
import argparse
import re
from argparse import RawTextHelpFormatter

help="""
Formats a MAFFT aligned fasta file

Performs the most basic formatting aligning 60 characters at a time
    formatmafft.py input.fasta -o output.rtf

Formats as alignment and colors regions based on information in transmembrane.tsv
    formatmafft.py input.fasta -o output.rft -m transmembrane.tsv

As above, but also capitalize the highlighted regions and lowercases everything else
    formatmafft.py input.fasta -o output.rft -m transmembrane.tsv -c transmembrane.tsv

Formats aligned with bold showing which amino acids don't match the reference
    formatmafft.py input.fasta -o output.rft -n bold

Formats aligned with colored regions, non-matching amino acids within the regions are given a white background
    formatmafft.py input.fasta -o output.rft -m transmembrane.tsv -n bg

As above, but non-matching amino acids are also bolded
    formatmafft.py input.fasta -o output.rft -m transmembrane.tsv -n bold,bg

Greys out matching amino acids, colors similar amino acids black, and non-similar ones bold
    formatmafft.py input.fasta -o output.rft -s

Underlines the repeats stored in repeats.tsv
    formatmafft.py input.fasta -o output.rft -u repeats.tsv

Formats as a fasta file with colored sections according to the regions file
    formatmafft.py input.fasta -o output.rft -m transmembrane.tsv -f

Formats aligned but displays a color based on the individual amino acid
    formatmafft.py input.fasta -o output.rft -m transmembrane.tsv -a

Provides a custom color pallet for coloring the regions, in the example red, green, blue
    formatmafft.py input.fasta -o output.rft -m transmembrane.tsv --colors '255,0,0 0,255,0 0,0,255'

Provides a custom color pallet that is loaded from the file colors.txt (can be separeted by newlines and commented with #s)
    formatmafft.py input.fasta -o output.rft -m trasnmembrane.tsv --colors colors.txt

Colors regions from stdin instead of from a file
    fetchuniprot.py P35500 -m transmem | formatmafft.py input.fasta -o output.rft -m -

Combines many features, output in fasta format and is colored and capitalized by section, repeats are underlined, amino acids colored and bolded based on matches
    formatmafft.py input.fasta -o output.rft -m transmembrane.tsv -c transmembrane.tsv -u repeats -a -f -n bold
"""

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('fasta', type=argparse.FileType('r'), nargs="?", default=sys.stdin, help='the MAFFT aligned FASTA sequence(s)')
parser.add_argument('-o', "--output", type=argparse.FileType('w'), nargs="?", default=sys.stdout, help='Output file')
parser.add_argument('-f', "--fastaout", action="store_const", const=True, default=False, help="Output in a fasta format")

parser.add_argument("--wrap",dest="wrap", type=int, default="60",help='Wrap the fasta to this width')
parser.add_argument('-C', "--colors", default=None, help='Provide custom colors in the format "255,0,255 0,255,255"')
parser.add_argument('-A', "--aminocolors", default=None,help='Provide custom colors for amino acids, same format as -C')

parser.add_argument('-n', dest="nomatch", default="",help='Show nonmatching proteins (accepts "bold", "bg", or "amino". values can combined by comma')
parser.add_argument('-s', dest="similar", action="store_const", const=True, default=False,help='Show similar proteins')
parser.add_argument('-a', dest="coloraminoacids", action="store_const", const=True, default=False,help='Color amino acids according to their group')

parser.add_argument('-m', "--motifs", type=argparse.FileType('r'), nargs="?", help='Provide a regions file for coloration')
parser.add_argument('-c', "--capitalization", type=argparse.FileType('r'), nargs="?", help='Provides a region file for capitalization')
parser.add_argument('-b', "--bold", type=argparse.FileType('r'), nargs="?", help='Provides a region file for bolding')
parser.add_argument('-u', "--underline", type=argparse.FileType('r'), nargs="?", help='Provides a region file for underline')
parser.add_argument('-i', "--italics", type=argparse.FileType('r'), nargs="?", help='Provides a region file for italics')

def rtf_encode_char(unichar):
    code = ord(unichar)
    if code < 128:
        return str(unichar)
    return '\\u' + str(code if code <= 32767 else code-65536) + '?'

def rtf_encode(unistr):
    return ''.join(rtf_encode_char(c) for c in unistr)

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
    elif (c != orig and ss.nonmatch is not None):
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
                if "amino" in style["special"]:
                    style = appendAminoAcidColoring(style,c,ss);
                    style["special"].remove("amino")
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
                if "amino" in style["special"]:
                    style = appendAminoAcidColoring(style,c,ss);
                    style["special"].remove("amino")
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
        self["special"] = [];

        self.update(dic)
        self["special"] = self["special"][:]

        if (isinstance(self["special"],str)):
            self["special"] = [self["special"]];
    def add(self,addstyle):
        filtered = {k: v for k, v in addstyle.items() if v != None}
        del filtered["special"]
        self.update(filtered)
        if addstyle["special"]:
            self["special"] += addstyle["special"];
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

    a = a.replace("{","\\{")
    a = a.replace("}","\\}")
    a = rtf_encode(a);

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

def loadColorsOrDefault(argument,default):
    colorString = default
    if (argument is not None):
        if (os.path.exists(argument)):
            with open(argument,"r") as f:
                colorString = f.read()
        else:
            colorString = argument;

    return parseColorString(colorString);
        

def generateStyle(rtf,args):
    style = Bunch()
    style.default = Style();
    style.similar = None;
    style.nonmatch = None;

    if (args.coloraminoacids):
        style.default = style.default.add(Style({"special":"amino"}));

    if args.nomatch:
        modeStyles = {
            "bold":Style({"bold":True}),
            "bg":Style({"background":0}),
            "amino":Style({"special":"amino"}),
        }
        style.nonmatch = Style();
        for mode in args.nomatch.split(","):
            style.nonmatch = style.nonmatch.add(modeStyles[mode]);

    if args.capitalization: style.default["capitalize"]=False;
    if args.similar:
        style.default["foreground"]=2
        style.similar = Style({"foreground":0})
        style.nonmatch = Style({"bold":True,"foreground":0})

    if (rtf):

        reservedColors = parseColorString("255,255,255 150,150,150")
        aminoColors = loadColorsOrDefault(args.aminocolors,"0,255,0 255,255,0 255,0,0 0,0,255")
        colorList = loadColorsOrDefault(args.colors,"255,255,0 0,255,255 70,228,70 255,0,255 255,0,0 100,100,255")
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
        for i in range(6):
            style.region = [];
            style.region.append(Style({"background":i}));

            style.amino = [];
            style.amino.append(Style({"foreground":i}));
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
        import colorama

        coloramaBgColors = [colorama.Back.YELLOW, colorama.Back.CYAN, colorama.Back.GREEN, colorama.Back.MAGENTA,colorama.back.RED,colorama.Back.BLUE]
        coloramaFgColors = [colorama.Fore.YELLOW, colorama.Fore.CYAN, colorama.Fore.GREEN, colorama.Fore.MAGENTA,colorama.back.RED,colorama.Fore.BLUE]
        def write(a,fg="black",bg="white"):
            print(coloramaFgColorMapping[fg]+coloramaBgColorMapping[bg]+a+colorama.Fore.RESET+colorama.Back.RESET)

        if (args.fastaout):
            doFastaOutput(sequences,regions,ss,write,args)
        else:
            doAlignmentOutput(sequences,regions,ss,write,args)

if __name__ == "__main__":
    main()
