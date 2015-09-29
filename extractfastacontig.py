#!/usr/local/bin/python3
import sys

i = 0;
accum = None;
for line in sys.stdin:
    i+=1;
    if (line.startswith(">")):
        if (accum is not None):
            break;
        if (sys.argv[1] in line):
            accum = "";
    if (accum is not None):
        accum += line;

print(accum,end="");
