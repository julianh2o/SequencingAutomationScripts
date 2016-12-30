#!/bin/bash

if [ -t 0 ]; then
    A=`curl $1 2> /dev/null`
else
    A=`cat`
fi

ACCESSION=`echo -n $A | perl -0777 -ne 'if(m/ACCESSION\s+(\S+)/g){print "$1\n";}'`
BODY=`echo -n $A | perl -0777 -ne 'if(m/ORIGIN([\w\W]+)\/\//g){print "$1\n";}' | perl -pe 's/[\s0-9]//g'`

BODY=`echo -n $BODY | awk '{print toupper($0)}'`

printf ">$ACCESSION\n$BODY\n"
