#! /usr/bin/env python

import sys

seq = 'PHPPHPHPPHH'
seq = 'HHPPHH'
col = {}
col['H'] = 'O   ALA'
col['P'] = 'N   ASP'

def writeAtom(i,x,y,at):
    print "ATOM%7i  %s%6i    %8.3f%8.3f%8.3f" % (i,col[at],i,x,y,0)

fo = open(sys.argv[1],'r')
c = 0
for l in fo:
    a = eval(l)
    print "MODEL %i" % c
    i = 0
    for x,y in a:
        writeAtom(i,x,y,seq[i])
        i += 1
    print "TER"
    print "END"
    c += 1

