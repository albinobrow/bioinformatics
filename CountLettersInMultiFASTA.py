#!/usr/bin/env python
# coding: utf-8

import re, sys, string

pattern = r"^>"
repattern = re.compile(pattern)

allfiles = sys.argv
n = len(allfiles)

for i in range(1, n):
    fp = open(allfiles[i], 'r')
    for j in fp:
        text = j.rstrip()
        patternmatching = repattern.match(text)
        if patternmatching:
            try:
                c
                print c
            except NameError:
                None
            print text
            c=0
        else:
            text=''.join(text.split())
            c+=len(text)
    print c
    del c
    fp.close()
