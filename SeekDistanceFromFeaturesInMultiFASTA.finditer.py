#!/usr/bin/env python
# coding: utf-8

# python ver. 2.4 or more higher
# [SYNONYM]$ python < this script > < reference genome fasta > < gtf or gff3 > < query sequence fasta >

import re, sys, string

pattern = r"^>"
repattern = re.compile(pattern)
restringnumber = re.compile(r'\w')

def ReverseComplement( seq ):
    return seq.translate(string.maketrans("ATGCRYSWKMBDHVatgcryswkmbdhv", "TACGYRWSMKVHDBtacgyrwsmkvhdb"))[::-1]

def ReadSingleFASTA( text, tag, seq ):
    while True:
        text = text.rstrip()
        QueryPatternMatch = repattern.match( text )
        if QueryPatternMatch or not text:
            if tag == "":
                tag = text
            else:
                return [ tag, seq[:10] ]
                if not text:
                    fpq.close()
                    break
                tag = text
                seq = ""
        else:
            seq += text


    if QueryPatternMatch:
        if tag == "":
            tag = text
        else:
            retrun [ tag,  "Query!!!", seq[:10] ]
            tag = ""
            tag = ""
    else:
        seq += text

# Index of first nucleotide of chromosome is 1 not 0
def PatternSearch( query, reference, queryname, referencename):

    modifiedquery = query.upper()
    modifiedquery = re.sub("R", r"[AG]", modifiedquery)
    modifiedquery = re.sub("Y", r"[CT]", modifiedquery)
    modifiedquery = re.sub("S", r"[CG]", modifiedquery)
    modifiedquery = re.sub("W", r"[AT]", modifiedquery)
    modifiedquery = re.sub("K", r"[GT]", modifiedquery)
    modifiedquery = re.sub("M", r"[AC]", modifiedquery)
    modifiedquery = re.sub("B", r"[CGT]", modifiedquery)
    modifiedquery = re.sub("D", r"[AGT]", modifiedquery)
    modifiedquery = re.sub("H", r"[ACT]", modifiedquery)
    modifiedquery = re.sub("V", r"[ACG]", modifiedquery)
    modifiedquery = re.sub("\.", "", modifiedquery) # "\" plays an important role to distinguish normal dot from regular expression.
    modifiedquery = re.sub("-", "", modifiedquery)
    modifiedquery = re.sub("N", ".", modifiedquery)
    MatchObject = re.finditer( modifiedquery, reference )
    if MatchObject:
        return [ queryname, referencename, MatchObject ]
    else:
        None

###
def PatternFindALL( query, reference, queryname, referencename):

    MatchObject = re.findall( query, reference )
    if MatchObject:
        return [ queryname, referencename, MatchObject ]
    else:
        None


ChromAccession = ""
QueryAccession = ""
ChromSequence = ""
QuerySequence = ""

AllArgumentFiles = sys.argv
n = len(AllArgumentFiles)

# Reference genome
print "Reference genome file is", AllArgumentFiles[ 1 ]
# Gff3 file
print "Concern Gff3 or GTF file is", AllArgumentFiles[ 2 ]
# Query fasta files
print "Query fasta files are", AllArgumentFiles[ 3: ]



for i in range(3, n):
    with open( AllArgumentFiles[i], 'r' ) as fpq:
        while True:
            text = fpq.readline().rstrip()
            QueryPatternMatch = repattern.match( text )
            if QueryPatternMatch or not text:
                if QueryAccession == "":
                    QueryAccession = text                 
                else:

                    ####
                    print QueryAccession
                    
                    with open( AllArgumentFiles[1], 'r' ) as fpr:
                        while True:
                            chromtext = fpr.readline().rstrip()
                            ChromPatternMatch = repattern.match( chromtext )
                            if ChromPatternMatch or not chromtext:
                                if ChromPatternMatch == "":
                                    ChromPatternMatch = chromtext
                                else:
                                    ####
                                    with open( AllArgumentFiles[2], 'r' ) as fpg:
                                        while True:
                                            gtftext = fpg.readline().rstrip()
                                            GtfPatternMatch = restringnumber.match( gtftext )
                                            if not gtftext:
                                                fpg.close()
                                                break
                                            elif GtfPatternMatch:
                                                gtf = gtftext.split('\t')
                                                GtfChrom = gtf[0]
                                                GtfAssembly = gtf[1]
                                                GtfFunction = gtf[2]
                                                GtfStart = int(gtf[3]) # Start position of the feature, with sequence numbering starting at 1
                                                GtfEnd = int(gtf[4]) # End position of the feature, with sequence numbering starting at 1
                                                GtfStrand = gtf[6] # defined as + (forward) or - (reverse) on the reference sequence
                                                #print ">"+GtfChrom, ChromAccession
                                                if re.match(">"+GtfChrom, ChromAccession) :# and GtfFunction == "gene":
                                                    if GtfStrand == "+":
                                                        x = PatternSearch(  QuerySequence, ChromSequence, QueryAccession, ChromAccession )
                                                        if x != None:
                                                            for j in x[2]:
                                                                print gtftext
                                                                print "absolute location cis\t", j.start()+1, j.end()
                                                                print "relative location cis\t",
                                                                #print j.start()+1 - GtfStart, j.end() - GtfStart
                                                                if GtfStart > j.end():
                                                                    print j.start()+1 - GtfStart, j.end() - GtfStart
                                                                elif GtfStart > j.start()+1 and GtfStart <= j.end():
                                                                    print j.start()+1 - GtfStart, j.end() - GtfStart + 1 # feature starts from +1 not 0
                                                                elif GtfStart <= j.start()+1:
                                                                    print j.start()+1 - GtfStart + 1, j.end() - GtfStart + 1 # feature starts from +1 not 0
                                                        y = PatternSearch(  ReverseComplement( QuerySequence ), ChromSequence, QueryAccession, ChromAccession )
                                                        if y != None:
                                                            for j in y[2]:
                                                                print gtftext
                                                                print "absolute location trans\t", j.start()+1, j.end()
                                                                print "relative location trans\t",
                                                                #print j.start()+1 - GtfStart, j.end() - GtfStart
                                                                if GtfStart > j.end():
                                                                    print j.start()+1 - GtfStart, j.end() - GtfStart
                                                                elif GtfStart > j.start()+1 and GtfStart <= j.end():
                                                                    print j.start()+1 - GtfStart, j.end() - GtfStart + 1 # feature starts from +1 not 0
                                                                elif GtfStart <= j.start()+1:
                                                                    print j.start()+1 - GtfStart + 1, j.end() - GtfStart + 1 # feature starts from +1 not 0

                                                    elif GtfStrand == "-":
                                                        x = PatternSearch(  ReverseComplement( QuerySequence ), ChromSequence, QueryAccession, ChromAccession )
                                                        if x != None:
                                                            for j in x[2]:
                                                                print gtftext
                                                                print "absolute location cis\t", j.start()+1, j.end()
                                                                print "relative location cis\t",
                                                                #print GtfEnd - j.end(), GtfEnd - j.start()+1
                                                                if GtfEnd < j.start()+1:
                                                                    print GtfEnd - j.end(), GtfEnd - ( j.start()+1 )
                                                                elif GtfEnd < j.end() and GtfEnd >= j.start()+1:
                                                                    print GtfEnd - j.end(), GtfEnd - ( j.start()+1 ) + 1 # feature starts from +1 not 0
                                                                elif GtfEnd >= j.end():
                                                                    print GtfEnd - j.end() + 1, GtfEnd - ( j.start()+1 ) + 1 # feature starts from +1 not 0
                                                        y = PatternSearch(  QuerySequence, ChromSequence, QueryAccession, ChromAccession )
                                                        if y != None:
                                                            for j in y[2]:
                                                                print gtftext
                                                                print "absolute location trans\t", j.start()+1, j.end()
                                                                print "relative location trans\t",
                                                                #print GtfEnd - j.end(), GtfEnd - j.start()+1
                                                                if GtfEnd < j.start()+1:
                                                                    print GtfEnd - j.end(), GtfEnd - (j.start()+1)
                                                                elif GtfEnd < j.end() and GtfEnd >= j.start()+1:
                                                                    print GtfEnd - j.end(), GtfEnd - (j.start()+1) + 1 # feature starts from +1 not 0
                                                                elif GtfEnd >= j.end():
                                                                    print GtfEnd - j.end() + 1, GtfEnd - (j.start()+1) + 1 # feature starts from +1 not 0
                                    ####
                                    if not chromtext:
                                        fpr.close()
                                        break
                                    ChromAccession = chromtext
                                    ChromSequence = ""
                            else:
                                ChromSequence += chromtext
                    ####
                    if not text:
                        fpq.close()
                        break
                    QueryAccession = text
                    QuerySequence = ""
            else:
                QuerySequence += text
