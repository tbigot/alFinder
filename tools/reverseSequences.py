# -*- coding: utf-8 -*-
# thomas.bigot@tmgo.net


import sys,string

if len(sys.argv) != 3:
    print("Usage :")
    print(sys.argv[0] + " inputFile outputFile")
    sys.exit()

complementaries = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

def reverseComplementary(seq):
        complementary = []
        for character in reversed(seq):
            complementary.append(complementaries[character])
        return "".join(complementary)
    

fhandle = open(sys.argv[1])
ofhandle = open(sys.argv[2],"w")

currSeq = []
for ligne in fhandle:
    if ligne.startswith('>') or not ligne.endswith('\n'):
        if not len(currSeq) == 0:
            ofhandle.write(currSeq[0]+"\n")
            ofhandle.write(reverseComplementary(string.join(currSeq[1:],""))+"\n")
            currSeq = []
    currSeq.append(ligne.strip().upper())