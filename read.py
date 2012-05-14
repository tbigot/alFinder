# -*- coding: utf-8 -*-
# thomas.bigot@tmgo.net

import string


class Read:
    
    wobbles = {'r': ('A','G'), 'w': ('A','T'), 'm': ('A','C'), 'y': ('C','T'), 's': ('C','G'), 'k': ('G','T'), 'b': ('C','T','G'), 'd': ('A','G','T'), 'h': ('A','C','T'), 'v': ('A','C','G'), 'n': ('A','C','G','T')}
    
    
    _reads = []
    complementaries = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    
    
    @staticmethod
    def unWobble(inString,strand):
        allowedList = list(inString)
        for it in range(len(allowedList)):
            if allowedList[it] in Read.wobbles.keys():
                allowedList[it] = Read.wobbles[allowedList[it]]
            else:
                allowedList[it] = [allowedList[it]]
        if(strand == "R"):
            finalList = []
            for it in allowedList[::-1]:
                sl = []
                for el in it:
                    sl.append(Read.complementaries[el])
                finalList.append(sl)
            return(finalList)
        else:
            return(allowedList)
                
    
    
    def __init__(self,seq):
	self.name = seq[0][1:].split()[0]
	self.individual = None
	self.locus = None
	self.seq = string.join(seq[1:],'')
	self.rev = False
	self.allele = None
	
    def __str__(self):
	return 'Sequence ' + self.name + "; Individual " + str(self.individualName) + "; Locus " + str(self.locus)
	
    def getName(self):
	return self.name
	
    def getSeq(self):
	return self.seq
	
    def getTags(self):
	return((self.seq[0:6],self.seq[len(self.seq)-6:len(self.seq)]))
	
    def reverse(self):
	self.seq = Read.reverseComplementary(self.seq)
	self.rev = not self.rev
	
    @staticmethod
    def getReads():
	return Read._reads
	
    @staticmethod
    def loadFromFile(path):
	fhandle = open(path)
	currSeq = []
	for ligne in fhandle:
	    if ligne.startswith('>'):
		if not len(currSeq) == 0:
		    Read._reads.append(Read(currSeq))
		    currSeq = []
	    currSeq.append(ligne.strip().upper())
	    
    @staticmethod
    def reverseComplementary(seq):
	complementary = []
	for character in reversed(seq):
	    complementary.append(Read.complementaries[character])
	return "".join(complementary)
	
    @staticmethod
    def writeTo(path,showUnidentified=0):
	o = open(path,'w')
	o.write("seqName,individual,locus,allele")
	for currRead in Read._reads:
	    currRead.oneWriteTo(o,showUnidentified)
	    
    def oneWriteTo(self,output,showUnidentified):
	if (showUnidentified == 2 or (showUnidentified == 1 and self.individual != None) or (showUnidentified == 0 and self.allele != None)):
            if self.individual != None:
                indivName = self.individual.getName()
            else:
                indivName = "<unidentified>"
            if self.locus != None:
                locusName = str(self.locus+1)
            else:
                locusName = "<unidentified>"
            if self.allele != None:
                alleleName = self.allele
            else:
                alleleName = '<unidentified>'
	    output.write("\n" + self.name + "," + indivName + "," + locusName + "," + alleleName)
    
    @staticmethod
    def identify(IndividualClass):
	for currRead in Read._reads:
	    currRead.oneIdentify(IndividualClass)
    
    def oneIdentify(self,IndividualClass):
        
        # step 1 which locus is it?
        for currLocusRE in range(len(IndividualClass._lociRE)):
            if IndividualClass._lociRE[currLocusRE].match(self.seq):
                self.locus = currLocusRE
                break
        if self.locus == None and not self.rev:
            self.reverse()
            self.oneIdentify(IndividualClass)
        
        
        # step 2: knowing locus, which individual is it?
        if self.locus == None or self.individual != None:
            return
            
	tags = self.getTags()
	if (self.locus,tags) in IndividualClass._tags:
	    self.individual = IndividualClass._tags[self.locus,tags]
	    
	    
    @staticmethod
    def match(alleles):
	for currRead in Read._reads:
	    if currRead.individual != None:
                if currRead.allele == None:
                    currRead.oneMatch(alleles)
    
    def oneMatch(self,alleles):
	currAlleles = alleles[self.locus]
	for currAlleleName in currAlleles.keys():
	    if currAlleles[currAlleleName] in self.seq:
		self.allele = currAlleleName
		break
	self.individual.raiseSeqNr()
	if(self.allele == None):
	    self.individual.addUnknownAllele(self.locus,self.seq[6:-6])

	 
   
		