# -*- coding: utf-8 -*-
# thomas.bigot@tmgo.net

import string


class Read:
    
    _reads = []
    complementaries = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    
    
    def __init__(self,seq):
	self.name = seq[0][1:].split()[0]
	self.individual = None
	self.individualName = '<unidentified>'
	self.locus = '<unidentified>'
	self.seq = string.join(seq[1:],'')
	self.rev = False
	self.identified = False
	self.allelle = '<unknown>'
	self.matched = False
	
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
	if (showUnidentified == 2 or (showUnidentified == 1 and self.identified) or (showUnidentified == 0 and self.matched)):
	    output.write("\n" + self.name + "," + self.individualName + "," + str(self.locus) + "," + self.allelle)
    
    @staticmethod
    def identify(indivList):
	for currRead in Read._reads:
	    currRead.oneIdentify(indivList)
    
    def oneIdentify(self,indivList):
	tags = self.getTags()
	#Debug
	#print("Looking for " + str(tags) + " in the list.")
	if tags in indivList:
	    found = indivList[self.getTags()]
	    self.individual = found[0]
	    self.individualName = self.individual.getName()
	    self.locus = found[1]
	    self.identified = True
	    # debug
	    #print("Séquence " + self.name + " identifiée. Individual = "+ self.individualName.getName() + ". Locus = " + str(self.locus))
	elif not self.rev:
	    self.reverse()
	    Read.oneIdentify(self,indivList)
	#else:
	    #print("Impossible d’identifier la séquence "+ self.name +". Désolé.")
	    
    @staticmethod
    def match(alleles):
	for currRead in Read._reads:
	    if currRead.identified:
                if not currRead.matched:
                    currRead.oneMatch(alleles)
    
    def oneMatch(self,alleles):
	currAllelles = alleles[self.locus]
	for currAllelleName in currAllelles.keys():
	    if currAllelles[currAllelleName] in self.seq:
		self.allelle = currAllelleName
		self.matched = True
		break
	self.individual.raiseSeqNr()
	if(not self.matched):
	    self.individual.addUnknownAllele(self.locus,self.seq[6:-6])

	 
   
		