# -*- coding: utf-8 -*-
# Copyright or © or Copr. Université Claude Bernard Lyon 1
# contributor : Thomas Bigot (2012-2014)
#
# thomas.bigot@univ-lyon1.fr
#
# This software is a bioinformatics computer program whose purpose is
# recognize and post-process sequence data (assosiate a sequence read to
# an individual, a locus, and one peculiar variant of this locus).
# It is described in this article:
# "Large-scale genotyping by next generation sequencing:
# how to overcome the challenges to reliably genotype individuals?"
# Ferrandiz-Rovira M, Bigot T, Allainé D, Callait-Cardinal M-P,
# Radwan J, Cohas A.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

import string

## Class Individual
#
#  The object Read is designed to manage the reads in the NGS meaning:
#  one read = one sequence in a fasta file
#  Also, other static methods have been included in this class,
#  mainly ones related to sequence management.
class Read:
    
    # Static, for reading only.
    # is used to convert wobbles to a list of corresponding nucleotides
    wobbles = {'r': ('A','G'), 'w': ('A','T'), 'm': ('A','C'), 'y': ('C','T'), 's': ('C','G'), 'k': ('G','T'), 'b': ('C','T','G'), 'd': ('A','G','T'), 'h': ('A','C','T'), 'v': ('A','C','G'), 'n': ('A','C','G','T')}
    
    # Static: the list of all Read objects
    _reads = []
    
    # Static, for reading only.
    # is used to convert a nucleotide to its complementary
    complementaries = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    
    
    # Static: convert a sequence containing wobbles to a list of possible
    # nucleotides of each position
    # @param inString the input string (eg: rwAT)
    # @param strand F or R
    # @return a list of possible nucleotide for each position (for the above example:
    #                   [('A','G'),('A','T'),['A'],['T']])
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
                
    
    # Constructor
    # @param a raw sequence string (including name and several lines)
    def __init__(self,seq):
	self.name = seq[0][1:].split()[0]
	self.individual = None
	self.locus = None
	self.seq = string.join(seq[1:],'')
	self.rev = False
	self.allele = None
	
    # For debug purpose
    # printing a sequence
    def __str__(self):
	return 'Sequence ' + self.name + "; Individual " + str(self.individualName) + "; Locus " + str(self.locus)
	
    # Name accessor
    def getName(self):
	return self.name
	
    # Sequence accessor
    def getSeq(self):
	return self.seq
	
    # Extract tags of this read
    def getTags(self):
	return((self.seq[0:6],self.seq[len(self.seq)-6:len(self.seq)]))
	
    # Reverse in place the current sequence
    def reverse(self):
	self.seq = Read.reverseComplementary(self.seq)
	self.rev = not self.rev
	
    # True if this sequences is in the reverse sense
    def isReversed(self):
        return seq.rev
	
    # Static: accessor of the list of all Read objects
    @staticmethod
    def getReads():
	return Read._reads
    
    # Static: returns the number of all Read objects
    @staticmethod
    def getNumberOfReads():
        return len(Read._reads)
    
    # Static: load all the reads from a file
    # @param path the path of the file
    @staticmethod
    def loadFromFile(path):
	fhandle = open(path)
	currSeq = []
	for ligne in fhandle:
	    if ligne.startswith('>') or not ligne.endswith('\n'):
		if not len(currSeq) == 0:
		    Read._reads.append(Read(currSeq))
		    currSeq = []
	    currSeq.append(ligne.strip().upper())
	    
    # Static: returns the reverse complementary of a sequence.
    # @param seq an input sequence such as ATTG
    # @return the reverse complementary such as CAAT
    @staticmethod
    def reverseComplementary(seq):
	complementary = []
	for character in reversed(seq):
	    complementary.append(Read.complementaries[character])
	return "".join(complementary)
	
    # Static write the result of the identification to a file
    # @param path path of the file in which the result have to be written
    # @param showUnidentified an integer:
    #     0 : write only reads for which we have identified individual, locus and previously decribed allele
    #     1 : write only reads for which we have identified individual and locus
    #     2 : write all
    @staticmethod
    def writeTo(path,showUnidentified=0):
	o = open(path,'w')
	o.write("seqName,strand,individual,locus,allele")
	for currRead in Read._reads:
	    currRead.oneWriteTo(o,showUnidentified)
	 
    # Write this read identification
    # @param output a file handler
    # @param showUnidentified an integer (see writeTo doc)
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
            if self.rev:
                sens = "R"
            else:
                sens = "F"
	    output.write("\n" + self.name + "," + sens + "," + indivName + "," + locusName + "," + alleleName)
    
    # Static: identify all the reads
    # @param IndividualClass the Individual class to which link the identidified individuals
    @staticmethod
    def identify(IndividualClass):
        identifiedLoci = 0
        identifiedIndividuals = 0
	for currRead in Read._reads:
	    currRead.oneIdentify(IndividualClass)
	    if currRead.locus != None:
                identifiedLoci += 1
            if currRead.individual != None:
                identifiedIndividuals += 1
        return((identifiedLoci,identifiedIndividuals))
    
    # Identify the current sequence
    # IE: a read is associted to an individual and a locus
    # @param IndividualClass the Individual class to which link the identidified individuals
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
	    self.individual.raiseSeqNr()
	    
    # Static: associate each identified read to an allele
    # @param alleles the list of alleles (see Individual._alleles)
    # @param allelesSortedBySize the list of alleles sorted
    #                   by size(see Individual._allelesSortedBySize)
    # @param minNumberOfSeqsPerIndividual minimum number of sequences
    #                  per individual to the sequence to be kept
    #                  as a new allele
    @staticmethod
    def match(alleles,allelesSortedBySize,minNumberOfSeqsPerIndividual):
        numberMatching = 0
	for currRead in Read._reads:
	    if currRead.individual != None:
                if currRead.allele == None:
                    currRead.oneMatch(alleles,allelesSortedBySize,minNumberOfSeqsPerIndividual)
            if currRead.allele != None:
                numberMatching += 1
        return(numberMatching)
    
    # Associate this read to an allele
    # @param alleles the list of alleles (see Individual._alleles)
    # @param allelesSortedBySize the list of alleles sorted
    #                   by size(see Individual._allelesSortedBySize)
    # @param minNumberOfSeqsPerIndividual minimum number of sequences
    #                  per individual to the sequence to be kept
    #                  as a new allele
    def oneMatch(self,alleles,allelesSortedBySize,minNumberOfSeqsPerIndividual):
        if self.individual.getSeqNr() < minNumberOfSeqsPerIndividual:
            return
	currAlleles = alleles[self.locus]
	for currAlleleName in allelesSortedBySize[self.locus]:
	    if currAlleles[currAlleleName] in self.seq:
		self.allele = currAlleleName
		break
	if(self.allele == None):
	    self.individual.addUnknownAllele(self.locus,self.seq[6:-6])

	 
   
		