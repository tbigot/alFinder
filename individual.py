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


import read
import sys
import re
import operator


## Class Individual
#
#  The object Individual is designed to collect information about one individual.
#  One individual correspond to one individual line in the tag file.
#  In the static attributes, this class also manage general information about the process,
#  such as the list of the tags and the allele sequences.

class Individual:
    
    # static: list of all the individuals attached to this class
    _individuals = []
    
    # static: dictionnary of all the pair of tags,
    # associating them to an individual
    # (locus, (tag1, tag2)) -> individual
    _tags = {}
    
    # static: number of loci managed by the process
    _nrLoci = 0
    
    # static: list of the different alleles by locus in the form:
    # [alleles_of_locus1, alleles_of_locus2, ... alleles_of_locusN] with N = _nrLoci
    # with
    # allele_of_locus_x a dictionnary in the form:
    # sequence_name -> sequence
    _alleles = []
    
    # static: the same data, but sorted by decreasing size
    _allelesSortedBySize = []
    
    # static: list of the regular expression matching to the locus associated
    # to each primer in the form:
    # [re_locus_1, re_locus_2, ... re_of_locus_N] with N = _nrLoci
    # with re_locus_x a string written in RE syntax matching with 
    # init_tag_length + the degenerate starting primer + a sequence
    #     + the degenerate ending primer + terminal_tag_length
    _lociRE = []

    # static: dictionnary associating a list of new alleles per individual
    # in the form:
    # individual -> [list of new alleles]
    _newAllelesPerIndividual = {}
    
    # static: list containing all the new detected alleles
    _newAlleles = []
    
    # static: the same, sorted by size, in the decreasing order
    _newAllelesSortedBySize=[]
    
    # Constructor
    # 
    #   instanciate an individual
    #   @param name (a string) the name of the individual
    #   @tags the list of the tags in the form:
    #         [locus1_tag1, locus1_tag2, locus2_tag1, locus2_tag2, ... locusN_tag1, locusN_tag2]
     #        with N = _nrLoci
    def __init__(self,name,tags):
	self.unknownAlleles = {}
	self.seqNr = 0
	self.name = name
	self.tags = []
	if len(tags)%2 != 0:
	    print("Odd number of tags: problems will occur.")
	currAllele = 0
	while currAllele*2 != len(tags):
	    tag1 = tags[currAllele*2]
	    tag2 = read.Read.reverseComplementary(tags[currAllele*2+1])
	    #tag2 = tags[currAllele*2+1]
	    self.tags.append((tag1,tag2))
	    
	    if (currAllele,tag1,tag2) in Individual._tags.keys():
                print("Duplicate entry for allele " + str(currAllele+1) + str((tag1,tag2)))
            else :
                Individual._tags[(currAllele,(tag1,tag2))]=self
	    currAllele += 1
	    
	    
    # Tags accessor
    def getTags(self):
	return (self.tags)
    
    # Name accessor
    def getName(self):
	return self.name
	
    # Static: Accessor of the individual list
    @staticmethod
    def getIndividuals():
	return Individual._individuals
    
    # Static: converting a primer string to the corresponding Python regexp
    # notably using the "wobbles" codes (eg: w = A or T).
    # @param sequence the sequence to convert
    # @return the converted sequence
    @staticmethod
    def seqToRegExp(sequence):
        reLst = []
        for nt in sequence:
            if len(nt) == 1:
                reLst.append(nt[0])
            else:
                reLst.append("["+''.join(nt) + ']')
        reStr = ''.join(reLst)
        return(reStr)
        
                    
    # Static: add to _lociRE (the class list of RE) a pair of markers (= primers)
    # @param markers a pair of primers
    @staticmethod
    def setLociMarkers(markers):
        currPair = 0
        while currPair*2 != len(markers):
            markerF = read.Read.unWobble(markers[str(currPair*2)],"F")
            markerR = read.Read.unWobble(markers[str(currPair*2+1)],"R")
            Individual._lociRE.append(re.compile(".{6}" + Individual.seqToRegExp(markerF)+".+"+Individual.seqToRegExp(markerR) + ".{6}$"))
            currPair += 1
        
	
    # Static: accessor of the tags list
    @staticmethod
    def getTagsList():
	return Individual._tags
	
    # Raise the number of sequences of this individual
    def raiseSeqNr(self):
	self.seqNr += 1
	
    # Accessor of the number of sequences of this individual
    def getSeqNr(self):
        return(self.seqNr)
	
    # Adds a sequence to the list of unknown allele of this individual
    # must specify a peculiar locus
    # @param locus the number of the locus to which
    #           this sequence is associated
    # @param sequence the sequence itself                
    def addUnknownAllele(self,locus,sequence):
	# cropping sequence
	sequence[6:-6]
	if not locus in self.unknownAlleles.keys():
	    self.unknownAlleles[locus] = {}
	if not sequence in self.unknownAlleles[locus].keys():
	    self.unknownAlleles[locus][sequence] = 1
	else:
	    self.unknownAlleles[locus][sequence] += 1
	
    # Static: load individuals from the tags file
    # @param path the path of the tags file
    @staticmethod
    def loadFromFile(path):
	fhandle = open(path)
	for ligne in fhandle:
	    sligne = ligne.rstrip().split(",")
	    Individual._individuals.append(Individual(sligne[0],sligne[1:]))
	Individual._nrLoci = len(Individual._individuals[0].getTags())
	
    # Static: sort alleles and put the result in sortedCurrAlleles
    @staticmethod
    def sortAllelesBySize(alleleSet,sortedAlleleSet):
      for (currAlleles) in alleleSet:
	sortedCurrAlleles = sorted(currAlleles, key=lambda k: len(currAlleles[k]),reverse=True)
	sortedAlleleSet.append(sortedCurrAlleles)
	
    # Static: load the loci list from the files
    # @param files the list of the files in the form
    #        [file_for_allele1, file_for_allele2, ... file_for_alleleN] with N = _nrLoci
    @staticmethod
    def loadLociFromFiles(files):
	if(len(files) != Individual._nrLoci):
	    print("You must provide " + str(Individual._nrLoci) + " files. I cannot proceed with " + str(len(files)) +".")
	    return
	else:
            for currFile in files:
		currName = ''
		currSeq = ''
		result = {}
		for currLigne in currFile:
		    if currLigne.startswith(">"):
			if len(currSeq) != 0:
			    result[currName] = currSeq
			    currSeq = ''
			currName = currLigne.strip()[1:]
		    else:
			currSeq += currLigne.strip().upper()
		result[currName] = currSeq
		Individual._alleles.append(result)
		
    # Static: discover new alleles of all the individuals and write them to a file
    # @param files list of the file in which the new alleles will be collected
    #       in the form
    #       [file_for_allele1, file_for_allele2, ... file_for_alleleN] with N = _nrLoci
    # @param cropLengths a list integer correspondig to the real primers size (defined in the ini file)
    #       in the form
    #       [locus1_size1, locus1_size2, locus2_size1, locus2_size2, ... locusN_size1, locusN_size2]
    # @param threshold (deprecated, might be set to 1) an integer telling the minimum number sequences
    #                   of a new allele to be found to be considered as a real allele

    @staticmethod
    def discoverNewAlleles(files,cropLengths,threshold):
        
        for currIndiv in Individual._individuals:
            currIndiv.oneDiscoverNewAlleles()
	# _newAllelesPerIndividual is now filled
        
        # now, deleting sequences that are below the threshold
        for currLocus in Individual._newAllelesPerIndividual.keys():
            allelesToBeDeleted = []
            for currAllele in Individual._newAllelesPerIndividual[currLocus]:
                if (len(currAllele)-(int(cropLengths[currLocus*2])+int(cropLengths[currLocus*2+1])))%1 != 0 or Individual._newAllelesPerIndividual[currLocus][currAllele] < threshold:
                    allelesToBeDeleted.append(currAllele)
            for currATBD in allelesToBeDeleted:
                Individual._newAllelesPerIndividual[currLocus].pop(currATBD)
        
        
	for i in range(len(Individual._alleles)):
            Individual._newAlleles.append({})
	
	allelleCount = 1
	for locus in Individual._newAllelesPerIndividual.keys():
	    for sequence in Individual._newAllelesPerIndividual[locus].keys():
                Individual._newAlleles[locus]["NewAllele_" + str(allelleCount)] = sequence
		nrIndividuals = Individual._newAllelesPerIndividual[locus][sequence]
		files[locus].write(">NewAllele_" + str(allelleCount) + " count="+ str(nrIndividuals) + " threshold=" + str(threshold))
		
		# cropping sequence
		sequence = sequence[int(cropLengths[locus*2]):-int(cropLengths[locus*2+1])]
		
		for nbChar in range(len(sequence)):
		    if (nbChar % 60 == 0):
			files[locus].write("\n")
		    files[locus].write(sequence[nbChar])
		files[locus].write("\n")
		print("."),
		sys.stdout.flush()
		allelleCount += 1
	    
    # Discovering new alleles for this individual
    def oneDiscoverNewAlleles(self):
	for locus in self.unknownAlleles.keys():
	    for allele in self.unknownAlleles[locus].keys():
                if not locus in Individual._newAllelesPerIndividual.keys():
                    Individual._newAllelesPerIndividual[locus] = {}
                if not allele in Individual._newAllelesPerIndividual[locus].keys():
                    Individual._newAllelesPerIndividual[locus][allele] = self.unknownAlleles[locus][allele]
                    #print("\nIndiv = "+ self.name + " ; " + str(self.unknownAlleles[locus][allele] / float(self.seqNr)) + " - " + str(self.seqNr))
                else:
                    Individual._newAllelesPerIndividual[locus][allele] += self.unknownAlleles[locus][allele]
                    
		   
