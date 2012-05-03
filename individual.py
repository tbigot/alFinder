# thomas.bigot@tmgo.net

import read
import sys

class Individual:
    
    _individuals = []
    _tags = {}
    _nrLoci = 0
    _alleles = []
    
    _newAllelesPerIndividual = {}
    _newAlleles = []
    
    def __init__(self,name,tags):
	self.unknownAlleles = {}
	self.seqNr = 0
	self.name = name
	self.tags = []
	if len(tags)%2 != 0:
	    print("Odd number of tags: problems will occur.")
	currPair = 0
	while currPair != len(tags)/2:
	    tag1 = tags[currPair*2]
	    tag2 = read.Read.reverseComplementary(tags[currPair*2+1])
	    self.tags.append((tag1,tag2))
	    Individual._tags[(tag1,tag2)]=(self,currPair)
	    currPair += 1
	    
	    
	
    def getTags(self):
	return (self.tags)
    
    def getName(self):
	return self.name
	
    @staticmethod
    def getIndividuals():
	return Individual._individuals
	
    @staticmethod
    def getTagsList():
	return Individual._tags
	
    def raiseSeqNr(self):
	self.seqNr += 1
	
    def addUnknownAllele(self,locus,sequence):
	if not locus in self.unknownAlleles.keys():
	    self.unknownAlleles[locus] = {}
	if not sequence in self.unknownAlleles[locus].keys():
	    self.unknownAlleles[locus][sequence] = 1
	else:
	    self.unknownAlleles[locus][sequence] += 1
	
    @staticmethod
    def loadFromFile(path):
	fhandle = open(path)
	for ligne in fhandle:
	    sligne = ligne.rstrip().split(",")
	    Individual._individuals.append(Individual(sligne[0],sligne[1:]))
	Individual._nrLoci = len(Individual._individuals[0].getTags())
	
    @staticmethod
    def loadLociFromFiles(files):
	if(len(files) != Individual._nrLoci):
	    print("You must provide " + str(Individual._nrLoci) + " files. I cannot proceed with " + str(len(files)) +".")
	    return
	else:
	    for currPath in files:
		currFile = open(currPath)
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
		
    @staticmethod
    def discoverNewAlleles(files,threshold):
	openFiles = []
        for file in files:
            openFiles.append(open(file,'w'))           
        for currIndiv in Individual._individuals:
            currIndiv.oneDiscoverNewAlleles(threshold)
	# _newAllelesPerIndividual is now filled

	for i in range(len(Individual._alleles)):
            Individual._newAlleles.append({})
	
	allelleCount = 1
	for locus in Individual._newAllelesPerIndividual.keys():
	    for sequence in Individual._newAllelesPerIndividual[locus].keys():
                Individual._newAlleles[locus]["NewAllele_" + str(allelleCount)] = sequence
		nrIndividuals = Individual._newAllelesPerIndividual[locus][sequence]
		openFiles[locus].write(">NewAllele_" + str(allelleCount) + " nbIndividuals="+ str(nrIndividuals) + " threshold=" + str(threshold))
		for nbChar in range(len(sequence)):
		    if (nbChar % 60 == 0):
			openFiles[locus].write("\n")
		    openFiles[locus].write(sequence[nbChar])
		openFiles[locus].write("\n")
		print("."),
		sys.stdout.flush()
		allelleCount += 1
	    
    def oneDiscoverNewAlleles(self,threshold):
	for locus in self.unknownAlleles.keys():
	    for allele in self.unknownAlleles[locus].keys():
		if self.unknownAlleles[locus][allele] >= threshold:
		    if not locus in Individual._newAllelesPerIndividual.keys():
			Individual._newAllelesPerIndividual[locus] = {}
		    if not allele in Individual._newAllelesPerIndividual[locus].keys():
			Individual._newAllelesPerIndividual[locus][allele] = 1
			#print("\nIndiv = "+ self.name + " ; " + str(self.unknownAlleles[locus][allele] / float(self.seqNr)) + " - " + str(self.seqNr))
		    else:
			Individual._newAllelesPerIndividual[locus][allele] += 1
			
		   
