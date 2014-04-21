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

import read
import individual
import sys

print("       __               \n _  | |_  o __  _| _  __\n(_| | |   | | |(_|(/_ | \n")

INIVERSION = 5


# mini config parser inspired from http://www.decalage.info/fr/python/configparser
# (see the comment by “trinar”)

class ParseINI(dict):
  def __init__(self, f):
    self.f = f
    self.__read()

  def __read(self):
    with open(self.f, 'r') as f:
      soi = self
      for line in f:
        if not line.startswith("#") and not line.startswith(';') and line.strip() != "":
          line = line.replace('=', ':')
          line = line.replace(';', '#')
          index = line.find('#')
          line = line[:index]
          line = line.strip()
          if line.startswith("["):
            sections = line[1:-1].split('.')
            soi = self
            for section in sections:
              if section not in soi:
                soi[section] = {}
              soi = soi[section]
          else:
            if not self:
              soi['global'] = {}
              soi = soi['global']
            parts = line.split(":", 1)
            soi[parts[0].strip()] = parts[1].strip()

  def items(self, section):
    try:
      return self[section]
    except KeyError:
      return []
      
def minNumberOfSeqsPerIndividualctToList(pdict):
    rlist = []
    keys = pdict.keys()
    keys.sort()
    for currEntry in keys:
        rlist.append(pdict[currEntry])
    return(rlist)

try:
    ini = ParseINI('settings.ini')
except IOError, e:
    print("You have no configuration file. Please copy settings.sample.ini to settings.ini and update the settings.")
    sys.exit()
    
if int(ini['global']['iniVersion']) < INIVERSION: 
    print("Your version of settings is too old. Please copy settings.sample.ini to settings.ini and update the settings.")
    sys.exit()
    
    
#FIXME: for compat purpose, might be deleted.

minNumberOfSeqsPerIndividual = 1



# output
# ======
# Output file in which the results will be written

resultFile = ini['Files']['fastaSequences'].split('.')[0] +"_result.csv"

### SOME FUNCTIONS ###

def suffixFile(filename,suffix):
    filename = filename.split(".")
    begining = filename[:-1]
    for i in range(len(begining)-1):
        begining.insert(i*2+1,".")
    extension = filename[-1]
    return("".join(begining) + suffix + "." + extension)

def printAndLog(log,string):
    print(string)
    log.append(string+"\n")
    
def printTheLog(log):
    print(''.join(log))

### EXECUTING THE IDENTIFICATION PROCESS ###
### all these function are documented in the corresponding class file

log = list()

print("Loading sequences from file " + ini['Files']['fastaSequences'] + "…"),
sys.stdout.flush()
read.Read.loadFromFile(ini['Files']['fastaSequences'])
print("    [DONE]")
printAndLog(log,"*** Input data: " +  str(read.Read.getNumberOfReads()) + " sequences loaded.")


print("Loading individuals from file " + ini['Files']['tags'] + "…"),
sys.stdout.flush()
individual.Individual.loadFromFile(ini['Files']['tags'])
print("    [DONE]")

print("Loading loci markers…"),
sys.stdout.flush()
print(ini['Primers'])
individual.Individual.setLociMarkers(ini['Primers'])
print("    [DONE]")


print("Sequences are now being associated to loci/individual according to their markers/tags…"),
sys.stdout.flush()
(identifiedLoci,identifiedIndividuals) = read.Read.identify(individual.Individual)
print("    [DONE]")

printAndLog(log,"*** Identified Loci: " + str(identifiedLoci) + "/" + str(read.Read.getNumberOfReads()) + " (" + str(int(100*identifiedLoci / float(read.Read.getNumberOfReads()))) + "%)" )
printAndLog(log,"*** Identified Individuals: " + str(identifiedIndividuals) + "/" + str(read.Read.getNumberOfReads()) + " (" + str(int(100*identifiedIndividuals / float(read.Read.getNumberOfReads()))) + "%)" )


## getting alleles files


print("Loading alleles from " + str(ini['Files']['Alleles'].values()) +  " file…"),
sys.stdout.flush()

allelesFilesF = []

allelesFiles = minNumberOfSeqsPerIndividualctToList(ini['Files']['Alleles'])
for cFile in allelesFiles:
    allelesFilesF.append(open(cFile,"r"))
    
individual.Individual.loadLociFromFiles(allelesFilesF)

#now alleles are loaded, sorting them by decreasing size
individual.Individual.sortAllelesBySize(individual.Individual._alleles,individual.Individual._allelesSortedBySize)

for currFile in allelesFilesF:
    currFile.close()

print("    [DONE]")

print("Sequences are now being associated to allelles…"),
sys.stdout.flush()
numberMatching = read.Read.match(individual.Individual._alleles,individual.Individual._allelesSortedBySize,int(minNumberOfSeqsPerIndividual))
print("    [DONE]")
printAndLog(log,"*** Matching Read (sequences assiociated with a known allele): " + str(numberMatching) + "/" + str(read.Read.getNumberOfReads()) + " (" + str(int(100*numberMatching / float(read.Read.getNumberOfReads()))) + "%)" )


if ini['AlleleDiscovering']['discovering'].upper() == "TRUE":
    print("Discovering new alleles…"),
    sys.stdout.flush()
    
    newAllelesFilesF = []
    
    if ini['AlleleDiscovering']['toNewFiles'].upper() == "TRUE":  
        for cFile in allelesFiles:
            newAllelesFilesF.append(open(suffixFile(cFile,"_new"),"w"))
    else:
        for cFile in allelesFiles:
            newAllelesFilesF.append(open(cFile,"a"))
    
    
    individual.Individual.discoverNewAlleles(newAllelesFilesF,minNumberOfSeqsPerIndividualctToList(ini['AlleleDiscovering']['CropLengths']),int(ini['AlleleDiscovering']['threshold']))
    
    for currFile in newAllelesFilesF:
        currFile.close()
    
    print("    [DONE]")
    
    individual.Individual.sortAllelesBySize(individual.Individual._newAlleles,individual.Individual._newAllelesSortedBySize)
    
    
    print("Unidentified sequences are now being associated to new allelles…"),
    sys.stdout.flush()
    numberMatching = read.Read.match(individual.Individual._newAlleles,individual.Individual._newAllelesSortedBySize,int(minNumberOfSeqsPerIndividual))
    print("    [DONE]")
    print ("*** Matching Read : " + str(numberMatching) + "/" + str(read.Read.getNumberOfReads()) + " (" + str(int(100*numberMatching / float(read.Read.getNumberOfReads()))) + "%)" )
    

print("Writing result to file " + resultFile +"…"),
sys.stdout.flush()

read.Read.writeTo(resultFile,int(ini['Results']['showUnidentified']),)
print("    [DONE]")

print("\n--- STATS ---")

printTheLog(log)
