# -*- coding: utf-8 -*-

######################################################################
# SETTINGS HERE: Please modify according to your needs

## This example is based on the method described in the article: 
## Large-scale genotyping by next generation sequencing: how to overcome the challenges to reliably genotype individuals?
## Ferrandiz-Rovira M, Bigot T, Allainé D, Callait-Cardinal M-P, Radwan J, Cohas A.

## Output data file from alFinder

inputFile = "data_result.csv"

# Elimination of reads too short or too long
## For each locus:
##    1/specify the names of the file containing the alleles sequences (previously decribed and potential new alleles) (e.g. alleles_locus1.fas)
##    2/specify min: the minimum length of reads to be kept
##    3/specify max: the maximum length of reads to be kept

readSizes = (
  ("alleles_locus1.fas", 192,213),    # 95% and 105% of Mama-DRB1 alleles length
  ("alleles_locus2.fas", 218,241),    # 95% and 105% of Mama-DRB2 alleles length
  ("alleles_locus3.fas", 166,184),    # 95% and 105% of Mama-UB alleles length
  ("alleles_locus4.fas", 171,189)     # 95% and 105% of Mama-UD alleles length
  )

# Should new variants (e.g. indels) be renamed as previously described alleles 

stepRename = True
# If true provide a file with the correspondence between new variants (e.g. indels) names and names of previously decribed alleles
correspondenceFile = "correspondenceData.csv"


# Elimination of individuals with less than 12 reads

minSeqPerIndividualPerLocus = 12

# Within an individuals and a locus: elimination of variants with less than 3 reads per variant

minSeqPerIndividualPerLocusPerVariant = 3

# Within an individuals and a locus: elimination of variants present in less than 0.25% of the variant with the maximum number of reads

minProportionOfMostPopularVariantPerIndividual = 0.25

## Output file

outputFile = "data_filteredResults.csv"

########################################################################
# !!! DO NOT MODIFY BELOW THIS LINE !!!!
########################################################################

# Before anything, defining an output function to save the data
# (intermediate)

def writeData(filename):
  dataOutput = open(filename,"w")
  dataOutput.write("seqName,strand,individual,locus,allele\n")
  for line in data:
    dataOutput.write(",".join(line)+"\n")
  dataOutput.close()




# After alFinder (available codes at https://github.com/tbigot/alFinder)
###################
### Reading the file
dataFile = open(inputFile)

data = list()

for currLine in dataFile:
  currLineSplit = currLine.strip().split(",")
  data.append(currLineSplit)

#removing the headline
data.pop(0)


### Elimination of singletons
cleanedData = list()
for currData in data:
  if not currData[4].startswith("<unidentified>"):
    cleanedData.append(currData)

data = cleanedData

writeData("intermediateResults_no_singleton.csv")


##################
### Elimination of reads too short or too long
import sys,string

readsToBeKept = list()

for currReadSize in readSizes:
  fhandle = open(currReadSize[0])
  ofhandle = open(currReadSize[0]+"_filtered","w")
  currSeq = []
  for ligne in fhandle:
      if ligne.startswith('>') or not ligne.endswith('\n'):
	  if not len(currSeq) == 0:
	      currSeqJoined = string.join(currSeq[1:],"")
	      if len(currSeqJoined) <= currReadSize[2] and len(currSeqJoined) >= currReadSize[1]:
		  ofhandle.write(currSeq[0]+"\t"+ str(len(currSeqJoined)) +"\n")
		  for cfl in currSeq[1:]:
		      ofhandle.write(cfl+"\n")
		  readsToBeKept.append(currSeq[0].split()[0])
	      currSeq = []
      currSeq.append(ligne.strip()[1:])


#Now cleaning the data
cleanedData = list()
for currData in data:
  if currData[4].split()[0] in readsToBeKept:
    cleanedData.append(currData)
    
data = cleanedData

writeData("intermediateResults_size.csv")


#######################
### Rename new variants (e.g. indels) as previously described alleles
if stepRename:
  corresp = dict()
  correspFile = open(correspondenceFile)
  for currCorr in correspFile:
    splitCorr = currCorr.split(",")
    corresp[splitCorr[0]] = splitCorr[1]
  for currRead in data:
    variantId = currRead[4].split()
    if variantId[0] in corresp.keys():
      variantId[0] = corresp[variantId[0]]
      currRead[4] = ' '.join(variantId)
      

writeData("intermediateResults_indels.csv")


#######################
### Counting numbers of variant / locus / individual
count = dict()
# count {"indivName" −> {locus}}
# with locus {"locus" -> [nb_seq, variants]}
# with variants {"variant" -> nb_seq}

for currRead in data:
  indivName = currRead[2]
  locus = currRead[3]
  variant = currRead[4].split()[0]
  
  if not indivName in count.keys():
    count[indivName] = dict()
  currIndiv = count[indivName]
  if not locus in currIndiv.keys():
    currIndiv[locus] = [ 1, dict()]
  else:
    currIndiv[locus][0] += 1
  currLocusDict = currIndiv[locus][1]
  if not variant in currLocusDict:
    currLocusDict[variant] = 1
  else:
    currLocusDict[variant] += 1

# now counting is done, filtering with
# minSeqPerIndividualPerLocus
# AND
# minSeqPerIndividualPerLocusPerVariant

# WARNING: count will become :
# count {"indivName" −> {locus}}
# with locus {"locus" -> [nbOfmostPopularVariant, variants]}
# with variants {"variant" -> nb_seq}

for individual in count.keys():
  for locus in count[individual].keys():
    nbSeqForThisLocus = count[individual][locus][0]
    if nbSeqForThisLocus < minSeqPerIndividualPerLocus:
      count[individual].pop(locus)
    else:
      nbOfMostPopular = 0
      for variant in count[individual][locus][1].keys():
        nbSeqForThisVariant = count[individual][locus][1][variant]
        if nbSeqForThisVariant < minSeqPerIndividualPerLocusPerVariant:
          count[individual][locus][1].pop(variant)
        else:
          if nbSeqForThisVariant > nbOfMostPopular:
            nbOfMostPopular = nbSeqForThisVariant
      if len(count[individual][locus][1]) == 0:
        count[individual].pop(locus)
      else:
        count[individual][locus][0] = nbOfMostPopular
    if len(count[individual]) == 0:
      count.pop(individual)

writeData("intermediateResults_coverage.csv")



# second round, according to minProportionOfMostPopularVariantPerIndividual
for individual in count.keys():
  for locus in count[individual].keys():
    nbOfMostPopular = float(count[individual][locus][0])
    for variant in count[individual][locus][1].keys():
      if (float(count[individual][locus][1][variant]) / nbOfMostPopular) < float(minProportionOfMostPopularVariantPerIndividual):
        count[individual][locus][1].pop(variant)
    if len(count[individual][locus][1]) == 0:
      count[individual].pop(locus)
  if len(count[individual]) == 0:
    count.pop(individual)

# outputing the sequences files
writeData("intermediateResults_percentage.csv")



### outputing results, not sequences but counts
outputHandler = open(outputFile,"w")
outputHandler.write("individual,locus,variant,count\n")
for individual in count.keys():
  for locus in count[individual].keys():
    for variant in count[individual][locus][1].keys():
      outputHandler.write(individual + "," + locus + "," + variant + "," + str(count[individual][locus][1][variant])+"\n")
outputHandler.close()
