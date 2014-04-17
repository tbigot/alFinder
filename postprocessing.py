# -*- coding: utf-8 -*-

######################################################################
# SETTINGS HERE: Please modify according to your needs

## Output data file from alFinder

inputFile = "alFinder_data.csv"

# Elimination of reads too short or too long
## For each locus:
##    1/specify the names of the file containing the alleles sequences (previously decribed and potential new alleles) (e.g. alleles_locus1.fas)
##    2/specify min: the minimum length of reads to be kept
##    3/specify max: the maximum length of reads to be kept

readSizes = (
  ("alleles_locus1.fas", min,max),
  ("alleles_locus1.fas",min,max)
  )

# Should new variants (e.g. indels) be renamed as previously described alleles 

stepRename = True
# If true provide a file with the correspondence between new variants (e.g. indels) names and names of previously decribed alleles
correspondenceFile = "correspondenceData.csv"


# Elimination of individuals with less than xx reads

minSeqPerIndividualPerLocus = xx	

# Within an individuals and a locus: elimination of variants with less than yy reads per variant

minSeqPerIndividualPerLocusPerVariant = yy

# Within an individuals and a locus: elimination of variants present in less than zz% of the variant with the maximum number of reads

minProportionOfMostPopularVariantPerIndividual = zz

## Output file

outputFile = "filteredResults.csv"

########################################################################
# !!! DO NOT MODIFY BELOW THIS LINE !!!!
########################################################################

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

##################
### Elimination of reads too short or too long
import sys,string

readsToBeKept = list()

for currReadSize in readSizes:
  fhandle = open(currReadSize[0])
  ofHandle = open(currReadSize[0]+"_filtered","w")
  currSeq = []
  for ligne in fhandle:
      if ligne.startswith('>') or not ligne.endswith('\n'):
	  if not len(currSeq) == 0:
	      currSeqJoined = string.join(currSeq[1:],"")
	      if len(currSeqJoined) <= currReadSize[2] and len(currSeqJoined) >= currReadSize[1]:
		  ofhandle.write(currSeq[0]+"\t"+ str(len(currSeqJoined)) +"\n")
		  for cfl in currSeq[1:]:
		      ofhandle.write(cfl+"\n")
		  print("Kept sequence "+currSeq[0]+" : length = "+ str(len(currSeqJoined)))
		  readsToBeKept.append(currSeq[0])
	      else:
		  print("Rejected sequence "+currSeq[0]+" : length = "+ str(len(currSeqJoined)))
	      currSeq = []
      currSeq.append(ligne.strip())

#Now cleaning the data
cleanedData = list()
for currData in data:
  if currData[4].split()[0] in readsToBeKept:
    cleanedData.append(currData)
    
data = cleanedData

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

### outputing results
outputHandler = open(outputFile,"w")
for individual in count.keys():
  for locus in count[individual].keys():
    for variant in count[individual][locus][1].keys():
      outputHandler.write(individual + "," + locus + "," + variant + "," + str(count[individual][locus][1][variant]))
