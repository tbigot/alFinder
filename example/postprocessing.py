# -*- coding: utf-8 -*-

from math import factorial,pow,log


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

# errorRate: sequencing error rate per nucleotide

errorRate = 0.001

# Calculate all possible combinations of homozygous and heterozygous genotypes for a given amplicon based on the obtained variants per amplicon
unfilteredResults = "likelihood_all_beforeDecision.csv"

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
    if len(line) > 0:
      dataOutput.write(",".join(line)+"\n")
      if ",".join(line).count(",") < 1:
	print(line)
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
  ofhandle = open(currReadSize[0][:-4]+"_filtered.fas","w")
  currSeq = []
  for ligne in fhandle:
      if ligne.startswith('>') or not ligne.endswith('\n'):
	  if not len(currSeq) == 0:
	      currSeqJoined = string.join(currSeq[1:],"")
	      if len(currSeqJoined) <= currReadSize[2] and len(currSeqJoined) >= currReadSize[1]:
		  ofhandle.write(currSeq[0]+"\t"+ str(len(currSeqJoined)) +"\n")
		  for cfl in currSeq[1:]:
		      ofhandle.write(cfl+"\n")
		  readsToBeKept.append(currSeq[0].split()[0][1:])
		  # print(str(len(currSeqJoined))+ " accepted; needed " + str(currReadSize[2]) + " - " + str(currReadSize[1]))
	      #else:
		  # print(str(len(currSeqJoined))+ " rejected; needed " + str(currReadSize[2]) + " - " + str(currReadSize[1]))
	      currSeq = []
      currSeq.append(ligne.strip())


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
    splitCorr = currCorr.strip().split(",")
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
        #updating the nb of most popular nbOfMostPopular
        count[individual][locus][0] = nbOfMostPopular
    if len(count[individual]) == 0:
      count.pop(individual)

writeData("intermediateResults_coverage.csv")

lkfile = open(outputFile,"w")
lkfile.write("individual,locus,homozygote,isSignificant,LRT_value,allele1,allele2,loglk,totalCoverage,allele1Count,allele2Count,otherVariantsSum")
logFile = open(unfilteredResults,"w")
logFile.write("individual,locus,homozygousAllele,heterozygousAllele1,heterozygousAllele2,logLikelihood")

for individual in count.keys():
  # computing homozygosity
  for locus in count[individual].keys():
    bestLogLkHomo = -99999
    for variant in count[individual][locus][1].keys():
      otherVariants = count[individual][locus][1].keys()
      otherVariants.remove(variant)
      variantCount = count[individual][locus][1][variant]
      # getting otherVariantsCounts
      otherVariantsCounts = []
      for currOtherVariant in otherVariants:
        otherVariantsCounts.append(count[individual][locus][1][currOtherVariant])
      ## Computing some numbers
      # total numbers, n in Hohenlohe2010
      ntotal = sum(otherVariantsCounts)+variantCount
      # product of factorials (n1!n2!n3!n4! in Hohenlohe2010)
      prodFact = factorial(variantCount)
      for nx in otherVariantsCounts:
        prodFact *= factorial(nx)
      nbOtherVariants = len(otherVariantsCounts)
      nbVariants = nbOtherVariants+1
      likelihood = (factorial(ntotal)/(prodFact)) * pow((1-nbOtherVariants*errorRate/nbVariants),variantCount)*pow((errorRate/(nbVariants)),sum(otherVariantsCounts))
      if(likelihood != 0):
        print("\nHomozygote " + variant + " " + str(variantCount) + "\nOther variants "+ str(otherVariants) + " " + str(otherVariantsCounts))
        print("lik=" + str(log(likelihood)))
        print("ntotal = " + str(ntotal))
        if log(likelihood) > bestLogLkHomo:
          bestLogLkHomo = log(likelihood)
          bestVariant = variant
          bestVariantCount = variantCount
          bestOtherVariantsSum = sum(otherVariantsCounts)
      logFile.write("\n" + individual + "," + locus + "," + variant + ",,," + str(log(likelihood) if likelihood != 0 else -99999))
    
    # computing homozygosity
    bestLogLkHetero = -99999
    nbVariants = len(count[individual][locus][1])
    if(nbVariants >= 2):
      for i in range(nbVariants):
        for j in range(i+1,nbVariants):
          variant1 = count[individual][locus][1].keys()[i]
          variant2 = count[individual][locus][1].keys()[j]
          otherVariants = count[individual][locus][1].keys()
          otherVariants.remove(variant1)
          otherVariants.remove(variant2)
          variant1Count = count[individual][locus][1][variant1]
          variant2Count = count[individual][locus][1][variant2]
          # getting otherVariantsCounts
          otherVariantsCounts = []
          for currOtherVariant in otherVariants:
            otherVariantsCounts.append(count[individual][locus][1][currOtherVariant])
          
          ## Computing some numbers
          # total numbers, n in Hohenlohe2010
          ntotal = sum(otherVariantsCounts)+variant1Count+variant2Count
          # product of factorials (n1!n2!n3!n4! in Hohenlohe2010)
          prodFact = factorial(variant1Count) * factorial(variant2Count)
          for nx in otherVariantsCounts:
            prodFact *= factorial(nx)
          nbOtherVariants = len(otherVariantsCounts)
          nbVariants = nbOtherVariants+1
          likelihood = (factorial(ntotal)/(prodFact)) * pow((.5-errorRate/nbVariants),(variant1Count+variant2Count))*pow((errorRate/(nbVariants)),sum(otherVariantsCounts))
          if(likelihood != 0):
            print("\nHétérozygote 1 " + variant1 + " " + str(variant1Count) + "\n2 " + variant2 + " "+ str(variant2Count) + "\nOther V "+ str(otherVariants) + " -> " + str(otherVariantsCounts))
            print("lik=" + str(log(likelihood)))
            print("ntotal = " + str(ntotal))
            if log(likelihood) > bestLogLkHetero:
              bestLogLkHetero = log(likelihood)
              bestVariant1 = variant1
              bestVariant2 = variant2
              bestVariant1Count = variant1Count
              bestVariant2Count = variant2Count
              bestOtherVariantsSum = sum(otherVariantsCounts)
          logFile.write("\n" + individual + "," + locus + ",," + variant1 + "," + variant2 + "," + str(log(likelihood) if likelihood != 0 else -99999))
    LRT = 2*max(bestLogLkHomo,bestLogLkHetero) - 2*min(bestLogLkHomo,bestLogLkHetero)
    significance = LRT > 3.84
    homo = (bestLogLkHetero < bestLogLkHomo)
    lkfile.write("\n" + individual + "," + locus + "," + ("1" if homo else "0") + "," + str(significance) + "," + str(LRT)+ "," + (bestVariant if homo else bestVariant1) + "," + ("NA" if homo else bestVariant2)+ "," + str(bestLogLkHomo if homo else bestLogLkHetero) + "," + str(ntotal) + "," + str(bestVariantCount if homo else bestVariant1Count) + "," + str(0 if homo else bestVariant2Count) + "," + str(bestOtherVariantsSum))

