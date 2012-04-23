# -*- coding: utf-8 -*-
# thomas.bigot@tmgo.net


import read
import individual
import sys


### VARIABLES ###

# séquences Fasta
# ===============
# pour un run donné

fastaSequencesFile = "data/1.fas"


# tags d’identification des individus et des locus
# ================================================
# une ligne par individu, puis des paires de tags (brin droit, reverse)
# pour chaque locus, dans un ordre déterminé
# Exemple:
# individu1,tagLocus1Forward,tagLocus1Reverse,tagLocus2Forward,tagLocus2Reverse
# soit n le nombre de locus, on a ici 2n tags

tagsFile = "data/tagsRun1.csv"

# allèles
# =======
# liste de n fichiers (n étant le nombre de locus)
# TRÈS IMPORTANT:
#    1/ Il doit y avoir strictement n fichiers, n étant fixé par le fichier précédent (tags)
#    2/ Les fichiers doivent être dans l’ordre des locus de l’étape précédente (tags)
# dans chaque fichier, on doit avoir autant de séquences fasta que d’Allèles

allelesFiles = ["data/Mama-DRB1.fas","data/Mama-DRB2.fas","data/Mama-UB.fas","data/Mama-UD.fas"]

# sortie
# ======
# Fichier de sortie où seront écrits les résultats

resultFile = fastaSequencesFile.split('.')[0] +"_result.csv"

# doit-on afficher dans les résultats les allèles non identifiés
# [0,1,2]
# 0 : n’afficher que les séquences identifiées à l’allèle connu
# 1 : n’afficher que les séquences identifiées, allèle potentiellement inconnu
# 2 : tout afficher

showUnidentified = 1

# allele discovering
# ==================

# écrit les allèles dans le fichier d’allèles
alleleDiscovering = True

# seuil nb seq nouvelle de l’individu / nb seq total de l’individu
threshold = 0.2


### EXÉCUTION DES FONCTIONS ###

print("Loading sequences from file " + fastaSequencesFile + "…"),
sys.stdout.flush()
read.Read.loadFromFile(fastaSequencesFile)
print("    [DONE]")

print("Loading individuals from file " + tagsFile + "…"),
sys.stdout.flush()
individual.Individual.loadFromFile(tagsFile)
print("    [DONE]")

print("Sequences are now being associated to individuals/loci according to their tags…"),
sys.stdout.flush()
read.Read.identify(individual.Individual.getTagsList())
print("    [DONE]")

print("Loading alleles from " + str(len(allelesFiles)) +  " file…"),
sys.stdout.flush()
individual.Individual.loadLociFromFiles(allelesFiles)
print("    [DONE]")

print("Sequences are now being associated to allelles…"),
sys.stdout.flush()
read.Read.match(individual.Individual._alleles)
print("    [DONE]")

print("Writing result to file " + resultFile +"…"),
sys.stdout.flush()
read.Read.writeTo(resultFile,showUnidentified)
print("    [DONE]")

if alleleDiscovering:
    print("Discovering new alleles…"),
    sys.stdout.flush()
    individual.Individual.discoverNewAlleles(allelesFiles,threshold)
    print("    [DONE]")