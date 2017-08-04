#Clear workspace
rm(list = ls())

#Plyr must be loaded after dplyr (different definition of count() function)

library(dplyr)
library(plyr)
library(stringr)
library(zoo)
library(devtools)
library(PheWAS)

#Assign names of data files to be read
casesToRead <- 'HLCases02.csv'
controlsToRead <- 'HLControls02.csv'
ICD9CodesRaw <- 'ICD9codes02.txt'

source("loadShapeData.r")
source("PheWAS_General_Setup.r")
source("stratifyTemporal.r")
source("phenomeWideFisherTest.r")
source("FisherPheWASCorTest.r")

#***All cases and controls analysis***

#Unadjusted
runPheWAS(phenotypes, attributes, "AssociationsUnadjusted",FALSE)
#Adjusted
runPheWAS(phenotypes, attributes, "AssociationsAdjusted",TRUE)

#Fisher's Exact Test
pheomeWideFisherTest(phenotypes,attributes,"AssociationsFisher")

#Stratify attribute data
#Race
attributesRaceWhite <- attributes[which(attributes$RACE=='W'),]
attributesRaceBlack <- attributes[which(attributes$RACE=='B'),]

attributesRaceWhite$RACE <- NULL
attributesRaceBlack$RACE <- NULL

#Sex
attributesSexMale <- attributes[which(attributes$SEX=='M'),]
attributesSexFemale <- attributes[which(attributes$SEX=='F'),]

attributesSexMale$SEX <- NULL
attributesSexFemale$SEX <- NULL

#Race and Sex
attributesRaceWhiteSexMale <- attributesRaceWhite[which(attributesRaceWhite$SEX=='M'),]
attributesRaceBlackSexMale <- attributesRaceBlack[which(attributesRaceBlack$SEX=='M'),]
attributesRaceWhiteSexFemale <- attributesRaceWhite[which(attributesRaceWhite$SEX=='F'),]
attributesRaceBlackSexFemale <- attributesRaceBlack[which(attributesRaceBlack$SEX=='F'),]

attributesRaceWhiteSexMale$SEX <- NULL
attributesRaceBlackSexMale$SEX <- NULL
attributesRaceWhiteSexFemale$SEX <- NULL
attributesRaceBlackSexFemale$SEX <- NULL

#PheWAS
#Stratified by Race
runPheWAS(phenotypes,attributesRaceWhite,"AssociationsUnadjustedRaceWhite",FALSE)
runPheWAS(phenotypes,attributesRaceWhite,"AssociationsAdjustedRaceWhite",TRUE)
runPheWAS(phenotypes,attributesRaceBlack,"AssociationsUnadjustedRaceBlack",FALSE)
runPheWAS(phenotypes,attributesRaceBlack,"AssociationsAdjustedRaceBlack",TRUE)
#Stratified by Sex
runPheWAS(phenotypes,attributesSexMale,"AssociationsUnadjustedSexMale",FALSE)
runPheWAS(phenotypes,attributesSexMale,"AssociationsAdjustedSexMale",TRUE)
runPheWAS(phenotypes,attributesSexFemale,"AssociationsUnadjustedSexFemale",FALSE)
runPheWAS(phenotypes,attributesSexFemale,"AssociationsAdjustedSexFemale",TRUE)
#Stratifed by both Race and Sex
runPheWAS(phenotypes,attributesRaceWhiteSexMale,"AssociationsUnadjustedRaceWhiteSexMale",FALSE)
runPheWAS(phenotypes,attributesRaceWhiteSexMale,"AssociationsAdjustedRaceWhiteSexMale",TRUE)
runPheWAS(phenotypes,attributesRaceWhiteSexFemale,"AssociationsUnadjustedRaceWhiteSexFemale",FALSE)
runPheWAS(phenotypes,attributesRaceWhiteSexFemale,"AssociationsAdjustedRaceWhiteSexFemale",TRUE)
runPheWAS(phenotypes,attributesRaceBlackSexMale,"AssociationsUnadjustedRaceBlackSexMale",FALSE)
runPheWAS(phenotypes,attributesRaceBlackSexMale,"AssociationsAdjustedRaceBlackSexMale",TRUE)
runPheWAS(phenotypes,attributesRaceBlackSexFemale,"AssociationsUnadjustedRaceBlackSexFemale",FALSE)
runPheWAS(phenotypes,attributesRaceBlackSexFemale,"AssociationsAdjustedRaceBlackSexFemale",TRUE)

#***Before Temporal Stratification***
newPhenotypesBeforeHL <- stratifyTemporal('B',phenotypes,uniquePhecodeMappedSelected)

#Unadjusted
runPheWAS(newPhenotypesBeforeHL,attributes,"AssociationsUnadjustedTemporalBefore",FALSE)
#Adjusted
runPheWAS(newPhenotypesBeforeHL,attributes,"AssociationsAdjustedTemporalBefore",TRUE)
#Stratified by Race
runPheWAS(newPhenotypesBeforeHL,attributesRaceWhite,"AssociationsUnadjustedTemporalBeforeRaceWhite",FALSE)
runPheWAS(newPhenotypesBeforeHL,attributesRaceWhite,"AssociationsAdjustedTemporalBeforeRaceWhite",TRUE)
runPheWAS(newPhenotypesBeforeHL,attributesRaceBlack,"AssociationsUnadjustedTemporalBeforeRaceBlack",FALSE)
runPheWAS(newPhenotypesBeforeHL,attributesRaceBlack,"AssociationsAdjustedTemporalBeforeRaceBlack",TRUE)
#Stratified by Sex
runPheWAS(newPhenotypesBeforeHL,attributesSexMale,"AssociationsUnadjustedTemporalBeforeSexMale",FALSE)
runPheWAS(newPhenotypesBeforeHL,attributesSexMale,"AssociationsAdjustedTemporalBeforeSexMale",TRUE)
runPheWAS(newPhenotypesBeforeHL,attributesSexFemale,"AssociationsUnadjustedTemporalBeforeSexFemale",TRUE)
runPheWAS(newPhenotypesBeforeHL,attributesSexFemale,"AssociationsAdjustedTemporalBeforeSexFemale",TRUE)
#Stratifed by both Race and Sex
runPheWAS(newPhenotypesBeforeHL,attributesRaceWhiteSexMale,"AssociationsUnadjustedTemporalBeforeRaceWhiteSexMale",FALSE)
runPheWAS(newPhenotypesBeforeHL,attributesRaceWhiteSexMale,"AssociationsAdjustedTemporalBeforeRaceWhiteSexMale",TRUE)
runPheWAS(newPhenotypesBeforeHL,attributesRaceWhiteSexFemale,"AssociationsUnadjustedTemporalBeforeRaceWhiteSexFemale",FALSE)
runPheWAS(newPhenotypesBeforeHL,attributesRaceWhiteSexFemale,"AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale",TRUE)
runPheWAS(newPhenotypesBeforeHL,attributesRaceBlackSexMale,"AssociationsUnadjustedTemporalBeforeRaceBlackSexMale",FALSE)
runPheWAS(newPhenotypesBeforeHL,attributesRaceBlackSexMale,"AssociationsAdjustedTemporalBeforeRaceBlackSexMale",TRUE)
runPheWAS(newPhenotypesBeforeHL,attributesRaceBlackSexFemale,"AssociationsUnadjustedTemporalBeforeRaceBlackSexFemale",FALSE)
runPheWAS(newPhenotypesBeforeHL,attributesRaceBlackSexFemale,"AssociationsAdjustedTemporalBeforeRaceBlackSexFemale",TRUE)

#Fisher's Exact Test

#All
pheomeWideFisherTest(newPhenotypesBeforeHL,attributes,"AssociationsFisherTemporalBefore")

#Race
pheomeWideFisherTest(newPhenotypesBeforeHL,attributesRaceBlack,"AssociationsFisherTemporalBeforeRaceBlack")
pheomeWideFisherTest(newPhenotypesBeforeHL,attributesRaceWhite,"AssociationsFisherTemporalBeforeRaceWhite")

#Sex
pheomeWideFisherTest(newPhenotypesBeforeHL,attributesSexMale,"AssociationsFisherTemporalBeforeSexMale")
pheomeWideFisherTest(newPhenotypesBeforeHL,attributesSexFemale,"AssociationsFisherTemporalBeforeSexFemale")

#Race and Sex
pheomeWideFisherTest(newPhenotypesBeforeHL,attributesRaceWhiteSexMale,"AssociationsFisherTemporalBeforeRaceWhiteSexMale")
pheomeWideFisherTest(newPhenotypesBeforeHL,attributesRaceWhiteSexFemale,"AssociationsFisherTemporalBeforeRaceWhiteSexFemale")
pheomeWideFisherTest(newPhenotypesBeforeHL,attributesRaceBlackSexMale,"AssociationsFisherTemporalBeforeRaceBlackSexMale")
pheomeWideFisherTest(newPhenotypesBeforeHL,attributesRaceBlackSexFemale,"AssociationsFisherTemporalBeforeRaceBlackSexFemale")

#***After Temporal Stratification***
newPhenotypesAfterHL <- stratifyTemporal('A',phenotypes,uniquePhecodeMappedSelected)

#Unadjusted
runPheWAS(newPhenotypesAfterHL,attributes,"AssociationsUnadjustedTemporalAfter",FALSE)
#Adjusted
runPheWAS(newPhenotypesAfterHL,attributes,"AssociationsAdjustedTemporalAfter",TRUE)

#Stratified by Race
runPheWAS(newPhenotypesAfterHL,attributesRaceWhite,"AssociationsUnadjustedTemporalAfterRaceWhite",FALSE)
runPheWAS(newPhenotypesAfterHL,attributesRaceWhite,"AssociationsAdjustedTemporalAfterRaceWhite",TRUE)
runPheWAS(newPhenotypesAfterHL,attributesRaceBlack,"AssociationsUnadjustedTemporalAfterRaceBlack",FALSE)
runPheWAS(newPhenotypesAfterHL,attributesRaceBlack,"AssociationsAdjustedTemporalAfterRaceBlack",TRUE)
#Stratified by Sex
runPheWAS(newPhenotypesAfterHL,attributesSexMale,"AssociationsUnadjustedTemporalAfterSexMale",FALSE)
runPheWAS(newPhenotypesAfterHL,attributesSexMale,"AssociationsAdjustedTemporalAfterSexMale",TRUE)
runPheWAS(newPhenotypesAfterHL,attributesSexFemale,"AssociationsUnadjustedTemporalAfterSexFemale",TRUE)
runPheWAS(newPhenotypesAfterHL,attributesSexFemale,"AssociationsAdjustedTemporalAfterSexFemale",TRUE)
#Stratifed by both Race and Sex
runPheWAS(newPhenotypesAfterHL,attributesRaceWhiteSexMale,"AssociationsUnadjustedTemporalAfterRaceWhiteSexMale",FALSE)
runPheWAS(newPhenotypesAfterHL,attributesRaceWhiteSexMale,"AssociationsAdjustedTemporalAfterRaceWhiteSexMale",TRUE)
runPheWAS(newPhenotypesAfterHL,attributesRaceWhiteSexFemale,"AssociationsUnadjustedTemporalAfterRaceWhiteSexFemale",FALSE)
runPheWAS(newPhenotypesAfterHL,attributesRaceWhiteSexFemale,"AssociationsAdjustedTemporalAfterRaceWhiteSexFemale",TRUE)
runPheWAS(newPhenotypesAfterHL,attributesRaceBlackSexMale,"AssociationsUnadjustedTemporalAfterRaceBlackSexMale",FALSE)
runPheWAS(newPhenotypesAfterHL,attributesRaceBlackSexMale,"AssociationsAdjustedTemporalAfterRaceBlackSexMale",TRUE)
runPheWAS(newPhenotypesAfterHL,attributesRaceBlackSexFemale,"AssociationsUnadjustedTemporalAfterRaceBlackSexFemale",FALSE)
runPheWAS(newPhenotypesAfterHL,attributesRaceBlackSexFemale,"AssociationsAdjustedTemporalAfterRaceBlackSexFemale",TRUE)

#Fisher's Exact Test

#All
pheomeWideFisherTest(newPhenotypesAfterHL,attributes,"AssociationsFisherTemporalAfter")

#Race
pheomeWideFisherTest(newPhenotypesAfterHL,attributesRaceBlack,"AssociationsFisherTemporalAfterRaceBlack")
pheomeWideFisherTest(newPhenotypesAfterHL,attributesRaceWhite,"AssociationsFisherTemporalAfterRaceWhite")

#Sex
pheomeWideFisherTest(newPhenotypesAfterHL,attributesSexMale,"AssociationsFisherTemporalAfterSexMale")
pheomeWideFisherTest(newPhenotypesAfterHL,attributesSexFemale,"AssociationsFisherTemporalAfterSexFemale")

#Race and Sex
pheomeWideFisherTest(newPhenotypesAfterHL,attributesRaceWhiteSexMale,"AssociationsFisherTemporalAfterRaceWhiteSexMale")
pheomeWideFisherTest(newPhenotypesAfterHL,attributesRaceWhiteSexFemale,"AssociationsFisherTemporalAfterRaceWhiteSexFemale")
pheomeWideFisherTest(newPhenotypesAfterHL,attributesRaceBlackSexMale,"AssociationsFisherTemporalAfterRaceBlackSexMale")
pheomeWideFisherTest(newPhenotypesAfterHL,attributesRaceBlackSexFemale,"AssociationsFisherTemporalAfterRaceBlackSexFemale")

#General correlation test between Fisher's Exact Test odds ratios and PheWAS odds ratios
FisherversusLogisticCor("AssociationsFisher","AssociationsUnadjusted")
FisherversusLogisticCor("AssociationsFisher","AssociationsAdjusted")

#Specific correlation test between Fisher's Exact Test odds ratios and PheWAS odds ratios
FisherversusLogisticCor("AssociationsFisherTemporalBeforeRaceBlackSexMale","AssociationsUnadjustedTemporalBeforeRaceBlackSexMale")
FisherversusLogisticCor("AssociationsFisherTemporalBeforeRaceBlackSexFemale","AssociationsAdjustedTemporalBeforeRaceBlackSexFemale")