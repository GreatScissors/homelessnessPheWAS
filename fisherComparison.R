fisherSelection <- function(filename,bonferroniCorrection=FALSE,FDR=FALSE,minOR=1,pvalueThreshold=0.05,includeInf=FALSE,minCases=20) {
  fisherData <- data.frame(t(read.csv(paste0(filename,".csv"),row.names = 1)))
  
  fisherData <- fisherData[which(fisherData$CasesHL >=minCases),]
  
  #https://stackoverflow.com/questions/9704213/r-remove-part-of-string
  phecodes <- sapply(strsplit(rownames(fisherData), split='X', fixed=TRUE), function(x) (x[2]))
  
  rownames(fisherData) <- phecodes
  
  #FDR overides Bonferroni if both true
  passedFisherData <- NULL
  if(FDR == TRUE) {
    p.fdr <-p.adjust(p=fisherData$P.Value, method="fdr")
    passedFisherData <- fisherData[which(p.fdr < pvalueThreshold),]
  } else {
    if(bonferroniCorrection == TRUE) {
      p.bon <-p.adjust(p=fisherData$P.Value, method="bonferroni")
      passedFisherData <- fisherData[which(p.bon < pvalueThreshold),]
    } else {
      passedFisherData <- fisherData[which(fisherData$P.Value < pvalueThreshold),]
    }
  }

  #Use to map phecode to description
  AssociationsAdjustedTemporalBeforeSexFemale <- read.csv("AssociationsAdjustedTemporalBeforeSexFemale.csv",header = TRUE)
  
  matchedIndices <- match(as.numeric(row.names(passedFisherData)),AssociationsAdjustedTemporalBeforeSexFemale$phecode)
  passedFisherData <- cbind(AssociationsAdjustedTemporalBeforeSexFemale$description[matchedIndices],passedFisherData)
  
  selected <- passedFisherData[which(passedFisherData$Odds.Ratio>minOR & (is.finite(passedFisherData$Odds.Ratio) | includeInf) ),]
  
  colnames(selected)[1] <- "description"
  
  return(selected)
}

#White Male Temporal Comparison
BeforeRaceWhiteSexMale <- fisherSelection("AssociationsFisherTemporalBeforeRaceWhiteSexMale",FDR = TRUE)
AfterRaceWhiteSexMale <- fisherSelection("AssociationsFisherTemporalAfterRaceWhiteSexMale",FDR = TRUE)

x <- AfterRaceWhiteSexMale[!(AfterRaceWhiteSexMale$description %in% BeforeRaceWhiteSexMale$description),]
BeforeRaceWhiteSexMale$description[!(BeforeRaceWhiteSexMale$description %in% AfterRaceWhiteSexMale$description)]

#White Female Temporal Comparison
BeforeRaceWhiteSexFemale <- fisherSelection("AssociationsFisherTemporalBeforeRaceWhiteSexFemale",FDR = TRUE)
AfterRaceWhiteSexFemale <- fisherSelection("AssociationsFisherTemporalAfterRaceWhiteSexFemale",FDR = TRUE)

AfterRaceWhiteSexFemale$description[!(AfterRaceWhiteSexFemale$description %in% BeforeRaceWhiteSexFemale$description)]

#Black Male Temporal Comparison
BeforeRaceBlackSexMale <- fisherSelection("AssociationsFisherTemporalBeforeRaceBlackSexMale",FDR = TRUE)
AfterRaceBlackSexMale <- fisherSelection("AssociationsFisherTemporalAfterRaceBlackSexMale",FDR = TRUE)

AfterRaceBlackSexMale$description[!(AfterRaceBlackSexMale$description %in% BeforeRaceBlackSexMale$description)]

#Black Female Temporal Comparison
BeforeRaceBlackSexfemale <- fisherSelection("AssociationsFisherTemporalBeforeRaceBlackSexfemale",FDR = TRUE)
AfterRaceBlackSexfemale <- fisherSelection("AssociationsFisherTemporalAfterRaceBlackSexfemale",FDR = TRUE)

AfterRaceBlackSexfemale$description[!(AfterRaceBlackSexfemale$description %in% BeforeRaceBlackSexfemale$description)]