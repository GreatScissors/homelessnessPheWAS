#rm(list = ls())

numberofCases <- 40
oddsRatioUpperThreshold <- 500
oddsRatioLowerThreshold <- 0

oddsRatioComparison <- function(df1Raw,df2Raw,minCases=20) {
  df1 <- df1Raw[which(df1Raw$n_cases >= minCases & df1Raw$bonferroni==TRUE),]
  df2 <- df2Raw[which(df2Raw$n_cases >= minCases),]
  
  intersects <- intersect(df1$description,df2$description)
  
  df1Indices <- match(intersects,df1$description)
  df2Indices <- match(intersects,df2$description)
  
  ORsComparison <- df1[df1Indices,"OR"] / df2[df2Indices,"OR"]
  names(ORsComparison) <- intersects
  
  x <- abs(log(df1[df1Indices,"OR"]) - log(df2[df2Indices,"OR"]))
  SE <- sqrt(df1[df1Indices,"SE"]^2+df2[df2Indices,"SE"]^2)
  z <- x/SE
  Pvalue <- 2*(1-pnorm(z))
  
  ORsComparison <- rbind(df1[df1Indices,"n_cases"],df1[df1Indices,"p"],df1[df1Indices,"OR"],df2[df2Indices,"n_cases"],df2[df2Indices,"p"],df2[df2Indices,"OR"],ORsComparison,Pvalue)
  
  row.names(ORsComparison) <- c("n_cases_first","p_value_first","OR_first","n_cases_second","p_value_second","OR_second","OR_multiplier","Z_test_pvalue")
  
  return(ORsComparison)
}

#Load Sex Temporal Before
AssociationsAdjustedTemporalBeforeSexFemale <- read.csv("AssociationsAdjustedTemporalBeforeSexFemale.csv",header = TRUE)
AssociationsAdjustedTemporalBeforeSexMale <- read.csv("AssociationsAdjustedTemporalBeforeSexMale.csv",header = TRUE)

#Load Sex Temporal After
AssociationsAdjustedTemporalAfterSexFemale <- read.csv("AssociationsAdjustedTemporalAfterSexFemale.csv",header = TRUE)
AssociationsAdjustedTemporalAfterSexMale <- read.csv("AssociationsAdjustedTemporalAfterSexMale.csv",header = TRUE)

#Load Race Temporal Before
AssociationsAdjustedTemporalBeforeRaceWhite <- read.csv("AssociationsAdjustedTemporalBeforeRaceWhite.csv",header = TRUE)
AssociationsAdjustedTemporalBeforeRaceBlack <- read.csv("AssociationsAdjustedTemporalBeforeRaceBlack.csv",header = TRUE)

#Load Race Temporal After
AssociationsAdjustedTemporalAfterRaceWhite <- read.csv("AssociationsAdjustedTemporalAfterRaceWhite.csv",header = TRUE)
AssociationsAdjustedTemporalAfterRaceBlack <- read.csv("AssociationsAdjustedTemporalAfterRaceBlack.csv",header = TRUE)

#Load Temporal Before Data
AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale <- read.csv("AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale.csv",header = TRUE)
AssociationsAdjustedTemporalBeforeRaceWhiteSexMale <- read.csv("AssociationsAdjustedTemporalBeforeRaceWhiteSexMale.csv",header = TRUE)
AssociationsAdjustedTemporalBeforeRaceBlackSexFemale <- read.csv("AssociationsAdjustedTemporalBeforeRaceBlackSexFemale.csv",header = TRUE)
AssociationsAdjustedTemporalBeforeRaceBlackSexMale <- read.csv("AssociationsAdjustedTemporalBeforeRaceBlackSexMale.csv",header = TRUE)

#Load Temporal After 
AssociationsAdjustedTemporalAfterRaceWhiteSexFemale <- read.csv("AssociationsAdjustedTemporalAfterRaceWhiteSexFemale.csv",header = TRUE)
AssociationsAdjustedTemporalAfterRaceWhiteSexMale <- read.csv("AssociationsAdjustedTemporalAfterRaceWhiteSexMale.csv",header = TRUE)
AssociationsAdjustedTemporalAfterRaceBlackSexFemale <- read.csv("AssociationsAdjustedTemporalAfterRaceBlackSexFemale.csv",header = TRUE)
AssociationsAdjustedTemporalAfterRaceBlackSexMale <- read.csv("AssociationsAdjustedTemporalAfterRaceBlackSexMale.csv",header = TRUE)

#Temporal Before Race, Sex Combinations
data1 <- AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale[which(AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale$n_cases >= numberofCases & AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale$bonferroni==TRUE & AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale$OR >= oddsRatioLowerThreshold),]
data2 <- AssociationsAdjustedTemporalBeforeRaceWhiteSexMale[which(AssociationsAdjustedTemporalBeforeRaceWhiteSexMale$n_cases >= numberofCases & AssociationsAdjustedTemporalBeforeRaceWhiteSexMale$bonferroni==TRUE & AssociationsAdjustedTemporalBeforeRaceWhiteSexMale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalBeforeRaceWhiteSexMale$OR >= oddsRatioLowerThreshold),]
data3 <- AssociationsAdjustedTemporalBeforeRaceBlackSexFemale[which(AssociationsAdjustedTemporalBeforeRaceBlackSexFemale$n_cases >= numberofCases & AssociationsAdjustedTemporalBeforeRaceBlackSexFemale$bonferroni==TRUE & AssociationsAdjustedTemporalBeforeRaceBlackSexFemale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalBeforeRaceBlackSexFemale$OR >= oddsRatioLowerThreshold),]
data4 <- AssociationsAdjustedTemporalBeforeRaceBlackSexMale[which(AssociationsAdjustedTemporalBeforeRaceBlackSexMale$n_cases >= numberofCases & AssociationsAdjustedTemporalBeforeRaceBlackSexMale$bonferroni==TRUE & AssociationsAdjustedTemporalBeforeRaceBlackSexMale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalBeforeRaceBlackSexMale$OR >= oddsRatioLowerThreshold),]

#Black Male to White Male
BlackMaletoWhiteMaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalBeforeRaceBlackSexMale,AssociationsAdjustedTemporalBeforeRaceWhiteSexMale,minCases = 50)
#White Male to Black Male
WhiteMaletoBlackMaleORSComparison <- oddsRatioComparison(AssociationsAdjustedTemporalBeforeRaceWhiteSexMale,AssociationsAdjustedTemporalBeforeRaceBlackSexMale,minCases = 50)

#Black Female to White Female
BlackFemaletoWhiteFemaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalBeforeRaceBlackSexFemale,AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale,50)
#White Female to Black Female
WhiteFemaletoBlackFemaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale,AssociationsAdjustedTemporalBeforeRaceBlackSexFemale,50)

#Black Female to Black Male 
BlackFemaletoBlackMaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalBeforeRaceBlackSexFemale,AssociationsAdjustedTemporalBeforeRaceBlackSexMale,50)
#Black Male to Black Female
BlackMaletoBlackFemaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalBeforeRaceBlackSexMale,AssociationsAdjustedTemporalBeforeRaceBlackSexFemale,50)

#White Male to White Female
WhiteMaletoWhiteFemaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalBeforeRaceWhiteSexMale,AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale,minCases = 50)
#White Female to White Male
WhiteFemaletoWhiteMaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalBeforeRaceWhiteSexFemale,AssociationsAdjustedTemporalBeforeRaceWhiteSexMale,0)

Reduce(intersect,list(data1$description,data2$description,data3$description,data4$description))

#Temporal After Race, Sex Combinations
data1 <- AssociationsAdjustedTemporalAfterRaceWhiteSexFemale[which(AssociationsAdjustedTemporalAfterRaceWhiteSexFemale$n_cases >= numberofCases & AssociationsAdjustedTemporalAfterRaceWhiteSexFemale$bonferroni==TRUE & AssociationsAdjustedTemporalAfterRaceWhiteSexFemale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalAfterRaceWhiteSexFemale$OR >= oddsRatioLowerThreshold),]
data2 <- AssociationsAdjustedTemporalAfterRaceWhiteSexMale[which(AssociationsAdjustedTemporalAfterRaceWhiteSexMale$n_cases >= numberofCases & AssociationsAdjustedTemporalAfterRaceWhiteSexMale$bonferroni==TRUE & AssociationsAdjustedTemporalAfterRaceWhiteSexMale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalAfterRaceWhiteSexMale$OR >= oddsRatioLowerThreshold),]
data3 <- AssociationsAdjustedTemporalAfterRaceBlackSexFemale[which(AssociationsAdjustedTemporalAfterRaceBlackSexFemale$n_cases >= numberofCases & AssociationsAdjustedTemporalAfterRaceBlackSexFemale$bonferroni==TRUE & AssociationsAdjustedTemporalAfterRaceBlackSexFemale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalAfterRaceBlackSexFemale$OR >= oddsRatioLowerThreshold),]
data4 <- AssociationsAdjustedTemporalAfterRaceBlackSexMale[which(AssociationsAdjustedTemporalAfterRaceBlackSexMale$n_cases >= numberofCases & AssociationsAdjustedTemporalAfterRaceBlackSexMale$bonferroni==TRUE & AssociationsAdjustedTemporalAfterRaceBlackSexMale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalAfterRaceBlackSexMale$OR >= oddsRatioLowerThreshold),]

data1After <- data1
data2After <- data2
data3After <- data3
data4After <- data4

#Black Male to White Male
BlackMaletoWhiteMaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalAfterRaceBlackSexMale,AssociationsAdjustedTemporalAfterRaceWhiteSexMale,minCases = 50)
#White Male to Black Male
WhiteMaletoBlackMaleORSComparison <- oddsRatioComparison(AssociationsAdjustedTemporalAfterRaceWhiteSexMale,AssociationsAdjustedTemporalAfterRaceBlackSexMale,minCases = 0)

#Black Female to White Female
BlackFemaletoWhiteFemaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalAfterRaceBlackSexFemale,AssociationsAdjustedTemporalAfterRaceWhiteSexFemale,50)
#White Female to Black Female
WhiteFemaletoBlackFemaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalAfterRaceWhiteSexFemale,AssociationsAdjustedTemporalAfterRaceBlackSexFemale,50)

#Black Female to Black Male 
BlackFemaletoBlackMaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalAfterRaceBlackSexFemale,AssociationsAdjustedTemporalAfterRaceBlackSexMale,50)
#Black Male to Black Female
BlackMaletoBlackFemaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalAfterRaceBlackSexMale,AssociationsAdjustedTemporalAfterRaceBlackSexFemale,50)

#White Male to White Female
WhiteMaletoWhiteFemaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalAfterRaceWhiteSexMale,AssociationsAdjustedTemporalAfterRaceWhiteSexFemale,minCases = 50)
#White Female to White Male
WhiteFemaletoWhiteMaleORsComparison <- oddsRatioComparison(AssociationsAdjustedTemporalAfterRaceWhiteSexFemale,AssociationsAdjustedTemporalAfterRaceWhiteSexMale,50)

Reduce(intersect,list(data1$description,data2$description,data3$description,data4$description))

#Temporal Before Sex
data1 <- AssociationsAdjustedTemporalBeforeSexFemale[which(AssociationsAdjustedTemporalBeforeSexFemale$n_cases >= numberofCases & AssociationsAdjustedTemporalBeforeSexFemale$bonferroni==TRUE & AssociationsAdjustedTemporalBeforeSexFemale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalBeforeSexFemale$OR >= oddsRatioLowerThreshold),]
data2 <- AssociationsAdjustedTemporalBeforeSexMale[which(AssociationsAdjustedTemporalBeforeSexMale$n_cases >= numberofCases & AssociationsAdjustedTemporalBeforeSexMale$bonferroni==TRUE & AssociationsAdjustedTemporalBeforeSexMale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalBeforeSexMale$OR >= oddsRatioLowerThreshold),]

Reduce(intersect,list(data1$description,data2$description))

sexBefore <- oddsRatioComparison(AssociationsAdjustedTemporalBeforeSexMale,AssociationsAdjustedTemporalBeforeSexFemale,50)

#Temporal After Sex
data1 <- AssociationsAdjustedTemporalAfterSexFemale[which(AssociationsAdjustedTemporalAfterSexFemale$n_cases >= numberofCases & AssociationsAdjustedTemporalAfterSexFemale$bonferroni==TRUE & AssociationsAdjustedTemporalAfterSexFemale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalAfterSexFemale$OR >= oddsRatioLowerThreshold),]
data2 <- AssociationsAdjustedTemporalAfterSexMale[which(AssociationsAdjustedTemporalAfterSexMale$n_cases >= numberofCases & AssociationsAdjustedTemporalAfterSexMale$bonferroni==TRUE & AssociationsAdjustedTemporalAfterSexMale$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalAfterSexMale$OR >= oddsRatioLowerThreshold),]

Reduce(intersect,list(data1$description,data2$description))

sexAfter <- oddsRatioComparison(AssociationsAdjustedTemporalAfterSexMale,AssociationsAdjustedTemporalAfterSexFemale,50)

#Temporal Before Race
data1 <- AssociationsAdjustedTemporalBeforeRaceWhite[which(AssociationsAdjustedTemporalBeforeRaceWhite$n_cases >= numberofCases & AssociationsAdjustedTemporalBeforeRaceWhite$bonferroni==TRUE & AssociationsAdjustedTemporalBeforeRaceWhite$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalBeforeRaceWhite$OR >= oddsRatioLowerThreshold),]
data2 <- AssociationsAdjustedTemporalBeforeRaceBlack[which(AssociationsAdjustedTemporalBeforeRaceBlack$n_cases >= numberofCases & AssociationsAdjustedTemporalBeforeRaceBlack$bonferroni==TRUE & AssociationsAdjustedTemporalBeforeRaceBlack$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalBeforeRaceBlack$OR >= oddsRatioLowerThreshold),]

Reduce(intersect,list(data1$description,data2$description))

raceBefore <- oddsRatioComparison(AssociationsAdjustedTemporalBeforeRaceWhite,AssociationsAdjustedTemporalBeforeRaceBlack,20)

#Temporal After Race
data1 <- AssociationsAdjustedTemporalAfterRaceWhite[which(AssociationsAdjustedTemporalAfterRaceWhite$n_cases >= numberofCases & AssociationsAdjustedTemporalAfterRaceWhite$bonferroni==TRUE & AssociationsAdjustedTemporalAfterRaceWhite$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalAfterRaceWhite$OR >= oddsRatioLowerThreshold),]
data2 <- AssociationsAdjustedTemporalAfterRaceBlack[which(AssociationsAdjustedTemporalAfterRaceBlack$n_cases >= numberofCases & AssociationsAdjustedTemporalAfterRaceBlack$bonferroni==TRUE & AssociationsAdjustedTemporalAfterRaceBlack$OR <= oddsRatioUpperThreshold & AssociationsAdjustedTemporalAfterRaceBlack$OR >= oddsRatioLowerThreshold),]

Reduce(intersect,list(data1$description,data2$description))

raceAfter <- oddsRatioComparison(AssociationsAdjustedTemporalAfterRaceBlack,AssociationsAdjustedTemporalAfterRaceWhite,50)

greaterORs <- function(df) {
  test <- data.frame(t(df))
  return(test[which(test$OR_multiplier>1),])
}

#White Male After to Before
whiteMaleAftertoBefore <- oddsRatioComparison(AssociationsAdjustedTemporalAfterRaceWhiteSexMale,AssociationsAdjustedTemporalBeforeRaceWhiteSexMale)
greaterORs(WhiteMaleAftertoBefore)