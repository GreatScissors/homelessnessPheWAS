library(DescTools)

#This performs tests for homogeneity for HIV status


createContigencyTableCounts <- function(filename,phecode) {
  data <- read.csv(paste0(filename,".csv"))
  
  contigencyTable <- matrix(nrow = 2,ncol = 2)
  rownames(contigencyTable) <- c("TRUE","FALSE")
  colnames(contigencyTable) <- c("HL","NH")
  
  contigencyTable[1,1] <- data[1,paste0("X",phecode)]
  contigencyTable[2,1] <- data[2,paste0("X",phecode)]
  contigencyTable[1,2] <- data[4,paste0("X",phecode)]
  contigencyTable[2,2] <- data[5,paste0("X",phecode)]
  
  return(t(contigencyTable))
}

#Race

x <- createContigencyTableCounts("AssociationsFisherTemporalBeforeRaceBlackSexMale","071")
y <- createContigencyTableCounts("AssociationsFisherTemporalBeforeRaceWhiteSexMale","071")

data <- array(c(x,y),dim = c(2,2,2),dimnames=list(c("HL","NH"),c("TRUE","FALSE"),c("Black","White")))

mantelhaen.test(data, exact = TRUE)

BreslowDayTest(data, OR = NA, correct = FALSE)

#Sex

x <- createContigencyTableCounts("AssociationsFisherTemporalBeforeRaceBlackSexMale","071")
y <- createContigencyTableCounts("AssociationsFisherTemporalBeforeRaceBlackSexFemale","071")

data <- array(c(x,y),dim = c(2,2,2),dimnames=list(c("HL","NH"),c("TRUE","FALSE"),c("Male","Female")))

mantelhaen.test(data, exact = TRUE)

BreslowDayTest(data, OR = NA, correct = FALSE)



CTHomeless <- function(filename1, filename2, phecode) {
  data1 <- read.csv(paste0(filename1,".csv"))
  data2 <- read.csv(paste0(filename2,".csv"))
  
  contigencyTable <- matrix(nrow = 2,ncol = 2)
  rownames(contigencyTable) <- c("Black","White")
  colnames(contigencyTable) <- c("TRUE","FALSE")
  
  contigencyTable[1,1] <- data1[1,paste0("X",phecode)]
  contigencyTable[1,2] <- data1[2,paste0("X",phecode)]
  contigencyTable[2,1] <- data2[1,paste0("X",phecode)]
  contigencyTable[2,2] <- data2[2,paste0("X",phecode)]
  
  return(t(contigencyTable))
}

CTNeverHomeless <- function(filename1, filename2, phecode) {
  data1 <- read.csv(paste0(filename1,".csv"))
  data2 <- read.csv(paste0(filename2,".csv"))
  
  contigencyTable <- matrix(nrow = 2,ncol = 2)
  rownames(contigencyTable) <- c("Black","White")
  colnames(contigencyTable) <- c("TRUE","FALSE")
  
  contigencyTable[1,1] <- data1[4,paste0("X",phecode)]
  contigencyTable[1,2] <- data1[5,paste0("X",phecode)]
  contigencyTable[2,1] <- data2[4,paste0("X",phecode)]
  contigencyTable[2,2] <- data2[5,paste0("X",phecode)]
  
  return(t(contigencyTable))
}

#HL Status

x <- CTHomeless("AssociationsFisherTemporalBeforeRaceBlackSexMale","AssociationsFisherTemporalBeforeRaceWhiteSexMale","071")
y <- CTNeverHomeless("AssociationsFisherTemporalBeforeRaceBlackSexMale","AssociationsFisherTemporalBeforeRaceWhiteSexMale","071")

data <- array(c(x,y),dim = c(2,2,2),dimnames=list(c("True","False"),c("Black","White"),c("HL","NH")))

mantelhaen.test(data, exact = TRUE)

BreslowDayTest(data, OR = NA, correct = FALSE)
