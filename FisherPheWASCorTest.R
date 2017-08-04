FisherversusLogisticCor <- function(fisherDataFile, logisticDataFile) {
  
  FisherData <- read.csv(paste0(fisherDataFile,".csv"))
  LogisticData <- read.csv(paste0(logisticDataFile,".csv"))
  
  #remove first column names
  truncateData <- FisherData[2:ncol(FisherData)]
  #acquire phecode names and assign them to truncateData
  colnames(truncateData) <- names(phenotypes[2:ncol(phenotypes)])
  
  usableDataIndexes <- which(!is.na(LogisticData$OR))
  #include only phecodes without NA values from the PheWAS
  phecodes <- LogisticData$phecode[usableDataIndexes]
  FisherIndexes <- match(phecodes,colnames(truncateData))
  
  #Exlude phecodes that don't match
  ORLogistic <- LogisticData[usableDataIndexes[which(!is.na(FisherIndexes))],"OR"]
  ORFisher <- truncateData[7,FisherIndexes[which(!is.na(FisherIndexes))]]
  
  #include only finite odds ratios
  finite <- which(is.finite(as.numeric(ORFisher)))
  logistic <- ORLogistic[finite]
  fisher <- as.numeric(ORFisher[finite])
  
  #print phecodes
  print("Phecodes Utilized")
  print(names(ORFisher))
  print(cor.test(logistic,fisher))
  plot(logistic,fisher)
}