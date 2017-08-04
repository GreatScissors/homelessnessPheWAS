stratifyTemporal <- function(temporalSelection, currentPhenotypes, firstInstancePhecodeTemporal) {
  #Parameter: phecode name and PID selected uniquePhecodeMappedSelected data frame
  processPhecode <- function(phecodeName,PIDOnlyDF) {
    #Select inputed phecode only from PIDOnlyDF
    phecodeRow <- PIDOnlyDF[which(PIDOnlyDF$phecode == phecodeName),]
    
    allTemporal <- c("B","T","A")
    #Non-selected temporal value will remain in list
    temporalExclude <- setdiff(allTemporal,temporalSelection)
    
    #Check if phecode exists
    if(nrow(phecodeRow) > 0) {
      if(phecodeRow$temporal %in% temporalExclude) {
        #Omit phecode if temporal value is not the selected temporal value
        return(NA)
      } else {
        #True or 1: false*true = false, true*true=true, NA*true=NA (unchanged with multiplication)
        return(TRUE)
      }
    } else {
      return(TRUE)
    }
  }
  
  #Parameter: PheWAS table row
  excludePhecode <- function(row) {
    #First element is the PID
    PID <- as.numeric(row[1])
    #Print row number to monitor progress
    print(which(phenotypes$id == PID))
    
    #Retrieve all rows belonging to a PID number from uniquePhecodeMappedSelected data frame (passed as argument to stratifyTemporal())
    PIDOnlyDF <- firstInstancePhecodeTemporal[which(firstInstancePhecodeTemporal$PID == PID),]
    #Remove element containing PID
    phecodeValues <- row[-1]
    #Check if PIDOnly has rows.  If yes, then the PID is a case with at least one phecode.
    if(nrow(PIDOnlyDF) > 0){
      #Process each phecode name for temporal value
      phecodeLogical <- sapply(names(phecodeValues),processPhecode,PIDOnlyDF=PIDOnlyDF)
      return(as.logical(phecodeLogical*row[2:length(row)]))
    } else {
      #PID has no phecodes or is not a case.
      return(phecodeValues)
    }
  }
  
  #format data in same format as createPhewasTable()
  tester <- rbind(t(apply(currentPhenotypes,1,excludePhecode)))
  newPhenotypes <- as.data.frame(apply(tester,2,as.logical))
  newPhenotypes <- cbind(currentPhenotypes$id,newPhenotypes)
  colnames(newPhenotypes) <- colnames(phenotypes)
  
  return(newPhenotypes)
}