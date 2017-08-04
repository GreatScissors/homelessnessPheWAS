pheomeWideFisherTest <- function(phenotypesInput,attributesInput, filename) {
  HLAttributes <- attributesInput[which(attributesInput$HL_Status==1),]
  NHAttributes <- attributesInput[which(attributesInput$HL_Status==0),]
  
  HLPhenotypes <- phenotypesInput[match(HLAttributes$id,phenotypesInput$id),]
  NHPhenotypes <- phenotypesInput[match(NHAttributes$id,phenotypesInput$id),]
  
  trueCountCases <- colSums(HLPhenotypes,na.rm = TRUE)
  NACountCases <- colSums(is.na(HLPhenotypes))
  totalCountCases <- replicate(ncol(HLPhenotypes),nrow(HLPhenotypes))
  falseCountCases <- totalCountCases - trueCountCases - NACountCases
  
  trueCountControls <- colSums(NHPhenotypes,na.rm = TRUE)
  NACountControls <- colSums(is.na(NHPhenotypes))
  totalCountControls <- replicate(ncol(NHPhenotypes),nrow(NHPhenotypes))
  falseCountControls <- totalCountControls - trueCountControls - NACountControls
  
  CasesControls <- rbind(trueCountCases,falseCountCases,trueCountControls,falseCountControls)
  
  runFisher <- function(phenotype) {
    
    contigencyTable <- matrix(phenotype,nrow = 2,ncol = 2)
    rownames(contigencyTable) <- c("TRUE","FALSE")
    colnames(contigencyTable) <- c("Homeless","Never Homeless")
     
    return(c(contigencyTable["TRUE","Homeless"],contigencyTable["FALSE","Homeless"],100*(contigencyTable["TRUE","Homeless"]/(contigencyTable["TRUE","Homeless"]+contigencyTable["FALSE","Homeless"])),contigencyTable["TRUE","Never Homeless"],contigencyTable["FALSE","Never Homeless"],100*(contigencyTable["TRUE","Never Homeless"]/(contigencyTable["TRUE","Never Homeless"]+contigencyTable["FALSE","Never Homeless"])),fisher.test(contigencyTable)$estimate,fisher.test(contigencyTable)$p.value,fisher.test(contigencyTable)$conf.int))
  }
  
  results <- apply(CasesControls[,2:ncol(CasesControls)],2,runFisher)
  
  row.names(results) <- c("CasesHL","ControlsHL","%CasesHL","CasesNH","ControlsNH","%CasesNH","Odds.Ratio","P.Value","LowerConfidenceInterval","UpperConfidenceInterval")
  
  write.csv(results,paste0(filename,".csv"))
}
