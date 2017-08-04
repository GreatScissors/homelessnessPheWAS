countTemporal <- function(phenotypedf) {
  test <- replace(phenotypedf,is.na(phenotypedf),FALSE)
  
  truncTest <- test[,2:ncol(test)]
  
  convertedTruncTest <- apply(truncTest,2, function(x) as.numeric(x))
  
  return(sum(apply(convertedTruncTest,2,sum)))
}

countTemporal(phenotypes)
countTemporal(newPhenotypesBeforeHL)
countTemporal(newPhenotypesAfterHL)

length(which(uniquePhecodeMappedSelected$temporal == 'B'))
length(which(uniquePhecodeMappedSelected$temporal == 'A'))