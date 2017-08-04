calculateCI <- function(filename, phecode, alpha=0.05){
  data <- read.csv(paste0(filename,".csv"))
  entry <- data[which(data$phecode==phecode),]
  
  print(entry$OR)
  return(exp(log(entry$OR)+c(-1,1)*qnorm(1-alpha/2)*entry$SE))
}


calculateCI("AssociationsAdjusted","297.1")
calculateCI("AssociationsAdjusted","317.1")
calculateCI("AssociationsAdjusted","301")
calculateCI("AssociationsAdjusted","295.1")
calculateCI("AssociationsAdjusted","316")

calculateCI("AssociationsAdjustedTemporalBeforeRaceBlackSexMale","71")
calculateCI("AssociationsAdjustedTemporalBeforeRaceWhiteSexMale","71")
