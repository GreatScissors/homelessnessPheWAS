#vignette("PheWAS-package")

csv.phenotypes <- read.csv("id.icd9.count.csv",header = TRUE)
#rename to be compliant with example
colnames(csv.phenotypes) <- c("id","icd9","count")

attributes <- read.csv("attrib.csv",header = TRUE)
colnames(attributes)[1] <- "id"

set.seed(1)

phenotypes = createPhewasTable(csv.phenotypes)

runPheWAS <- function(phenotypesInput,attributesInput,filename,adjusted) {
  
  results <- NULL
  
  #Don't adjust for SEX if it's not contained in the attributes file
  if(adjusted == TRUE) {
    if('RACE' %in% colnames(attributesInput)) {
      if('SEX' %in% colnames(attributesInput)) {
        results = phewas(phenotypes = phenotypesInput,genotypes = attributesInput[,c("id","HL_Status")],covariates = attributesInput[,c("id",c("AGE","SEX","RACE"))],cores=4,significance.threshold = c("bonferroni"))
        } else {
        results = phewas(phenotypes = phenotypesInput,genotypes = attributesInput[,c("id","HL_Status")],covariates = attributesInput[,c("id",c("AGE","RACE"))],cores=4,significance.threshold = c("bonferroni"))
      }
    } else {
      if('SEX' %in% colnames(attributesInput)) {
        results = phewas(phenotypes = phenotypesInput,genotypes = attributesInput[,c("id","HL_Status")],covariates = attributesInput[,c("id",c("AGE","SEX"))],cores=4,significance.threshold = c("bonferroni"))
      } else {
        results = phewas(phenotypes = phenotypesInput,genotypes = attributesInput[,c("id","HL_Status")],covariates = attributesInput[,c("id",c("AGE"))],cores=4,significance.threshold = c("bonferroni"))
      }
    }
  } else {
    results = phewas(phenotypes = phenotypesInput,genotypes = attributesInput[,c("id","HL_Status")],cores=4,significance.threshold = c("bonferroni"))
  }
  
  results_d = addPhecodeInfo(results,descriptions = T,groups = T,groupnums = F,groupcolors = F)
  write.csv(results_d[order(-results_d$OR),],paste0(filename,".csv"))
}