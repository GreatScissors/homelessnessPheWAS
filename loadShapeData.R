#Load raw case and control data
cases <- read.csv(casesToRead)
controls <- read.csv(controlsToRead)

#Create attributeData df
casesFilt <- cases[,c(1,3,5,6)]
HL_Status <- 1:nrow(casesFilt)
#Integer 1 for homeless
HL_Status <- HL_Status/HL_Status
casesFilt <- cbind(casesFilt,HL_Status)

controlsFilt <- controls[,c(1,4:6)]
HL_Status <- 1:nrow(controlsFilt)
#Integer 0 for not homeless
HL_Status <- HL_Status*0
controlsFilt <- cbind(controlsFilt,HL_Status)

attributeData <- rbind(casesFilt,controlsFilt)

#Raw ICD9 codes
#First variable in each line is PID, followed by date ICD9 date ICD9 date ICD9 ...etc
lines <- readLines(ICD9CodesRaw)

ICD9Codes <- NULL

processLine <- function(line) {
  #Extract vector of strings of each variable in line
  elements <- unlist(str_split(line," "))
  
  #Extract and remove PID
  PID <- elements[1]
  elements <- elements[-1]
  
  #Create matrix with date and ICD9 column, filling by row
  matrixElements <- matrix(elements,ncol=2,nrow = length(elements)/2,byrow = TRUE)
  
  df <- as.data.frame(matrixElements,stringsAsFactors = FALSE)
  
  colnames(df) <- c("Date","ICD9")
  
  #count ICD9 numbers
  patientEntry <- count(df,"ICD9")
  
  #Repeat the PID for the number of distinct ICD9 numbers
  PID <- rep(PID,nrow(patientEntry))
  colnames(patientEntry) <- c("Code","Code_Count")
  patientDataFrame <- cbind(PID,patientEntry)
  
  return(patientDataFrame)
}

# returns PID Date ICD9 as columns for a data frame
processLineICD9Date <- function(line) {
  #Extract vector of strings of each variable in line
  elements <- unlist(str_split(line," "))
  
  #Extract and remove PID
  PID <- elements[1]
  elements <- elements[-1]
  
  #Create matrix with date and ICD9 column, filling by row
  matrixElements <- matrix(elements,ncol=2,nrow = length(elements)/2,byrow = TRUE)
  
  df <- as.data.frame(matrixElements,stringsAsFactors = FALSE)
  
  colnames(df) <- c("Date","ICD9")
  
  #keep only the first date of an ICD9 instance
  patientEntry <- df[order(df$ICD9,df$Date),]
  patientEntry <- patientEntry[!duplicated(patientEntry$ICD9),]
  
  #Repeat the PID for the number of distinct ICD9 numbers
  PID <- rep(PID,nrow(patientEntry))
  colnames(patientEntry) <- c("Date","ICD9")
  patientDataFrame <- cbind(PID,patientEntry)
  
  return(patientDataFrame)
}

ICD9.df <- bind_rows(lapply(lines, function(x) processLine(x)))

#Write data in CSV format to files
write.csv(attributeData,file="attrib.csv",row.names = FALSE)
write.csv(ICD9.df,file="id.icd9.count.csv",row.names = FALSE)

#***The following code produces the uniquePhecodeMappedSelected data frame which assigns temporal
#data to each unique phecode per PID

#retrive a data frame with columns PID Date ICD9
dateICD9 <- bind_rows(lapply(lines, function(x) processLineICD9Date(x)))

#merge dateICD9 and cases data frame by PID
mergedCasesDateICD9 <- merge(dateICD9,cases,by="PID")

#two PID in cases are not in dateICD9
notIn <- cases$PID[!is.element(cases$PID, mergedCasesDateICD9$PID)]

#Determine the patient's (PID) year at each of their ICD9s
yearsatICD9 <- as.numeric(difftime(mergedCasesDateICD9$Date,mergedCasesDateICD9$DOB, unit="weeks"))/52.25
mergedCasesDateICD9 <- cbind(mergedCasesDateICD9,yearsatICD9,mergedCasesDateICD9$ICD9)
#Rename for compatibility with mapICD9ToPhecodes()
colnames(mergedCasesDateICD9)[3] <- "icd9"

#map ICD9s to phecodes (many to many relationship)
phecodeMapped <- mapICD9ToPhecodes(mergedCasesDateICD9)

#some PIDs did not make it in.  Not all ICD9s map to a phecode
notIn <- mergedCasesDateICD9$PID[!is.element(mergedCasesDateICD9$PID,phecodeMapped$PID)]

temporal <- character(length = nrow(phecodeMapped))

#Assign to default "T" transition (phecode occurs within 1 year of initial homelessness date)
temporal[1:length(temporal)] <- "T"

#Assign temporal codes
temporal[which((phecodeMapped$AGE_FIRST_HLSTATUS-phecodeMapped$yearsatICD9) > 1)] <- "B"
temporal[which((phecodeMapped$AGE_FIRST_HLSTATUS-phecodeMapped$yearsatICD9) < -1)] <- "A"

phecodeMapped <- cbind(phecodeMapped,temporal)

#Order by PID, then phecode, then year at ICD9
phecodeMapped <- phecodeMapped[order(phecodeMapped$PID,phecodeMapped$phecode,phecodeMapped$yearsatICD9),]
#Remove duplicate PID, Phecode combination (Note: first phecode instance will remain)
uniquePhecodeMapped <- phecodeMapped[!duplicated(phecodeMapped[c(1,10)]),]

#Remove uneeded columns
uniquePhecodeMappedSelected <- uniquePhecodeMapped[,c(1,10,11)]