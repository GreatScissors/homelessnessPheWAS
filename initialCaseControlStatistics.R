#clear workspace
#rm(list = ls())

#load raw case and control data
cases <- read.csv("HLCases02.csv")
controls <- read.csv("HLControls02.csv")

#sanity test
length(intersect(cases$PID,controls$PID))

#***Cases***

#Raw Summary Statistics
summary(cases[,3:6])

#Age Statistics

#https://www.tutorialspoint.com/r/r_mean_median_mode.htm
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

median(cases$AGE)
getmode(cases$AGE)
mean(cases$AGE)
sd(cases$AGE)

hist(cases$AGE,probability = T)

lines(density(cases$AGE), # density plot
      lwd = 2, # thickness of line
      col = "aquamarine4")

abline(v = mean(cases$AGE),
       col = "palegreen3",
       lwd = 2)

abline(v = median(cases$AGE),
       col = "darkolivegreen3",
       lwd = 2)

legend(x = "topright", # location of legend within plot area
       c("Density plot", "Mean", "Median"),
       col = c("aquamarine4", "palegreen3", "darkolivegreen3"),
       lwd = c(2, 2, 2))

#Initial Age of Homelessness Status Statistics
median(cases$AGE_FIRST_HLSTATUS)
getmode(cases$AGE_FIRST_HLSTATUS)
mean(cases$AGE_FIRST_HLSTATUS)
sd(cases$AGE_FIRST_HLSTATUS)

hist(cases$AGE_FIRST_HLSTATUS,probability = T)

lines(density(cases$AGE_FIRST_HLSTATUS), # density plot
      lwd = 2, # thickness of line
      col = "aquamarine4")

abline(v = mean(cases$AGE_FIRST_HLSTATUS),
       col = "palegreen3",
       lwd = 2)

abline(v = median(cases$AGE_FIRST_HLSTATUS),
       col = "darkolivegreen3",
       lwd = 2)

legend(x = "topright", # location of legend within plot area
       c("Density plot", "Mean", "Median"),
       col = c("aquamarine4", "palegreen3", "darkolivegreen3"),
       lwd = c(2, 2, 2))

#Sex Statistics
summary(cases$SEX)/sum(summary(cases$SEX))
plot(cases$SEX)

#Exact Binomial Test for proportions
binom.test(length(which(cases$SEX == "F")),length(cases$SEX),0.519,alternative="two.sided")

#Racial Statistics
#Data (percents) derived from 2015 US Census for Davidson County 
racePerc <- c(0.655,0.28,0.10,0.005,0.035)
names(racePerc) <- c("W","B","H","I","A")
barplot(racePerc)

summary(cases$RACE)
summary(cases$RACE)/sum(summary(cases$RACE))*100
plot(cases$RACE)

#Exact Binomial Test for proportions
binom.test(length(which(cases$RACE == "W")),length(cases$RACE),0.655,alternative="two.sided")
binom.test(length(which(cases$RACE == "B")),length(cases$RACE),0.28,alternative="two.sided")
binom.test(length(which(cases$RACE == "H")),length(cases$RACE),0.10,alternative="two.sided")
binom.test(length(which(cases$RACE == "I")),length(cases$RACE),0.005,alternative="two.sided")
binom.test(length(which(cases$RACE == "A")),length(cases$RACE),0.035,alternative="two.sided")

#***Controls***

#Raw Summary Statistics
summary(controls[,4:6])

#Age Statistics
median(controls$AGE)
getmode(controls$AGE)
mean(controls$AGE)
sd(controls$AGE)

hist(controls$AGE,probability = T)

lines(density(controls$AGE), # density plot
      lwd = 2, # thickness of line
      col = "aquamarine4")

abline(v = mean(controls$AGE),
       col = "palegreen3",
       lwd = 2)

abline(v = median(controls$AGE),
       col = "darkolivegreen3",
       lwd = 2)

legend(x = "topright", # location of legend within plot area
       c("Density plot", "Mean", "Median"),
       col = c("aquamarine4", "palegreen3", "darkolivegreen3"),
       lwd = c(2, 2, 2))

#Sex Statistics
summary(controls$SEX)/sum(summary(controls$SEX))
plot(controls$SEX)

#Racial Statistics
summary(controls$RACE)
summary(controls$RACE)/sum(summary(controls$RACE))*100
plot(controls$RACE)


summary(cases[which(cases$RACE=='W'),]$SEX)
summary(cases[which(cases$RACE=='B'),]$SEX)
summary(controls[which(controls$RACE=='W'),]$SEX)
summary(controls[which(controls$RACE=='B'),]$SEX)