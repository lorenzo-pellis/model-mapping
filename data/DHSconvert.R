require(haven)
system('ls')
args <- commandArgs(trailingOnly = TRUE)
inName <- "ZAPR31FL.SAV" ##args[1]     ## name of output file.
outName <- "ZAPR31FL.CSV" ##args[2]    ## name of local data file.

data<-read_sav(inName)
data
d3<-data[c("HV001","HV002","HV005","HV009","HV105")]
write.csv(d3,outName,row.names=FALSE)

