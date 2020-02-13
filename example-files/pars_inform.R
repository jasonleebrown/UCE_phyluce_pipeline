#script to count parsimony informative sites and output to csv
#load phyloch
library(phyloch)

#setwd where the reads are
setwd("/home/bender/Desktop/tutorial/5_taxon-sets/all/muscle-nexus-clean-75p")

#get list of file names and number of files
listoffiles <- list.files(pattern="*.nex*")
nooffiles <- length(listoffiles)

#these are the column names
record <- c("locusname","pis","length")

#loop to calculate PIS and write the PIS info to a text file
for (j in 1:nooffiles) {
  write.table((gsub("?","N",(readLines(listoffiles[j])),fixed=TRUE)),"list_of_pis_by_locus.txt",sep="",quote=FALSE,row.names=FALSE,col.names=FALSE)
  tempfile <- read.nex("list_of_pis_by_locus.txt")
  templength <- dim(tempfile)[2]
  temppis <- pis(tempfile)
  temp <- cbind(listoffiles[j],temppis,templength)
  record <- rbind(record,temp)
}
#add column names to text file
write.table(record, "list_of_pis_by_locus.txt",quote=FALSE, row.names=FALSE,col.names=FALSE)

#find 200 most parsimony-informative loci and save their names
pistab <- read.delim("list_of_pis_by_locus.txt", header=TRUE, sep = " ")
pistab_sorted <- pistab[order(-pistab$pis),]
pistab_sorted_truncated <- pistab_sorted[1:200,]
topnames <- pistab_sorted_truncated$locusname
write.table(topnames, file = "top200names.txt")

###calculate summary data and visualize PIS variables###
par_data <- as.data.frame(read.table("list_of_pis_by_locus.txt", header = TRUE))
plot(par_data$pis) #visualize number of PIS
inform <- subset(par_data, pis >= 3) #subset data set for loci with PIS > 3
inform <- subset(inform, pis < 15) #remove outliers (change based on obvious outliers to data set)
plot(inform$pis) #visualize informative loci
length(inform$pis) #calculate the number of informative loci
sum(inform$pis) #total number of PIS for informative loci
fivenum(inform$pis) #summary stats for informative loci (min, lower quartile, median, upper quartile, max)
inform_names <- inform$locusname #get locus names of informative loci for locus filtering
#write locus names to a file
write.table(inform_names, file = "inform_names_PIS_15.txt")
