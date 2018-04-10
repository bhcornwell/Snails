# COMPARATIVE POPULATION GENETICS OF NUCELLA  SP. MICROSATELLITE DATA

# MERGING FILES

setwd("/Users/christineewers/Dropbox/Nucella/Microsatellite_development/Geneious Data Sept17")

# make function to count number of amplified loci
qc <- function(input) {
  for (i in 1:nrow(input)) {
    x[i] <- (length(which(input[i,2:(ncol(input))] > 0)))/2
  }
  return(as.numeric(x))
}

# read and clean data

#mp <- read.csv("NC_multiplexD_E_G_H.csv", sep=",", header=T)
#mp <- read.csv("NE_multiplexD_E_G_H.csv", sep=",", header=T)
#mp <- read.csv("NF_Nlw1_2.csv", sep=",", header=T)
#mp <- read.csv("NH_Nlw1_2.csv", sep=",", header=T)
#mp <- read.csv("NI_multiplexD_E_G_H.csv", sep=",", header=T)
#mp <- read.csv("NL_multiplexG_H.csv", sep=",", header=T)
#mp <- read.csv("NO_multiplexD_E.csv", sep=",", header=T)
#mp <- read.csv("NP_Nlw1_2_3.csv", sep=",", header=T)

head(mp)
t(apply(data.frame(mp$Name, mp$Name.1, mp$Name.2, mp$Name.3), 1, FUN=duplicated))
t(apply(data.frame(mp$Name, mp$Name.1), 1, FUN=duplicated))
# this shows that the same individuals are in each row, and we can remove the three duplicated Name columns
mp <- mp[,!colnames(mp) == "Name.1" & !colnames(mp) == "Name.2" & !colnames(mp) == "Name.3" & !colnames(mp) == "PET...1"]

# remove well labels (first 4 characters) and sequencer code (last 4 characters)
x <- as.character(mp$Name)
mp$Name <- substr(x, 5, nchar(x)-4)
mp$Name

#remove EMPTY columns
empty <- mp[mp$Name == "EMPTY",]
mp <- mp[!mp$Name == "EMPTY",]
head(mp)

# replace "No peaks" and "No peaks in locus" with NA
mp[mp == "No peaks"] <- NA
mp[mp == "No peaks in locus"] <- NA
mp[mp == "Unbinned peaks in locus"] <- NA
mp[mp == "Unbinned peaks"] <- NA
mp[mp == "Peaks outside loci"] <- NA
mp[mp == "No peaks in bins"] <- NA
mp[mp == "Multiple peaks"] <- NA

mp[,2:ncol(mp)] <- apply(mp[,2:ncol(mp)], 2, FUN=as.character)
mp[,2:ncol(mp)] <- apply(mp[,2:ncol(mp)], 2, FUN=as.numeric)
str(mp)


# how many loci amplified per ind?

QC <- qc(mp)
mp$qc <- QC[!is.na(QC)]

# distribution of missing loci
hist(QC)
(ncol(mp)-2)/2


# save file

#write.table(mp, "NC_multiplexD_E_G_H.csv", sep=",")
#write.table(mp, "NE_multiplexD_E_G_H.csv", sep=",")
#write.table(mp, "NF_Nlw1_2.csv", sep=",")
#write.table(mp, "NH_Nlw1_2.csv", sep=",")
#write.table(mp, "NI_multiplexD_E_G_H.csv", sep=",")
#write.table(mp, "NL_multiplexG_H.csv", sep=",")
#write.table(mp, "NO_multiplexD_E.csv", sep=",")
#write.table(mp, "NP_Nlw1_2_3.csv", sep=",")
