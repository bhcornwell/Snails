#Script to take cleaned files and start parentage analysis
library(related)
setwd("/Users/Lab4349/Dropbox/2018 Nucella Data/Geneious Data Sept17/")
genos<-read.table("NL_multiplexG_H.csv",sep=',',header=T)
species<-as.character("NL")
#ID pops and offspring

#relatedness needs individuals to be represented by two columns, with each row a new marker:
n<-ncol(genos)-1
  #Remove QC Column:
ReformGenos<-genos[,1:n]
ReformGenos<-ReformGenos[which(substr(ReformGenos$Name,5,5)!="l"),]#Excludes RBC samples (historical from museum for COI purposes)
#Count number of unique Clutches (each can contain multiple capsules, and each of those can contain multiple eggs)
x <- as.character(ReformGenos$Name)
x1<-substr(x,6,6)
ReformGenos$Class <- paste(x1)
ReformGenos$Class
Eggs<-ReformGenos[which(ReformGenos$Class=="E"),]
Adults<-ReformGenos[which(ReformGenos$Class!="E"),]
Pops<-unique(substr(Adults$Name,5,5))

#Calculate allele frequencies uniquely for each population

Frequencies=list()
for(i in 1:length(Pops)){
  SubPop<-Adults[which(substr(Adults$Name,5,5)==Pops[i]),]
  t<-list()
  for(j in 1:((ncol(Adults)/2)-1)){
    Loc<-unlist(SubPop[,(j*2):((j*2)+1)])
    Alleles<-as.data.frame(table(Loc))
    AlName<-unlist(Alleles$Loc)
    Sum<-sum(Alleles$Freq)
    Freq<-Alleles$Freq#/Sum #Took out /Sum to get counts
    name<-colnames(SubPop[j*2])
    t[[j]]<-data.frame(Freq,AlName)#t[[name]]<-data.frame(Freq,AlName)
  }
  nameTwo<-as.character(Pops[i])
  Frequencies[[i]]<-t#Frequencies[[nameTwo]]<-t
}

#Need to separate out eggs into their respective populations

n<-unique(substr(as.character(Eggs$Name),7,8)) #Clutch IDs
m<-unique(substr(as.character(Eggs$Name),7,11)) #Unique Capsules
numberCap<-0
numberClu<-0
#Extract number of genotyped eggs in clutches
for(j in 1:length(n)){
  Temp<-Eggs[which(substr(as.character(Eggs$Name),7,8)==n[j]),]
  numberClu[j]<-nrow(Temp)
}
#Extract number of genotyped eggs in capsules
for(j in 1:length(m)){
   Temp<-Eggs[which(substr(as.character(Eggs$Name),7,11)==m[j]),]
    numberCap[j]<-nrow(Temp)
  }

Names<-names(ReformGenos)
Names<-Names[c(F,T)]
Names<-Names[1:length(Names)-1]

###need to calculate allele frequencies from adults, declare no "Candidate" males or females
###needs to be run on a per population basis

#######
for(i in 1:length(Pops)){
TempFreqs<-Frequencies[[i]]
TempEggs<-Eggs[which(substr(as.character(Eggs$Name),5,5)==Pops[i]),]
#Adds a dummy locus in adult population to ensure all alleles detected in the offspring are "found" in the population
AlCount<-0
for(p in 1:((ncol(TempEggs)/2)-1)){
  EggAlleles<-unique(unlist(TempEggs[,((p*2)):((p*2)+1)]))
  EggAlleles<-EggAlleles[!is.na(EggAlleles)]
  AlName<-as.factor(setdiff(EggAlleles,TempFreqs[[p]]$AlName))
  for(q in 1:length(Miss)){
    Len<-length(TempFreqs[[p]]$AlName)
    Freq<-rep(1,length.out=length(AlName))
    New<-data.frame(Freq,AlName)
  }
  TempFreqs[[p]]<-rbind(TempFreqs[[p]],New)
  AlCount[p]<-length(TempFreqs[[p]]$Freq)
}


if(nrow(TempEggs)!=0){
TempEggs[is.na(TempEggs)]<-0
Colony<-matrix(nrow=(21+nrow(TempEggs)+2+1+1+1+nrow(TempEggs)+numberClu+1+1+1+1),ncol=sum(numberClu)+1)#to account for rows needed in Colony2.dat file #Add more rows for the allele frequency lines
Colony[1,1]<-species #Project Name
Colony[2,1]<-paste(species,Pops[1],as.character("Colony"),sep="") #output file name
Colony[3,1]<-nrow(TempEggs) #number of offspring in sample
Colony[4,1]<-(ncol(Adults)/2)-1 #number of loci
Colony[5,1]<-sample(1:10000,1) #seed for random number generator
Colony[6,1]<-1 #not updating/updating allele freq 0/1
Colony[7,1]<-2 #dioecious/monoecious species 2/1
Colony[8,1]<-0 #0/1=no inbreeding/inbreeding
Colony[9,1]<-0 #0/1=Diploid species/HaploDiploid species
Colony[10,1]<-0 #0/1=Polygamy/Monogamy for males & females
Colony[10,2]<-0 #0/1=Polygamy/Monogamy for males & females
Colony[11,1]<-0 #0/1=Clone inference =No/Yes
Colony[12,1]<-1 #0/1=Scale full sibship=No/Yes   ######### Not sure what is correct choice here 4.11.18
Colony[13,1]<-0 #The program manual has two things on line 13, so using 1 to give weak prior and indication to use allele freqs. 0/1/2/3=No/Weak/Medium/Strong sibship prior; 4=optimal sibship prior ######### Not sure what is correct choice here 4.11.18
Colony[14,1]<-1 #0/1=Unknown/Known population allele frequency ################# Should this be given based on adults?
Colony[15,1:((ncol(Adults)/2)-1)]<-as.matrix(AlCount)
#Subscripts in Colony[x...y...z,1] are going to get messed up, fix after for loop works
#Insert allele freqs: 1 line of allele names (microsat #) next line their frequencies in population
count<-16
for(k in 1:((ncol(Adults)/2)-1)){
  y<-as.character(TempFreqs[[k]]$AlName)
  z<-as.numeric(TempFreqs[[k]]$Freq)
  for(j in 1:length(y)){
    Colony[count,j]<-y[j]
    Colony[count+1,j]<-z[j]/sum(z)
  }
  count<-count+2
}

Colony[count,1]<-1 #Number of runs
Colony[count+1,1]<- 2 #1/2/3/4=short/medium/long/very long run
Colony[count+2,1]<-0 #0/1=Monitor method by Iterate#/Time in second must be 0 when not running with Windows GUI
Colony[count+3,1]<-100 #Monitor interval in Iterate# / in seconds must be 0 when not running with Windows GUI
Colony[count+4,1]<-0 #0/1=No/Yes for run with Windows GUI
Colony[count+5,1]<-2 #0/1/2=PairLikelihood score/Fulllikelihood/FPLS ######### Not sure what is correct choice here 4.11.18
Colony[count+6]<-2 #0/1/2/3=Low/Medium/High/Very high precision with Fulllikelihood

Colony[count+7,1:length(Names)]<-Names
Colony[count+8,1]<-as.character("0@") #denotes codominant alleles
Colony[count+9,1]<-as.character("0.0001@") #allelic dropout rate, maybe should just be 0?
Colony[count+10,1]<-as.character("0.0001@") #allelic dropout rate, maybe should just be 0?
Colony[(count+11):(count+10+nrow(TempEggs)),1:(ncol(TempEggs)-1)]<-as.matrix(TempEggs[,1:(ncol(TempEggs)-1)]) #Offspring Genotypes
count<-count+10+nrow(TempEggs)
whoa<-matrix(0,nrow=nrow(TempEggs),ncol=2)
Colony[(count+1):(count+nrow(TempEggs)),1:2]<-as.matrix(whoa)
count<-count+nrow(TempEggs)
Colony[count+1,1:2]<-0 #Says no candidate males or females
Colony[count+2,1]<-0 #number of known paternity 
Colony[count+2,2]<-0 #exclusion threshold from known paternity
Colony[count+3,1]<-0 #number of known maternity 
Colony[count+3,2]<-0 #exclusion threshold from known maternity
Colony[count+4,1]<-0 #known paternal sibships
#Colony[count+5,1]<-length(unique(substr(as.character(TempEggs$Name),7,8))) #known maternal sibships, should this be by capsule or clutch? Currently by clutch.
for(q in 1:length(unique(substr(as.character(TempEggs$Name),7,11)))){
  o<-unique(substr(as.character(TempEggs$Name),7,11))
  s<-TempEggs[which(substr(as.character(TempEggs$Name),5,5)==Pops[i]),]
  l<-s[which(substr(as.character(s$Name),7,11)==o[q]),]
  Colony[count+4+q,1]<-nrow(l)
  Colony[count+4+q,1:nrow(s[which(substr(as.character(s$Name),7,11)==o[q]),])+1]<-as.character(l[1:nrow(l),1])
      }
count<-count+4+q
Colony[count+1,1]<-0 #known excluded paternal sibships
Colony[count+2,1]<-0 #known excluded maternal sibships
write.table(Colony,paste(species,Pops[i],as.character("Colony"),sep="_"),row.names = F, col.names = F, na= "",quote=F)
  }
} 


