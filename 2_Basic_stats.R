# COMPARATIVE POPULATION GENETICS OF NUCELLA  SP. MICROSATELLITE DATA

# BASIC POP GEN

setwd("/Users/christineewers/Dropbox/Nucella/Microsatellite_evolution/Comparative_microsats")
setwd("/Users/Rick/Dropbox/Nucella/Microsatellite_evolution/Comparative_microsats")
setwd("/Users/christineewers/Dropbox/Nucella/Microsatellite_development/Geneious Data Sept17")

#mp <- read.csv("NC_multiplexD_E_G_H.csv", sep=",", header=T)
#mp <- read.csv("NE_multiplexD_E_G_H.csv", sep=",", header=T)
#mp <- read.csv("NF_Nlw1_2.csv", sep=",", header=T)
#mp <- read.csv("NH_Nlw1_2.csv", sep=",", header=T)
#mp <- read.csv("NI_multiplexD_E_G_H.csv", sep=",", header=T)
#mp <- read.csv("NL_multiplexG_H.csv", sep=",", header=T)
#mp <- read.csv("NO_multiplexD_E.csv", sep=",", header=T)
mp <- read.csv("NP_Nlw1_2_3.csv", sep=",", header=T)


### make species and pop columns ####

# ID example: 16NCEX01
# meaning: year (2 characters), species (2 characters), location (1 character), gender/eggs (1 character), ind ID (2 characters)

# rename mis-identified sp
x <- as.character(mp$Name)

# COI showed that Nucella from Petersburg close to the harbor are N. lima, not N. canaliculata. Nucella from Hungry Point are still N. canaliculata
NIE.old <- c("16NCEX09", "16NCEX10", "16NCEX11", "16NCEX12", "16NCEX13", "16NCEX14", "16NCEX15", "16NCEX16", "16NCEX17", "16NCEX18")
NIE.new <- c("16NIEX09", "16NIEX10", "16NIEX11", "16NIEX12", "16NIEX13", "16NIEX14", "16NIEX15", "16NIEX16", "16NIEX17", "16NIEX18")
x[x %in% NIE.old] <- NIE.new

# COI showed that Nucella from Prince Rupert are N. ostrina
NOP.old <- c("16NCPF01", "16NCPF02" ,"16NCPF03","16NCPF04", "16NCPF05", "16NCPF06", "16NCPF07", "16NCPF08", 
             "16NCPF09", "16NCPF10", "16NCPF11", "16NCPF12", "16NCPF13", "16NCPM01", "16NCPM02", "16NCPM03", 
             "16NCPM04", "16NCPM05", "16NCPM06", "16NCPM07", "16NCPM08", "16NCPM09", "16NCPM10", "16NCPM11",
             "16NCPM12", "16NCPM13", "16NCPM14")
NOP.new <- c("16NOPF01", "16NOPF02" ,"16NOPF03","16NOPF04", "16NOPF05", "16NOPF06", "16NOPF07", "16NOPF08", 
             "16NOPF09", "16NOPF10", "16NOPF11", "16NOPF12", "16NOPF13", "16NOPM01", "16NOPM02", "16NOPM03", 
             "16NOPM04", "16NOPM05", "16NOPM06", "16NOPM07", "16NOPM08", "16NOPM09", "16NOPM10", "16NOPM11",
             "16NOPM12", "16NOPM13", "16NOPM14")
x[x %in% NOP.old] <- NOP.new

x2 <- x[order(x)]
mp <- mp[order(x2),]

# define species for each ind
sp <- substr(x2, 3, 4)

# define pop for each ind
pop <- substr(x2, 5, 5)

### what data do we have? ####

# number of ind per pop and sp
table(pop, sp)

# loci names
colnames(mp)
x <- colnames(mp)[2:(ncol(mp)-1)] # get names for alleles only (not for ind ID, pop or sp)
#x <- colnames(mp)[2:(ncol(mp))] # get names for alleles only (not for ind ID, pop or sp)
x
lociNames <- unique(substr(x,1,nchar(x)-4))

# number of loci
noLoci <- length(lociNames)
noLoci 
# 27 loci for NC
# 36 loci for NE
#  6 loci for NF
#  6 loci for NH
# 39 loci for NI
# 22 loci for NL
# 17 loci for NO
# 13 loci for NP

### bring table in correct formats for genind ####

# duplicate each loci name
lociNamesDup <- substr(x,1,nchar(x)-4)

# allele calls only without ID names
mp2 <- mp[,2:(ncol(mp)-1)]
#mp2 <- mp[,2:(ncol(mp))]
noPop <- length(unique(pop))

# concatenate allele calls into same column, separated by "/"
loci <- matrix(ncol=noLoci, nrow=nrow(mp2))
colnames(loci) <- lociNames
rownames(loci) <- mp$Name

for (i in 1:noLoci) {
  for (j in 1:nrow(loci)) {
    z <- as.character(mp2[j,lociNamesDup == lociNames[i]])
    loci[j,i] <- paste0(z, collapse="/")
  }
}

head(loci)
nrow(loci)
library(stringr)
colnames(loci) <- str_replace_all(colnames(loci), "\\.", "")

write.table(loci, "GenindInput_NP.csv")

### allelic diversity ####


library(adegenet)
g <- df2genind(loci, sep="/", ncode=3, ind.names = mp$Name, NA.char = "NA/NA",
               ploidy=2, type="codom", pop=rep("pop", length=nrow(loci)))
#m@other$xy <- data.frame(mh$lat, mh$long)
loc.names <- unique(g@loc.fac)
pop.names <- unique(pop2)

# number of alleles per locus
write.table(g@loc.n.all, "noAlleles_NP.csv", sep=",")

library(PopGenReport)
#par(mfcol=c(8,5))
pdf("allele_freq_heatmap_NC.pdf", height=4, width=4)
awesome <- allele.dist(g)
dev.off()


Hs(g)
#      NCE       NCL       NIE       NLN       NLS       NLW       NLZ       NOC       NOK       NON 
#0.6469253 0.7400626 0.6595394 0.4912855 0.5108072 0.4977905 0.5242984 0.3858710 0.3797123 0.4218933 
#NOP       NPR 
#0.6658375 0.4822187 

library(PopGenReport)
no.chr <- allel.rich(g)$pop.sizes
x <- allel.rich(g)
x
# doesnt calculate allelic div because all loci are completely NA for some pop


allel.rich.fake <- allel.rich(g)$all.richness
write.table(allel.rich.fake, "AllelicRichness_NP.csv", sep=",")

### number and range of alleles per locus ####

range <- sapply(g@all.names, range)
write.table(range, "AlleleRange_NP.csv", sep=",")
# not sure why there are NA in max
