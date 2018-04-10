# Snails
These files are the analysis pipeline for the Nucella project. Trying to keep naming conventions the same as when Christine started them.

File 1: Cleaning files
  -Counts number of amplified loci
  -Reads in the csv files of individuals and loci (these will probably change to some master file as more loci are added)
  -Removes duplicated individual names
  -Replaces Geneious calles of "No Peaks", "Unbinned Peaks", etc... with NA
  -Writes files as .csv when cleaned
File 2: Basic Stats
  -Reads in files
  -Replaces mid-IDd individual names with new names to reflect the correct species ID (4th character of individual ID)
  -Characterizes the species and population for each individual (extracted from the individual names)
  -Finds # of loci
  -Reformats to be compatible with genind
  -Calculates allelic richness, Hs (heterozygosity?), using PopGenReport
  
