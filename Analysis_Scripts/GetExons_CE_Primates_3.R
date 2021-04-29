


#This was added on September 22 of 2020:

#To get conserved circRNAs of primates we are going to filter potential exons that can form the circRNA
#For this we are going to use TWO filters, all CE with a PSI >=5 and another dataframe with a PSI >= 10

#Script to get conserved circRNAs of primates we are going to filter POTENTIAL EXONS that can form the circRNA
#For this we are going to use TWO filters, all CE with a PSI >=5 and another dataframe with a PSI >= 10

library(dplyr)

##First we read filter .psi.gz file with orthologues genes from each primate

path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/"

Human_Shared = read.delim(file=paste(path,"HumanThatAreOrth_PSI.txt", sep=""), 
                        header = TRUE, as.is = TRUE) #94,149

Chimp_Shared = read.delim(file=paste(path,"ChimpThatAreOrth_PSI.txt", sep=""), 
                          header = TRUE, as.is = TRUE) #98,792

Lemur_Shared = read.delim(file=paste(path,"LemurThatAreOrth_PSI.txt", sep=""), 
                          header = TRUE, as.is = TRUE) #96,760

Macaque_Shared = read.delim(file=paste(path,"MacaqueThatAreOrth_PSI.txt", sep=""), 
                            header = TRUE, as.is = TRUE) #97,986

Macaca_Shared = read.delim(file=paste(path,"MacacaThatAreOrth_PSI.txt", sep=""), 
                           header = TRUE, as.is = TRUE) #98,282

Marmoset_Shared = read.delim(file=paste(path,"MarmosetThatAreOrth_PSI.txt", sep=""), 
            header = TRUE, as.is = TRUE) #108,426

SquirrelMonkey_Shared = read.delim(file=paste(path,"SquirrelMonkeyThatAreOrth_PSI.txt", sep=""), 
            header = TRUE, as.is = TRUE) #97,100

Baboon_Shared = read.delim(file=paste(path,"BaboonThatAreOrth_PSI.txt", sep=""), 
            header = TRUE, as.is = TRUE) #97,612


###Make a function to read all .psi.gz files from each primate and filter according to its respective <Primate>_Share

Files_Shared <- function(Primate_Shared, path_primate) {
  
  files = list.files(path = path_primate, pattern = ".psi.gz")
  
  FilesDF = data.frame(Files = files, Name = sub("(_[^_]+)_.*", "\\1", files),stringsAsFactors = FALSE)
  
  Primate_PSI = vector(mode="list" ,length = nrow(FilesDF))
  names(Primate_PSI) = FilesDF$Name
  
  for (t in names(Primate_PSI)) {
    
    Primate_PSI[[t]] = read.delim(gzfile( paste(path_primate, FilesDF[grep(t, FilesDF$Name), "Files"], sep = "/") ),  
                                  header = TRUE, as.is = TRUE)
    
    #Some primates have in their Gene Id the Ensmebl_geneID<.>number
    #Remove the <.>number part
    Primate_PSI[[t]][["Gene"]] = gsub("\\..*","", Primate_PSI[[t]][["Gene"]])
    
    #Subset those that are "CE"
    Primate_PSI[[t]] = Primate_PSI[[t]][grep("CE", Primate_PSI[[t]][["Type"]]),]
    
    #Subset those that are in Shared file
    shared_genes = intersect(Primate_Shared$Gene, Primate_PSI[[t]][["Gene"]])
    
    #
    Primate_PSI[[t]] = Primate_PSI[[t]][Primate_PSI[[t]][["Gene"]] %in%  shared_genes,]
    
    
  }
  
  return(Primate_PSI)
  
}

#This takes a time

Human_PSI = Files_Shared(Human_Shared, "/home/gaby/lab_Garvan/Primates/Human_SamplesData_CircAtlas/WhippetOut")

Lemur_PSI = Files_Shared(Lemur_Shared, "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/Lemur")

Macaque_PSI = Files_Shared(Macaque_Shared, "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/Macaque")
Macaca_PSI = Files_Shared(Macaca_Shared, "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/Macaca")
Marmoset_PSI = Files_Shared(Marmoset_Shared, "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/Marmoset")
SquiMonkey_PSI = Files_Shared(SquirrelMonkey_Shared, "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/SquirrelMonkey/")
Baboon_PSI = Files_Shared(Baboon_Shared, "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/Baboon/")

Chimp_PSI = Files_Shared(Chimp_Shared, "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/Chimp/")

###Make a function to filter each tissue sample of each primate according to PSI filter in each tissue (sample)
PSI_filter <- function(Primate_PSI, psi_cutoff) {
  
  for (t in names(Primate_PSI)) {
    
    Primate_PSI[[t]] = Primate_PSI[[t]][which(Primate_PSI[[t]][["Psi"]] >= psi_cutoff),]
    
  }
  
  return(Primate_PSI)
  
}


Human_PSI_10 = PSI_filter(Human_PSI, 0.1)
Human_PSI_5 = PSI_filter(Human_PSI, 0.05)

Lemur_PSI_10 = PSI_filter(Lemur_PSI, 0.1)
Lemur_PSI_5 = PSI_filter(Lemur_PSI, 0.05)

Macaque_PSI_10 = PSI_filter(Macaque_PSI, 0.1)
Macaque_PSI_5 = PSI_filter(Macaque_PSI, 0.05)

Macaca_PSI_10 = PSI_filter(Macaca_PSI, 0.1)
Macaca_PSI_5 = PSI_filter(Macaca_PSI, 0.05)

Marmoset_PSI_10 = PSI_filter(Marmoset_PSI, 0.1)
Marmoset_PSI_5 = PSI_filter(Marmoset_PSI, 0.05)

SquiMonkey_PSI_10 = PSI_filter(SquiMonkey_PSI, 0.1)
SquiMonkey_PSI_5 = PSI_filter(SquiMonkey_PSI, 0.05)

Baboon_PSI_10 = PSI_filter(Baboon_PSI, 0.1)
Baboon_PSI_5 = PSI_filter(Baboon_PSI, 0.05)

Chimp_PSI_10 = PSI_filter(Chimp_PSI, 0.1)
Chimp_PSI_5 = PSI_filter(Chimp_PSI,0.05)

#Make a single dataframe with the "Gene", "Node", "Coord" of all "CE" events that pass the PSI filter
JoinPSI_DF <-function(Primate_PSI_val, cols) {
  
  #Filter columns of interest
  ColsDF = lapply(Primate_PSI_val, "[", , cols)
  
  ##Add and ID to each CE ("Gene"_"Coord")
  for (t in names(ColsDF)) {
   
    ColsDF[[t]][["ID"]] = paste(ColsDF[[t]][["Gene"]],ColsDF[[t]][["Coord"]], sep = "_" ) 
    
    
  }
  
  #Make a dataframe with all values
  PrimateCoords = bind_rows(ColsDF, .id = "ID") #The ID will become the sample/tissue where it comes
  
  #Remove duplicated rows 
  PrimateCoords = unique( PrimateCoords[ , c("Gene","Node","Coord") ] )
  
  ##Add ID again
  PrimateCoords$ID = paste(PrimateCoords$Gene, PrimateCoords$Coord, sep="_")
  
  return(PrimateCoords)
  
}


col = c("Gene", "Node", "Coord")

HumanCoord10 = JoinPSI_DF(Human_PSI_10, col) #91,766
HumanCoord5 = JoinPSI_DF(Human_PSI_5, col) #91,978

LemurCoord10 = JoinPSI_DF(Lemur_PSI_10, col) #90,770
LemurCoord5 = JoinPSI_DF(Lemur_PSI_5, col) #90,825

MacaqueCoord10 = JoinPSI_DF(Macaque_PSI_10, col) #93,547
MacaqueCoord5 = JoinPSI_DF(Macaque_PSI_5, col) #93,662

MacacaCoord10 = JoinPSI_DF(Macaca_PSI_10, col) #94,182
MacacaCoord5 = JoinPSI_DF(Macaca_PSI_5, col) #94,273

MarmosetCoord10 = JoinPSI_DF(Marmoset_PSI_10, col) #96,019
MarmosetCoord5 = JoinPSI_DF(Marmoset_PSI_5, col) #96,143

SquiMonkeyCoord10 = JoinPSI_DF(SquiMonkey_PSI_10, col) #92,204
SquiMonkeyCoord5 = JoinPSI_DF(SquiMonkey_PSI_5, col) #92,255

BaboonCoord10 = JoinPSI_DF(Baboon_PSI_10, col) #92,935
BaboonCoord5 = JoinPSI_DF(Baboon_PSI_5, col) #93,004

ChimpCoord10 = JoinPSI_DF(Chimp_PSI_10, col) #94,270
ChimpCoord5 = JoinPSI_DF(Chimp_PSI_5, col) #94,375,



##SAVE info <Primate>Coord<PSIcutoff> for retriving the coords to use on bedintersect and LiftOver
path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE"

#######################For PSI 10
write.table(HumanCoord10, file=paste(path, "PSI_10","HumanThatAreOrth_PSI_10.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(ChimpCoord10, file=paste(path, "PSI_10", "ChimpThatAreOrth_PSI_10.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(LemurCoord10, file=paste(path, "PSI_10","LemurThatAreOrth_PSI_10.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(MacaqueCoord10, file=paste(path, "PSI_10","MacaqueThatAreOrth_PSI_10.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(MacacaCoord10, file=paste(path, "PSI_10","MacacaThatAreOrth_PSI_10.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(MarmosetCoord10, file=paste(path, "PSI_10","MarmosetThatAreOrth_PSI_10.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(SquiMonkeyCoord10, file=paste(path, "PSI_10","SquirrelMonkeyThatAreOrth_PSI_10.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(BaboonCoord10, file=paste(path, "PSI_10","BaboonThatAreOrth_PSI_10.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")



##############For PSI 5
write.table(HumanCoord5, file=paste(path, "PSI_5","HumanThatAreOrth_PSI_5.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(ChimpCoord5, file=paste(path, "PSI_5","ChimpThatAreOrth_PSI_5.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(LemurCoord5, file=paste(path, "PSI_5","LemurThatAreOrth_PSI_5.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(MacaqueCoord5, file=paste(path, "PSI_5","MacaqueThatAreOrth_PSI_5.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(MacacaCoord5, file=paste(path, "PSI_5","MacacaThatAreOrth_PSI_5.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(MarmosetCoord5, file=paste(path, "PSI_5","MarmosetThatAreOrth_PSI_5.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(SquiMonkeyCoord5, file=paste(path, "PSI_5","SquirrelMonkeyThatAreOrth_PSI_5.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(BaboonCoord5, file=paste(path, "PSI_5","BaboonThatAreOrth_PSI_5.txt", sep="/"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")


###SAVE List of whippet output of each Primate to later check if the circRNA and exon have accordingly PSI values 
#to be considered EXPRESS in such tissue

######PSI 10
save(Human_PSI_10, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_10/HumanList_PSI10")
save(Baboon_PSI_10, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_10/BaboonList_PSI10")
save(Lemur_PSI_10, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_10/LemurList_PSI10")
save(Macaca_PSI_10, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_10/MacacaList_PSI10")
save(Macaque_PSI_10, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_10/MacaqueList_PSI10")
save(Marmoset_PSI_10, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_10/MarmosetList_PSI10")
save(SquiMonkey_PSI_10, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_10/SquiMonkeyList_PSI10")

save(Chimp_PSI_10, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_10/ChimpList_PSI10")

######PSI 5
save(Human_PSI_5, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_5/HumanList_PSI5")
save(Baboon_PSI_5, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_5/BaboonList_PSI5")
save(Lemur_PSI_5, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_5/LemurList_PSI5")
save(Macaca_PSI_5, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_5/MacacaList_PSI5")
save(Macaque_PSI_5, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_5/MacaqueList_PSI5")
save(Marmoset_PSI_5, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_5/MarmosetList_PSI5")
save(SquiMonkey_PSI_5, file = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_5/SquiMonkeyList_PSI5")

save(Chimp_PSI_5, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_5/ChimpList_PSI5")




