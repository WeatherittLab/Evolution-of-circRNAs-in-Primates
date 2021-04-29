
##Script to make heatmap of all circRNAs
library(pheatmap)

#Get conserved and non-conserved circRNAs of all primates, get PSI values across tissues


InfoConservedCircRNAs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ConservedCircRNAs_Bob1.txt",
                                   header = TRUE, as.is = TRUE) #11,974


InfoTissueConservedCircRNAs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/TissueConservedCircRNAs_Bob2.txt",
                                         header = TRUE, as.is = TRUE) #11,735


#Get all primates conserved and non conserved circRNAs


HumanCircRNAs10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/HumanInfoCircRNAs10.txt",
            header = TRUE, as.is = TRUE) #50,002


BaboonHumanCircs10Consv = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/BaboonHumanCirc10.txt",
                                     header = TRUE, as.is = TRUE)

BaboonSpeciesCircs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/BaboonHumanCirc10ExonsConvsNoBSJ.txt",
                                header = TRUE, as.is = TRUE)
###Add column of specie
BaboonSpeciesCircs$Primate = "Baboon"
BaboonSpeciesCircs$SpeciesID = 9555

LemurHumanCircs10Consv = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/LemurHumanCirc10.txt",
            header = TRUE, as.is = TRUE)  

LemurSpeciesCircs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/LemurHumanCirc10ExonsConsvNoBSJ.txt",
                               header = TRUE, as.is = TRUE) 
##Add column of specie
LemurSpeciesCircs$Primate = "Lemur"
LemurSpeciesCircs$SpeciesID = 30608

MacaqueHumanCircs10Consv = read.delim( file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MacaqueHumanCirc10.txt",
                                       header = TRUE, as.is = TRUE)

MacaqueSpeciesCircs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MacaqueHumanCirc10ExonsConsvNoBSJ.txt",
                                 header = TRUE, as.is = TRUE)

##Add column of specie
MacaqueSpeciesCircs$Primate = "Macaque"
MacaqueSpeciesCircs$SpeciesID = 9544


MacacaHumanCircs10Consv = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MacacaHumanCirc10.txt",
                                     header = TRUE, as.is = TRUE)  

MacacaSpeciesCircs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MacacaHumanCirc10ExonsConsvNoBSJ.txt",
                                header = TRUE, as.is = TRUE) 

##Add column of specie
MacacaSpeciesCircs$Primate = "Macaca"
MacacaSpeciesCircs$SpeciesID = 9541


MarmosetHumanCircs10Consv = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MarmosetHumanCirc10.txt",
                                       header = TRUE, as.is = TRUE)

MarmosetSpeciesCircs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MarmosetHumanCirc10ExonsConsvNoBSJ.txt",
                                  header = TRUE, as.is = TRUE)

##Add column of specie
MarmosetSpeciesCircs$Primate = "Marmoset"
MarmosetSpeciesCircs$SpeciesID = 9483

SquiMonkeyHumanCircs10Consv = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/SquiMonkeyHumanCirc10.txt",
                                         header = TRUE, as.is = TRUE)

SquiMonkeySpeciesCircs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/SquiMonkeyHumanCirc10ExonsConsvNoBSJ.txt",
                                    header = TRUE, as.is = TRUE)

##Add column of specie
SquiMonkeySpeciesCircs$Primate = "Squirrel Monkey"
SquiMonkeySpeciesCircs$SpeciesID = 9521

ChimpHumanCircs10Consv = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ChimpHumanCirc10.txt",
                                    header = TRUE, as.is = TRUE)


ChimpSpeciesCircs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ChimpHumanCirc10ExonsConsvNoBSJ.txt",
                               header = TRUE, as.is = TRUE)

##Add column of specie
ChimpSpeciesCircs$Primate = "Chimpanzee"
ChimpSpeciesCircs$SpeciesID = 9598
  
#Make a single file of Species-Specific circRNAs

PrimatesSpecific = rbind(ChimpSpeciesCircs, BaboonSpeciesCircs,
                         MacaqueSpeciesCircs, MacacaSpeciesCircs,
                         MarmosetSpeciesCircs, SquiMonkeySpeciesCircs,
                         LemurSpeciesCircs)

#write.csv(PrimatesSpecific, file = "~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/Primates_SpeciesSpecific_circRNA.csv",
#           quote=FALSE, row.names = FALSE )



###Read PSI values of circRNAs

Human_PSI = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Human_PSI.txt", 
                       header = TRUE, as.is = TRUE) #43,681
#Add colnames with primate name at the begining
colnames(Human_PSI) = paste("Human", colnames(Human_PSI), sep="_")


Baboon_PSI = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Baboon_PSI.txt", 
                        header = TRUE, as.is = TRUE) #16,123

colnames(Baboon_PSI) = paste("Baboon", colnames(Baboon_PSI), sep="_")


Lemur_PSI = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Lemur_PSI.txt", 
                       header = TRUE, as.is = TRUE) #10,173
colnames(Lemur_PSI) = paste("Lemur", colnames(Lemur_PSI), sep="_")


Macaca_PSI = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Macaca_PSI.txt", 
                        header = TRUE, as.is = TRUE) #16,262
colnames(Macaca_PSI) = paste("Macaca", colnames(Macaca_PSI), sep="_")


Macaque_PSI = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Macaque_PSI.txt", 
                         header = TRUE, as.is = TRUE) #17,580
colnames(Macaque_PSI) = paste("Macaque", colnames(Macaque_PSI), sep="_")


Marmoset_PSI = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Marmoset_PSI.txt",
                          header = TRUE, as.is = TRUE) #13,897
colnames(Marmoset_PSI) = paste("Marmoset", colnames(Marmoset_PSI), sep="_")


SquirrelMonkey_PSI = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/SquirrelMonkey_PSI.txt", 
                                header = TRUE, as.is = TRUE) #15,546
colnames(SquirrelMonkey_PSI) = paste("SquirrelMonkey", colnames(SquirrelMonkey_PSI), sep="_")


Chimp_PSI = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Chimp_PSI.txt",
                       header = TRUE, as.is = TRUE) #18,017
colnames(Chimp_PSI) = paste("Chimp", colnames(Chimp_PSI), sep="_")


Human_PSI = Human_PSI*100
Baboon_PSI = Baboon_PSI*100
Lemur_PSI = Lemur_PSI*100
Macaca_PSI = Macaca_PSI*100
Macaque_PSI = Macaque_PSI*100
Marmoset_PSI = Marmoset_PSI*100
SquirrelMonkey_PSI = SquirrelMonkey_PSI*100
Chimp_PSI = Chimp_PSI*100

##Subset Human_PSI with InfoConsvIDs
Human_PSI_filt = Human_PSI[rownames(InfoConservedCircRNAs),]



#Get primate IDs of consv circRNAs 
BaboonConsvIDs = BaboonHumanCircs10Consv$Primate_CircID  #1,952
LemurConsvIDs = LemurHumanCircs10Consv$Primate_CircID #397
MacacaConsvIDs = MacacaHumanCircs10Consv$Primate_CircID #1586
MacaqueConsvIDs = MacaqueHumanCircs10Consv$Primate_CircID #1661
MarmosetConsvIDs = MarmosetHumanCircs10Consv$Primate_CircID #1082
SquirrelMonkConsvIDs = SquiMonkeyHumanCircs10Consv$Primate_CircID #839
ChimpConsvIDs = ChimpHumanCircs10Consv$Primate_CircID #2,452

#Subset PSI of consv circRNAs, and change to human circRNA ID 
Baboon_ConsvPSI = Baboon_PSI[BaboonConsvIDs,]
Lemur_ConsvPSI = Lemur_PSI[LemurConsvIDs,]
Macaca_ConsvPSI = Macaca_PSI[MacacaConsvIDs,]
Macaque_ConsvPSI = Macaque_PSI[MacaqueConsvIDs,]
Marmoset_ConsvPSI = Marmoset_PSI[MarmosetConsvIDs,]
SquiMonk_ConsvPSI = SquirrelMonkey_PSI[SquirrelMonkConsvIDs,]
Chimp_ConsvPSI = Chimp_PSI[ChimpConsvIDs,]
      
####To ConsvPSI add a column with human IDs
Baboon_ConsvPSI$HumanID = BaboonHumanCircs10Consv$Human_CircID
Lemur_ConsvPSI$HumanID = LemurHumanCircs10Consv$Human_CircID
Macaca_ConsvPSI$HumanID = MacacaHumanCircs10Consv$Human_CircID
Macaque_ConsvPSI$HumanID = MacaqueHumanCircs10Consv$Human_CircID
Marmoset_ConsvPSI$HumanID = MarmosetHumanCircs10Consv$Human_CircID
SquiMonk_ConsvPSI$HumanID = SquiMonkeyHumanCircs10Consv$Human_CircID
Chimp_ConsvPSI$HumanID = ChimpHumanCircs10Consv$Human_CircID

####Check for duplicated human IDs on ConsvPSI
#Make a function to Order _ConsvPSI according to HumanID column, find duplicated HumanIDs
#Go through duplicated Human IDs and keep the row with the duplicated value that has less NAs across samples


#Baboon_ConsvPSI = Baboon_ConsvPSI[order(Baboon_ConsvPSI$HumanID),]


#rowSums(is.na(Baboon_ConsvPSI[grep("ENSG00000082269_17-8", Baboon_ConsvPSI$HumanID),]))

#temp[which.max(temp)]
#Baboon_ConsvPSI[!rownames(Baboon_ConsvPSI) %in% "ENSPANG00000000816_11-4",]

RemDupHumIDs <-function(Primate_ConsvPSI) {
  
  Primate_ConsvPSI = Primate_ConsvPSI[order(Primate_ConsvPSI$HumanID),]
  
  Ids_dup = unique(Primate_ConsvPSI$HumanID[duplicated(Primate_ConsvPSI$HumanID)])
  
  print(length(Ids_dup))
  
  id2rem = c()
  
  for (d in Ids_dup) {
    
    temp = rowSums(is.na(Primate_ConsvPSI[grep(d, Primate_ConsvPSI$HumanID),]))
    
    id2rem = c(id2rem, names(temp[which.max(temp)]))
    
  }
  
  Primate_ConsvPSI = Primate_ConsvPSI[!rownames(Primate_ConsvPSI) %in% id2rem,]
  
  #Check if there are still duplicated values
  if (sum(duplicated(Primate_ConsvPSI$HumanID)) > 0) {
    
    Ids_dup = unique(Primate_ConsvPSI$HumanID[duplicated(Primate_ConsvPSI$HumanID)])
    
    print(length(Ids_dup))
    
    id2rem = c()
    
    for (d in Ids_dup) {
      
      temp = rowSums(is.na(Primate_ConsvPSI[grep(d, Primate_ConsvPSI$HumanID),]))
      
      id2rem = c(id2rem, names(temp[which.max(temp)]))
      
    }
    
    Primate_ConsvPSI = Primate_ConsvPSI[!rownames(Primate_ConsvPSI) %in% id2rem,]
    
    #Rename PSI_filt with HumanIDs and remove HumanID column
    rownames(Primate_ConsvPSI) = Primate_ConsvPSI$HumanID
    Primate_ConsvPSI$HumanID = NULL
    
    
  } else {
  
  #Rename PSI_filt with HumanIDs and remove HumanID column
  rownames(Primate_ConsvPSI) = Primate_ConsvPSI$HumanID
  Primate_ConsvPSI$HumanID = NULL
  
  }
  
  return(Primate_ConsvPSI)
  
}

Baboon_ConsvPSI_filt = RemDupHumIDs(Baboon_ConsvPSI)
Lemur_ConsvPSI_filt = RemDupHumIDs(Lemur_ConsvPSI)
Macaca_ConsvPSI_filt = RemDupHumIDs(Macaca_ConsvPSI)
Macaque_ConsvPSI_filt = RemDupHumIDs(Macaque_ConsvPSI)
Marmoset_ConsvPSI_filt = RemDupHumIDs(Marmoset_ConsvPSI)
SquiMonk_ConsvPSI_filt = RemDupHumIDs(SquiMonk_ConsvPSI)
Chimp_ConsvPSI_filt = RemDupHumIDs(Chimp_ConsvPSI)

      
#Get primate IDs of non-conserv circRNAs for each primate
BaboonNonConsvIDs = BaboonSpeciesCircs$Primate_CircID #474
LemurNonConsvIDs =  LemurSpeciesCircs$Primate_CircID #133
MacacaNonConsvIDs = MacacaSpeciesCircs$Primate_CircID #399
MacaqueNonConsvIDs = MacaqueSpeciesCircs$Primate_CircID #412
MarmosetNonConsvIDs = MarmosetSpeciesCircs$Primate_CircID #335
SquirrelMonkNonConsvIDs = SquiMonkeySpeciesCircs$Primate_CircID #268
ChimpNonConsvIDs = ChimpSpeciesCircs$Primate_CircID #431
        
#With IDs subset <Primate>_PSI data
Baboon_NonConsvPSI = Baboon_PSI[BaboonNonConsvIDs,]
Lemur_NonConsvPSI =  Lemur_PSI[LemurNonConsvIDs,]
Macaca_NonConsvPSI = Macaca_PSI[MacacaNonConsvIDs,]
Macaque_NonConsvPSI = Macaque_PSI[MacaqueNonConsvIDs,]
Marmoset_NonConsvPSI = Marmoset_PSI[MarmosetNonConsvIDs,]
SquiMonk_NonConsvPSI = SquirrelMonkey_PSI[SquirrelMonkNonConsvIDs,]
Chimp_NonConsvPSI = Chimp_PSI[ChimpNonConsvIDs,]
      

  
#Make a single DF with the filtered consv and non-consv PSI data of each primate

#Joing all HumanIDs
HumanIDs = unique(c(rownames(Human_PSI_filt), rownames(Baboon_ConsvPSI_filt), 
                    rownames(Macaca_ConsvPSI_filt), rownames(Macaque_ConsvPSI_filt),
                    rownames(Marmoset_ConsvPSI_filt), rownames(SquiMonk_ConsvPSI_filt),
                    rownames(Chimp_ConsvPSI_filt), rownames(Lemur_ConsvPSI_filt))) #18,005
#Make PrimatesIDs
PrimatesIDs = c(rownames(Baboon_NonConsvPSI), rownames(Macaca_NonConsvPSI),
                rownames(Macaque_NonConsvPSI), rownames(Chimp_NonConsvPSI),
                rownames(Marmoset_NonConsvPSI), rownames(SquiMonk_NonConsvPSI),
                rownames(Lemur_NonConsvPSI))

numcols = sum(ncol(Baboon_ConsvPSI_filt), ncol(Macaca_ConsvPSI_filt), ncol(Macaque_ConsvPSI_filt),
              ncol(Marmoset_ConsvPSI_filt), ncol(SquiMonk_ConsvPSI_filt), ncol(Chimp_ConsvPSI_filt),
              ncol(Lemur_ConsvPSI_filt), ncol(Human_PSI_filt))

colnamesPSI = c(colnames(Baboon_ConsvPSI_filt), colnames(Macaca_ConsvPSI_filt), colnames(Macaque_ConsvPSI_filt),
                      colnames(Marmoset_ConsvPSI_filt), colnames(SquiMonk_ConsvPSI_filt), colnames(Chimp_ConsvPSI_filt),
                      colnames(Lemur_ConsvPSI_filt), colnames(Human_PSI_filt))

#Remove "_PSI" to all colnamesPSI
colnamesPSI = gsub("_PSI","",colnamesPSI)
#Remove "HumanID"

colnamesPSI = colnamesPSI[!colnamesPSI %in% "HumanID"]


PSI_all = as.data.frame(matrix(nrow = sum(length(HumanIDs), length(PrimatesIDs)) , ncol = length(colnamesPSI) ))
rownames(PSI_all) = c(HumanIDs, PrimatesIDs)
colnames(PSI_all) = colnamesPSI

#Add values to PSI_all

#BAboon
PSI_all[rownames(Baboon_ConsvPSI_filt), grep("Baboon", colnames(PSI_all))] = Baboon_ConsvPSI_filt
PSI_all[rownames(Baboon_NonConsvPSI), grep("Baboon", colnames(PSI_all))] = Baboon_NonConsvPSI


#LEmur
PSI_all[rownames(Lemur_ConsvPSI_filt), grep("Lemur", colnames(PSI_all))] = Lemur_ConsvPSI_filt
PSI_all[rownames(Lemur_NonConsvPSI), grep("Lemur", colnames(PSI_all))] = Lemur_NonConsvPSI

#Macaca
PSI_all[rownames(Macaca_ConsvPSI_filt), grep("Macaca", colnames(PSI_all))] = Macaca_ConsvPSI_filt
PSI_all[rownames(Macaca_NonConsvPSI), grep("Macaca", colnames(PSI_all))] = Macaca_NonConsvPSI

#Macaque
PSI_all[rownames(Macaque_ConsvPSI_filt), grep("Macaque", colnames(PSI_all))] = Macaque_ConsvPSI_filt
PSI_all[rownames(Macaque_NonConsvPSI), grep("Macaque", colnames(PSI_all))] = Macaque_NonConsvPSI

#Marmoset
PSI_all[rownames(Marmoset_ConsvPSI_filt), grep("Marmoset", colnames(PSI_all))] = Marmoset_ConsvPSI_filt
PSI_all[rownames(Marmoset_NonConsvPSI), grep("Marmoset", colnames(PSI_all))] = Marmoset_NonConsvPSI

#Squi Monkey
PSI_all[rownames(SquiMonk_ConsvPSI_filt), grep("SquirrelMonkey", colnames(PSI_all))] = SquiMonk_ConsvPSI_filt
PSI_all[rownames(SquiMonk_NonConsvPSI), grep("SquirrelMonkey", colnames(PSI_all))] = SquiMonk_NonConsvPSI

#Chimp
PSI_all[rownames(Chimp_ConsvPSI_filt), grep("Chimp", colnames(PSI_all))] = Chimp_ConsvPSI_filt
PSI_all[rownames(Chimp_NonConsvPSI), grep("Chimp", colnames(PSI_all))] = Chimp_NonConsvPSI

#Human
PSI_all[rownames(Human_PSI_filt), grep("Human", colnames(PSI_all))] = Human_PSI_filt




#Make NAs as zero
PSI_all[is.na(PSI_all)] <- 0

##Calculate correlation

Corr_PSI = as.data.frame(matrix(ncol=ncol(PSI_all), nrow = ncol(PSI_all)))
colnames(Corr_PSI) = colnames(PSI_all)
rownames(Corr_PSI) = colnames(PSI_all)


for (c in colnames(Corr_PSI)) {
  
  for (d in rownames(Corr_PSI)) {
    
    cor = cor.test(PSI_all[,c], PSI_all[,d],
                   method = "pearson", alternative = "greater")
    
    
    #print(c(c, d))
    
    Corr_PSI[d,c] = cor$estimate
    
  }
  
}


pheatmap(data.matrix(Corr_PSI),
         fontsize_row = 8, fontsize_col = 8,
         cellwidth = 8, cellheight = 8,
         treeheight_row = 0, treeheight_col = 0, main = "circRNAs")



#Make heatmap





