

#Script to retrieve conserved circRNAs between human and primates

#The main aim for this script is to retrieve the conserved circRNAs between primate_A and Human

#Also to get those circRNAs that have conserved exons but are not conserved between human

library(plyr)

#####Read files of LO Primate exons coords that matched with Human Exons coords (of circRNAs (BSJ) both, primate and human)

##PSI_10

Baboon_LO10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Baboon_Human_ExonsCircs_GTF-CE10.txt",
                         header = FALSE, as.is = TRUE) #523,726

Lemur_LO10= read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Lemur_Human_ExonsCircs_GTF-CE10.txt",
                       header = FALSE, as.is = TRUE) #314,959

Macaca_LO10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Macaca_Human_ExonsCircs_GTF-CE10.txt",
                         header = FALSE, as.is = TRUE) #499,709


Macaque_LO10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Macaque_Human_ExonsCircs_GTF-CE10.txt",
                          header = FALSE, as.is = TRUE) #608,872

Marmoset_LO10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Marmoset_Human_ExonsCircs_GTF-CE10.txt",
                           header = FALSE, as.is = TRUE) #68,696

SquiMonkey_LO10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/SquirrelMonkey_Human_ExonsCircs_GTF-CE10.txt",
                             header = FALSE, as.is = TRUE) #309,571

Chimp_LO10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Chimp_Human_ExonsCircs_GTF-CE10.txt",
                        header = FALSE, as.is = TRUE) #594,821


##PSI_5


#Baboon_LO5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_5/GTF-CE/Baboon_Human_ExonsCircs_GTF-CE5.txt",
#                        header = FALSE, as.is = TRUE) #523,813

#Lemur_LO5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_5/GTF-CE/Lemur_Human_ExonsCircs_GTF-CE5.txt",
#                       header = FALSE, as.is = TRUE) #315,005

#Macaca_LO5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_5/GTF-CE/Macaca_Human_ExonsCircs_GTF-CE5.txt",
#                        header = FALSE, as.is = TRUE) #499,861

#Macaque_LO5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_5/GTF-CE/Macaque_Human_ExonsCircs_GTF-CE5.txt",
#                         header = FALSE, as.is = TRUE) #609,044

#Marmoset_LO5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_5/GTF-CE/Marmoset_Human_ExonsCircs_GTF-CE5.txt",
#                          header = FALSE, as.is = TRUE) #68,785

#SquiMonkey_LO5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_5/GTF-CE/SquiMonkey_Human_ExonsCircs_GTF-CE5.txt",
#                            header = FALSE, as.is = TRUE) #462,435

#Chimp_LO5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_5/GTF-CE/Chimp_Human_ExonsCircs_GTF-CE5.txt",
#                       header = FALSE, as.is = TRUE) #595,038

#Make a function to:
###Add colnames to LO files and remove lines that do not have human coords (Exon data)

Filt_LO <- function(Primate_LO_cutoffPSI) {
  
  #Add colnames
  colnames(Primate_LO_cutoffPSI) = c("Chr_PrimateLO", "Start_PrimateLO", "End_PrimateLO", "ID_PrimateExon",
                                     "Chr_Human", "Start_Human", "End_Human", "ID_HumanExon")
  
  #Remove rows that dont have a matching coordinate with human exons
  Primate_LO_cutoffPSI = Primate_LO_cutoffPSI[grep("\\.", Primate_LO_cutoffPSI$ID_HumanExon, invert = TRUE),]
  
  
  ###Remove duplicated rows
  no_dup = which(!duplicated(Primate_LO_cutoffPSI))
  
  Primate_LO_cutoffPSI = Primate_LO_cutoffPSI[no_dup,]
  
  return(Primate_LO_cutoffPSI)
  
}

Baboon_LO10 = Filt_LO(Baboon_LO10) #33,765
Lemur_LO10 = Filt_LO(Lemur_LO10) #24,789
Macaca_LO10 = Filt_LO(Macaca_LO10) #33,803
Macaque_LO10 = Filt_LO(Macaque_LO10) #36,373
Marmoset_LO10 = Filt_LO(Marmoset_LO10) #33,019
SquiMonkey_LO10 = Filt_LO(SquiMonkey_LO10) #31,145
Chimp_LO10 = Filt_LO(Chimp_LO10) #34,586

###Get number of unique exons between primates that mapped to human
ExonsMapped = unique(c(Baboon_LO10$ID_HumanExon, Lemur_LO10$ID_HumanExon, Macaca_LO10$ID_HumanExon, 
                       Macaque_LO10$ID_HumanExon, Marmoset_LO10$ID_HumanExon, SquiMonkey_LO10$ID_HumanExon,
                       Chimp_LO10$ID_HumanExon))


#Baboon_LO5 = Filt_LO(Baboon_LO5) #33,778
#Lemur_LO5 = Filt_LO(Lemur_LO5) #24,791
#Macaca_LO5 = Filt_LO(Macaca_LO5) #33,815
#Macaque_LO5 = Filt_LO(Macaque_LO5) #36,398
#Marmoset_LO5 = Filt_LO(Marmoset_LO5) #554
#SquiMonkey_LO5 = Filt_LO(SquiMonkey_LO5) #31,158
#Chimp_LO5 = Filt_LO(Chimp_LO5) #34,610


##Load objects of CE IDs to removed as they are specific to samples that we removed from the analysis 
#according to weird clustering in GE

#HUMAN
# load(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/HumanRemove_CE10_IDs") #14
#load(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/HumanRemove_CE5_IDs") #9

#MACAQUE
#load(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/MacaqueRemove_CE10_IDs") #1,830
#load(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/MacaqueRemove_CE5_IDs") #1833

####################Remove such IDs

#for Macaque you only need to remove the IDs in corresponding Macaque_LO<PSI>
#Macaque_LO10 = Macaque_LO10[!Macaque_LO10$ID_PrimateExon %in% MacaqueRemove_CE10_IDs,] ##36,373 to 36,313
#Macaque_LO5 = Macaque_LO5[!Macaque_LO5$ID_PrimateExon %in% MacaqueRemove_CE5_IDs,] ##36,398 to 36,335 

#For Human you have to remove the IDs to ALL primates accoridng to their <PSI>

#RemoveHumanCE <-function(Primate_LO, IDs_Human) {
  
  #Primate_LO = Primate_LO[!Primate_LO$ID_HumanExon %in% IDs_Human,]
 
  #print(dim(Primate_LO))
   
  #return(Primate_LO)
  
#}



#Baboon_LO10 = RemoveHumanCE(Baboon_LO10, HumanRemove_CE10_IDs) #33,764
#Lemur_LO10 = RemoveHumanCE(Lemur_LO10, HumanRemove_CE10_IDs) #24,788
#Macaca_LO10 = RemoveHumanCE(Macaca_LO10, HumanRemove_CE10_IDs) #33,802
#Macaque_LO10 = RemoveHumanCE(Macaque_LO10, HumanRemove_CE10_IDs) #36,312
#Marmoset_LO10 = RemoveHumanCE(Marmoset_LO10, HumanRemove_CE10_IDs) #33,019
#SquiMonkey_LO10 = RemoveHumanCE(SquiMonkey_LO10, HumanRemove_CE10_IDs) #31,144
#Chimp_LO10 = RemoveHumanCE(Chimp_LO10, HumanRemove_CE10_IDs) #34,585


#Baboon_LO5 = RemoveHumanCE(Baboon_LO5, HumanRemove_CE5_IDs) #33,778
#Lemur_LO5 = RemoveHumanCE(Lemur_LO5, HumanRemove_CE5_IDs) #24,791
#Macaca_LO5 = RemoveHumanCE(Macaca_LO5, HumanRemove_CE5_IDs) #33,815
#Macaque_LO5 = RemoveHumanCE(Macaque_LO5, HumanRemove_CE5_IDs) #36,398
#Marmoset_LO5 =RemoveHumanCE(Marmoset_LO5, HumanRemove_CE5_IDs) #33,043
#SquiMonkey_LO5 = RemoveHumanCE(SquiMonkey_LO5, HumanRemove_CE5_IDs) #31,158
#Chimp_LO5 = RemoveHumanCE(Chimp_LO5, HumanRemove_CE5_IDs) #34,610


###Read files of Exons Coords (CE) that intersect with circRNAs of THE SAME PRIMATE
##ExonCirc info
#Since 7 Oct the Exon Coords are Whippet CE (merged) intersected with GTF

###PSI 10

Human_Circ10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Human_Int_ExonCirc_GTF-CE10.txt", 
                          header = TRUE, as.is = TRUE) #238,795; 36,739 circRNAs

Baboon_Circ10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Baboon_Int_ExonCirc_GTF-CE10.txt", 
                           header = TRUE, as.is = TRUE) #62,959

Lemur_Circ10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Lemur_Int_ExonCirc_GTF-CE10.txt", 
                          header = TRUE, as.is = TRUE) #43,258

Macaca_Circ10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Macaca_Int_ExonCirc_GTF-CE10.txt", 
                           header = TRUE, as.is = TRUE) #61,910

Macaque_Circ10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Macaque_Int_ExonCirc_GTF-CE10.txt", 
                            header = TRUE, as.is = TRUE) #75,514

Marmoset_Circ10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Marmoset_Int_ExonCirc_GTF-CE10.txt", 
                             header = TRUE, as.is = TRUE) #65,068

SquiMonkey_Circ10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/SquiMonkey_Int_ExonCirc_GTF-CE10.txt", 
                               header = TRUE, as.is = TRUE) #62,005

Chimp_Circ10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Chimp_Int_ExonCirc_GTF-CE10.txt",
                          header = TRUE, as.is = TRUE) #67,844


##PSI 5

#Human_Circ5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Human_Int_ExonCirc_GTF-CE5.txt", 
#                         header = TRUE, as.is = TRUE) #239,129

#Baboon_Circ5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Baboon_Int_ExonCirc_GTF-CE5.txt", 
#                          header = TRUE, as.is = TRUE) #62,982

#Lemur_Circ5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Lemur_Int_ExonCirc_GTF-CE5.txt", 
#                         header = TRUE, as.is = TRUE) #43,280

#Macaca_Circ5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Macaca_Int_ExonCirc_GTF-CE5.txt", 
#                          header = TRUE, as.is = TRUE) #61,957

#Macaque_Circ5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Macaque_Int_ExonCirc_GTF-CE5.txt", 
#                           header = TRUE, as.is = TRUE) #75,587

#Marmoset_Circ5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Marmoset_Int_ExonCirc_GTF-CE5.txt", 
#                            header = TRUE, as.is = TRUE) #65,162

#SquiMonkey_Circ5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/SquiMonkey_Int_ExonCirc_GTF-CE5.txt", 
#                              header = TRUE, as.is = TRUE) #62,047

#Chimp_Circ5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Chimp_Int_ExonCirc_GTF-CE5.txt",
#                         header = TRUE, as.is = TRUE)


####Load IDs of circRNAs/BSJ that are in samples that we removed according to GE data visualization

#load(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/HumanRemove_BSJ_IDs") #51

#load(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/MacaqueRemove_BSJ_IDs") #2,089


##############################################################Remove such IDs
#For macaque IDs only do it in respective PSI macaque data
#Macaque_Circ10 = Macaque_Circ10[!Macaque_Circ10$ID_Circ %in% MacaqueRemove_BSJ_IDs,] #75514 to 65,821
#Macaque_Circ5 = Macaque_Circ5[!Macaque_Circ5$ID_Circ %in% MacaqueRemove_BSJ_IDs,] #75587 to 65,883

#For human IDs only do it in respective PSI human data
#Human_Circ10 = Human_Circ10[!Human_Circ10$ID_Circ %in% HumanRemove_BSJ_IDs,] #238,795 to 238,460
#Human_Circ5 = Human_Circ5[!Human_Circ5$ID_Circ %in% HumanRemove_BSJ_IDs,] #239,129 to 238,794
  


#Make a table with the total number of circRNAs in each primate

AllCirc10 = data.frame("Human"= length(unique(Human_Circ10$ID_Circ)),
                       "Chimp" = length(unique(Chimp_Circ10$ID_Circ)),
                       "Baboon" = length(unique(Baboon_Circ10$ID_Circ)),
                       "Macaque" = length(unique(Macaque_Circ10$ID_Circ)),
                       "Macaca" = length(unique(Macaca_Circ10$ID_Circ)),
                       "Lemur" = length(unique(Lemur_Circ10$ID_Circ)),
                       "Marmoset" = length(unique(Marmoset_Circ10$ID_Circ)) ,
                       "Squirrel_Monkey" = length(unique(SquiMonkey_Circ10$ID_Circ)))


#AllCirc5 = data.frame("Human" = length(unique(Human_Circ5$ID_Circ)),
#                      "Chimp" = length(unique(Chimp_Circ5$ID_Circ)),
#                      "Baboon" = length(unique(Baboon_Circ5$ID_Circ)),
#                      "Macaque" = length(unique(Macaque_Circ5$ID_Circ)),
#                      "Macaca" = length(unique(Macaca_Circ5$ID_Circ)),
#                      "Lemur" = length(unique(Lemur_Circ5$ID_Circ)),
#                      "Marmoset" = length(unique(Marmoset_Circ5$ID_Circ)) ,
#                      "Squirrel_Monkey" = length(unique(SquiMonkey_Circ5$ID_Circ)))



#Save table(s) to later barplot
write.table(AllCirc10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/PrimatesNumOfCircRNAs.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

#write.table(AllCirc5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/PrimatesNumOfCircRNAs.txt",
#            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)






##Make function to order <Primate>_Circ<PSIcutoff> according to ID_Circ
#and then make a column with coordinates such as (Chr_Exon:Start_Exon-End_Exon)

CircOrder_Coord <-function(Primate_Circ_PSIcutoff) {
  
  #order according to ID_Circ
  Primate_Circ_PSIcutoff = Primate_Circ_PSIcutoff[order(Primate_Circ_PSIcutoff$ID_Circ),]
  
  #Add coordinate info
  Primate_Circ_PSIcutoff$Coord_Exon  = paste(Primate_Circ_PSIcutoff$Chr_Exon,
                                             paste(Primate_Circ_PSIcutoff$Start_Exon, Primate_Circ_PSIcutoff$End_Exon, sep = "-"), sep = ":")
  
  return(Primate_Circ_PSIcutoff)
  
}

###############PSI_10

Human_Circ10 = CircOrder_Coord(Human_Circ10)
Baboon_Circ10 = CircOrder_Coord(Baboon_Circ10)
Lemur_Circ10 = CircOrder_Coord(Lemur_Circ10)
Macaca_Circ10 = CircOrder_Coord(Macaca_Circ10)
Macaque_Circ10 = CircOrder_Coord(Macaque_Circ10)
Marmoset_Circ10 = CircOrder_Coord(Marmoset_Circ10)
SquiMonkey_Circ10 = CircOrder_Coord(SquiMonkey_Circ10)
Chimp_Circ10 = CircOrder_Coord(Chimp_Circ10)

##############PSI_5

#Human_Circ5 = CircOrder_Coord(Human_Circ5)
#Baboon_Circ5 = CircOrder_Coord(Baboon_Circ5)
#Lemur_Circ5 = CircOrder_Coord(Lemur_Circ5)
#Macaca_Circ5 = CircOrder_Coord(Macaca_Circ5)
#Macaque_Circ5 = CircOrder_Coord(Macaque_Circ5)
#Marmoset_Circ5 = CircOrder_Coord(Marmoset_Circ5)
#SquiMonkey_Circ5 = CircOrder_Coord(SquiMonkey_Circ5)
#Chimp_Circ5 = CircOrder_Coord(Chimp_Circ5)


###Make a function to count all exons associated to a circRNA (unique exons) (number of exons within the BSJ coordinate)
NumExons <- function(Primate_Circ_PSIcutoff) {
  
  #Get number of exons without counting duplicated exons coordinates
  #by each circRNA
  
  no_dups = which(!duplicated(Primate_Circ_PSIcutoff[,c("ID_Circ", "Coord_Exon")]))
  
  #How many times the ID_Circ is repited (how many exons are associated to that circRNA)
  Num = as.data.frame(table(Primate_Circ_PSIcutoff[no_dups, "ID_Circ"]))
  
  #Get the total of ID_Circ to repeat the number of each ID Circ according to Num result
  
  Tot = as.data.frame(table(Primate_Circ_PSIcutoff$ID_Circ))
  
  Num = Num[rep(row.names(Num), Tot$Freq), 1:2]
  
  Primate_Circ_PSIcutoff$NumExonsPerCirc = Num$Freq
  
  return(Primate_Circ_PSIcutoff)
  
}


#PSI_10

Human_Circ10 = NumExons(Human_Circ10)
Baboon_Circ10 = NumExons(Baboon_Circ10)
Lemur_Circ10 = NumExons(Lemur_Circ10)
Macaca_Circ10 = NumExons(Macaca_Circ10)
Macaque_Circ10 = NumExons(Macaque_Circ10)
Marmoset_Circ10 = NumExons(Marmoset_Circ10) #tot number of exons: 29,125  #Tot number of circRNAs: 11,368
SquiMonkey_Circ10 = NumExons(SquiMonkey_Circ10)
Chimp_Circ10 = NumExons(Chimp_Circ10)

#PSI_5

#Human_Circ5 = NumExons(Human_Circ5)
#Baboon_Circ5 = NumExons(Baboon_Circ5)
#Lemur_Circ5 = NumExons(Lemur_Circ5)
#Macaca_Circ5 = NumExons(Macaca_Circ5)
#Macaque_Circ5 = NumExons(Macaque_Circ5)
#Marmoset_Circ5 = NumExons(Marmoset_Circ5)
#SquiMonkey_Circ5 = NumExons(SquiMonkey_Circ5)
#Chimp_Circ5 = NumExons(Chimp_Circ5)

#Make a function that outputs a dataframe with unique circIDs as one column and the other column the number of exons (for each primate)
FilterNumExonsPerCirc = function(Primate_CircPSI) {
  
  PrimateNumExonsCirc = Primate_CircPSI[!duplicated(Primate_CircPSI[,c("ID_Circ","NumExonsPerCirc")]),c("ID_Circ","NumExonsPerCirc")]
  
  return(PrimateNumExonsCirc)
  
} 

Human_NumExonsPerCirc10 = FilterNumExonsPerCirc(Human_Circ10)
Baboon_NumExonsPerCirc10 = FilterNumExonsPerCirc(Baboon_Circ10)
Lemur_NumExonsPerCirc10 = FilterNumExonsPerCirc(Lemur_Circ10)
Macaca_NumExonsPerCirc10 = FilterNumExonsPerCirc(Macaca_Circ10)
Macaque_NumExonsPerCirc10 = FilterNumExonsPerCirc(Macaque_Circ10)
Marmoset_NumExonsPerCirc10 = FilterNumExonsPerCirc(Marmoset_Circ10)
SquiMonkey_NumExonsPerCirc10 = FilterNumExonsPerCirc(SquiMonkey_Circ10)
Chimp_NumExonsPerCirc10 = FilterNumExonsPerCirc(Chimp_Circ10)


#Human_NumExonsPerCirc5 = FilterNumExonsPerCirc(Human_Circ5)
#Baboon_NumExonsPerCirc5 = FilterNumExonsPerCirc(Baboon_Circ5)
#Lemur_NumExonsPerCirc5 = FilterNumExonsPerCirc(Lemur_Circ5)
#Macaca_NumExonsPerCirc5 = FilterNumExonsPerCirc(Macaca_Circ5)
#Macaque_NumExonsPerCirc5 = FilterNumExonsPerCirc(Macaque_Circ5)
#Marmoset_NumExonsPerCirc5 = FilterNumExonsPerCirc(Marmoset_Circ5)
#SquiMonkey_NumExonsPerCirc5 = FilterNumExonsPerCirc(SquiMonkey_Circ5)
#Chimp_NumExonsPerCirc5 = FilterNumExonsPerCirc(Chimp_Circ5)

#################################################################################################################################

################# Saving files of NumExonsPerCirc outputs #######################################################################

#PSI10
write.table(Human_NumExonsPerCirc10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/HumanCircNumExons10.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Baboon_NumExonsPerCirc10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/BaboonCircNumExons10.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Lemur_NumExonsPerCirc10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/LemurCircNumExons10.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Macaca_NumExonsPerCirc10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/MacacaCircNumExons10.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Macaque_NumExonsPerCirc10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/MacaqueCircNumExons10.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Marmoset_NumExonsPerCirc10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/MarmosetCircNumExons10.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(SquiMonkey_NumExonsPerCirc10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/SquiMonkeyCircNumExons10.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(Chimp_NumExonsPerCirc10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/ChimpCircNumExons10.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


#PSI 5
#write.table(Human_NumExonsPerCirc5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/HumanCircNumExons5.txt",
#            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(Baboon_NumExonsPerCirc5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/BaboonCircNumExons5.txt",
#            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(Lemur_NumExonsPerCirc5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/LemurCircNumExons5.txt",
#            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(Macaca_NumExonsPerCirc5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/MacacaCircNumExons5.txt",
#            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(Macaque_NumExonsPerCirc5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/MacaqueCircNumExons5.txt",
#            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(Marmoset_NumExonsPerCirc5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/MarmosetCircNumExons5.txt",
#            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(SquiMonkey_NumExonsPerCirc5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/SquiMonkeyCircNumExons5.txt",
#            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(Chimp_NumExonsPerCirc5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/ChimpCircNumExons5.txt",
#            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


#################################################################################################################################

#Make List with Primate (Info of Lift Over (orthologue exons)) and Primate_Circ (Info of circRNAs and Exons from the primate) info

List_Baboon10 = list("Primate" = Baboon_LO10, "Primate_Circ" = Baboon_Circ10)
List_Lemur10 = list("Primate" = Lemur_LO10, "Primate_Circ" = Lemur_Circ10)
List_Macaca10 = list("Primate" = Macaca_LO10, "Primate_Circ" = Macaca_Circ10)
List_Macaque10 = list("Primate" = Macaque_LO10, "Primate_Circ" = Macaque_Circ10)
List_Marmoset10 = list("Primate" = Marmoset_LO10, "Primate_Circ" = Marmoset_Circ10)
List_SquiMonkey10 = list("Primate" = SquiMonkey_LO10, "Primate_Circ" = SquiMonkey_Circ10)
List_Chimp10 = list("Primate" = Chimp_LO10, "Primate_Circ" = Chimp_Circ10)

#List_Baboon5 = list("Primate" = Baboon_LO5, "Primate_Circ" = Baboon_Circ5)
#List_Lemur5 = list("Primate" = Lemur_LO5, "Primate_Circ" = Lemur_Circ5)
#List_Macaca5 = list("Primate" = Macaca_LO5, "Primate_Circ" = Macaca_Circ5)
#List_Macaque5 = list("Primate" = Macaque_LO5, "Primate_Circ" = Macaque_Circ5)
#List_Marmoset5 = list("Primate" = Marmoset_LO5, "Primate_Circ" = Marmoset_Circ5)
#List_SquiMonkey5 = list("Primate" = SquiMonkey_LO5, "Primate_Circ" = SquiMonkey_Circ5)
#List_Chimp5 = list("Primate" = Chimp_LO5, "Primate_Circ" = Chimp_Circ5)


#Make a function to get the info of each CircRNA (ID_Circ, Num_Exons, Tot_Length, Coord_Circ, Exons_IDs)

Info_Circs <- function(List_Primate) {
  
  Circs = unique(List_Primate$Primate_Circ$ID_Circ)
  
  #Make a dataframe with all the circRNAs of that primate and columns with above info
  Info_PrimateCircs = as.data.frame(matrix(nrow = length(Circs), ncol = 5))
  colnames(Info_PrimateCircs) = c("ID_Circ", "Num_Exons", "Tot_Length", "Coord_Circ", "Exons_IDs")
  
  Info_PrimateCircs$ID_Circ = Circs
  Info_PrimateCircs$Num_Exons = List_Primate$Primate_Circ[match(Circs, List_Primate$Primate_Circ$ID_Circ), "NumExonsPerCirc"]
  Info_PrimateCircs$Coord_Circ = paste(List_Primate$Primate_Circ[match(Circs, List_Primate$Primate_Circ$ID_Circ), "Chr_Circ"], 
                                       paste(List_Primate$Primate_Circ[match(Circs, List_Primate$Primate_Circ$ID_Circ), "Start_Circ"],
                                             List_Primate$Primate_Circ[match(Circs, List_Primate$Primate_Circ$ID_Circ), "End_Circ"],sep = "-"), sep = ":")
  
  
  #To add Tot_Length and Exons_IDs
  
  ###Total Length is calculated: according to the ID_Circ sum the of the Length_Exon of the exons associated to such ID_Circ
  TEMP = ddply(List_Primate$Primate_Circ,.(ID_Circ), summarize, Tot_Length = sum(Length_Exon), IDs_Exons = paste(ID_Exon, collapse = " | "))
  
  Info_PrimateCircs$Tot_Length = TEMP$Tot_Length
  Info_PrimateCircs$Exons_IDs = TEMP$IDs_Exons
  
  return(Info_PrimateCircs)
  
}

#Below objects have the circRNAID, number of exons in that circRNA, 
#the total length, the coordinate of the circRNA and the ExonsIDs in a single cell of the dataframe

#PSI 10

Info_BaboonCircs10 = Info_Circs(List_Baboon10) #12,513
Info_LemurCircs10 = Info_Circs(List_Lemur10) #8,572
Info_MacaqueCircs10 = Info_Circs(List_Macaque10) #14,122
Info_MacacaCircs10 = Info_Circs(List_Macaca10) #12,439
Info_MarmosetCircs10 = Info_Circs(List_Marmoset10) #11,368
Info_SquiMonkeyCircs10 = Info_Circs(List_SquiMonkey10) #11,986
Info_ChimpCircs10 = Info_Circs(List_Chimp10) #12,641

#PSI 5  
#Info_BaboonCircs5 = Info_Circs(List_Baboon5) #12,514
#Info_LemurCircs5 = Info_Circs(List_Lemur5) #8,572
#Info_MacaqueCircs5 = Info_Circs(List_Macaque5) #14,124
#Info_MacacaCircs5 = Info_Circs(List_Macaca5) #12,439
#Info_MarmosetCircs5 = Info_Circs(List_Marmoset5) #11,370
#Info_SquiMonkeyCircs5 = Info_Circs(List_SquiMonkey5) #11,988
#Info_ChimpCircs5 = Info_Circs(List_Chimp5) #12,642


###Make Info_HumanCircs<PSI>
Info_HumanCirc10 = ddply(Human_Circ10,.(ID_Circ, NumExonsPerCirc), summarize, 
                         Tot_Length = sum(Length_Exon), IDs_Exons = paste(ID_Exon, collapse = " | "))

Info_HumanCirc10$Coord_Circ = paste(Human_Circ10[match(Info_HumanCirc10$ID_Circ, Human_Circ10$ID_Circ), "Chr_Circ"], 
                                    paste(Human_Circ10[match(Info_HumanCirc10$ID_Circ, Human_Circ10$ID_Circ), "Start_Circ"], 
                                          Human_Circ10[match(Info_HumanCirc10$ID_Circ, Human_Circ10$ID_Circ), "End_Circ"], sep = "-") , sep=":") #36,739

#Info_HumanCirc5 = ddply(Human_Circ5,.(ID_Circ, NumExonsPerCirc), summarize, 
#                        Tot_Length = sum(Length_Exon), IDs_Exons = paste(ID_Exon, collapse = " | "))

#Info_HumanCirc5$Coord_Circ = paste(Human_Circ5[match(Info_HumanCirc5$ID_Circ, Human_Circ5$ID_Circ), "Chr_Circ"], 
#                                   paste(Human_Circ5[match(Info_HumanCirc5$ID_Circ, Human_Circ5$ID_Circ), "Start_Circ"], 
#                                         Human_Circ5[match(Info_HumanCirc5$ID_Circ, Human_Circ5$ID_Circ), "End_Circ"], sep = "-") , sep=":")


#The Info_<Primate>Circ<PSI> has the same number of circRNAs as the dataframe AllCirc<PSI>. Therefore at this point
#We dont save the data


##Function to Add info (length of LO exons) to LO data (List_Primate$Primate) 

#Make another dataframe with the info of List_Primate$Primate_Circ where the
#exons are not LiftOver(LO). Remove this from the List_Primate$Primate_Circ the circID which exones are not LiftOver (LO)

#I added the parameter: Output
#This parameter is to chose if I want as an output the #Filtered List_Primate or 
#if I want the number of circRNAs
#That I'm losing in each step of the filtering



Add_Info2LO <- function(List_Primate, Output) {
  
  
  #Make a list of potential
  FilterOutput = list("Total" = c(),
                      "LO_cutoff" = c(),
                      "Length_cutoff" = c())
  
  FilterOutput$Total = unique(List_Primate$Primate_Circ$ID_Circ)
  
  #Calculate length of exons in List_Primate<PSI>$Primate (LO exon length)
  List_Primate$Primate$LengthExon_LO = List_Primate$Primate$End_PrimateLO - List_Primate$Primate$Start_PrimateLO
  
  #Calculate length of Human exons in List_Primate<PSI>$Primate
  List_Primate$Primate$LengthExon_Human = List_Primate$Primate$End_Human - List_Primate$Primate$Start_Human
  
  ####Exons ID with LiftOver (are in primate and in human) and in circRNA (are in a circRNA)
  
  print("Number of Exons that are LiftedOver (orthologue to Human)")
  print(length(unique(List_Primate$Primate$ID_PrimateExon)))
  
  print("Number of Exons associated to a circRNA in the Primate")
  print(length(unique(List_Primate$Primate_Circ$ID_Exon)))
  
  
  commonIDexon = intersect(List_Primate$Primate$ID_PrimateExon, List_Primate$Primate_Circ$ID_Exon)
  
  print("Number of exons that are LiftedOver (orthologue to Human) and are associated to a circRNA in the Primate")
  print(length(commonIDexon))
  
  #Marmoset_PSI10 case: you start with 404 commonIDexons
  #later when you filter the NotLO (the ones that ARE NOT Lifted Over)
  #and you removed all the circRNAs that are formed by NotLO and make the
  #NewCommon then you are left with only 50 exons, from which later their length wont match
  
  #Exons that form a circRNA in the primate but are not LO
  NotLO = setdiff(List_Primate$Primate_Circ$ID_Exon, commonIDexon)
  
  ################Esto es nuevo
  
  #Using the NotLO subset the List_Primate$Primate_Circ 
  circRNAsWithNotLO = List_Primate$Primate_Circ[List_Primate$Primate_Circ$ID_Exon %in% NotLO,]
  
  
  #For each circRNA that is associated to a NotLO check if the circRNA is made of ONLY NotLO
  #For this you add a column to circRNAsWithNotLO names (NumExonsNotLO)
  
  circRNAsWithNotLO$NumExonsNotLO  = NumExons(circRNAsWithNotLO)[,"NumExonsPerCirc"]
  
  #If its made of ONLY NotLO Remove it from List_Primate$Primate_Circ and keep it in another object (TO BE DEFINED)
  CircRNAs_OnlyWithNotLO = circRNAsWithNotLO[which(circRNAsWithNotLO$NumExonsPerCirc == circRNAsWithNotLO$NumExonsNotLO),]
  
  #ID_Circ of  circRNAsWithNotLO but that ARE NOT in CircRNAs_OnlyWithNotLO
  circRNAsIDs_SomeWithLO = setdiff(circRNAsWithNotLO$ID_Circ, CircRNAs_OnlyWithNotLO$ID_Circ)
  
  #Keep the circRNAsIDs_SomeWithLO in the circRNAsWithNotLO (e.i: The dataframe now only has the circRNAs info of circRNAs that have at least an exon 
  #that is not LO but is made of some other exons that are LO)
  circRNAsWithNotLO = circRNAsWithNotLO[circRNAsWithNotLO$ID_Circ %in% circRNAsIDs_SomeWithLO,]
  
  
  print("Number of total circRNAs")
  print(length(unique(List_Primate$Primate_Circ$ID_Circ)))
  
  print("Number of circRNAs that at least have an exon that is NOT ORTHOLOGUE (Lifted Over)")
  print(length(unique(List_Primate$Primate_Circ[List_Primate$Primate_Circ$ID_Exon %in% NotLO,"ID_Circ"])))
  
  #Get the proportion of circRNAs that have NotLO and LO and those that only have LO exons
  
  print("Number of circRNAs that only have NOT ORTHOLOGUE exons")
  print(length(unique(CircRNAs_OnlyWithNotLO$ID_Circ)))
  
  print("Number of circRNAs that at least have ONE NOT ORTHOLOGUE exon")
  print(length(unique(circRNAsWithNotLO$ID_Circ)))
  
  
  #Remove from List_Primate$Primate_Circ ID_Circ that are in CircRNAs_OnlyWithNotLO 
  List_Primate$Primate_Circ = List_Primate$Primate_Circ[!List_Primate$Primate_Circ$ID_Circ %in% CircRNAs_OnlyWithNotLO$ID_Circ,]
  
  
  ###ADD: Number of circRNAs (circIDs) that are NOT in CircRNAs_OnlyWithNotLO 
  FilterOutput$LO_cutoff = unique(List_Primate$Primate_Circ$ID_Circ)
  
  ####################################################################################################################################
  ##############################IMPORTANTE##################################################################################
  
  ####Exon Length filtering##### No usaste este filtro el 3 de marzo de 2021
  #Este dia estabas pensando en ya publicar tus resultados, List_Primate fue probado con
  #los datos de Macaque10 y despues de remover los circRNAs q solo estaban formados por exones que no tienen LiftOVER
  #tenias 12,799 de un total de 14,122.
  #Cuando probaste ese dataset con el filtro del 90% de la longitud entre exones the quedaste con
  #4,414 (muy pocos)
  
  #Solo agreastr la información (en el objeto FilterOutput$Length_cutoff) de cuántos circRNAs quedarian si proseguias con el filtrado 
  
  
  #Calculate % of length of LengthExon_LO if LengthExon_Human is 100%
  List_Primate$Primate$PercLength_ExonLO = (List_Primate$Primate$LengthExon_LO *100)/List_Primate$Primate$LengthExon_Human
  
  #Calculate % of length of LengthExon_Humanif LengthExon_LO is 100%
  List_Primate$Primate$PercLength_ExonHuman = (List_Primate$Primate$LengthExon_Human *100)/List_Primate$Primate$LengthExon_LO
  
  
  #Filter those that are less than 90% PercLength_ExonHuman or PercLength_ExonLO
  rows_90 = List_Primate$Primate$PercLength_ExonLO >= 90 & List_Primate$Primate$PercLength_ExonHuman >= 90
  
  #Which rows have less than 90% of Length of ExonLO and ExonHuman
  rows_rm = !(List_Primate$Primate$PercLength_ExonLO >= 90 & List_Primate$Primate$PercLength_ExonHuman >= 90)
  
  Exons_rm = unique(List_Primate$Primate[rows_rm,"ID_PrimateExon"])
  
  #Keep exons that fullfill the length % cutoff
  #List_Primate$Primate = List_Primate$Primate[rows_90,]
  
  #Remove those ID_PrimateExon that have different length between LengthExon_LO and LengthExon_Human
  #List_Primate$Primate = List_Primate$Primate[List_Primate$Primate$LengthExon_LO == List_Primate$Primate$LengthExon_Human,]
  
  #Add another element to List_Primate (all circRNAs that have at least one exon that didn't pass the length cutoff)
  rmCircs = unique(List_Primate$Primate_Circ[List_Primate$Primate_Circ$ID_Exon %in% Exons_rm,"ID_Circ"])
  
  #List_Primate$Primate_CircNotLength = List_Primate$Primate_Circ[List_Primate$Primate_Circ$ID_Circ %in% rmCircs,]
  
  #Get the CircRNAs that are not in the rmCircs
  #List_Primate$Primate_Circ = List_Primate$Primate_Circ[!List_Primate$Primate_Circ$ID_Circ %in% rmCircs,]
  
  FilterOutput$Length_cutoff = unique(List_Primate$Primate_Circ[!List_Primate$Primate_Circ$ID_Circ %in% rmCircs,"ID_Circ"])
  
  #Get same exons in List_Primate$Primate_Circ
  #List_Primate$Primate = List_Primate$Primate[List_Primate$Primate$ID_PrimateExon %in% List_Primate$Primate_Circ$ID_Exon,]
  
  #There could be the case that there are duplicated rows (exact info in the rows)
  #Remove such rows
  #List_Primate$Primate = unique(List_Primate$Primate)
  
  print("Number of circRNAs after removing circRNAs that are only made of NO orthologue Exons (No LiftedOver)")
  print(length(unique(List_Primate$Primate_Circ$ID_Circ)))
  
  if (Output == "List") {
    
    return(List_Primate)
    
  } 
  
  if (Output == "Cutoffs" ) {
    
    return(FilterOutput)
    
  }
  
}


#Add_Info2LO function removes from List_<Primate><PSI> $ Primate_Circ those
#circRNAs that are formed by at least one exon that is not LiftedOver, these are saved
#in the List_<Primate><PSI> $ Primate_CircNotLO element of the list

#Also it Remove those ID_PrimateExon that have different length between LengthExon_LO and LengthExon_Human,
#Then also remove the circRNas that do not pass the length of exons cutoff in the List_Primate$Primate 
#those circRNas that did not pass the cutoff of exon length are saved in a new list element named List_Primate$Primate_CircNotLength

#####PSI_10
List_Baboon10_filt = Add_Info2LO(List_Baboon10, "List")
List_Lemur10_filt = Add_Info2LO(List_Lemur10, "List")
List_Macaca10_filt = Add_Info2LO(List_Macaca10, "List")
List_Macaque10_filt = Add_Info2LO(List_Macaque10, "List")

List_Marmoset10_filt = Add_Info2LO(List_Marmoset10, "List")
List_SquiMonkey10_filt = Add_Info2LO(List_SquiMonkey10, "List")
List_Chimp10_filt = Add_Info2LO(List_Chimp10, "List") ###ALGO ESTA PASANDO

##Loss by cutoff
Baboon10_cutoff = Add_Info2LO(List_Baboon10, "Cutoffs")
Lemur10_cutoff = Add_Info2LO(List_Lemur10, "Cutoffs")
Macaca10_cutoff = Add_Info2LO(List_Macaca10, "Cutoffs")
Macaque10_cutoff = Add_Info2LO(List_Macaque10, "Cutoffs")
Marmoset10_cutoff = Add_Info2LO(List_Marmoset10, "Cutoffs")
SquiMonkey10_cutoff = Add_Info2LO(List_SquiMonkey10, "Cutoffs")
Chimp10_cutoff = Add_Info2LO(List_Chimp10, "Cutoffs")


###Make table of how many circRNAs (in each primate I lose after each cutoff)
##PSI 10
UniqueCircRNAs_AfterLOcutoff_10 = data.frame("Primate" = c("Macaque", "Chimp", 
                                                           "Macaca", "Marmoset", 
                                                           "Lemur", "Squirrel Monkey","Baboon"),
                                             "Total circRNAs" = c(length(Macaque10_cutoff$LO_cutoff), length(Chimp10_cutoff$LO_cutoff),
                                                                  length(Macaca10_cutoff$LO_cutoff), length(Marmoset10_cutoff$LO_cutoff),
                                                                  length(Lemur10_cutoff$LO_cutoff), length(SquiMonkey10_cutoff$LO_cutoff),
                                                                  length(Baboon10_cutoff$LO_cutoff)))


write.table(UniqueCircRNAs_AfterLOcutoff_10, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/BS/NumberOfCircRNAs/TotalNumberUniqueCircRNAs_LOcutoffCE10.txt",
col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")


##PSI_5
#List_Baboon5_filt = Add_Info2LO(List_Baboon5, "List")
#List_Lemur5_filt = Add_Info2LO(List_Lemur5, "List")
#List_Macaca5_filt = Add_Info2LO(List_Macaca5, "List")
#List_Macaque5_filt = Add_Info2LO(List_Macaque5, "List")
#List_Marmoset5_filt = Add_Info2LO(List_Marmoset5, "List")
#List_SquiMonkey5_filt = Add_Info2LO(List_SquiMonkey5, "List")
#List_Chimp5_filt = Add_Info2LO(List_Chimp5, "List")

##Loss by cutoff
#Baboon5_cutoff = Add_Info2LO(List_Baboon5, "Cutoffs")
#Lemur5_cutoff = Add_Info2LO(List_Lemur5, "Cutoffs")
#Macaca5_cutoff = Add_Info2LO(List_Macaca5, "Cutoffs")
#Macaque5_cutoff = Add_Info2LO(List_Macaque5, "Cutoffs")
#Marmoset5_cutoff = Add_Info2LO(List_Marmoset5, "Cutoffs")
#SquiMonkey5_cutoff = Add_Info2LO(List_SquiMonkey5, "Cutoffs")
#Chimp5_cutoff = Add_Info2LO(List_Chimp5, "Cutoffs")

#Table
#UniqueCircRNAs_AfterLOcutoff_5 = data.frame("Primate" = c("Macaque", "Chimp", 
#                                                          "Macaca", "Marmoset", 
#                                                          "Lemur", "Squirrel Monkey","Baboon"),
#                                            "Total circRNAs" = c(length(Macaque5_cutoff$LO_cutoff), length(Chimp5_cutoff$LO_cutoff),
#                                                                 length(Macaca5_cutoff$LO_cutoff), length(Marmoset5_cutoff$LO_cutoff),
#                                                                 length(Lemur5_cutoff$LO_cutoff), length(SquiMonkey5_cutoff$LO_cutoff),
#                                                                 length(Baboon5_cutoff$LO_cutoff)))


###Save data 
#write.table(UniqueCircRNAs_AfterLOcutoff_5, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/BS/NumberOfCircRNAs/TotalNumberUniqueCircRNAs_LOcutoffCE5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")



##Make table of how many circRNAs (in each primate I lose after each cutoff)
#After Length and LO cutoff
UniqueCircRNAs_AfterLengthCutoff_10 = data.frame("Primate" = c("Macaque", "Chimp", 
                                                               "Macaca", "Marmoset", 
                                                               "Lemur", "Squirrel Monkey","Baboon"),
                                                 "Total circRNAs" = c(length(Macaque10_cutoff$Length_cutoff), length(Chimp10_cutoff$Length_cutoff),
                                                                      length(Macaca10_cutoff$Length_cutoff), length(Marmoset10_cutoff$Length_cutoff),
                                                                      length(Lemur10_cutoff$Length_cutoff), length(SquiMonkey10_cutoff$Length_cutoff),
                                                                      length(Baboon10_cutoff$Length_cutoff)))

write.table(UniqueCircRNAs_AfterLengthCutoff_10, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/BS/NumberOfCircRNAs/TotalNumberUniqueCircRNAs_LengthCutoffCE10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

#UniqueCircRNAs_AfterLengthCutoff_5 = data.frame("Primate" = c("Macaque", "Chimp", 
#                                                              "Macaca", "Marmoset", 
#                                                              "Lemur", "Squirrel Monkey","Baboon"),
#                                                "Total circRNAs" = c(length(Macaque5_cutoff$Length_cutoff), length(Chimp5_cutoff$Length_cutoff),
 #                                                                    length(Macaca5_cutoff$Length_cutoff), length(Marmoset5_cutoff$Length_cutoff),
#                                                                     length(Lemur5_cutoff$Length_cutoff), length(SquiMonkey5_cutoff$Length_cutoff),
#                                                                     length(Baboon5_cutoff$Length_cutoff)))


###Save data 
#write.table(UniqueCircRNAs_AfterLengthCutoff_5, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/BS/NumberOfCircRNAs/TotalNumberUniqueCircRNAs_LengthCutoffCE5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")


#Define which exonID is the ones that are part of the BSJ (this will later help us to defined conserved circRNAs according to BSJ)

#Remember that List_Primate_filt$Primate_Circ has LiftedOver exons and NonLiftedOverExons per circRNA

DefineExonsInBSJ <-function(List_Primate_filt) {
  
  #Get ID_Exon that Start coord matches with Start_Circ Coord
  List_Primate_filt$Primate_Circ$StartOrEndExon_BSJ = NA
  
  starts = which(List_Primate_filt$Primate_Circ$Start_Exon == List_Primate_filt$Primate_Circ$Start_Circ)
  ends = which(List_Primate_filt$Primate_Circ$End_Exon == List_Primate_filt$Primate_Circ$End_Circ)
  
  starts_end = intersect(starts, ends)
  
  List_Primate_filt$Primate_Circ[starts, "StartOrEndExon_BSJ"] = rep("Start", length(starts))
  List_Primate_filt$Primate_Circ[ends, "StartOrEndExon_BSJ"] = rep("End", length(ends))
  
  List_Primate_filt$Primate_Circ[starts_end, "StartOrEndExon_BSJ"] = rep("Start_End", length(starts_end))
  
  
  #Which circRNAs do not have Start Coord
  IDsCoords = unique(List_Primate_filt$Primate_Circ[complete.cases(List_Primate_filt$Primate_Circ$StartOrEndExon_BSJ),"ID_Circ"])
  IDsNoCoords = setdiff(List_Primate_filt$Primate_Circ$ID_Circ, IDsCoords)
  
  TempNoCoord = List_Primate_filt$Primate_Circ[List_Primate_filt$Primate_Circ$ID_Circ %in% IDsNoCoords,]
  
  #Calculate the difference between Start_Exon and Start_Circ, and End_Circ and End_Exon
  TempNoCoord$DiffStart = abs(TempNoCoord$Start_Exon - TempNoCoord$Start_Circ)
  TempNoCoord$DiffEnd = abs(TempNoCoord$End_Exon - TempNoCoord$End_Circ)
  
  #which Start_Diff is <= 100
  starts_nocoord = which(TempNoCoord$DiffStart <= 100)
  
  #which End_Diff is <= 100
  end_nocoord = which(TempNoCoord$DiffEnd <= 100)
  
  #which Start_Diff And End_Diff <= 100 (should be CircRNAs made of only one Exon)
  starts_end_nocoord = intersect(starts_nocoord, end_nocoord)
  
  ##Add to TempNoCoord Start and End
  TempNoCoord[starts_nocoord, "StartOrEndExon_BSJ"] = "Start"
  TempNoCoord[end_nocoord, "StartOrEndExon_BSJ"] = "End"
  TempNoCoord[starts_end_nocoord, "StartOrEndExon_BSJ"] = "Start_End"
  
  ###From TempNoCoord Get ID_Circ that dont Start or End defined
  NoStartCirc = setdiff(TempNoCoord$ID_Circ, TempNoCoord[complete.cases(TempNoCoord$StartOrEndExon_BSJ),"ID_Circ"])
  
  #Remove from TempNoCoord the NoStartCirc
  TempNoCoord = TempNoCoord[!TempNoCoord$ID_Circ %in% NoStartCirc,]
  
  
  ##From List_Primate_filt$Primate_Circ remove those circRNAs used to make TempNoCoord
  
  List_Primate_filt$Primate_Circ = List_Primate_filt$Primate_Circ[!List_Primate_filt$Primate_Circ$ID_Circ %in% IDsNoCoords,]
  
  
  ###Add to List_Primate_filt the TempNoCoord (now has Start and End Coord)
  List_Primate_filt$Primate_Circ = rbind(List_Primate_filt$Primate_Circ, TempNoCoord[,grep("Diff", colnames(TempNoCoord), invert = TRUE)])
  
  
  ###Filter List_PRimate_filt$Primate keeping the ID_PrimateExon that are in List_Primate_filt$Primate_Circ$ID_Exon
  List_Primate_filt$Primate = List_Primate_filt$Primate[List_Primate_filt$Primate$ID_PrimateExon %in% List_Primate_filt$Primate_Circ$ID_Exon,]
  
  
  
  return(List_Primate_filt)
  
}


#########ADD Start and End Exons of BSJ
#####PSI_10
List_Baboon10_filt = DefineExonsInBSJ(List_Baboon10_filt) #circRNAs:6791
List_Lemur10_filt = DefineExonsInBSJ(List_Lemur10_filt) # 4642
List_Macaca10_filt = DefineExonsInBSJ(List_Macaca10_filt) #6877
List_Macaque10_filt = DefineExonsInBSJ(List_Macaque10_filt) #7645
List_Marmoset10_filt = DefineExonsInBSJ(List_Marmoset10_filt) #6425
List_SquiMonkey10_filt = DefineExonsInBSJ(List_SquiMonkey10_filt) #6285
List_Chimp10_filt = DefineExonsInBSJ(List_Chimp10_filt) #7051


###PSI_5
#List_Baboon5_filt = DefineExonsInBSJ(List_Baboon5_filt)
#List_Lemur5_filt = DefineExonsInBSJ(List_Lemur5_filt)
#List_Macaca5_filt = DefineExonsInBSJ(List_Macaca5_filt)
#List_Macaque5_filt = DefineExonsInBSJ(List_Macaque5_filt)
#List_Marmoset5_filt = DefineExonsInBSJ(List_Marmoset5_filt)
#List_SquiMonkey5_filt = DefineExonsInBSJ(List_SquiMonkey5_filt)
#List_Chimp5_filt = DefineExonsInBSJ(List_Chimp5_filt)

####Remove circRNAs in List_PrimatPSI_filt which Start or End do not have a LiftedOver Coord in Human

#Make a function that checks if ID_Exon with Start,End or Start_End have a human orthologue (those exons are conserved)
#If either the ID_Exon with Start,End or Start_End dont share a human orthologue, remove the comple circRNA (BSJ is not conserved)

RemoveNotConservedBSJ <-function(List_Primate_filt) {
  
  #Get the Exons that are the Start,End or Start_End of the circRNA (BSJ)
  ExonsInBSJ = unique(List_Primate_filt$Primate_Circ[complete.cases(List_Primate_filt$Primate_Circ$StartOrEndExon_BSJ),"ID_Exon"])
  
  #Get the Exons that are the Start,End or Start_End of the circRNA (BSJ) but dont have a human orthologue
  NotConsv_ExonsInBSJ = setdiff(ExonsInBSJ, List_Primate_filt$Primate$ID_PrimateExon)
  
  #Get the circRNAsIDs that are the  Start,End or Start_End of the circRNA (BSJ) but are NotConsv_ExonsInBSJ
  TempCircRNAs_BSJinfo = List_Primate_filt$Primate_Circ[complete.cases(List_Primate_filt$Primate_Circ$StartOrEndExon_BSJ),]
  
  #From the TempCircRNAs_BSJinfo get the ID_Circ of the NotConsv_ExonsInBSJ
  NotConsv_CircRNAs = unique(TempCircRNAs_BSJinfo[TempCircRNAs_BSJinfo$ID_Exon %in% NotConsv_ExonsInBSJ, "ID_Circ"])
  
  #Remove the NotConsv_CircRNAs from List_Primate_filt$Primate_Circ
  List_Primate_filt$Primate_Circ = List_Primate_filt$Primate_Circ[!List_Primate_filt$Primate_Circ$ID_Circ %in% NotConsv_CircRNAs,]
  
  
  return(List_Primate_filt)
  
}


#########Remove Not conserved BSJ
#####PSI_10
List_Baboon10_filt = RemoveNotConservedBSJ(List_Baboon10_filt) #circRNAs:5310
List_Lemur10_filt = RemoveNotConservedBSJ(List_Lemur10_filt) # 3323
List_Macaca10_filt = RemoveNotConservedBSJ(List_Macaca10_filt) #5296
List_Macaque10_filt = RemoveNotConservedBSJ(List_Macaque10_filt) #5953
List_Marmoset10_filt = RemoveNotConservedBSJ(List_Marmoset10_filt) #4857
List_SquiMonkey10_filt = RemoveNotConservedBSJ(List_SquiMonkey10_filt) # 4374
List_Chimp10_filt = RemoveNotConservedBSJ(List_Chimp10_filt) #5819


###PSI_5
#List_Baboon5_filt = RemoveNotConservedBSJ(List_Baboon5_filt)
#List_Lemur5_filt = RemoveNotConservedBSJ(List_Lemur5_filt)
#List_Macaca5_filt = RemoveNotConservedBSJ(List_Macaca5_filt)
#List_Macaque5_filt = RemoveNotConservedBSJ(List_Macaque5_filt)
#List_Marmoset5_filt =RemoveNotConservedBSJ(List_Marmoset5_filt)
#List_SquiMonkey5_filt = RemoveNotConservedBSJ(List_SquiMonkey5_filt)
#List_Chimp5_filt = RemoveNotConservedBSJ(List_Chimp5_filt)



#The Info_<Primate><PSI> will only be filtered according to the relaxed cutoff of
#Lifted Over exons, where I only removed the circRNAs of the primates that were ONLY FORMED by NON ORTOLOGUE (no Lifted Over) exons
#(e.i, we only kept the circRNAs that are made of at least one exon that is LiftedOver)
#And that the BSJ is made of conserved Exons
#Filter Info_<Primate>Circs according to LiftOver results of List_Primate$Primate_Circ
#

##PSI 10

#Number of circRNAs according to 3 of march 2021 filters 

Info_BaboonCircs10_filt = Info_BaboonCircs10[Info_BaboonCircs10$ID_Circ %in% List_Baboon10_filt$Primate_Circ$ID_Circ,] 
Info_LemurCircs10_filt  = Info_LemurCircs10[Info_LemurCircs10$ID_Circ %in% List_Lemur10_filt$Primate_Circ$ID_Circ,] #
Info_MacaqueCircs10_filt = Info_MacaqueCircs10[Info_MacaqueCircs10$ID_Circ %in% List_Macaque10_filt$Primate_Circ$ID_Circ,] #
Info_MacacaCircs10_filt = Info_MacacaCircs10[Info_MacacaCircs10$ID_Circ %in% List_Macaca10_filt$Primate_Circ$ID_Circ,] #
Info_MarmosetCircs10_filt = Info_MarmosetCircs10[Info_MarmosetCircs10$ID_Circ %in% List_Marmoset10_filt$Primate_Circ$ID_Circ,] #
Info_SquiMonkeyCircs10_filt = Info_SquiMonkeyCircs10[Info_SquiMonkeyCircs10$ID_Circ %in% List_SquiMonkey10_filt$Primate_Circ$ID_Circ,] #
Info_ChimpCircs10_filt = Info_ChimpCircs10[Info_ChimpCircs10$ID_Circ %in% List_Chimp10_filt$Primate_Circ$ID_Circ,] #


## PSI 5
#Info_BaboonCircs5_filt = Info_BaboonCircs5[Info_BaboonCircs5$ID_Circ %in% List_Baboon5_filt$Primate_Circ$ID_Circ,]
#Info_LemurCircs5_filt = Info_LemurCircs5[Info_LemurCircs5$ID_Circ %in% List_Lemur5_filt$Primate_Circ$ID_Circ,]
#Info_MacaqueCircs5_filt = Info_MacaqueCircs5[Info_MacaqueCircs5$ID_Circ %in% List_Macaque5_filt$Primate_Circ$ID_Circ,]
#Info_MacacaCircs5_filt = Info_MacacaCircs5[Info_MacacaCircs5$ID_Circ %in% List_Macaca5_filt$Primate_Circ$ID_Circ,]
#Info_MarmosetCircs5_filt = Info_MarmosetCircs5[Info_MarmosetCircs5$ID_Circ %in% List_Marmoset5_filt$Primate_Circ$ID_Circ,] 
#Info_SquiMonkeyCircs5_filt = Info_SquiMonkeyCircs5[Info_SquiMonkeyCircs5$ID_Circ %in% List_SquiMonkey5_filt$Primate_Circ$ID_Circ,]
#Info_ChimpCircs5_filt = Info_ChimpCircs5[Info_ChimpCircs5$ID_Circ %in% List_Chimp5_filt$Primate_Circ$ID_Circ,]

  
#The objects Info_<Primate>Circs<PSI>_filt have the circRNA info (circRNAID and Exons_IDs OF PRIMATE data)

#With the below function we make an Info_PrimatePSI dataframe but with LiftedOver data
#Meaning the same ID_Circ_Primate as in Info_<Primate>Circs<PSI> will be in the resulted dataframe
#but with info according to LiftOver results (human coords and human exons IDs)

###Function to make a Info_PrimateCirc<PSI> type DF but with LO data: ID_Circ_Primate, ID_HumanExon (matched according to ID_PrimateExon), 
#Tot_Length_LO (According to each length of exons from LengthExon_LO in the $Primate element from List_Primate)


Make_InfoPrimateLO <- function(Info_PrimatePSI_filt, List_Primate_filt) {
  
  #Defined the final DF that you want to output:
  #The columns we want: ID_Circ_Primate, NumExons_LO, Tot_Length_LO, HumanExons_IDs
  Info_PrimateLO = as.data.frame(matrix(nrow = nrow(Info_PrimatePSI_filt) , ncol = 9))
  colnames(Info_PrimateLO) = c("ID_Circ_Primate", "NumExons_LO", "Tot_Length_LO", 
                               "Tot_Length_Human","HumanExons_IDs",
                               "StartExonCircRNA_Primate", "StartExonCircRNA_Primate_LO",
                               "EndExonCircRNA_Primate", "EndExonCircRNA_Primate_LO")
  
  #Add info we already have
  Info_PrimateLO$ID_Circ_Primate = Info_PrimatePSI_filt$ID_Circ
  
  
  #Make a list with ID_Circ and identifiers and Exons _IDs (as vector) in each element of the list
  TempList_CircExons = split(List_Primate_filt$Primate_Circ$ID_Exon, List_Primate_filt$Primate_Circ$ID_Circ)
  #The List is: name of element in List (circRNAID), values in name are the Primate_ExonID
  
  
  #Make a loop? to: 1) take the ID_HumanExon according to the Exons_IDs from above list, 
  #2) calculate the total length (Total length of the circRNA) according to the LengthExon_LO
  
  for (circ in names(TempList_CircExons)) {
    
    #print(circ)
    
    #Take Exons IDs 
    IDs = TempList_CircExons[[circ]]
    
    Info_PrimateLO[match(circ, Info_PrimateLO$ID_Circ_Primate), "NumExons_LO"] = length(List_Primate_filt$Primate[match(IDs, List_Primate_filt$Primate$ID_PrimateExon, nomatch = FALSE), 
                                                                                                                  "ID_HumanExon"])
    
    Info_PrimateLO[match(circ, Info_PrimateLO$ID_Circ_Primate),"Tot_Length_LO"] = sum(List_Primate_filt$Primate[match(IDs, List_Primate_filt$Primate$ID_PrimateExon, nomatch = FALSE), 
                                                                                                                "LengthExon_LO"])
    
    Info_PrimateLO[match(circ, Info_PrimateLO$ID_Circ_Primate),"Tot_Length_Human"] = sum(List_Primate_filt$Primate[match(IDs, List_Primate_filt$Primate$ID_PrimateExon, nomatch = FALSE), 
                                                                                                                   "LengthExon_Human"])
    
    Info_PrimateLO[match(circ, Info_PrimateLO$ID_Circ_Primate),"HumanExons_IDs"] = paste(List_Primate_filt$Primate[match(IDs, List_Primate_filt$Primate$ID_PrimateExon, nomatch = FALSE), 
                                                                                                                   "ID_HumanExon"], collapse = " | ")
    
  } 
  
  #Add to Info_PrimateLO the ID_Exon that starts the BSJ
  #Subset List_Primate_filt$Primate_Circ with Start info
  StartTemp = List_Primate_filt$Primate_Circ[grep("Start", List_Primate_filt$Primate_Circ$StartOrEndExon_BSJ),]
  Info_PrimateLO[,"StartExonCircRNA_Primate"] = StartTemp[match(Info_PrimateLO$ID_Circ_Primate, StartTemp$ID_Circ),"ID_Exon"]
  
  EndTemp = List_Primate_filt$Primate_Circ[grep("End", List_Primate_filt$Primate_Circ$StartOrEndExon_BSJ),]
  Info_PrimateLO[,"EndExonCircRNA_Primate"] = EndTemp[match(Info_PrimateLO$ID_Circ_Primate, EndTemp$ID_Circ),"ID_Exon"]
  
  
  #Add orthologue Exon ID for Stat and End Exon
  Info_PrimateLO[, "StartExonCircRNA_Primate_LO"] = List_Primate_filt$Primate[match(Info_PrimateLO$StartExonCircRNA_Primate, List_Primate_filt$Primate$ID_PrimateExon),
                                                                              "ID_HumanExon"]
  
  Info_PrimateLO[, "EndExonCircRNA_Primate_LO"] =  List_Primate_filt$Primate[match(Info_PrimateLO$EndExonCircRNA_Primate, List_Primate_filt$Primate$ID_PrimateExon),
                                                                             "ID_HumanExon"]
  

  #Note: Tot_Length_LO and Tot_Length_Human will differ in some cases as the we didnt filter according
  #to exon length as example: the ID_Primate exon: ENSMMUG00000000021_1:221040728-221041220 Start LO is 244552322 and End 244552814
  #the associated Exon coordinate of human of such exon is Start: 244552322 and End: 244552475.
  #The end exon of the LO exon is longer therefore its length is going to be higher
  
  
  return(Info_PrimateLO)
  
}


##PSI 10


Info_BaboonCircs10_LO = Make_InfoPrimateLO(Info_BaboonCircs10_filt, List_Baboon10_filt) #8,067 #Num change:
Info_LemurCircs10_LO  = Make_InfoPrimateLO(Info_LemurCircs10_filt, List_Lemur10_filt) #4,842
Info_MacaqueCircs10_LO = Make_InfoPrimateLO(Info_MacaqueCircs10_filt, List_Macaque10_filt) #8,855
Info_MacacaCircs10_LO = Make_InfoPrimateLO(Info_MacacaCircs10_filt, List_Macaca10_filt) #7,959
Info_MarmosetCircs10_LO = Make_InfoPrimateLO(Info_MarmosetCircs10_filt, List_Marmoset10_filt) #Empty
Info_SquiMonkeyCircs10_LO = Make_InfoPrimateLO(Info_SquiMonkeyCircs10_filt, List_SquiMonkey10_filt) # 7,169
Info_ChimpCircs10_LO = Make_InfoPrimateLO(Info_ChimpCircs10_filt, List_Chimp10_filt)

## PSI 5
#Info_BaboonCircs5_LO = Make_InfoPrimateLO(Info_BaboonCircs5_filt, List_Baboon5_filt) 
#Info_LemurCircs5_LO = Make_InfoPrimateLO(Info_LemurCircs5_filt, List_Lemur5_filt)
#Info_MacaqueCircs5_LO = Make_InfoPrimateLO(Info_MacaqueCircs5_filt, List_Macaque5_filt)
#Info_MacacaCircs5_LO = Make_InfoPrimateLO(Info_MacacaCircs5_filt, List_Macaca5_filt)
#Info_MarmosetCircs5_LO = Make_InfoPrimateLO(Info_MarmosetCircs5_filt, List_Marmoset5_filt)
#Info_SquiMonkeyCircs5_LO = Make_InfoPrimateLO(Info_SquiMonkeyCircs5_filt, List_SquiMonkey5_filt)
#Info_ChimpCircs5_LO = Make_InfoPrimateLO(Info_ChimpCircs5_filt, List_Chimp5_filt)

####Remove NAs from Info_PrimatePSI_LO (e.i.: Start/EndExons that dont have a LiftedOver Start/End)
#PSI 10
Info_BaboonCircs10_LO = Info_BaboonCircs10_LO[complete.cases(Info_BaboonCircs10_LO),]
Info_LemurCircs10_LO  = Info_LemurCircs10_LO[complete.cases(Info_LemurCircs10_LO),]
Info_MacaqueCircs10_LO = Info_MacaqueCircs10_LO[complete.cases(Info_MacaqueCircs10_LO),]
Info_MacacaCircs10_LO = Info_MacacaCircs10_LO[complete.cases(Info_MacacaCircs10_LO),]
Info_MarmosetCircs10_LO = Info_MarmosetCircs10_LO[complete.cases(Info_MarmosetCircs10_LO),]
Info_SquiMonkeyCircs10_LO = Info_SquiMonkeyCircs10_LO[complete.cases(Info_SquiMonkeyCircs10_LO),]
Info_ChimpCircs10_LO = Info_ChimpCircs10_LO[complete.cases(Info_ChimpCircs10_LO),]


#PSI 5
#Info_BaboonCircs5_LO = Info_BaboonCircs5_LO[complete.cases(Info_BaboonCircs5_LO),]
#Info_LemurCircs5_LO  = Info_LemurCircs5_LO[complete.cases(Info_LemurCircs5_LO),]
#Info_MacaqueCircs5_LO = Info_MacaqueCircs5_LO[complete.cases(Info_MacaqueCircs5_LO),]
#Info_MacacaCircs5_LO = Info_MacacaCircs5_LO[complete.cases(Info_MacacaCircs5_LO),]
#Info_MarmosetCircs5_LO = Info_MarmosetCircs5_LO[complete.cases(Info_MarmosetCircs5_LO),]
#Info_SquiMonkeyCircs5_LO = Info_SquiMonkeyCircs5_LO[complete.cases(Info_SquiMonkeyCircs5_LO),]
#Info_ChimpCircs5_LO = Info_ChimpCircs5_LO[complete.cases(Info_ChimpCircs5_LO),]




##Make a function to add to Info_Primate<PSI>_LO the Human_CircID, the NumExons_Human and the Tot_Length_Human

Add_HumanInfo <- function(Info_PrimatePSI_LO, Info_HumanPSI, output) {
  
  print("Number of circRNAs in Human")
  print(length(unique(Info_HumanPSI$ID_Circ)))
  
  ##Coord_Circ in Info_HumanPSI, separate with start and End coords
  Info_HumanPSI$StartCirc = as.numeric(gsub(".*:","",gsub("-.*","", Info_HumanPSI$Coord_Circ)))
  
  Info_HumanPSI$EndCirc =  as.numeric(gsub(".*-","",gsub(".*:","", Info_HumanPSI$Coord_Circ)))
  
  
  ###Add columns to Info_Primate<PSI>  
  Info_PrimatePSI_LO$Human_CircID = NA
  Info_PrimatePSI_LO$NumExons_Human = NA
  Info_PrimatePSI_LO$Tot_Length_InfoHuman = NA
  
  
  
  #List of Primate circIDs with HUMAN exons IDs
  TempIDs_Primate = strsplit(Info_PrimatePSI_LO$HumanExons_IDs, split =  " | ", fixed = TRUE)
  names(TempIDs_Primate) = Info_PrimatePSI_LO$ID_Circ_Primate
  
  #Transform it into a DF
  TempIDs_Primate_DF = data.frame(ID = rep(names(TempIDs_Primate), sapply(TempIDs_Primate, length)),
                                  NumExons = rep(sapply(TempIDs_Primate, length), sapply(TempIDs_Primate, length)),
                                  ExonsIDs = unlist(TempIDs_Primate), stringsAsFactors = FALSE)
  rownames(TempIDs_Primate_DF) =NULL
  
  
  ####Add start and End info of Primate and LO
  TempIDs_Primate_DF$StartExonCircRNA_Primate = Info_PrimatePSI_LO[match(TempIDs_Primate_DF$ID, Info_PrimatePSI_LO$ID_Circ_Primate), 
                                                                   "StartExonCircRNA_Primate"]
  
  TempIDs_Primate_DF$StartExonCircRNA_Primate_LO = Info_PrimatePSI_LO[match(TempIDs_Primate_DF$ID, Info_PrimatePSI_LO$ID_Circ_Primate), 
                                                                      "StartExonCircRNA_Primate_LO"]
  
  TempIDs_Primate_DF$EndExonCircRNA_Primate =  Info_PrimatePSI_LO[match(TempIDs_Primate_DF$ID, Info_PrimatePSI_LO$ID_Circ_Primate), 
                                                                  "EndExonCircRNA_Primate"]
  
  TempIDs_Primate_DF$EndExonCircRNA_Primate_LO =  Info_PrimatePSI_LO[match(TempIDs_Primate_DF$ID, Info_PrimatePSI_LO$ID_Circ_Primate), 
                                                                     "EndExonCircRNA_Primate_LO"]
  
  #List of Human circIDs with HUMAN exons IDs
  TempIDs_Human = strsplit(Info_HumanPSI$IDs_Exons, split = " | ", fixed = TRUE)
  names(TempIDs_Human) = Info_HumanPSI$ID_Circ
  
  #Transform it into a DF
  TempIDs_Human_DF = data.frame(ID = rep(names(TempIDs_Human), sapply(TempIDs_Human, length)),
                                NumExons = rep(sapply(TempIDs_Human, length), sapply(TempIDs_Human, length)),
                                ExonsIDs = unlist(TempIDs_Human), stringsAsFactors = FALSE)
  
  #Add to TempIDs_Human_DF:  Start ExonCoord, End ExonCoord,  StartCirc_Coord, EndCircCord
  TempIDs_Human_DF$StartExon_Coord = as.numeric(gsub(".*:","",gsub("-.*","", TempIDs_Human_DF$ExonsIDs)))
  TempIDs_Human_DF$EndExon_Coord = as.numeric(gsub(".*-","",gsub(".*:","", TempIDs_Human_DF$ExonsIDs)))
  
  TempIDs_Human_DF$StartCirc_Coord = Info_HumanPSI[match(TempIDs_Human_DF$ID, Info_HumanPSI$ID_Circ),"StartCirc"]
  TempIDs_Human_DF$EndCirc_Coord =  Info_HumanPSI[match(TempIDs_Human_DF$ID, Info_HumanPSI$ID_Circ),"EndCirc"]
  
  
  #Define which exon is the Start, End, Start-End of Human data  
  TempIDs_Human_DF$StartOrEndExon_BSJ = NA
  TempIDs_Human_DF$StartOrEndExon_BSJ = NA
  TempIDs_Human_DF$StartOrEndExon_BSJ = NA
  
  starts = which(TempIDs_Human_DF$StartExon_Coord == TempIDs_Human_DF$StartCirc_Coord)
  ends = which(TempIDs_Human_DF$EndExon_Coord == TempIDs_Human_DF$EndCirc_Coord)
  
  starts_end = intersect(starts, ends)
  
  TempIDs_Human_DF[starts, "StartOrEndExon_BSJ"] = rep("Start", length(starts))
  TempIDs_Human_DF[ends, "StartOrEndExon_BSJ"] = rep("End", length(ends))
  
  TempIDs_Human_DF[starts_end, "StartOrEndExon_BSJ"] = rep("Start_End", length(starts_end))
  
  #Those circRNAs that dont have perfect Start/End Match
  
  
  #Which circRNAs do not have Start Coord
  IDsCoords = unique(TempIDs_Human_DF[complete.cases(TempIDs_Human_DF$StartOrEndExon_BSJ),"ID"])
  IDsNoCoords = setdiff(TempIDs_Human_DF$ID, IDsCoords)
  
  
  TempNoCoord = TempIDs_Human_DF[TempIDs_Human_DF$ID %in% IDsNoCoords,]
  
  #Calculate the difference between Start_Exon and Start_Circ, and End_Circ and End_Exon
  TempNoCoord$DiffStart = abs(TempNoCoord$StartExon_Coord - TempNoCoord$StartCirc_Coord)
  TempNoCoord$DiffEnd = abs(TempNoCoord$EndExon_Coord - TempNoCoord$EndCirc_Coord)
  
  #which Start_Diff is <= 100
  starts_nocoord = which(TempNoCoord$DiffStart <= 100)
  
  #which End_Diff is <= 100
  end_nocoord = which(TempNoCoord$DiffEnd <= 100)
  
  #which Start_Diff And End_Diff <= 100 (should be CircRNAs made of only one Exon)
  starts_end_nocoord = intersect(starts_nocoord, end_nocoord)
  
  ##Add to TempNoCoord Start and End
  TempNoCoord[starts_nocoord, "StartOrEndExon_BSJ"] = "Start"
  TempNoCoord[end_nocoord, "StartOrEndExon_BSJ"] = "End"
  TempNoCoord[starts_end_nocoord, "StartOrEndExon_BSJ"] = "Start_End"
  
  
  ####From TempNoCoord Get ID_Circ that dont Start or End defined
  NoStartCirc = setdiff(TempNoCoord$ID_Circ, TempNoCoord[complete.cases(TempNoCoord$StartOrEndExon_BSJ),"ID"])
  
  #Remove from TempNoCoord the NoStartCirc
  #TempNoCoord = TempNoCoord[!TempNoCoord$ID_Circ %in% NoStartCirc,]
  
  
  ##From List_Primate_filt$Primate_Circ remove those circRNAs used to make TempNoCoord
  
  TempIDs_Human_DF = TempIDs_Human_DF[!TempIDs_Human_DF$ID %in% IDsNoCoords,]
  
  
  ###Add to List_Primate_filt the TempNoCoord (now has Start and End Coord)
  TempIDs_Human_DF = rbind(TempIDs_Human_DF, TempNoCoord[,grep("Diff", colnames(TempNoCoord), invert = TRUE)])
  

  #########################Comparing TempIDs_Human_DF and TempIDs_Primate_DF###########################
  
  #Save Also Primates and Human CircRNAs that are made of conserved Exons but are not circRNAs also
  #In the primate or in the human (not shared between human and primate)
  
  
  ####################################################################################################
 
  ####Get here how many circRNAs how many circRNAs IDs are lost when complete.cases(TempIDs_Human_DF$StartOrEndExon_BSJ)
  MissingHumanCircRNAs = unique(TempIDs_Human_DF[is.na(TempIDs_Human_DF$StartOrEndExon_BSJ),"ID"])
   
  #For TempIDs_Human_DF subset those that onlye have START,End,Start-End Annotation (exon matches to BSJ)
  
  Annot_Human_DF = TempIDs_Human_DF[complete.cases(TempIDs_Human_DF$StartOrEndExon_BSJ),]
  
  ###How many circRNAs are lost when looking for start and end exon
  MissingCircRNAsIDs = setdiff(MissingHumanCircRNAs, unique(Annot_Human_DF$ID))
  
  
  #From Annot_Human_DF Get ExonIDs of those that are START
  
  HumanStart_ID = unique(Annot_Human_DF[grep("Start", Annot_Human_DF$StartOrEndExon_BSJ), "ExonsIDs"])
  
  HumanEnd_ID = unique(Annot_Human_DF[grep("End", Annot_Human_DF$StartOrEndExon_BSJ), "ExonsIDs"]) 
  
  
  ##Which HumanStart_ID match to TempIDs_Primate_DF$StartExonCircRNA_Primate_LO
  TempIDs_Primate_DF_MatchStartHuman = TempIDs_Primate_DF[TempIDs_Primate_DF$StartExonCircRNA_Primate_LO %in% HumanStart_ID,]
  
  #Using The StartExonCircRNA_Primate_LO in TempIDs_Primate_DF_MatchStartHuman subset Annot_Human_DF
  TempIDs_Primate_DF_MatchEndHuman = TempIDs_Primate_DF[TempIDs_Primate_DF$EndExonCircRNA_Primate_LO %in% HumanEnd_ID,]
  
  # make a single dataframe with the circRNAID, numExons,startExonPrimate, startExon_LO, endExonPrimate, endExon_LO,per circRNA
  
  SinglePrimateHumanDF = data.frame("Primate_CircID" = intersect(TempIDs_Primate_DF_MatchStartHuman$ID, TempIDs_Primate_DF_MatchEndHuman$ID),
                                    "Primate_numExons" = NA,
                                    "Primate_StartExon" = NA,
                                    "Primate_EndExon" = NA,
                                    "LO_StartExon" = NA,
                                    "LO_EndExon" =NA,
                                    "Human_CircID" = NA)
  
  SinglePrimateHumanDF$Primate_numExons = TempIDs_Primate_DF_MatchStartHuman[match(SinglePrimateHumanDF$Primate_CircID, TempIDs_Primate_DF_MatchStartHuman$ID), 
                                                                             "NumExons"]
  SinglePrimateHumanDF$Primate_StartExon = TempIDs_Primate_DF_MatchStartHuman[match(SinglePrimateHumanDF$Primate_CircID, TempIDs_Primate_DF_MatchStartHuman$ID), 
                                                                              "StartExonCircRNA_Primate"]
  SinglePrimateHumanDF$Primate_EndExon = TempIDs_Primate_DF_MatchEndHuman[match(SinglePrimateHumanDF$Primate_CircID, TempIDs_Primate_DF_MatchEndHuman$ID), 
                                                                          "EndExonCircRNA_Primate"]
  
  SinglePrimateHumanDF$LO_StartExon = TempIDs_Primate_DF_MatchStartHuman[match(SinglePrimateHumanDF$Primate_CircID, TempIDs_Primate_DF_MatchStartHuman$ID), 
                                                                         "StartExonCircRNA_Primate_LO"]
  
  SinglePrimateHumanDF$LO_EndExon = TempIDs_Primate_DF_MatchEndHuman[match(SinglePrimateHumanDF$Primate_CircID, TempIDs_Primate_DF_MatchEndHuman$ID), 
                                                                     "EndExonCircRNA_Primate_LO"]
  
  #Add the circRNA according to both the start and end ExonLO
  #HumanAnnotStart = Annot_Human_DF[grep("Start", Annot_Human_DF$StartOrEndExon_BSJ),]
  #HumanAnnotEnd = Annot_Human_DF[grep("End", Annot_Human_DF$StartOrEndExon_BSJ),]
  
  for (r in 1:nrow(SinglePrimateHumanDF) ) {
    
    IDs_start = Annot_Human_DF[grep(SinglePrimateHumanDF[r,"LO_StartExon"], Annot_Human_DF$ExonsIDs),"ID"]
    
    IDs_end = Annot_Human_DF[grep(SinglePrimateHumanDF[r,"LO_EndExon"], Annot_Human_DF$ExonsIDs),"ID"]
    
    ID = intersect(IDs_start, IDs_end)
    
    if (length(ID) == 1) {
      
      SinglePrimateHumanDF[r,"Human_CircID"] = ID 
      
    } 
    
    else {
      
      next
      
    }
  }
  
  ###This function returns 2 possible outputs
  
  PrimateWithOrthExons_DF = SinglePrimateHumanDF[is.na(SinglePrimateHumanDF$Human_CircID),]
  
  SinglePrimateHumanDF = SinglePrimateHumanDF[complete.cases(SinglePrimateHumanDF$Human_CircID),]
  
  
  #If ouput is human (which should be the same for all primates)
  #We get the info of each circRNA in human with its circRNAsIDs,NumExons, ExonIDs, Coordinates, StartCirc_Coord, EndCirc_Coord
  #StartOrEndExon_BSJ
  
  if (output == "Human") {
    
    return(Annot_Human_DF)
    
  }
  
  if (output == "Primate") {
    
    List_circRNAs = list("CircRNAsWithOrthExons" = PrimateWithOrthExons_DF,
                         "CircRNAsOrth2Human" = SinglePrimateHumanDF)
    
    return(List_circRNAs)
    
  }

  
}


###AQUI
HumanCircRNAs10 = Add_HumanInfo(Info_BaboonCircs10_LO, Info_HumanCirc10, "Human")

##Save HumanData

write.table(HumanCircRNAs10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/HumanInfoCircRNAs10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")



BaboonHumanCircs10 = Add_HumanInfo(Info_BaboonCircs10_LO, Info_HumanCirc10, "Primate")

LemurHumanCircs10_LO = Add_HumanInfo(Info_LemurCircs10_LO, Info_HumanCirc10,"Primate")

MacaqueHumanCircs10_LO = Add_HumanInfo(Info_MacaqueCircs10_LO, Info_HumanCirc10, "Primate")

MacacaHumanCircs10_LO  = Add_HumanInfo(Info_MacacaCircs10_LO, Info_HumanCirc10, "Primate")

MarmosetHumanCircs10_LO = Add_HumanInfo(Info_MarmosetCircs10_LO, Info_HumanCirc10, "Primate")

SquiMonkeyHumanCircs10_LO = Add_HumanInfo(Info_SquiMonkeyCircs10_LO, Info_HumanCirc10, "Primate")

ChimpHumanCircs10_LO = Add_HumanInfo(Info_ChimpCircs10_LO, Info_HumanCirc10, "Primate")

###Save data!!!!

write.table(BaboonHumanCircs10$CircRNAsOrth2Human, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/BaboonHumanCirc10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(BaboonHumanCircs10$CircRNAsWithOrthExons, 
            file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/BaboonHumanCirc10ExonsConvsNoBSJ.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")



write.table(LemurHumanCircs10_LO$CircRNAsOrth2Human, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/LemurHumanCirc10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")  

write.table(LemurHumanCircs10_LO$CircRNAsWithOrthExons, 
            file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/LemurHumanCirc10ExonsConsvNoBSJ.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t") 



  write.table(MacaqueHumanCircs10_LO$CircRNAsOrth2Human, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MacaqueHumanCirc10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(MacaqueHumanCircs10_LO$CircRNAsWithOrthExons, 
            file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MacaqueHumanCirc10ExonsConsvNoBSJ.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")


write.table(MacacaHumanCircs10_LO$CircRNAsOrth2Human, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MacacaHumanCirc10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")  

write.table(MacacaHumanCircs10_LO$CircRNAsWithOrthExons, 
            file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MacacaHumanCirc10ExonsConsvNoBSJ.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t") 



write.table(MarmosetHumanCircs10_LO$CircRNAsOrth2Human, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MarmosetHumanCirc10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(MarmosetHumanCircs10_LO$CircRNAsWithOrthExons, 
            file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/MarmosetHumanCirc10ExonsConsvNoBSJ.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")



write.table(SquiMonkeyHumanCircs10_LO$CircRNAsOrth2Human, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/SquiMonkeyHumanCirc10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(SquiMonkeyHumanCircs10_LO$CircRNAsWithOrthExons, 
            file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/SquiMonkeyHumanCirc10ExonsConsvNoBSJ.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

################################

write.table(ChimpHumanCircs10_LO$CircRNAsOrth2Human, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ChimpHumanCirc10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")


write.table(ChimpHumanCircs10_LO$CircRNAsWithOrthExons, 
            file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ChimpHumanCirc10ExonsConsvNoBSJ.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
























