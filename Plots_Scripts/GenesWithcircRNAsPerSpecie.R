
#Script to count from orthologue genes how many of these orthologues have circRNAs in different primates

#For example: if Gene A in human has circRNA (at least one) count it one, if the same gene has circRNA in Chimp you count 2, and that..


#
#All orth
All_orth=read.delim(file="/home/gaby/lab_Garvan/Primates/One2One_Orthologues/AllOrth_AllPrimates.txt",
                    header = TRUE, as.is = TRUE) #11,649
colnames(All_orth) = c("Human", "Chimp", "Lemur","Macaque", "Macaca", "Marmoset", "SquirrelMonkey", "Baboon")


#Gene Expression

PrimatesExpre = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/ExpressionOrthPrimates.txt",
                           header = TRUE, as.is = TRUE)


InfoConservedCircRNAs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ConservedCircRNAs_Bob1.txt",
                                   header = TRUE, as.is = TRUE)


load( file ="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/Primates_circRNA_BSJconsv")


#SAve Coords for liftover
path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver"

#<file>_p1 is for files filtered as PSI >= 0.1
#<file>_p05 is files filtered as PSI >= 0.05 and 5 reads

Macaque_Coord = read.delim(file=paste(path, "MacaqueCirc_Coord_p05.txt", sep="/"),
                           header = FALSE, as.is = TRUE)
colnames(Macaque_Coord) = c("Chr", "Start", "End", "ID")

Chimp_Coord = read.delim(file=paste(path, "ChimpCirc_Coord_p05.txt", sep="/"),
                         header = FALSE, as.is = TRUE)
colnames(Chimp_Coord) = c("Chr", "Start", "End", "ID")


Human_Coord = read.delim(file=paste(path, "HumanCirc_Coord_p05.txt", sep="/"),
                         header = FALSE, as.is = TRUE)
colnames(Human_Coord) = c("Chr", "Start", "End", "ID")


Baboon_Coord = read.delim(file=paste(path, "BaboonCirc_Coord_p05.txt", sep="/"),
                          header = FALSE, as.is = TRUE)
colnames(Baboon_Coord) = c("Chr", "Start", "End", "ID")


Lemur_Coord = read.delim(file=paste(path, "LemurCirc_Coord_p05.txt", sep="/"),
                         header = FALSE, as.is = TRUE)
colnames(Lemur_Coord) = c("Chr", "Start", "End", "ID")


Marmoset_Coord = read.delim(file=paste(path, "MarmosetCirc_Coord_p05.txt", sep="/"),
                            header = FALSE, as.is = TRUE)
colnames(Marmoset_Coord) = c("Chr", "Start", "End", "ID")


SquiMonkey_Coord = read.delim(file=paste(path, "SquiMonkeyCirc_Coord_p05.txt", sep="/"),
                              header = FALSE, as.is = TRUE)
colnames(SquiMonkey_Coord) = c("Chr", "Start", "End", "ID")

Macaca_Coord = read.delim(file=paste(path, "MacacaCirc_Coord_p05.txt", sep="/"),
                          header = FALSE, as.is = TRUE)
colnames(Macaca_Coord) = c("Chr", "Start", "End", "ID")


#Add the gene ID to all dataframes
Macaque_Coord$PrimateGeneID = gsub("_.*", "", Macaque_Coord$ID)
Chimp_Coord$PrimateGeneID = gsub("_.*", "", Chimp_Coord$ID)
Human_Coord$PrimateGeneID = gsub("_.*", "", Human_Coord$ID)
Baboon_Coord$PrimateGeneID = gsub("_.*", "", Baboon_Coord$ID)
Lemur_Coord$PrimateGeneID = gsub("_.*", "", Lemur_Coord$ID)
Marmoset_Coord$PrimateGeneID = gsub("_.*", "", Marmoset_Coord$ID)
SquiMonkey_Coord$PrimateGeneID = gsub("_.*", "", SquiMonkey_Coord$ID)
Macaca_Coord$PrimateGeneID = gsub("_.*", "", Macaca_Coord$ID)
  
  
#For the non human data frames add the human GeneID

AddHumanID <- function(Primate_Coord, All_orth, primate) {
  
  Primate_Coord$HumanGeneID = All_orth[match(Primate_Coord$PrimateGeneID, All_orth[,primate]), "Human"]
  
  return(Primate_Coord)
  
}

Macaque_Coord = AddHumanID(Macaque_Coord, All_orth, "Macaque")
Chimp_Coord = AddHumanID(Chimp_Coord, All_orth, "Chimp")
Baboon_Coord = AddHumanID(Baboon_Coord, All_orth, "Baboon")
Lemur_Coord = AddHumanID(Lemur_Coord, All_orth, "Lemur")
Marmoset_Coord = AddHumanID(Marmoset_Coord, All_orth, "Marmoset")
SquiMonkey_Coord = AddHumanID(SquiMonkey_Coord, All_orth, "SquirrelMonkey")
Macaca_Coord = AddHumanID(Macaca_Coord, All_orth, "Macaca")

#Subset for each primate only the Human Gene IDs

HumanIDs = unique(Human_Coord$PrimateGeneID) #6,318
MacaqueIDs = unique(Macaque_Coord$HumanGeneID) #3,426
ChimpIDs = unique(Chimp_Coord$HumanGeneID) #4,028
BaboonIDs = unique(Baboon_Coord$HumanGeneID) #4,000
LemurIDs = unique(Lemur_Coord$HumanGeneID) #2,808
MarmosetIDs = unique(Marmoset_Coord$HumanGeneID) #3,282
SquiMonkeyIDs = unique(SquiMonkey_Coord$HumanGeneID) #3,896
MacacaIDs = unique(Macaca_Coord$HumanGeneID) #4,095


#Total number of orthologue genes with circRNAs (all genes)
AllOrthWithCirc = unique(c(HumanIDs, MacaqueIDs, ChimpIDs,BaboonIDs, LemurIDs, MarmosetIDs, SquiMonkeyIDs, MacacaIDs)) #7,285 (passes psi and reads cutoff)



#Join all the IDs and count how many the IDs are repeated (how many times the gene is in different primates)
AllIDs = c(HumanIDs, MacaqueIDs, ChimpIDs, BaboonIDs, LemurIDs, MarmosetIDs, SquiMonkeyIDs, MacacaIDs)
TimesGenes = table(AllIDs)

#Make a Datafrmae NumSpeciesShareGene (1-8) and how many genes are in such categories 

TimesGenesShared = data.frame("NumSpeciesShareGene" = c(1:8),
                              "Num" = NA)

TimesGenesShared[grep(1, TimesGenesShared$NumSpeciesShareGene), "Num"] = length(TimesGenes[TimesGenes == 1])
TimesGenesShared[grep(2, TimesGenesShared$NumSpeciesShareGene), "Num"] = length(TimesGenes[TimesGenes == 2])
TimesGenesShared[grep(3, TimesGenesShared$NumSpeciesShareGene), "Num"] = length(TimesGenes[TimesGenes == 3])
TimesGenesShared[grep(4, TimesGenesShared$NumSpeciesShareGene), "Num"] = length(TimesGenes[TimesGenes == 4])
TimesGenesShared[grep(5, TimesGenesShared$NumSpeciesShareGene), "Num"] = length(TimesGenes[TimesGenes == 5])
TimesGenesShared[grep(6, TimesGenesShared$NumSpeciesShareGene), "Num"] = length(TimesGenes[TimesGenes == 6])
TimesGenesShared[grep(7, TimesGenesShared$NumSpeciesShareGene), "Num"] = length(TimesGenes[TimesGenes == 7])
TimesGenesShared[grep(8, TimesGenesShared$NumSpeciesShareGene), "Num"] = length(TimesGenes[TimesGenes == 8])

###Total circRNAs without any cutoff
GenesWithcircRNAperSpecie = ggplot(data=TimesGenesShared, aes(x = NumSpeciesShareGene, y=Num)) +
  geom_bar(stat="identity", fill="deeppink4") +
  geom_text(aes(label = Num), vjust= 1.6, color="white", size= 3.5)+
  #scale_x_continuous(limits = c(1,8))+
  theme_bw()

#to remove axis as scientific notation
#require(scales)
#FiltPSIReads_circRNAs + scale_y_continuous(labels = comma) + labs(x = "Primate", y = "Number of circRNAs")





#Make GenesWithcircRNAperSpecie as % and add the proportion of circRNAs that are conserved in 1-2-3 organisms and the human specif



#####Make function to find circRNAs shared between specific primates
#Get intersections as in upsetplot

fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}


BinTable = fromList(Primates_circRNA_BSJconsv)

#Sum the times a circRNA is present
RowSumsConsv = rowSums(BinTable)

###Get IDs of all circRNAs that are shared between 2 species
SharedBy2 = RowSumsConsv[RowSumsConsv ==2] #2855

  
##Get IDs of all circRNAs that are shared between 3 species
SharedBy3 = RowSumsConsv[RowSumsConsv == 3] #1243


##Get IDs of all circRNAs that are shared between 4 species
SharedBy4 = RowSumsConsv[RowSumsConsv == 4] #679

##Get IDs of all circRNAs that are shared between 5 species
SharedBy5 = RowSumsConsv[RowSumsConsv == 5] #366


###Get IDs of all circRNAs thate are shared between 6 species
SharedBy6 = RowSumsConsv[RowSumsConsv == 6] #148


###Get IDs of all circRNAs that are shared between 7 species
SharedBy7 = RowSumsConsv[RowSumsConsv == 7] #55


###Get IDs of all circRNAs that are shared between  8 species
SharedBy8 = RowSumsConsv[RowSumsConsv == 8] #6



#Add to TimesGenesShare the values of circRNAs shared between species and the human specific ones (specie 1)
#Calculate % of genes shared and circRNAs shares (and human specific)
colnames(TimesGenesShared) = c("NumSpeciesShareGene", "Num_Genes")
TimesGenesShared$Percent_Genes = NA
TimesGenesShared$Num_circRNAs = NA
TimesGenesShared$Percent_circRNAs = NA

#Add values of Num_circRNAs =
TimesGenesShared[1,"Num_circRNAs"] = nrow(InfoConservedCircRNAs[grep("Non Conserved", InfoConservedCircRNAs$Category),])
TimesGenesShared[2,"Num_circRNAs"] = length(SharedBy2)
TimesGenesShared[3,"Num_circRNAs"] = length(SharedBy3)
TimesGenesShared[4,"Num_circRNAs"] = length(SharedBy4)
TimesGenesShared[5, "Num_circRNAs"] = length(SharedBy5)
TimesGenesShared[6,"Num_circRNAs"] = length(SharedBy6)
TimesGenesShared[7, "Num_circRNAs"] = length(SharedBy7)
TimesGenesShared[8, "Num_circRNAs"] = length(SharedBy8)

#calculate percentages 
tot_genes = sum(TimesGenesShared$Num_Genes)

for (i in 1:nrow(TimesGenesShared) ) {
  
  TimesGenesShared[i,"Percent_Genes"] = (TimesGenesShared[i,"Num_Genes"]*100)/tot_genes
  
}

tot_consvcir = sum(TimesGenesShared[1:nrow(TimesGenesShared),"Num_circRNAs"])

for (j in 1:nrow(TimesGenesShared) ) {
  
  TimesGenesShared[j, "Percent_circRNAs"] = (TimesGenesShared[j,"Num_circRNAs"]*100)/tot_consvcir
  
}

#Human circRNAs (only shared with one specie are its 100%)
#TimesGenesShared[1,"Percent_circRNAs"] = 100

#Arrange TimesGenesShare data to make grouped barplots 
GenesCircRNA_Shared_data = data.frame("Number_of_Species" = rep(1:8, 2),
                                      "Type" = c(rep("orthologue genes",8), rep("circRNAs", 8)),
                                      "Percentage" = c(TimesGenesShared$Percent_Genes, TimesGenesShared$Percent_circRNAs))


ggplot(GenesCircRNA_Shared_data, aes(fill = Type, y = Percentage, x = Number_of_Species)) +
  geom_bar(position = "dodge", stat="identity", colour="black") +
  theme_bw() +
  scale_fill_manual(values = c('#F2F3F4','#4D5656'))



