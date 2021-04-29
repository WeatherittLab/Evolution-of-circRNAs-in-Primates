

#Script to compare the number of exons and position in transcript between groups (conserved and non-conserved)
library(ggplot2)


##Upload matrices of conserved and tissue conserved

InfoConservedCircRNAs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ConservedCircRNAs_Bob1.txt",
                                   header = TRUE, as.is = TRUE)


InfoTissueConservedCircRNAs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/TissueConservedCircRNAs_Bob2.txt",
                                         header = TRUE, as.is = TRUE)



#Get the info of circRNAs (in human)
HumanCircRNAs10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/HumanInfoCircRNAs10.txt",
                             header = TRUE, as.is = TRUE) #50,002; 22,667 circRNAs


################3Consv Data
ConsvNotCons_NumExons = as.data.frame(matrix(ncol = 2, nrow=nrow(InfoConservedCircRNAs)))
colnames(ConsvNotCons_NumExons) = c("Number_of_Exons", "Category")
rownames(ConsvNotCons_NumExons) = rownames(InfoConservedCircRNAs)

ConsvNotCons_NumExons$Number_of_Exons = HumanCircRNAs10[match(rownames(ConsvNotCons_NumExons), HumanCircRNAs10$ID), "NumExons"]
ConsvNotCons_NumExons$Category = InfoConservedCircRNAs[rownames(ConsvNotCons_NumExons),"Category"]


FewrExonsPval = wilcox.test(x = ConsvNotCons_NumExons[grep("Non Conserved", ConsvNotCons_NumExons$Category, invert = TRUE), "Number_of_Exons"],
                            y = ConsvNotCons_NumExons[grep("Non Conserved", ConsvNotCons_NumExons$Category), "Number_of_Exons"],
                            alternative = "less")

FewrExonsPval$p.value # 2.230608e-20


#Table the num of exons
Table_ConsvNotCons_Consv = table(ConsvNotCons_NumExons[grep("Non Conserved", ConsvNotCons_NumExons$Category, invert = TRUE), "Number_of_Exons"])

Table_ConsvNotCons_NoConsv = table(ConsvNotCons_NumExons[grep("Non Conserved", ConsvNotCons_NumExons$Category), "Number_of_Exons"])

###Make dataframe for barplot
NumExons_ConsNotCons = data.frame("Category" = c(rep("Conserved", 21),rep("Non Conserved", 21)), 
                                  "Number_of_Exons" = c(1:20,">=21"),
                                  "Number" = NA)

#Add info to the NumExons_ConsvNotCons
NumExons_ConsNotCons[grep("Non Conserved", NumExons_ConsNotCons$Category, invert = TRUE), "Number"] = as.vector(c(Table_ConsvNotCons_Consv[1:20], 
                                                                                                                  sum(Table_ConsvNotCons_Consv[21:length(Table_ConsvNotCons_Consv)])))

NumExons_ConsNotCons[grep("Non Conserved", NumExons_ConsNotCons$Category), "Number"] = as.vector(c(Table_ConsvNotCons_NoConsv[1:20], 
                                                                                                   sum(Table_ConsvNotCons_NoConsv[21:length(Table_ConsvNotCons_NoConsv)])))

NumExons_ConsNotCons$Number_of_Exons = factor(NumExons_ConsNotCons$Number_of_Exons, levels = unique(NumExons_ConsNotCons$Number_of_Exons))  

#Calculate %
NumExons_ConsNotCons$Percent = NA

tot_cons = sum(NumExons_ConsNotCons[grep("Non Conserved", NumExons_ConsNotCons$Category, invert = TRUE), "Number"])

for (i in grep("Non Conserved", NumExons_ConsNotCons$Category, invert = TRUE)) {
  
  NumExons_ConsNotCons[i, "Percent"] = (NumExons_ConsNotCons[i, "Number"]*100)/tot_cons
  
}  

tot_nocons = sum(NumExons_ConsNotCons[grep("Non Conserved", NumExons_ConsNotCons$Category), "Number"])

for (j in grep("Non Conserved", NumExons_ConsNotCons$Category)) {
  
  NumExons_ConsNotCons[j, "Percent"] = (NumExons_ConsNotCons[j, "Number"]*100)/tot_nocons
  
}
    

###################Tissue Consv
TissueConsvNotCons_NumExons = as.data.frame(matrix(ncol = 2, nrow=nrow(InfoTissueConservedCircRNAs)))
colnames(TissueConsvNotCons_NumExons) = c("Number_of_Exons", "Category")
rownames(TissueConsvNotCons_NumExons) = rownames(InfoTissueConservedCircRNAs)

TissueConsvNotCons_NumExons$Number_of_Exons = HumanCircRNAs10[match(rownames(TissueConsvNotCons_NumExons), HumanCircRNAs10$ID), "NumExons"]
TissueConsvNotCons_NumExons$Category = InfoTissueConservedCircRNAs[rownames(TissueConsvNotCons_NumExons), "Category"]

Table_TissueConsvNotCons_Consv = table(TissueConsvNotCons_NumExons[grep("Non Conserved", TissueConsvNotCons_NumExons$Category, invert = TRUE), "Number_of_Exons"])
  
Table_TissueConsvNotCons_NoConsv = table(TissueConsvNotCons_NumExons[grep("Non Conserved", TissueConsvNotCons_NumExons$Category), "Number_of_Exons"])


###Make dataframe for barplot
NumExons_TissueConsNotCons = data.frame("Category" = c(rep("Tissue Conserved", 21),rep("Non Conserved", 21)), 
                                        "Number_of_Exons" = c(1:20,">=21"),
                                        "Number" = NA)

NumExons_TissueConsNotCons[grep("Non Conserved", NumExons_TissueConsNotCons$Category, invert = TRUE), "Number"] = as.vector(c(Table_TissueConsvNotCons_Consv[1:20], 
                                                                                                                              sum(Table_TissueConsvNotCons_Consv[21:length(Table_TissueConsvNotCons_Consv)])))

NumExons_TissueConsNotCons[grep("Non Conserved", NumExons_TissueConsNotCons$Category), "Number"] = as.vector(c(Table_TissueConsvNotCons_NoConsv[1:20], 
                                                                                                               sum(Table_TissueConsvNotCons_NoConsv[21:length(Table_TissueConsvNotCons_NoConsv)])))


NumExons_TissueConsNotCons$Number_of_Exons = factor(NumExons_TissueConsNotCons$Number_of_Exons, levels = unique(NumExons_TissueConsNotCons$Number_of_Exons))


#Calculate %
NumExons_TissueConsNotCons$Percent = NA

tot_tisscons = sum(NumExons_TissueConsNotCons[grep("Non Conserved", NumExons_TissueConsNotCons$Category, invert = TRUE), "Number"])

for (i in grep("Non Conserved", NumExons_TissueConsNotCons$Category, invert = TRUE)) {
  
  NumExons_TissueConsNotCons[i, "Percent"] = (NumExons_TissueConsNotCons[i, "Number"]*100)/tot_tisscons
  
}  

tot_tissnocons = sum(NumExons_TissueConsNotCons[grep("Non Conserved", NumExons_TissueConsNotCons$Category), "Number"])

for (j in grep("Non Conserved", NumExons_TissueConsNotCons$Category)) {
  
  NumExons_TissueConsNotCons[j, "Percent"] = (NumExons_TissueConsNotCons[j, "Number"]*100)/tot_tissnocons
  
}



#######Make barplots

ConsvNotConsv = ggplot(NumExons_ConsNotCons, aes(x= Number_of_Exons, y = Percent, fill = Category )) +
  geom_bar(stat="identity", color="black", position = position_dodge())+
  theme_bw()+
  labs(title = "Number of Exons in Conserved and Non Conserved circRNAs" )

ConsvNotConsv + scale_fill_manual(values=c('#999999','#E69F00'))





TissueConsvNotConsv = ggplot(NumExons_TissueConsNotCons, aes(x= Number_of_Exons, y = Percent, fill = Category )) +
  geom_bar(stat="identity", color="black", position = position_dodge())+
  theme_bw()+
  labs(title = "Number of Exons in Tissue Conserved and Non Conserved circRNAs" )

TissueConsvNotConsv + scale_fill_manual(values=c('#999999','#E69F00'))




#################################
#Comparison of position in transcript between groups

#Info both InfoConservedCircRNAs and InfoTissueConservedCircRNAs make two new columns. 
#The first will have the  starting node and the second the ending node

InfoConservedCircRNAs$StartNode = gsub("-.*", "", rownames(InfoConservedCircRNAs))
InfoConservedCircRNAs$EndNode = paste(gsub("_.*", "", rownames(InfoConservedCircRNAs)) , gsub(".*-", "", gsub(".*_", "", rownames(InfoConservedCircRNAs))), sep="_")
  
InfoTissueConservedCircRNAs$StartNode = gsub("-.*", "", rownames(InfoTissueConservedCircRNAs))
InfoTissueConservedCircRNAs$EndNode = paste(gsub("_.*", "", rownames(InfoTissueConservedCircRNAs)), gsub(".*-", "", gsub(".*_", "", rownames(InfoTissueConservedCircRNAs))), sep="_" )
  

#Make function that given a table of Info<Tissue/Consv>CircRNAs I get the number of :
#circRNAs with unique Start and unique End
#circRNAs with Repeat Start
#circRNAs with Repeat End
#And if there are circRNAs with repeat Start and repeat End (?)

NumberPosCircRNAs <- function(InfoConsv, type_consv) {
  
  #Subset Conserved from Non Conserved
  Conserved = InfoConsv[grep("Non Conserved", InfoConsv$Category, invert = TRUE),]
  
  NonConserved = InfoConsv[grep("Non Conserved", InfoConsv$Category),]
  
  ##count starts and ends in Consv and NonCons
  TableConsvStart = table(Conserved$StartNode)
  TableConsvEnd = table(Conserved$EndNode)
  
  #########################Conserved analysis
  #Get which circRNAs have unique Start
  Circ_ConsvStartUniq = rownames(Conserved[ Conserved$StartNode %in% names(TableConsvStart[which(TableConsvStart == 1)]), ]) #576
  
  #Get which circRNAs have repeated Start
  Circ_ConsvRepeatStart = rownames(Conserved[Conserved$StartNode %in% names(TableConsvStart[which(TableConsvStart > 1)]), ]) #197
  
  #Get which circRNas have unique End
  Circ_ConsvEndUniq = rownames(Conserved[Conserved$EndNode %in% names(TableConsvEnd[which(TableConsvEnd == 1)]),]) #566
  
  #get which circRNAs have repeated End
  Circ_ConsvRepeatEnd = rownames(Conserved[Conserved$EndNode %in% names(TableConsvEnd[which(TableConsvEnd >1)]),]) #207
  
  
  #Which circRNAs are unique Start and Unique End
  UniqConsv = intersect(Circ_ConsvStartUniq, Circ_ConsvEndUniq) #424
  
  #Which circRNAs are repeated Start and Unique End
  RepStartUniqEnd = intersect(Circ_ConsvRepeatStart, Circ_ConsvEndUniq) #142
  
  #Which circRNAs are repeated End and Unique Start
  RepEndUniqStart = intersect(Circ_ConsvRepeatEnd, Circ_ConsvStartUniq) #152
  
  
  #Which circRNAs are repeated End and Repeated Start ?
  RepStartRepEnd = intersect(Circ_ConsvRepeatStart, Circ_ConsvRepeatEnd) #55
  
  
  #Make Dataframe to later use for barplot
  Pos_circRNAs = data.frame("Type" = rep(c("Unique_Start-End", "Repeated_Start",
                                       "Repeated_End", "Repeated_Start-End"), 2),
                            "Category"= c(rep(type_consv, 4), rep("Non Conserved", 4)),
                            "Number" = c(length(UniqConsv), length(RepStartUniqEnd), length(RepEndUniqStart), length(RepStartRepEnd), rep(NA,4)),
                            "Percentage" = NA)
  
  ########################Non Conserved analysis
  
  TableNonConsvStart = table(NonConserved$StartNode)
  TableNonConsvEnd = table(NonConserved$EndNode)
  
  ##Get which have unique Start
  Circ_NonConsvStartUniq = rownames(NonConserved[NonConserved$StartNode %in% names(TableNonConsvStart[which(TableNonConsvStart ==1)]), ]) #5,558
  
  #Get which circRNAs have repeated Start
  Circ_NonConsvRepeatStart = rownames(NonConserved[NonConserved$StartNode %in% names(TableNonConsvStart[which(TableNonConsvStart > 1)]), ]) #5,643
  
  #Get which circRNAs have unique End
  Circ_NonConsvEndUniq = rownames(NonConserved[NonConserved$EndNode %in% names(TableNonConsvEnd[which(TableNonConsvEnd == 1)]), ]) #5,346
  
  #Get which circRNAs have repeated End
  Circ_NonConsvRepeatEnd = rownames(NonConserved[NonConserved$EndNode %in% names(TableNonConsvEnd[which(TableNonConsvEnd > 1)]), ]) #5,855
  
  ###Which circRNAs are uniq Start and Uniq End
  UniqNonConsv = intersect(Circ_NonConsvStartUniq, Circ_NonConsvEndUniq) #2,813
  
  #"Which circRNAs are repeated Start and unique End
  NonConsv_RepStartUniqEnd = intersect(Circ_NonConsvRepeatStart, Circ_NonConsvEndUniq) #2,533
  
  #Which circRNAs are repeated End and Unique Start
  NonConsv_RepEndUniqStart = intersect(Circ_NonConsvRepeatEnd, Circ_NonConsvStartUniq) #2,745
  
  #Which circRNAs are repeated End and Repeated Start
  NonConsv_RepStartRepEnd = intersect(Circ_NonConsvRepeatStart, Circ_NonConsvRepeatEnd) #3,110
  
  #Add info to Pos_circRNAs data.frame
  Pos_circRNAs[grep("Non Conserved", Pos_circRNAs$Category),"Number"] = c(length(UniqNonConsv), length(NonConsv_RepStartUniqEnd), 
                                                                          length(NonConsv_RepEndUniqStart), length(NonConsv_RepStartRepEnd))
  
  ##Calculate %
  #for consverved
  rows_consv = grep("Non Conserved", Pos_circRNAs$Category, invert = TRUE)
  tot_consv = sum(Pos_circRNAs[rows_consv,"Number"])
  
  for (i in rows_consv) {
    
    Pos_circRNAs[i,"Percentage"] = (Pos_circRNAs[i,"Number"]*100)/tot_consv
    
  }
  
  #For non consv
  rows_nonconsv = grep("Non Conserved", Pos_circRNAs$Category)
  tot_nonconsv = sum(Pos_circRNAs[rows_nonconsv, "Number"])
  
  for (j in rows_nonconsv) {
    
    Pos_circRNAs[j, "Percentage"] = (Pos_circRNAs[j,"Number"]*100)/tot_nonconsv
    
  }
  
  return(Pos_circRNAs)
  
}

PosConsv_circRNAs = NumberPosCircRNAs(InfoConservedCircRNAs, "Conserved")

PosTissue_circRNAs = NumberPosCircRNAs(InfoTissueConservedCircRNAs, "Tissue Conserved")

###Calculate Chi-square

#UniqStartEnd_Pval =

#RepeatStart_Pval =
  
#RepeatEnd_Pval =
  
#Repeat_StartEnd_Pval =  

##Calculate Fisher exact test
MakeMatrixForFisherTest <- function(PosTable) {
  
  DF_fisher = matrix(nrow =2, ncol = 2)
  colnames(DF_fisher) = c("ToTest", "Others")
  rownames(DF_fisher) = c("Conserved", "Non Conserved")
  
  ToTestIDs = unique(PosTable$Type)
  
  PVal = vector(mode="list", length = length(ToTestIDs))
  names(PVal) = ToTestIDs
  
  for (t in ToTestIDs) {
  
  print(t)
    
  rows_typeconsv = intersect(grep("^Conserved$",PosTable$Category),  grep(paste("^", t, "$", sep=""), PosTable$Type))  
  rows_othersconsv = intersect(grep("^Conserved$",PosTable$Category),  grep(paste("^", t, "$", sep=""),PosTable$Type, invert = TRUE)) 
  
  DF_fisher["Conserved","ToTest"] = PosTable[rows_typeconsv,"Number"]
  DF_fisher["Conserved","Others"] = sum(PosTable[rows_othersconsv,"Number"])
  
  rows_typenonconsv = intersect(grep("^Non Conserved$",PosTable$Category),  grep(paste("^",t,"$", sep=""),PosTable$Type))  
  rows_othernonconsv = intersect(grep("^Non Conserved$",PosTable$Category),  grep(paste("^",t,"$", sep=""),PosTable$Type, invert = TRUE)) 
  
  DF_fisher["Non Conserved","ToTest"] = PosTable[rows_typenonconsv, "Number"]
  DF_fisher["Non Conserved","Others"] = sum(PosTable[rows_othernonconsv, "Number"])
   
  PVal[[t]] = fisher.test(DF_fisher)
  
  
   
  }
  
  return(PVal)
  
}
  
  
PValConsv = MakeMatrixForFisherTest(PosConsv_circRNAs)






  
##############Make barplots

PosConsv_plot = ggplot(PosConsv_circRNAs, aes(x= Type, y = Percentage, fill = Category )) +
  geom_bar(stat="identity", color="black", position = position_dodge())+
  theme_bw()+
  labs(title = "Types of Position in Transcript of circRNAs" )

PosConsv_plot + scale_fill_manual(values=c('#999999','#E69F00'))



PosTissueConsv_plot = ggplot(PosTissue_circRNAs, aes(x= Type, y = Percentage, fill = Category )) +
  geom_bar(stat="identity", color="black", position = position_dodge())+
  theme_bw()+
  labs(title = "Types of Position in Transcript of Tissue Conserved circRNAs" )

PosTissueConsv_plot + scale_fill_manual(values=c('#999999','#E69F00'))
