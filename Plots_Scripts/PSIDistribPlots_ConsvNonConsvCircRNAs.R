
#Script to compare PSI values between conserved and non conserved circRNAs

library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)


InfoConservedCircRNAs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ConservedCircRNAs_Bob1.txt",
                                   header = TRUE, as.is = TRUE) #11,974


InfoTissueConservedCircRNAs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/TissueConservedCircRNAs_Bob2.txt",
                                         header = TRUE, as.is = TRUE) #11,735

#Load list of conserved circRNAs
load(file ="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/Primates_circRNA_BSJconsv")





#Get PSI values of Human

Human_PSI = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Human_PSI.txt", 
                       header = TRUE, as.is = TRUE) #43,681



#Subset Human_PSI according to InfoConserved and InfoTissueConserved

InfoConsv_PSI = Human_PSI[rownames(InfoConservedCircRNAs),]
  
InfoTissueConsv_PSI = Human_PSI[rownames(InfoTissueConservedCircRNAs), ]

###Make a function to calculate Mean PSI value per tissue
GetMean_perTissue <-function(PSI_DF) {
  
  tissues = unique(gsub("_.*", "" ,colnames(PSI_DF)))
  
  PSI_Mean = as.data.frame(matrix(nrow = nrow(PSI_DF), ncol = length(tissues)))
  rownames(PSI_Mean) = rownames(PSI_DF)
  colnames(PSI_Mean) = tissues
  
  for (t in tissues) {
    
    PSI_Mean[,t] = rowMeans(PSI_DF[,grep(t, colnames(PSI_DF))], na.rm = TRUE)
    
  }
  
  return(PSI_Mean)
  
}

Consv_PSIMean = GetMean_perTissue(InfoConsv_PSI)

TissueConsv_PSIMean = GetMean_perTissue(InfoTissueConsv_PSI)


##Make dataframe used for ggplot violin plots according to category 

MakeDF_4plot <- function(PSI_DF, InfoConsv) {
  
  tissues = colnames(PSI_DF)
  
  DF = data.frame("PSI" = NA, 
                  "Category" = NA,
                  "Tissue" = NA,
                  "ID" = rep(rownames(PSI_DF), length(tissues)) )
  
  #Add Tissue 
  DF$Tissue = rep(tissues, nrow(PSI_DF))[order(rep(tissues, nrow(PSI_DF)))]
  
  #Add category
  DF$Category = InfoConsv[match(DF$ID, rownames(InfoConsv)), "Category"]
  
  #Add PSI values per tissue
  for (t in tissues) {
    
    DF[grep(t, DF$Tissue),"PSI"] = PSI_DF[,t]
    
  }
  
  return(DF)
  
}

Consv_DF = MakeDF_4plot(Consv_PSIMean, InfoConservedCircRNAs)

TissueConsv_DF = MakeDF_4plot(TissueConsv_PSIMean, InfoTissueConservedCircRNAs)

##Change NA to zero
#Consv_DF[is.na(Consv_DF)] <-0

#TissueConsv_DF[is.na(TissueConsv_DF)] <- 0

#Multiply by 100
Consv_DF$PSI = Consv_DF$PSI*100

TissueConsv_DF$PSI = TissueConsv_DF$PSI*100

###Make violin plots

ConsvPSIplot = ggplot(Consv_DF, aes(fill=Category, y=PSI, x=Tissue)) + 
    geom_violin(position="dodge", alpha=0.5) +
    scale_fill_viridis(discrete=T, name="") +
    theme_bw()  +
    xlab("") +
    labs(title = "PSI distributions across tissues in human for Conserved and Non Conserved circRNAs" ) +
    ylab("PSI")



TissueConsvPSIplot = ggplot(TissueConsv_DF, aes(fill=Category, y=PSI, x=Tissue)) + 
  geom_violin(position="dodge", alpha=0.5) +
  scale_fill_viridis(discrete=T, name="") +
  theme_bw()  +
  xlab("") +
  labs(title = "PSI distributions across tissues in human for Tissue Conserved and Non Conserved circRNAs" ) +
  ylab("PSI")


#Make function to calculate wilcox test between tissues comparing consv vs non consv

CompConsvNonConsv <-function(Consv) {
  
  #Make a table with Tissue and P-value to add after comparing
  #Consv vs non conserved
  DF_pval = data.frame("Tissue" = unique(Consv$Tissue),
                       "P-Value" = NA)
  
  for (t in DF_pval$Tissue) {
    
  tempConsv = Consv[grep(t, Consv$Tissue),]  
  tempWilcox = wilcox.test(tempConsv[grep("Non Conserved", tempConsv$Category, invert = TRUE), "PSI"],
                           tempConsv[grep("Non Conserved", tempConsv$Category), "PSI"], alternative = "greater")
    
   DF_pval[grep(t,DF_pval$Tissue),"P.Value"] = tempWilcox$p.value
    
  }
 
  return(DF_pval) 
  
}

Consv_Pval = CompConsvNonConsv(Consv_DF)

TissueConsv_Pval = CompConsvNonConsv(TissueConsv_DF)


#Make df with Mean values of circRNAs 
MedianDF <-function(ConsvPSI, InfoConsv) {
  
  DF_Median = data.frame("PSI_Median" = NA,
                       "Category" = InfoConsv$Category)
  
  rownames(DF_Median) = rownames(ConsvPSI)
  
  DF_Median$PSI_Median = apply(ConsvPSI, 1, median, na.rm = T)
  
  return(DF_Median)
  
}

ConsvAllTissues_MedianPSI = MedianDF(InfoConsv_PSI, InfoConservedCircRNAs)

TissueConsvAllTissue_MedianPSI = MedianDF(InfoTissueConsv_PSI, InfoTissueConservedCircRNAs)


#Make ecdf plot with Mean valu of all Consv vs Non Consv circRNAs 

Consv_ECDF = ggplot(ConsvAllTissues_MedianPSI, aes(PSI_Median, colour= Category) ) + 
            stat_ecdf() + 
            theme_bw()  +
            labs(title = "Cumulative plot of PSI distribution of Consverved and Non Conserved circRNAs" ) 

TissueConsv_ECDF = ggplot(TissueConsvAllTissue_MedianPSI, aes(PSI_Median, colour= Category) ) + 
                  stat_ecdf() +
                  theme_bw() +
                  labs(title="Cumulative plot of PSI distribution of Tissue Conserved and Non Conserved circRNAs")


#Wilcox text
ConsvPval = wilcox.test(x = ConsvAllTissues_MedianPSI[grep("Non Conserved", ConsvAllTissues_MedianPSI$Category , invert = TRUE), "PSI_Median"],
            y = ConsvAllTissues_MedianPSI[grep("Non Conserved", ConsvAllTissues_MedianPSI$Category), "PSI_Median"], alternative = "greater")


####Get set of interest of conserved circRNAs
#With those IDs filter the circRNA_PSI dataframe

GetIDs <-function(List_ConsvCircRNAs, primates) {
  
  #Get IDs of primates of interest
  IDs_primates = Reduce(intersect, List_ConsvCircRNAs[primates])
  
  #Get IDs of the other primates
  IDs_OtherPrimates = Reduce(union, List_ConsvCircRNAs[grep(paste(primates, collapse = "|"), names(List_ConsvCircRNAs), invert = TRUE)])
  
  #Get the IDs that are only in the primates of interest
  IDs_primates = setdiff(IDs_primates, IDs_OtherPrimates)
  
  ###With IDs_primates subset PSI_DF
  #primates= gsub(" ", "", primates) #In the case of Squirrel Monkey, List Name has a blank space but colnames in PSI_DF dont,
  #PSI_filt = PSI_DF[IDs_primates, grep(paste(primates, collapse = "|"), colnames(PSI_DF))]
  
  return(IDs_primates)
  
}

HumanIDs = rownames(InfoConservedCircRNAs[grep("Non Conserved", InfoConservedCircRNAs$Category),]) #11,201
HumanChimpIDs = GetIDs(Primates_circRNA_BSJconsv, c("Human","Chimp")) #995
HCB_IDs = GetIDs(Primates_circRNA_BSJconsv, c("Human","Chimp", "Baboon")) #141
HCBMcqsIDs = GetIDs(Primates_circRNA_BSJconsv, c("Human","Chimp","Baboon", "Macaque", "Macaca")) #109
HCBMcqsMarmoIDs =  GetIDs(Primates_circRNA_BSJconsv, c("Human","Chimp","Baboon", "Macaque", "Macaca", "Marmoset")) #40
UntilNewWorldIDs =  GetIDs(Primates_circRNA_BSJconsv, c("Human","Chimp","Baboon", "Macaque", "Macaca", "Marmoset", "Squirrel Monkey")) #32
#AllIDs =  GetIDs(Primates_circRNA_BSJconsv, c("Human","Chimp","Baboon", "Macaque", "Macaca", "Marmoset", "Squirrel Monkey", "Lemur")) #6


#After getting all the IDs of the intersections of interest subset the Human_PSI, calculate median PSI of each intersection of circRNAs
HumanPSI = apply(Human_PSI[HumanIDs,], 1, median, na.rm=T)
HumanChimpPSI = apply(Human_PSI[HumanChimpIDs,], 1, median, na.rm=T)
HCB_PSI = apply(Human_PSI[HCB_IDs,], 1, median, na.rm=T)
HCBMcqs_PSI = apply(Human_PSI[HCBMcqsIDs,], 1, median, na.rm=T)
HCBMcqsMarmo_PSI = apply(Human_PSI[HCBMcqsMarmoIDs,], 1, median, na.rm=T)
UntilNewWorld_PSI = apply(Human_PSI[UntilNewWorldIDs,], 1, median, na.rm=T)
#All_PSI = apply(Human_PSI[AllIDs,], 1, median, na.rm=T)

#To plot the ecdf of PSI distributions across intersections

#Make DF for ggplot
PSI_ok_consv =data.frame("PSI_Median" = c(HumanPSI, HumanChimpPSI, HCB_PSI, HCBMcqs_PSI, HCBMcqsMarmo_PSI, UntilNewWorld_PSI) ,
                         "ID" = c(names(HumanPSI),names(HumanChimpPSI), names(HCB_PSI),
                                  names(HCBMcqs_PSI), names(HCBMcqsMarmo_PSI), names(UntilNewWorld_PSI)),
                         "Category" = c(rep("Human", length(HumanPSI)),
                                        rep("Hominoids", length(HumanChimpPSI)),
                                        rep("Hominoids And Baboon", length(HCB_PSI)),
                                        rep("Hominoids And Old World Monkeys", length(HCBMcqs_PSI)),
                                        rep("Hominoids, Old World Monkeys and Marmoset", length(HCBMcqsMarmo_PSI)),
                                        rep("Hominoids, Old World Monkeys and New World Monkeys", length(UntilNewWorld_PSI))))



ggplot(PSI_ok_consv, aes(PSI_Median, colour= Category)) + 
  stat_ecdf(size=2) + 
  theme_bw()  +
  scale_color_brewer(palette = "Paired", direction =  -1) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size = 15)) +
  labs(title = "Cumulative plot of PSI distribution" ) 


#Make wilcox.text of PSI_of_consv between sets
  
#Human < Hominoids
Hominoids_Pval = wilcox.test(x = PSI_ok_consv[grep("^Hominoids$", PSI_ok_consv$Category),"PSI_Median"], 
            y = PSI_ok_consv[grep("Human", PSI_ok_consv$Category), "PSI_Median"], alternative = "greater")

Hominoids_Pval$p.value #8.065577e-19

#Human < HominoidsAnd Baboon
HomiBaboon_Pval = wilcox.test(x = PSI_ok_consv[grep("Hominoids And Baboon", PSI_ok_consv$Category),"PSI_Median"], 
                              y = PSI_ok_consv[grep("Human", PSI_ok_consv$Category), "PSI_Median"], alternative = "greater")

HomiBaboon_Pval$p.value # 1.519452e-09

#Hominoids < Old World
HominoidsOldWorld_Pval = wilcox.test(x = PSI_ok_consv[grep("^Hominoids And Old World Monkeys$", PSI_ok_consv$Category),"PSI_Median"], 
                             y = PSI_ok_consv[grep("^Hominoids$", PSI_ok_consv$Category), "PSI_Median"], alternative = "greater")

HominoidsOldWorld_Pval$p.value #1.029656e-08

#Hominods < Old World and Marmoset
HominoidsOldWorldMarmo_Pval = wilcox.test(x = PSI_ok_consv[grep("^Hominoids, Old World Monkeys and Marmoset$", PSI_ok_consv$Category),"PSI_Median"], 
                                     y = PSI_ok_consv[grep("^Hominoids$", PSI_ok_consv$Category), "PSI_Median"], alternative = "greater")

HominoidsOldWorldMarmo_Pval$p.value #1.983584e-08


#Old World < New World
OldWorldNewWorld_Pval =  wilcox.test(x = PSI_ok_consv[grep("^Hominoids, Old World Monkeys and New World Monkeys$", PSI_ok_consv$Category),"PSI_Median"], 
                                     y = PSI_ok_consv[grep("^Hominoids And Old World Monkeys$", PSI_ok_consv$Category), "PSI_Median"], alternative = "greater")

OldWorldNewWorld_Pval$p.value #  0.0947979

#Hominoids < New World
HominoidsNewWorld_Pval = wilcox.test(x = PSI_ok_consv[grep("^Hominoids, Old World Monkeys and New World Monkeys$", PSI_ok_consv$Category),"PSI_Median"], 
                                     y = PSI_ok_consv[grep("^Hominoids$", PSI_ok_consv$Category), "PSI_Median"], alternative = "greater")

HominoidsNewWorld_Pval$p.value # 2.136245e-06
