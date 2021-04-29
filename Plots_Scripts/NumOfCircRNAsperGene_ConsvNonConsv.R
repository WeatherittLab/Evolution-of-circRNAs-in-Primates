


#Script to plot how many circRNAs are made in genes with consv. circRNAs vs how many circRNAs are made in non-conserv sir
library(ggplot2)


#Load Info of conservd circRNAs

InfoConservedCircRNAs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ConservedCircRNAs_Bob1.txt",
                                   header = TRUE, as.is = TRUE) #2,382


#All orth
All_orth=read.delim(file="/home/gaby/lab_Garvan/Primates/One2One_Orthologues/AllOrth_AllPrimates.txt",
                    header = TRUE, as.is = TRUE) #11,649
colnames(All_orth) = c("Human", "Chimp", "Lemur","Macaque", "Macaca", "Marmoset", "SquirrelMonkey", "Baboon")






#From conserved list of circRNAs, get the genes same from conserved

GenesConsv = InfoConservedCircRNAs[grep("Non Conserved", InfoConservedCircRNAs$Category, invert = TRUE),]
#TableGenesConsv = table(GenesConsv$Ensembl_ID)


GenesNonConsv = InfoConservedCircRNAs[grep("Non Conserved", InfoConservedCircRNAs$Category),]
#TableGenesNonConsv = table(GenesNonConsv$Ensembl_ID)

#Table the number of EnsemblIDs in InfoConsvCircRNA
TableGenes = table(InfoConservedCircRNAs$Ensembl_ID)

#For shared Genes sum the number of circRNas
#shareGenes = intersect(names(TableGenesConsv), names(TableGenesNonConsv))



#From the dataframes of each primate subset those circRNAs that comes from genes of consv vs non-consv

#Make dataframe for plot

NumOfCircRNAs = data.frame("CircRNA_ID" = c(rownames(GenesConsv), rownames(GenesNonConsv)),
                           "GeneID" = c(GenesConsv$Ensembl_ID, GenesNonConsv$Ensembl_ID), 
                           "Num_circRNAs_in_same_gene" = NA,
                           "Category" = c(GenesConsv$Category, GenesNonConsv$Category))

#Add number of circRNAs
NumOfCircRNAs$Num_circRNAs_in_same_gene = TableGenes[match(NumOfCircRNAs$GeneID, names(TableGenes))]

ggplot(NumOfCircRNAs, aes(Num_circRNAs_in_same_gene, colour= Category) ) + 
  stat_ecdf(size=2) + 
  theme_bw()  +
  theme(text = element_text(size = 20)) +
  scale_color_manual(values = c("darkgoldenrod2", "azure4")) +
  labs(title = "Number of circRNAs" ) +
  xlab("Number of circRNAs in the same gene") +
  ylab("Fraction of Data (F(x))")



