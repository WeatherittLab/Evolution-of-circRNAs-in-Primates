
##Script to get circRNAs with tissue conservation profile

library(ggplot2)
library(gtable)

#PSI file
Primates_circRNA_PSI = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/PSI_Primates_ConsvBSJ.txt",
                                  header = TRUE, as.is = TRUE)

InfoConservedCircRNAs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ConservedCircRNAs_Bob1.txt",
                                   header = TRUE, as.is = TRUE) #11,974


load(file ="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/Primates_circRNA_BSJconsv")


#GeneNames files
GeneNames = read.delim(file="/home/gaby/lab_Garvan/Human_Data/Human_Ensembl_CodingGenes.txt",
                       header = TRUE, as.is = TRUE)

#Add gene Names to InfoConservedCircRNAs
InfoConservedCircRNAs$GeneName = GeneNames[match(InfoConservedCircRNAs$Ensembl_ID, GeneNames$ensembl_gene_id), "external_gene_name"]

#Make a function that gets the intersection od IDs of Human and each of the other primates 
#(Human And Chimp, Human And Baboon, ...)
#With the IDs of such intersections subset the PSI_DF and the samples of human and such primates
#(e.i if IDs are between Human And Chimp, I get samples of Human And Chimp)


GetIDs_PSI_ofInterest <-function(List_ConsvCircRNAs, PSI_DF, primates) {
  
  IDs_primatesPSI = vector(mode="list", length = length(primates))
  names(IDs_primatesPSI) = primates
  
  temp_IDs_primates = c()
  temp_IDs_OtherPrimates = c()
  
  for (p in primates) {
    
    #Get IDs of primates of interest
    temp_IDs_primates = Reduce(intersect, List_ConsvCircRNAs[c("Human",p)])
    
    #Get IDs of the other primates
    #temp_IDs_OtherPrimates = Reduce(union, List_ConsvCircRNAs[grep(paste(c("Human",p), collapse = "|"), names(List_ConsvCircRNAs), invert = TRUE)])
    
    #Get the IDs that are only in the primates of interest
    #temp_IDs_primates = setdiff(temp_IDs_primates, temp_IDs_OtherPrimates)
    
    
    ###With temp_IDs_primates subset PSI_DF
    temp_p =gsub(" ", "", p) #In the case of Squirrel Monkey, List Name has a blank space but colnames in PSI_DF dont,
    
    #Subset the PSI_DF accordin to shared IDs between Human and p and only get PSI samples of Human and p
    IDs_primatesPSI[[p]] = PSI_DF[temp_IDs_primates, grep(paste(c("Human",temp_p), collapse = "|"), colnames(PSI_DF))]

    
  }
  
  return(IDs_primatesPSI)
  
}

primates = names(Primates_circRNA_BSJconsv)[names(Primates_circRNA_BSJconsv) != "Human"]

HumanIntEachOtherPrimate_PSI = GetIDs_PSI_ofInterest(Primates_circRNA_BSJconsv, Primates_circRNA_PSI, primates)




#Make function to find if a same circRNA is expressed in the same tissue between species of the list of PSI dataframes filterd above

ConsvExp <- function(Primates_PSI) {
  
  #Which Primates are you comparing
  primates = unique(gsub("_.*", "", colnames(Primates_PSI)))
  
  #Get samples of each primate
  List_samples = vector(mode="list", length = length(primates))
  names(List_samples) = primates
  
  for (p in names(List_samples)) {
    
    List_samples[[p]] = colnames(Primates_PSI)[grep(p, colnames(Primates_PSI))]
    
  }
  
  
  #From the List_samples get only the tissues in those samples
  
  for (p in names(List_samples)) {
    
    List_samples[[p]] = unique(gsub("^(?:[^_]+_){1}([^_]+).*", "\\1", List_samples[[p]]))
    
  }
    
  #Which are the tissues in common
  common_tissues = Reduce(intersect, List_samples)
  
  if( !"Brain" %in% common_tissues ) {
    
    common_tissues = c(common_tissues, "Brain")
    
  }
  
  
  
  #Subset Primates_PSI according to common_tissues
  Primates_PSI_filt = Primates_PSI[,grep(paste(common_tissues, collapse = "|"), colnames(Primates_PSI))]
  
  #####To find if a circRNA is conserved in the same tissue
  
  #Make a list per tissue in common_tissue to keep such IDs that are expressed in the tissues in at least 2 samples
  List_ConsvTissues = vector(mode="list", length = length(common_tissues))
  names(List_ConsvTissues) = common_tissues
  
  #Loop through all the tissues in common_tissues
  
  for (t in common_tissues) {
    
    #Get the colnames of t tisue in th Primates_PSI_filt
    colnames = colnames(Primates_PSI_filt[,grep(t, colnames(Primates_PSI_filt))])
    
    ##double check that at all the primates of interest have columns of that tissue
    primates_dbcheck = unique(gsub("_.*", "", colnames))
    
    ###AQUI
    
    if (sum(primates %in% primates_dbcheck) == length(primates) ) {
      
      #Remove replicat number of the samples colnames and subset those samples that have more than one replicate
      samples = table(sub("_[^_]+$", "", colnames))
      samplesMoreThanOne = names(samples[samples >1])
      
      ###Check if samplesMoreThanOne are from both primates (human and the other primate)
      primates_samplesMoreThanOne = gsub("_.*", "", samplesMoreThanOne)
      
      #As we want to check if for that tissue, where both human and primate have replicates,
      #their circRNAs are expressed in the same tissue we have to get the IDs where
      #there is circRNA expression independently by primate (check expression in only human samples
      #and check expression in only primate samples)
      
      if (length(primates_samplesMoreThanOne) > 1) {
        
        tempIDs_Human =c()
        HumanSample = samplesMoreThanOne[grep("Human", samplesMoreThanOne)]
        
        for (c in colnames[grep(HumanSample, colnames)] ) {
        
          tempIDs_Human = unique(c(tempIDs_Human, rownames(Primates_PSI_filt[!is.na(Primates_PSI_filt[,c]),]) ))
        
        }  
        
        NotHumanSample = samplesMoreThanOne[grep("Human", samplesMoreThanOne, invert = TRUE)]
          
        tempIDs_OtherPrimate = c()
        
        for (c in colnames[grep(NotHumanSample, colnames)]) {
          
          tempIDs_OtherPrimate = unique(c(tempIDs_OtherPrimate, rownames(Primates_PSI_filt[!is.na(Primates_PSI_filt[,c]),]) ))
          
        }
        
        List_ConsvTissues[[t]] = intersect(tempIDs_Human, tempIDs_OtherPrimate)
        
          
      } else {
      
        TempIDs_MoreThanOne = c()
      
        #Loop through the samples that have more than one replicate
        #And check if in any of such replicates there is circRNA expression, IDs of such circRNAs
        #in the vector TempIDs_MoreThanOne
        for (c in colnames[grep(samplesMoreThanOne, colnames)]) {
        
          TempIDs_MoreThanOne = unique(c(TempIDs_MoreThanOne, rownames(Primates_PSI_filt[!is.na(Primates_PSI_filt[,c]),]) ))
        
        }
      
        #For samples with only one replicate
        samplesOnlyOne = names(samples[samples == 1])
        
        TempIDs_OnlyOneSample = c()
        
        #Get the IDs of circRNAs where in samples with only one replicate there is circRNA expression
        for (c in colnames[grep(samplesOnlyOne, colnames)]) {
       
          TempIDs_OnlyOneSample = unique(c(TempIDs_OnlyOneSample, rownames(Primates_PSI_filt[!is.na(Primates_PSI_filt[,c]),]))) 
        
        }
        
        List_ConsvTissues[[t]] = intersect(TempIDs_MoreThanOne, TempIDs_OnlyOneSample)
        
      }
      
      
    }
    
  }
  
    return(List_ConsvTissues)
  
}




HumanChimp_ConsvExpID = ConsvExp(HumanIntEachOtherPrimate_PSI$Chimp)

HumanBaboon_ConsvExpID = ConsvExp(HumanIntEachOtherPrimate_PSI$Baboon)

HumanMacaque_ConsvExpID = ConsvExp(HumanIntEachOtherPrimate_PSI$Macaque)

HumanMacaca_ConsvExpID =  ConsvExp(HumanIntEachOtherPrimate_PSI$Macaca)

HumanMarmoset_ConsvExpID = ConsvExp(HumanIntEachOtherPrimate_PSI$Marmoset)

HumanSquiMonkey_ConsvExpID = ConsvExp(HumanIntEachOtherPrimate_PSI$`Squirrel Monkey`)

HumanLemur_ConsvExpID = ConsvExp(HumanIntEachOtherPrimate_PSI$Lemur)


###Save neuronal shared circRNAs at least between human, chimp, baboon, (a macaque) and maybe marmoset if we are lucky :D

#Subset Human<Primate>_ConsvExpID of above primates of neuronal tissues
Neuronal = c("Brain","Cerebellum","FrontalCortex")

HumanChimp_ConsvNeuronal = HumanChimp_ConsvExpID[Neuronal]
HumanBaboon_ConsvNeuronal = HumanBaboon_ConsvExpID[Neuronal]
HumanMacaque_ConsvNeuronal = HumanMacaque_ConsvExpID[Neuronal]
HumanMacaca_ConsvNeuronal =  HumanMacaca_ConsvExpID[Neuronal]
HumanMarmoset_ConsvNeuronal = HumanMarmoset_ConsvExpID[Neuronal]


#Unlist IDs for neuronal circRNAs and get the union of them
HumanChimp_NeuronalIDs = unique(as.vector(unlist(HumanChimp_ConsvNeuronal))) #1,263
HumanBaboon_NeuronalIDs = unique(as.vector(unlist(HumanBaboon_ConsvNeuronal))) #948
HumanMacaque_NeuronalIDs = unique(as.vector(unlist(HumanMacaque_ConsvNeuronal))) #887
HumanMacaca_NeuronalIDs = unique(as.vector(unlist(HumanMacaca_ConsvNeuronal))) #548
HumanMarmoset_NeuronalIDs = unique(as.vector(unlist(HumanMarmoset_ConsvNeuronal))) #282


########Side interest (Get interesting conserved circRNAs for plot)
#

#Make the intersections
#HumanChimpBaboonMacaque = Reduce(intersect, list(HumanChimp_NeuronalIDs, HumanBaboon_NeuronalIDs, HumanMacaque_NeuronalIDs, HumanMacaca_NeuronalIDs))

#From tissue conserved, we got the ones in neuronal tissues from below primates
UntilMarmosetIDs = Reduce(intersect, list(HumanChimp_NeuronalIDs, HumanBaboon_NeuronalIDs, HumanMacaque_NeuronalIDs, HumanMacaca_NeuronalIDs, HumanMarmoset_NeuronalIDs))


#Retrived PSI values
NeuronalPSI = Primates_circRNA_PSI[UntilMarmosetIDs, grep(paste(Neuronal, collapse = "|"), colnames(Primates_circRNA_PSI))]

#Add gene names to NeuronalPSI
NeuronalPSI$GeneName = GeneNames[match(gsub("_.*", "", rownames(NeuronalPSI)), GeneNames$ensembl_gene_id),"external_gene_name"]


#From Gokool et al paper we found matching genes of circRNAs
IDs_Irina = c("CACNA1C", "ERC1", "CASC4", "ATRNL1", "CUL5", "ARGHGAP26", "KLHDC1")
circRNAIDs_Irina = rownames(NeuronalPSI[NeuronalPSI$GeneName %in% IDs_Irina,])


#Retrieved those circRNAIDs_Irina from all tissue, all primates samples
IntCircPSI = Primates_circRNA_PSI[circRNAIDs_Irina,]

#Add gene Name and rename rows with GeneName_Nodes
IntCircPSI$GeneName = GeneNames[match(gsub("_.*", "", rownames(IntCircPSI)), GeneNames$ensembl_gene_id), "external_gene_name"]

IntCircPSI$ID = paste(IntCircPSI$GeneName, gsub(".*_","",rownames(IntCircPSI)), sep="_")
rownames(IntCircPSI) = IntCircPSI$ID

#Remove GeneName and ID
IntCircPSI$GeneName = NULL
IntCircPSI$ID = NULL

#Do a function that for each interesting circRNA do a group barplot according to tissue (x) and color code according to primate

PlotPSI_intCirc <- function(CircPSI) {
  
  #Make dataframe per circPSI
  DF_circ = data.frame("Tissue" = gsub(".*_(.*)\\_.*", "\\1", colnames(CircPSI)),
                       "Primate" = gsub("_.*", "", colnames(CircPSI)),
                       "Samples" = colnames(CircPSI),
                       "PSI" = NA)
  
  #Order the Primate column phylogenetically
  phylo = c("Human", "Chimp", "Baboon", "Macaque", "Macaca", "Marmoset", "SquirrelMonkey", "Lemur")
  
  DF_circ = DF_circ[order(factor(DF_circ$Primate, levels = phylo)),]
  
  DF_circ$Primate =factor(DF_circ$Primate, levels = phylo )
  
  #Also order CircPSI columns as DF_circ$Samples
  CircPSI = CircPSI[,DF_circ$Samples]
  
  Plots = vector(mode="list", length = nrow(CircPSI))
  names(Plots) = rownames(CircPSI)
  
  for (circ in names(Plots)) {
   
    DF_circ$PSI = unlist(CircPSI[circ,])
    
    #For NA values change to 0
    DF_circ[is.na(DF_circ)] <- 0
     
    Plots[[circ]] = ggplot(DF_circ, aes(fill = Primate, y = PSI, x = Tissue)) +
      geom_bar(position = "dodge", stat="identity", colour="black") +
      scale_fill_brewer(palette = "YlOrBr", direction =  -1) +
      theme_bw() +
      theme(legend.title = element_text(size=16), legend.text = element_text(size = 15),
            axis.text = element_text(size = 12)) +
      ggtitle(paste("PSI values across primates and tissues of", circ, sep=" "))
    
    
  }
  
  return(Plots)
  
}

###AQUI

plots = PlotPSI_intCirc(IntCircPSI)

#Grid plot

#To only have one legende in the grid plot
#defined function to get legend from plots
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#Remove the legend from  ERC1_24-19
#plots$`ERC1_24-19` = plots$`ERC1_24-19` + theme(legend.position = "none")

#Save legend from #CACNA1C_11-9 and then remove it
#legend = get_legend(plots$`CACNA1C_11-9`)

#plots$`CACNA1C_11-9` = plots$`CACNA1C_11-9` + theme(legend.position = "none")

#For 2 fav circRNAs


g_CACNA1C = ggplotGrob(plots$`CACNA1C_11-9`) #chr12:2504436-2512984

g_ERC1 = ggplotGrob(plots$`ERC1_24-19`) #chr12:1180540-1204512


#For supplementary

g_ATRNL1_1 = ggplotGrob(plots$`ATRNL1_15-5`) #chr10:115120185-115215880
g_ATRNL1_2 = ggplotGrob(plots$`ATRNL1_30-28`)#chr10:115461941-115469329
g_CUL5 = ggplotGrob(plots$`CUL5_10-5`) #chr11:108046270-108054955
g_KLHDC1 = ggplotGrob(plots$`KLHDC1_14-3`) #chr14:49709159-49743805
g_CASC4 = ggplotGrob(plots$`CASC4_12-6`) #chr15:44328685-44380976


grid.arrange(g_CACNA1C, g_ERC1, nrow=2)


G_top2 = rbind(g_ERC1, g_CACNA1C, size = "first")

grid.arrange(g_ATRNL1_1, g_ATRNL1_2, g_CUL5, g_KLHDC1, g_CASC4, nrow=5)


#G_top2 = unit.(g_ERC1$widths, g_CACNA1C$widths)

grid.newpage()

#grid.draw(G_top2)










#Make a single vector of IDs per primate

ChimpIDs = unique(as.vector(unlist(HumanChimp_ConsvExpID)))
BaboonIDs = unique(as.vector(unlist(HumanBaboon_ConsvExpID)))
MacaqueIDs = unique(as.vector(unlist(HumanMacaque_ConsvExpID)))
MacacaIDs =  unique(as.vector(unlist(HumanMacaca_ConsvExpID)))
MarmosetIDs = unique(as.vector(unlist(HumanMarmoset_ConsvExpID)))
SquiMonkeyIDs = unique(as.vector(unlist(HumanSquiMonkey_ConsvExpID)))
LemurIDs = unique(as.vector(unlist(HumanLemur_ConsvExpID)))


#Make a List as Primates_circRNA_BSJconsv, where each primate is the vector of IDs that have expression same tissue as human (vector above)

Primates_circRNA_BSJ_ExprCons = list("Human" = unique(c(ChimpIDs, BaboonIDs, MacaqueIDs, MacacaIDs, MarmosetIDs, SquiMonkeyIDs, LemurIDs)),
                                     "Baboon" = BaboonIDs,
                                     "Chimp" = ChimpIDs,
                                     "Macaque" = MacaqueIDs,
                                     "Macaca" = MacacaIDs,
                                     "Marmoset" = MarmosetIDs,
                                     "Squirrel Monkey" = SquiMonkeyIDs,
                                     "Lemur" = LemurIDs)

#Save list of Primates_circRNA_BSJ_ExprCons
save(Primates_circRNA_BSJ_ExprCons, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/Primates_circRNA_BSJ_ExprCons")


#upset(fromList(Primates_circRNA_BSJ_ExprCons), order.by = "freq",  sets = names(Primates_circRNA_BSJ_ExprCons),
#      sets.bar.color = "purple", matrix.color = "purple", 
#      text.scale =2, shade.color = "purple", main.bar.color = "purple")







###################################33ABAJO NO

###SAVE THIS AGAIN

###Save Objects to use for heatmap
path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/TissueConsv"

save(HumanChimp_ConsvExpID, file=paste(path, "HumanChimp_ConsvID", sep="/"))
save(HumanBaboon_ConsvExpID, file=paste(path, "HumanBaboon_ConsvID", sep="/"))
save(HumanMacaque_ConsvExpID, file=paste(path, "HumanMacaque_ConsvID", sep="/"))
save(HumanMacaca_ConsvExpID, file=paste(path, "HumanMacaca_ConsvID", sep="/"))

save(HumanUntilOldWorld_ConsvExpID, file=paste(path, "HumanOldWorld_ConsvID", sep="/"))
save(HumanUntilMarmoset_ConsvExpID, file=paste(path, "HumanUntilMarmoset_ConsvID", sep="/"))
save(HumanUntilNewWorld_ConsvExpID, file=paste(path, "HumanUntilNewWorld_ConsvID", sep="/"))
save(HumanMarmoset_ConsvExpID, file=paste(path, "HumanMarmoset_ConsvID", sep="/"))
save(HumanSquiMonkey_ConsvExpID, file=paste(path, "HumanSquiMonkey_ConsvID", sep="/"))
save(HumanLemur_ConsvExpID, file=paste(path, "HumanLemur_ConsvID", sep="/"))

  
  
    