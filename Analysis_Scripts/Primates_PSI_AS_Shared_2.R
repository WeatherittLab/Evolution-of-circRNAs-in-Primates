
#Scriipt to get orhtologue genes of Primates that have splicing events

#read orthologue primates data

All_orth=read.delim(file="/home/gaby/lab_Garvan/Primates/One2One_Orthologues/AllOrth_AllPrimates.txt",
                    header = TRUE, as.is = TRUE)


#Read one file of Whippet PSI primates results
path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs"

HumanDF = read.delim(gzfile("/home/gaby/lab_Garvan/Primates/Human_SamplesData_CircAtlas/WhippetOut/Brain_1_SRR787271_1.psi.gz"),
                   header = TRUE, as.is = TRUE) # 386,203

BaboonDF = read.delim(gzfile(paste(path,"Baboon", "Baboon_Brain_SRR1758903SRR1758904_1.psi.gz", sep="/")),
                    header = TRUE, as.is = TRUE) # 259,994

ChimpDF = read.delim(gzfile(paste(path,"Chimp", "ChimpanzeeCerebellum_1_SRR1758916.psi.gz", sep="/")), #Check chimp samples
                   header = TRUE, as.is = TRUE) #  277,886

LemurDF = read.delim(gzfile(paste(path,"Lemur", "Lemur_Cerebellum_SRR1758989SRR1758990_1.psi.gz", sep="/")),
                   header = TRUE, as.is = TRUE) # 239,025
LemurDF$Gene = gsub("\\..*","", LemurDF$Gene)

#Check Lemur files


MacacaDF = read.delim(gzfile(paste(path,"Macaca", "Macaca_Cerebellum_SRR1758948SRR1758949_1.psi.gz", sep="/")),
                    header = TRUE, as.is = TRUE) #260,029

MacaqueDF = read.delim(gzfile(paste(path,"Macaque", "MacaqueCerebellum_1_314.psi.gz", sep="/")),
                     header = TRUE, as.is = TRUE)  #276,902

MarmosetDF = read.delim(gzfile(paste(path,"Marmoset", "Marmoset_Brain_SRR1758977SRR1758979SRR1758978_1.psi.gz", sep="/")),
                      header = TRUE, as.is = TRUE) #307,707
MarmosetDF$Gene = gsub("\\..*","", MarmosetDF$Gene)

SquirrelMonkeyDF = read.delim(gzfile(paste(path,"SquirrelMonkey", "SquirrelMonkeyCerebellum_1_SRR1759034.psi.gz", sep="/")),
                            header = TRUE, as.is = TRUE)
SquirrelMonkeyDF$Gene = gsub("\\..*","", SquirrelMonkeyDF$Gene) #253264



###Make Orthologues one2one

Human = unique(HumanDF$Gene) #28,161

Baboon = unique(BaboonDF$Gene) #20,037

Chimp = unique(ChimpDF$Gene) #21,693

Lemur = unique(LemurDF$Gene) #18,921

Macaca = unique(MacacaDF$Gene) #20,160

Macaque = unique(MacaqueDF$Gene) #22,064

Marmoset = unique(MarmosetDF$Gene) #21,609

SquirrelMonkey = unique(SquirrelMonkeyDF$Gene) #19,222

###Subset genes that have orthologues

All_orth_hum = All_orth[match(intersect(Human, All_orth$ensembl_gene_id), All_orth$ensembl_gene_id),] 

#Get intersection of orthologues
hum_baboon = All_orth_hum[match(intersect(All_orth_hum$panubis_homolog_ensembl_gene, Baboon), All_orth_hum$panubis_homolog_ensembl_gene),]
  
hum_baboon_chimp = hum_baboon[match(intersect(hum_baboon$ptroglodytes_homolog_ensembl_gene, Chimp), hum_baboon$ptroglodytes_homolog_ensembl_gene),]
  
hum_baboon_chimp_lemur = hum_baboon_chimp[match(intersect(hum_baboon_chimp$mmurinus_homolog_ensembl_gene, Lemur),
                                                hum_baboon_chimp$mmurinus_homolog_ensembl_gene),]
  
hum_baboon_chimp_lemur_maca = hum_baboon_chimp_lemur[match(intersect(hum_baboon_chimp_lemur$mfascicularis_homolog_ensembl_gene, Macaca),
                                                     hum_baboon_chimp_lemur$mfascicularis_homolog_ensembl_gene),]
  
hum_baboon_chimp_lemur_maca_macq = hum_baboon_chimp_lemur_maca[match(intersect(hum_baboon_chimp_lemur_maca$mmulatta_homolog_ensembl_gene, Macaque),
                                                                          hum_baboon_chimp_lemur_maca$mmulatta_homolog_ensembl_gene),]
  
hum_baboon_chimp_lemur_maca_macq_marmo = hum_baboon_chimp_lemur_maca_macq[match(intersect(hum_baboon_chimp_lemur_maca_macq$cjacchus_homolog_ensembl_gene, Marmoset),
                                                                                hum_baboon_chimp_lemur_maca_macq$cjacchus_homolog_ensembl_gene),]
  
#8,547
Shared_Orth = hum_baboon_chimp_lemur_maca_macq_marmo[match(intersect(hum_baboon_chimp_lemur_maca_macq_marmo$sbboliviensis_homolog_ensembl_gene,
                                                                     SquirrelMonkey),hum_baboon_chimp_lemur_maca_macq_marmo$sbboliviensis_homolog_ensembl_gene),]


#####Subset the Shared Orth
Orgs_Names = colnames(Shared_Orth)
names(Orgs_Names) = c("Human","Chimp","Lemur","Macaque","Macaca","Marmoset","SquirrelMonkey","Baboon")

Subset_Shared <- function(OrgDF, Orgs_Names, organism, SharedDF) {
  
  common = intersect(OrgDF$Gene, SharedDF[[ Orgs_Names[[organism]] ]])
  
  not_comm = setdiff(OrgDF$Gene, common)

  OrgShared  = subset(OrgDF, !OrgDF[,"Gene"] %in% not_comm)
  
  ###Subset only CE
  OrgShared = OrgShared[grep("CE", OrgShared$Type),]
  
  return(OrgShared)
    
}


Human_Shared = Subset_Shared(HumanDF, Orgs_Names, "Human", Shared_Orth) # 94,149
Chimp_Shared = Subset_Shared(ChimpDF, Orgs_Names, "Chimp", Shared_Orth) # 98,792 
Lemur_Shared = Subset_Shared(LemurDF, Orgs_Names, "Lemur", Shared_Orth)  # 96,760
Macaque_Shared = Subset_Shared(MacaqueDF, Orgs_Names, "Macaque", Shared_Orth) # 97,986
Macaca_Shared = Subset_Shared(MacacaDF, Orgs_Names, "Macaca", Shared_Orth) # 98,282
Marmoset_Shared = Subset_Shared(MarmosetDF, Orgs_Names, "Marmoset", Shared_Orth) # 108,426
SquirrelMonkey_Shared = Subset_Shared(SquirrelMonkeyDF, Orgs_Names, "SquirrelMonkey", Shared_Orth) #97,100
Baboon_Shared = Subset_Shared(BaboonDF, Orgs_Names, "Baboon", Shared_Orth) #97,612

###Save tables
path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/"


write.table(Human_Shared, file=paste(path,"HumanThatAreOrth_PSI.txt", sep=""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
write.table(Chimp_Shared, file=paste(path,"ChimpThatAreOrth_PSI.txt", sep=""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
write.table(Lemur_Shared, file=paste(path,"LemurThatAreOrth_PSI.txt", sep=""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
write.table(Macaque_Shared, file=paste(path,"MacaqueThatAreOrth_PSI.txt", sep=""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
write.table(Macaca_Shared, file=paste(path,"MacacaThatAreOrth_PSI.txt", sep=""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
write.table(Marmoset_Shared, file=paste(path,"MarmosetThatAreOrth_PSI.txt", sep=""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
write.table(SquirrelMonkey_Shared, file=paste(path,"SquirrelMonkeyThatAreOrth_PSI.txt", sep=""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
write.table(Baboon_Shared, file=paste(path,"BaboonThatAreOrth_PSI.txt", sep=""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")


