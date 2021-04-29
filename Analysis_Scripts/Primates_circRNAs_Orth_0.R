
# Read all PSI and gene expression files of primates and make a single
#Data frame for each tissue (one for PSI and one for Gene expression)

library(plyr)
#library(ggplot2)

#All orth
All_orth=read.delim(file="/home/gaby/lab_Garvan/Primates/One2One_Orthologues/AllOrth_AllPrimates.txt",
                    header = TRUE, as.is = TRUE)

Primates = c("Human", "Chimp", 
             "Macaque", "Baboon", 
             "Lemur", "Marmoset",
             "SquirrelMonkey", "Macaca")


#Read whippet 
path_whippet = "/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/BS/"


###Make a function to read all PSI/ gene.tpm files of a specific organisms and assign a variablename
#suffix_file is either ".psi.gz" or ".gene.tpm.gz"

read.files <- function(path, primate, suffix_file) {
  
  path = paste(path, primate, sep = "/")
  
  Files = list.files(path = path, pattern = suffix_file)
  
  ###Make a df with the tissue and the file of that tissue
  
  if (primate ==  "Macaque") {
    
    Files_DF = data.frame(Var_Name = gsub("_.*", "", gsub(".*ian_","", Files)),
                             File = paste(path, Files, sep = "/"), stringsAsFactors = FALSE)
    
    Files_DF = Files_DF[order(Files_DF$Var_Name),]
    
    
  } 
  
  
  if ( primate == "Chimp") {
    
    Files_DF = data.frame(Var_Name = gsub("_.*", "", gsub(".*ee_","", Files)),
                             File = paste(path, Files, sep = "/"), stringsAsFactors = FALSE)
    
    
    Files_DF = Files_DF[order(Files_DF$Var_Name),]
    
  }
  
  if ( primate == "Human" ){
    
    Files_DF = data.frame(Var_Name = gsub("_.*", "", Files),
                             File = paste(path, Files, sep = "/"), stringsAsFactors = FALSE)
    
    
    Files_DF = Files_DF[order(Files_DF$Var_Name),]
    
    #Remove Astrocyte and PrimaryNeuronStemCell
    
    Files_DF = Files_DF[grep("Astrocyte", Files_DF$Var_Name, invert = TRUE),]
    
    Files_DF = Files_DF[grep("PrimaryNeuronStemCell", Files_DF$Var_Name, invert = TRUE),]
    
  }
  
  if ( primate == "Baboon") {
    
    Files_DF = data.frame(Var_Name = gsub("_.*", "", gsub(".*oon_","", Files)),
                             File = paste(path, Files, sep = "/"), stringsAsFactors = FALSE)
    
    
    Files_DF = Files_DF[order(Files_DF$Var_Name),]
    
  }
  
  if ( primate == "Lemur") {
    
    Files_DF = data.frame(Var_Name = gsub("_.*", "", gsub(".*ur_","", Files)),
                             File = paste(path, Files, sep = "/"), stringsAsFactors = FALSE)
    
    
    Files_DF = Files_DF[order(Files_DF$Var_Name),]
    
  }
  
  
  if ( primate == "Marmoset") {
    
    Files_DF = data.frame(Var_Name = gsub("_.*", "", gsub(".*set_","", Files)),
                            File = paste(path, Files, sep = "/"), stringsAsFactors = FALSE)
    
    
    Files_DF = Files_DF[order(Files_DF$Var_Name),]
    
  }
  
  
  if ( primate == "SquirrelMonkey") {
    
    Files_DF = data.frame(Var_Name = gsub("_.*", "", gsub(".*key_","", Files)),
                             File = paste(path, Files, sep = "/"), stringsAsFactors = FALSE)
    
    
    Files_DF = Files_DF[order(Files_DF$Var_Name),]
    
  }
  
  
  if ( primate == "Macaca") {
    
    Files_DF = data.frame(Var_Name = gsub("_.*", "", gsub(".*an_","", Files)),
                             File = paste(path, Files, sep = "/"), stringsAsFactors = FALSE)
    
    
    Files_DF = Files_DF[order(Files_DF$Var_Name),]
    
  }
  
  ###Add the number of replicates to the Var_Name
  
  temp_table = c()
  
  for (v in unique(Files_DF$Var_Name)) {
    
    temp_table = table(Files_DF[grep(v, Files_DF$Var_Name), "Var_Name"])
    
    Files_DF[grep(v, Files_DF$Var_Name), "Var_Name"] = paste(Files_DF[grep(v, Files_DF$Var_Name), "Var_Name"],
                                                                   seq(1:temp_table), sep = "_")
    
  }
  
  #Read files and assign variable name according to Files_DF and save them in a LIST
  temp_file = data.frame()
  
  DF = vector(mode="list", length = length(Files_DF$File))
  
  names(DF) = Files_DF$Var_Name
  
  for (f in Files_DF$File) {
    
    temp_file = read.delim(gzfile(f), header = FALSE, as.is = TRUE)
    
    DF[[ Files_DF[match(f, Files_DF$File), "Var_Name"]  ]] = temp_file
    
    
    
  }
  
  return(DF)
  
}





#Read PSI files
Macaque_PSI_DFs = read.files(path_whippet, "Macaque", ".psi.gz_BS")
Chimp_PSI_DFs = read.files(path_whippet, "Chimp", ".psi.gz_BS")
Human_PSI_DFs = read.files(path_whippet, "Human", ".psi.gz_BS")
Baboon_PSI_DFs = read.files(path_whippet, "Baboon", ".psi.gz_BS")
Lemur_PSI_DFs = read.files(path_whippet, "Lemur", ".psi.gz_BS") 
Marmoset_PSI_DFs = read.files(path_whippet, "Marmoset", ".psi.gz_BS")
SquiMonkey_PSI_DFs = read.files(path_whippet, "SquirrelMonkey", ".psi.gz_BS")
Macaca_PSI_DFs = read.files(path_whippet, "Macaca", ".psi.gz_BS")

#Rename columns and Filter

cols = c("Gene", "Node", "Coord", "Strand", "Type", "Psi", "CI_Width", 
         "CI_Lo,Hi", "Total_Reads", "Complexity", "Entropy", "Inc_Paths","Exc_Paths", "Edges")

Macaque_PSI_DFs = lapply(Macaque_PSI_DFs, setNames, nm=cols)
Chimp_PSI_DFs = lapply(Chimp_PSI_DFs, setNames, nm=cols)
Human_PSI_DFs = lapply(Human_PSI_DFs, setNames, nm=cols)
Baboon_PSI_DFs = lapply(Baboon_PSI_DFs, setNames, nm=cols)
Lemur_PSI_DFs = lapply(Lemur_PSI_DFs, setNames, nm=cols)
Marmoset_PSI_DFs = lapply(Marmoset_PSI_DFs, setNames, nm=cols)
SquiMonkey_PSI_DFs = lapply(SquiMonkey_PSI_DFs, setNames, nm=cols)
Macaca_PSI_DFs = lapply(Macaca_PSI_DFs, setNames, nm=cols)

#Filter columns
col_filt = c("Gene", "Node", "Coord", "Strand", "Type", "Psi","Total_Reads")

Macaque_PSI_DFs = lapply(Macaque_PSI_DFs, "[", col_filt)
Chimp_PSI_DFs = lapply(Chimp_PSI_DFs, "[", col_filt)
Human_PSI_DFs = lapply(Human_PSI_DFs, "[", col_filt)
Baboon_PSI_DFs = lapply(Baboon_PSI_DFs,"[", col_filt)
Lemur_PSI_DFs = lapply(Lemur_PSI_DFs,"[", col_filt)
Marmoset_PSI_DFs = lapply(Marmoset_PSI_DFs,"[", col_filt)
SquiMonkey_PSI_DFs = lapply(SquiMonkey_PSI_DFs,"[", col_filt)
Macaca_PSI_DFs = lapply(Macaca_PSI_DFs, "[", col_filt)

###Lemur_PSI_DFs has gene IDs as ENS###.2, remove ".2"
for (t in names(Lemur_PSI_DFs)) {
  
  Lemur_PSI_DFs[[t]][["Gene"]] = gsub("\\..*","",Lemur_PSI_DFs[[t]][["Gene"]])
  
}

###Marmoset also..
for (t in names(Marmoset_PSI_DFs)) {
  
  Marmoset_PSI_DFs[[t]][["Gene"]] = gsub("\\..*","",Marmoset_PSI_DFs[[t]][["Gene"]])
  
}

##SquirrelMonkey to...
for (t in names(SquiMonkey_PSI_DFs)) {
  
  SquiMonkey_PSI_DFs[[t]][["Gene"]] = gsub("\\..*","",SquiMonkey_PSI_DFs[[t]][["Gene"]])
  
}

#check size of dataframes in tissues
lapply(Macaque_PSI_DFs, dim)
#Braincerebellum:16,304
#Brainfrontalcortex:16,082 
#Liver:5,960 #Lung:7,116 
#SkeletalMusc: 3538
#Spleen:7,912 

#

lapply(Chimp_PSI_DFs, dim)
#Brain 19,926
#Cerebellum:13,566
#Colon:12,243 #FrontalCortex:17,360
#Heart:7,325 
#Liver_1:1,722
#Lung_1:6,833 #Lung_2:7,760
#SkeletalMs_1:4,527
#Spleen_1:3,578 #Spleen_2:3,913




lapply(Human_PSI_DFs, dim)
#Brain_1:10,121 #Brain_2:10,584
#Brain_3:14,889
#Cerebellum_1:32,947  #Cerebellum_2:37,505
#Colon_1:2,947, Colon_2:775, Colon_3: 669, Colon_4: 8,774
#FrontalCortex_1: 57,737 #FrontalCortex_2:55,626
#Heart_1:17,124, #Heart_2: 25,921, 
#Heart_3: 27,333, #Heart_4:4,615
#Kidney:1,869 
#Liver_1:23,040, #Liver_2: 4,168, #Liver_3: 12,650 
#Lung_1:11,853, #Lung_2: 28,136, #Lung_3: 1,499
#SkeletalMuscle_1: 20,885, #SkeletalMuscle_2: 21,905
#SkeletalMuscle_3: 928
#Spleen_1: 6,664, #Spleen_2: 7,341, #Spleen_3:6,908


lapply(Baboon_PSI_DFs, dim)
#Cerebellum:22,764
#Colon:6,175 #FrontalCortex:17,948
#Heart:7,009 #Kidney:9,778 #Liver:6,810
#Lung:12,203 #LymphNode:10,848
#Pituitary:9,848 #SkeletalMus_1:6,032
#SkMus_2:4,553 #Spleen:12,812 #TempLobe:20,656

lapply(Lemur_PSI_DFs, dim)
#Cerebellum_1:7,926 #Cerebell_2: 6,959
#Colon:4,939 #FrontCortx_1:5,589 #FrontCortx_2:8,850
#Kidney_1:2,226 #Kidney_2:3,752 #Liver:4,269
#Lung:9,170 #SkMus_1:1,192 #SkuMus_2:2,708
#SkuMusc_3:1,155 #Spleen_1:1,162 #Spleen_2:2,111
#TempLobe:10,829

lapply(Marmoset_PSI_DFs, dim)
#BoneMarrow:10,365 #BrainLeftHemi:12,608 #BrainRightHemi:16,036
#Colon:7,157 #Heart_1:7,390 #Heart_2:11,670 #Kidney_1:8,817
#Liver_1:6,788 #Lung:11,613 #LymphNode:10,084 #Pituitary:9,441
#SkeletMus:4,527 #Spleen:9,915

lapply(SquiMonkey_PSI_DFs, dim)
#BoneMarrow:12,565 #Cereb:14,248 #Colon:7,845 #FrontCort:14,421
#Heart:7,132 #Kidney:9,206 #Liver:5,315 #Lung:9,381
#LymphNode:10,198 #Pituit:9,632 #SkeMus:4,896 #Spleen_1:6,285
#Spleen_2:5,440 #TempLob:15,099

lapply(Macaca_PSI_DFs, dim)
#Cereb_1:6,179 #Cereb_2:7,628 #Colon:10,650 
#FrontCort_1:9,583 #Heart:8,079 #Kidney:9,488
#Liver:11,903 #Lung:12,521 #LympNode:14,689
#Pituit:12,339 #SkeMus:6,109 #Spleen:12,814
#TempLob:8,423 #Tym:11,062



#Count how many unique BS you get for each primate

NumberUniqueBS <- function(Primate_PSI_DFs) {
  
  all_IDs = c()
  
  #loop through each sample tissue 
  for (n in names(Primate_PSI_DFs)) {
    
    #
    Primate_PSI_DFs[[n]][["ID_BS"]] = paste(Primate_PSI_DFs[[n]][["Gene"]], Primate_PSI_DFs[[n]][["Coord"]], sep="_" )
    
    all_IDs = c(all_IDs, Primate_PSI_DFs[[n]][["ID_BS"]] )
    
  }
  
  all_IDs = unique(all_IDs)
  
  return(all_IDs)
  
}

IDs_Macaque = NumberUniqueBS(Macaque_PSI_DFs)
IDs_Chimp = NumberUniqueBS(Chimp_PSI_DFs)
IDs_Macaca = NumberUniqueBS(Macaca_PSI_DFs)
IDs_Marmoset = NumberUniqueBS(Marmoset_PSI_DFs)
IDs_Lemur = NumberUniqueBS(Lemur_PSI_DFs)
IDs_SquiMonkey = NumberUniqueBS(SquiMonkey_PSI_DFs)
IDs_Human = NumberUniqueBS(Human_PSI_DFs)
IDs_Baboon = NumberUniqueBS(Baboon_PSI_DFs)

###Make a table

UniqueCircRNAs = data.frame("Primate" = c("Macaque", "Chimp", "Macaca", "Marmoset", "Lemur", "Squirrel Monkey",
                                          "Human", "Baboon"),
                            "Total circRNAs" = c(length(IDs_Macaque), length(IDs_Chimp), length(IDs_Macaca), length(IDs_Marmoset),
                                                 length(IDs_Lemur),length(IDs_SquiMonkey), length(IDs_Human), length(IDs_Baboon)))

#save table
write.table(UniqueCircRNAs, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/BS/NumberOfCircRNAs/TotalNumberUniqueCircRNAs.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

##Make a dataframe with the number of BJ for each organisms for each tissue
#before and after filtering accordint to different parameters


dataframe4filt <- function(Org_PSI_DFs) {
  
  DF = as.data.frame(matrix(nrow = length(Org_PSI_DFs)*3, 
                            ncol = length(Org_PSI_DFs)+1))
  
  colnames(DF) = c("Tissues", names(Org_PSI_DFs))
  
  DF[,"Tissues"] = c(names(Org_PSI_DFs),
                     paste(names(Org_PSI_DFs), "orth", sep="_"),
                     paste(names(Org_PSI_DFs), "exp_psi", sep = "_"))
  
  DF = DF[order(DF[,"Tissues"]),]
  
  for (c in colnames(DF)[2:ncol(DF)]) {
    
    DF[match(c, DF$Tissues),c] = nrow(Org_PSI_DFs[[grep(c, names(Org_PSI_DFs))]])
    
  }
  
  return(DF)
}

DF4filt_Macaque = dataframe4filt(Macaque_PSI_DFs)
DF4filt_Macaca = dataframe4filt(Macaca_PSI_DFs)
DF4filt_Chimp = dataframe4filt(Chimp_PSI_DFs)
DF4filt_Baboon = dataframe4filt(Baboon_PSI_DFs)
DF4filt_Lemur = dataframe4filt(Lemur_PSI_DFs)
DF4filt_Human = dataframe4filt(Human_PSI_DFs)
DF4filt_Marmoset = dataframe4filt(Marmoset_PSI_DFs)
DF4filt_SquiMonkey = dataframe4filt(SquiMonkey_PSI_DFs)
    
###Filter according to All_orth

filter_rows <- function(PSI_DF, primate, All_orth) {
  
  if ( primate == "Macaque") {
    
    for (t in names(PSI_DF)) {
    
      rows = intersect(All_orth[,"mmulatta_homolog_ensembl_gene"], PSI_DF[[t]][,"Gene"])  
      
      PSI_DF[[t]] = PSI_DF[[t]][PSI_DF[[t]][,"Gene"] %in% rows,]

    
    }
    
  }
  
  #####Chimp
  if ( primate == "Chimp") {
    
    for (t in names(PSI_DF)) {
      
      rows = intersect(All_orth[,"ptroglodytes_homolog_ensembl_gene"], PSI_DF[[t]][,"Gene"])  
      
      PSI_DF[[t]] = PSI_DF[[t]][PSI_DF[[t]][,"Gene"] %in% rows,]
      
      
    }
    
  }
  
  #####Human
 
  if ( primate == "Human") {
    
    for (t in names(PSI_DF)) {
      
      rows = intersect(All_orth[,"ensembl_gene_id"], PSI_DF[[t]][,"Gene"])  
      
      PSI_DF[[t]] = PSI_DF[[t]][PSI_DF[[t]][,"Gene"] %in% rows,]
      
      
    }
    
  }
  
  ###Baboon
  if ( primate == "Baboon") {
    
    for (t in names(PSI_DF)) {
      
      rows = intersect(All_orth[,"panubis_homolog_ensembl_gene"], PSI_DF[[t]][,"Gene"])  
      
      PSI_DF[[t]] = PSI_DF[[t]][PSI_DF[[t]][,"Gene"] %in% rows,]
      
      
    }
    
  }
  
  ####Lemur
  if ( primate == "Lemur") {
    
    for (t in names(PSI_DF)) {
      
      rows = intersect(All_orth[,"mmurinus_homolog_ensembl_gene"], PSI_DF[[t]][,"Gene"])  
      
      PSI_DF[[t]] = PSI_DF[[t]][PSI_DF[[t]][,"Gene"] %in% rows,]
      
      
    }
    
  }
  
  ####Marmoset
  if ( primate == "Marmoset") {
    
    for (t in names(PSI_DF)) {
      
      rows = intersect(All_orth[,"cjacchus_homolog_ensembl_gene"], PSI_DF[[t]][,"Gene"])  
      
      PSI_DF[[t]] = PSI_DF[[t]][PSI_DF[[t]][,"Gene"] %in% rows,]
      
      
    }
    
  }
  
  ####Squirrel Monkey
  if ( primate == "SquirrelMonkey") {
    
    for (t in names(PSI_DF)) {
      
      rows = intersect(All_orth[,"sbboliviensis_homolog_ensembl_gene"], PSI_DF[[t]][,"Gene"])  
      
      PSI_DF[[t]] = PSI_DF[[t]][PSI_DF[[t]][,"Gene"] %in% rows,]
      
      
    }
    
  }
  
  ####Macaca
  if ( primate == "Macaca") {
    
    for (t in names(PSI_DF)) {
      
      rows = intersect(All_orth[,"mfascicularis_homolog_ensembl_gene"], PSI_DF[[t]][,"Gene"])  
      
      PSI_DF[[t]] = PSI_DF[[t]][PSI_DF[[t]][,"Gene"] %in% rows,]
      
      
    }
    
  }
  
  ####
  return(PSI_DF)
  
}


###Rows Filter (According to orthologues)
Macaque_PSI_filt = filter_rows(Macaque_PSI_DFs, "Macaque", All_orth)
Chimp_PSI_filt = filter_rows(Chimp_PSI_DFs, "Chimp", All_orth)
Human_PSI_filt =filter_rows(Human_PSI_DFs, "Human", All_orth)
Baboon_PSI_filt = filter_rows(Baboon_PSI_DFs,"Baboon", All_orth)
Lemur_PSI_filt = filter_rows(Lemur_PSI_DFs,"Lemur", All_orth)
Marmoset_PSI_filt = filter_rows(Marmoset_PSI_DFs,"Marmoset", All_orth)
SquiMonkey_PSI_filt = filter_rows(SquiMonkey_PSI_DFs,"SquirrelMonkey", All_orth)
Macaca_PSI_filt = filter_rows(Macaca_PSI_DFs, "Macaca", All_orth)

###Calculate unique BSJ/circRNAs from orthologue genes
MacaqueIDs_orth = NumberUniqueBS(Macaque_PSI_filt)
ChimpIDs_orth = NumberUniqueBS(Chimp_PSI_filt)
HumanIDs_orth = NumberUniqueBS(Human_PSI_filt)
BaboonIDs_orth = NumberUniqueBS(Baboon_PSI_filt)
LemurIDs_orth =  NumberUniqueBS(Lemur_PSI_filt)
MarmosetIDs_orth = NumberUniqueBS(Marmoset_PSI_filt)
SquiMonkeyIDs_orth = NumberUniqueBS(SquiMonkey_PSI_filt)
MacacaIDs_orth = NumberUniqueBS(Macaca_PSI_filt)

###Make a table

UniqueCircRNAs_orth = data.frame("Primate" = c("Macaque", "Chimp", "Macaca", "Marmoset", "Lemur", "Squirrel Monkey",
                                          "Human", "Baboon"),
                            "Total circRNAs" = c(length(MacaqueIDs_orth), length(ChimpIDs_orth), length(MacacaIDs_orth), length(MarmosetIDs_orth),
                                                 length(LemurIDs_orth),length(SquiMonkeyIDs_orth), length(HumanIDs_orth), length(BaboonIDs_orth)))

#save table
write.table(UniqueCircRNAs_orth, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/BS/NumberOfCircRNAs/TotalNumberUniqueCircRNAs_orthGenes.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")



#####AQUI tienes q considerar Sequencing Depth of Organism - Tissue  


###Filter according to PSI and Reads
filt_psi_exp <- function(PSI_filt, psi_cutoff, read_cutoff) {
  
  
  for (t in names(PSI_filt)) {
    
    #Filter according to read_cutoff
    PSI_filt[[t]] = PSI_filt[[t]][PSI_filt[[t]][["Total_Reads"]] >= read_cutoff,]
    
    #Filter according to PSI
    PSI_filt[[t]] = PSI_filt[[t]][PSI_filt[[t]][["Psi"]] >= psi_cutoff, ]
    
    
  }
  
  return(PSI_filt)
  
}

#Filt PSI and exp
#Macaque_filt = filt_psi_exp(Macaque_PSI_filt, psi_cutoff = 0.2, read_cutoff = 10)
#Chimp_filt = filt_psi_exp(Chimp_PSI_filt, psi_cutoff = 0.2, read_cutoff = 10)
#Human_filt = filt_psi_exp(Human_PSI_filt, psi_cutoff = 0.2, read_cutoff = 10)
#Baboon_filt = filt_psi_exp(Baboon_PSI_filt, psi_cutoff = 0.2, read_cutoff = 10)
#Lemur_filt = filt_psi_exp(Lemur_PSI_filt, psi_cutoff = 0.2, read_cutoff = 10)
#Marmoset_filt = filt_psi_exp(Marmoset_PSI_filt, psi_cutoff = 0.2, read_cutoff = 10)
#SquiMonkey_filt = filt_psi_exp(SquiMonkey_PSI_filt, psi_cutoff = 0.2, read_cutoff = 10)
#Macaca_filt = filt_psi_exp(Macaca_PSI_filt,psi_cutoff = 0.2, read_cutoff = 10)

#######With lower PSI cutoff (0.1) and 10 reads
#Macaque_filt = filt_psi_exp(Macaque_PSI_filt, psi_cutoff = 0.1, read_cutoff = 10)
#Chimp_filt = filt_psi_exp(Chimp_PSI_filt, psi_cutoff = 0.1, read_cutoff = 10)
#Human_filt = filt_psi_exp(Human_PSI_filt, psi_cutoff = 0.1, read_cutoff = 10)
#Baboon_filt = filt_psi_exp(Baboon_PSI_filt, psi_cutoff = 0.1, read_cutoff = 10)
#Lemur_filt = filt_psi_exp(Lemur_PSI_filt, psi_cutoff = 0.1, read_cutoff = 10)
#Marmoset_filt = filt_psi_exp(Marmoset_PSI_filt, psi_cutoff = 0.1, read_cutoff = 10)
#SquiMonkey_filt = filt_psi_exp(SquiMonkey_PSI_filt, psi_cutoff = 0.1, read_cutoff = 10)
#Macaca_filt = filt_psi_exp(Macaca_PSI_filt,psi_cutoff = 0.1, read_cutoff = 10)

#####Even lower PSI cutoff (0.05)
Macaque_filt = filt_psi_exp(Macaque_PSI_filt, psi_cutoff = 0.05, read_cutoff = 5)
Chimp_filt = filt_psi_exp(Chimp_PSI_filt, psi_cutoff = 0.05, read_cutoff = 5)
Human_filt = filt_psi_exp(Human_PSI_filt, psi_cutoff = 0.05, read_cutoff = 5)
Baboon_filt = filt_psi_exp(Baboon_PSI_filt, psi_cutoff = 0.05, read_cutoff = 5)
Lemur_filt = filt_psi_exp(Lemur_PSI_filt, psi_cutoff = 0.05, read_cutoff = 5)
Marmoset_filt = filt_psi_exp(Marmoset_PSI_filt, psi_cutoff = 0.05, read_cutoff = 5)
SquiMonkey_filt = filt_psi_exp(SquiMonkey_PSI_filt, psi_cutoff = 0.05, read_cutoff = 5)
Macaca_filt = filt_psi_exp(Macaca_PSI_filt,psi_cutoff = 0.05, read_cutoff = 5)

###Get number of IDs and make a dataframe of the number of circRNAs IDs that passed the cutoff
MacaqueIDs_ExpFilt = NumberUniqueBS(Macaque_filt)
ChimpIDs_ExpFilt = NumberUniqueBS(Chimp_filt)
HumanIDs_ExpFilt = NumberUniqueBS(Human_filt)
BaboonIDs_ExpFilt = NumberUniqueBS(Baboon_filt)
LemurIDs_ExpFilt =  NumberUniqueBS(Lemur_filt)
MarmosetIDs_ExpFilt = NumberUniqueBS(Marmoset_filt)
SquiMonkeyIDs_ExpFilt = NumberUniqueBS(SquiMonkey_filt)
MacacaIDs_ExpFilt = NumberUniqueBS(Macaca_filt)


###Make a table

UniqueCircRNAs_ExpFilt = data.frame("Primate" = c("Macaque", "Chimp", "Macaca", "Marmoset", "Lemur", "Squirrel Monkey",
                                               "Human", "Baboon"),
                                 "Total circRNAs" = c(length(MacaqueIDs_ExpFilt), length(ChimpIDs_ExpFilt), length(MacacaIDs_ExpFilt), 
                                                      length(MarmosetIDs_ExpFilt), length(LemurIDs_ExpFilt),length(SquiMonkeyIDs_ExpFilt), 
                                                      length(HumanIDs_ExpFilt), length(BaboonIDs_ExpFilt)))

#save table
write.table(UniqueCircRNAs_ExpFilt, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/BS/NumberOfCircRNAs/TotalNumberUniqueCircRNAs_orthExpPSIFilt.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")


###Save info PSI of filtered BSJ
save(Macaque_filt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/BS/PSI_filt/Macaque_filt")
save(Chimp_filt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/BS/PSI_filt/Chimp_filt")
save(Human_filt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/BS/PSI_filt/Human_filt")
save(Baboon_filt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/BS/PSI_filt/Baboon_filt")
save(Lemur_filt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/BS/PSI_filt/Lemur_filt")
save(Marmoset_filt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/BS/PSI_filt/Marmoset_filt")
save(SquiMonkey_filt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/BS/PSI_filt/SquiMonkey_filt")
save(Macaca_filt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/BS/PSI_filt/Macaca_filt")

#Check size
lapply(Macaque_filt, dim)
####With PSI = 0.05 and 5 reads
#Cerebellum:5,174
#FrontalCortex:5,248 
#Liver:1,944 #Lung:2,137
#SkeletalMusc: 1,106
#Spleen:2,506

lapply(Chimp_filt, dim)
####With PSI = 0.05 and 5 reads
#Brain:4112 #Cerebellum:3825
#Colon:2568 #FrontalCortex:5136
#Heart:1863 
#Liver_1:694 ? 
#Lung_1:2116 #Lung_2:2224
#SkeletalMs_1:1527
#Spleen_1:1371 #Spleen_2:1348


lapply(Human_filt, dim)

####With PSI = 0.05 and 5 reads
#Brain_1:5,273 #Brain_2: 4,993
#Brain_3: 6,012 #Cerebellum_1:8,032
#Cerebellum_2: 10,875, #Colon_1:1,547
#Colon_2: 326, #Colon_3:270, #Colon_4: 2,841
#FrontalCortex_1: 16,043, #FrontalCortex_2:16,703
#Heart_1:5,009 #Heart_2: 6,226 #Heart_3: 4,486 #Heart_4:2,313
#Liver_1:3,504 #Liver_2:1,595 #Lung_1:2,897
#Lung_2: 6,713 #Lung_3: 3,126 #Lung_4:679
#SkeletalMuscle_1: 9,178 #SkeletalMuscle_2: 10,772
#SkeletalMuscle_3: 387, #Spleen_1: 2,622, #Spleen_2:2,454
#Spleen_3: 2,469


lapply(Baboon_filt, dim)

####With PSI = 0.05 and 5 reads
#Brain:4208 #Cerebellum:5451
#Colon:2021 #FrontalCortex:3,961
#Heart_1:1,702
#Liver_1:1,982
#Lung:2,701 
#SkeletalMus_1:1,465
#SkMus_2:1,273 #Spleen:2,745

lapply(Lemur_filt, dim)

####With PSI = 0.05 and 5 reads
#Cerebellum_1:2,609 #Cerebell_2:2,518
#Colon:1,725 #FrontCortx_1:2,616#FrontCortx_2:3,107
#Liver:875
#Lung:1,554 #SkMus_1:357 #SkuMus_2:506
#SkuMusc_3:348 #Spleen_1:528 #Spleen_2:779
#Brain:2,101

lapply(Marmoset_filt, dim)

####With PSI = 0.05 and 5 reads
#Brain:3,511
#Colon:2,065 #Heart_1:1,920 #Heart_2:2,643
#Liver_1:1,993 #Lung:2,458 #LymphNode:2,792
#SkeletMus:1,133 #Spleen:2,483


lapply(SquiMonkey_filt, dim)

####With PSI = 0.05 and 5 reads
#BoneMarrow:2,366 #Cereb:3,871 #Colon:1,518 #FrontCort:5,604
#Heart:1,679 #Kidney:1,796 #Liver:1,306 #Lung:2,077
#LymphNode:1837 #Pituit:2,868 #SkeMus:1,152 #Spleen_1:1,933
#Spleen_2:1,892 #TempLob:4,765

lapply(Macaca_filt, dim)

####With PSI = 0.05 and 5 reads
#Cereb_1:2,400 #Cereb_2:2,476 #Colon: 2,503
#FrontCort_1:2,513 #Heart:2,383 #Kidney:2,355
#Liver:3,578 #Lung:2,686 #LympNode:3,646
#Pituit:3,839 #SkeMus:1,542 #Spleen:3,229
#TempLob:2,168 #Tym:2,805



###Add to the DF4filt_Org the rows size of filtered by Orth and by Exp_PSI
AddInfo_DF4filt <- function(DF4filt_Org, Org_PSI_filt, Org_filt) {
  
  for (t in names(Org_PSI_filt)) {
  
      DF4filt_Org[match(paste(t, "orth", sep="_"), DF4filt_Org[,"Tissues"]), t] = nrow(Org_PSI_filt[[t]])
  
      DF4filt_Org[match(paste(t, "exp_psi", sep="_"),  DF4filt_Org[,"Tissues"]), t] = nrow(Org_filt[[t]])
      
  }
    
  return(DF4filt_Org)
  
}

DF4filt_Macaque = AddInfo_DF4filt(DF4filt_Macaque, Macaque_PSI_filt, Macaque_filt)
DF4filt_Macaca = AddInfo_DF4filt(DF4filt_Macaca, Macaca_PSI_filt, Macaca_filt)
DF4filt_Chimp = AddInfo_DF4filt(DF4filt_Chimp, Chimp_PSI_filt, Chimp_filt)
DF4filt_Baboon = AddInfo_DF4filt(DF4filt_Baboon, Baboon_PSI_filt, Baboon_filt)
DF4filt_Lemur = AddInfo_DF4filt(DF4filt_Lemur, Lemur_PSI_filt, Lemur_filt)
DF4filt_Human = AddInfo_DF4filt(DF4filt_Human, Human_PSI_filt, Human_filt)
DF4filt_Marmoset = AddInfo_DF4filt(DF4filt_Marmoset, Marmoset_PSI_filt, Marmoset_filt)
DF4filt_SquiMonkey = AddInfo_DF4filt(DF4filt_SquiMonkey, SquiMonkey_PSI_filt, SquiMonkey_filt)


#Make unique ID according to Gene_Node and make a single data.frame
#For each organism
#Also separate The Coord into columns Chr, Start, End

unique_ids = function(Org_filt) {
  
  for (t in names(Org_filt)) {
    
    Org_filt[[t]][,"ID"] = paste(Org_filt[[t]][["Gene"]], Org_filt[[t]][["Node"]], sep="_")
    
    Org_filt[[t]][,"Chr"] = gsub(":.*","", Org_filt[[t]][["Coord"]])
    Org_filt[[t]][,"Start"] = as.numeric(gsub(".*:","", gsub("-.*", "", Org_filt[[t]][["Coord"]])))
    Org_filt[[t]][,"End"] = as.numeric(gsub(".*-", "", Org_filt[[t]][["Coord"]]))
      
    rownames(Org_filt[[t]]) = Org_filt[[t]][["ID"]]
    
  }
  
  return(Org_filt)
  
}



###Get unique Ids and coords sep

Macaque_IDs = unique_ids(Macaque_filt)
Chimp_IDs = unique_ids(Chimp_filt)
Human_IDs = unique_ids(Human_filt)
Baboon_IDs = unique_ids(Baboon_filt)
Lemur_IDs = unique_ids(Lemur_filt)
Marmoset_IDs = unique_ids(Marmoset_filt)
SquiMonkey_IDs = unique_ids(SquiMonkey_filt)
Macaca_IDs = unique_ids(Macaca_filt)



###AQUI
#Make a single DF for each organisms with PSI and Total_Reads with each tissue info

DF_PSI = function(Org_IDs) {
  
  cols = c("Psi", "ID")

  
  Org_IDs = lapply(Org_IDs, "[", cols)
  
  #Make a single dataframe with all the IDs
  IDs = unique(as.vector(unlist(lapply(Org_IDs, "[", "ID"))))
  
  DF_PSI = as.data.frame(matrix(nrow=length(IDs), ncol = length(Org_IDs)))
  colnames(DF_PSI) = paste(names(Org_IDs), "PSI", sep = "_")
  rownames(DF_PSI) = IDs
  
  for (sample in names(Org_IDs)) {
    
    rows = rownames(Org_IDs[[sample]])
    DF_PSI[rows,paste(sample, "PSI", sep="_")] = Org_IDs[[sample]][,"Psi"]
    
  }
  
  return(DF_PSI)
  
}

Macaque_PSI = DF_PSI(Macaque_IDs)
Chimp_PSI = DF_PSI(Chimp_IDs)
Human_PSI = DF_PSI(Human_IDs)
Baboon_PSI = DF_PSI(Baboon_IDs)
Lemur_PSI = DF_PSI(Lemur_IDs)
Marmoset_PSI = DF_PSI(Marmoset_IDs)
SquiMonkey_PSI = DF_PSI(SquiMonkey_IDs)
Macaca_PSI = DF_PSI(Macaca_IDs)

###Save files
#/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific

write.table(Chimp_PSI, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Chimp_PSI.txt",
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")
write.table(Lemur_PSI, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Lemur_PSI.txt",
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")
write.table(Macaca_PSI, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Macaca_PSI.txt",
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")
write.table(Macaque_PSI, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Macaque_PSI.txt",
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")
write.table(Marmoset_PSI,file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Marmoset_PSI.txt",
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")
write.table(SquiMonkey_PSI, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/SquirrelMonkey_PSI.txt",
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")
write.table(Baboon_PSI, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Baboon_PSI.txt",
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")
write.table(Human_PSI, file="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver/ForTissueSpecific/Human_PSI.txt",
            col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")



#Make a single dataframe for each organism with the Chr, Start, End columns
#Check if Start is < End, if not switch the columns

DF_Coord = function(Org_IDs) {
  
  cols = c("Chr", "Start", "End", "ID")
  
  Org_IDs = lapply(Org_IDs, "[", cols)

  Coord = ldply(Org_IDs, rbind)
  
  Coord$.id = NULL

  Coord = Coord[!duplicated(Coord),]
  Coord = Coord[!duplicated(Coord$ID),]
  rownames(Coord) = Coord$ID
    
  Coord$ID = NULL
  
  #Order Start and End columns with Star < End
  rows = rownames(Coord[Coord$Start > Coord$End,])
  
  rows_switch = Coord[rows,c("Start","End")]
  
  Coord[rows, "Start"] = rows_switch[rows,"End"]  
  Coord[rows,"End"] = rows_switch[rows,"Start"]
  
  return(Coord)
  
}


Macaque_Coord = DF_Coord(Macaque_IDs)#1,941     #0.1 psi: 5,464 #0.05 psi and 5 reads: 
Chimp_Coord = DF_Coord(Chimp_IDs) #1,755     #0.1 psi: 5,504 #0.05 psi and 5 reads: 

Human_Coord = DF_Coord(Human_IDs) #4,052     #0.1 psi: 5,929? #0.05 psi and 5 reads: 


Baboon_Coord = DF_Coord(Baboon_IDs) #1,889     #0.1 psi: 5,622      #0.05 psi and 5 reads:

Lemur_Coord = DF_Coord(Lemur_IDs) #1,026     #0.1 psi: 3,148      #0.05 psi and 5 reads: 

Marmoset_Coord = DF_Coord(Marmoset_IDs) #1,397     #0.1 psi: 4,137    #0.05 psi and 5 reads: 
SquiMonkey_Coord = DF_Coord(SquiMonkey_IDs)#1,592     #0.1 psi: 4,773     #0.05 psi and 5 reads:
Macaca_Coord = DF_Coord(Macaca_IDs) #1,598     #0.1 psi: 4,829        #0.05 psi and 5 reads:



#Add IDs
Macaque_Coord$ID = rownames(Macaque_Coord)
Chimp_Coord$ID = rownames(Chimp_Coord)
Human_Coord$ID = rownames(Human_Coord)
Baboon_Coord$ID = rownames(Baboon_Coord)
Lemur_Coord$ID = rownames(Lemur_Coord)
Marmoset_Coord$ID = rownames(Marmoset_Coord)
SquiMonkey_Coord$ID = rownames(SquiMonkey_Coord)
Macaca_Coord$ID = rownames(Macaca_Coord)



#SAve Coords for liftover
path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/ForLiftOver"

#<file>_p1 is for files filtered as PSI >= 0.1
#<file>_p05 is files filtered as PSI >= 0.05 and 5 reads

write.table(Macaque_Coord, file=paste(path, "MacaqueCirc_Coord_p05.txt", sep="/"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Chimp_Coord, file=paste(path, "ChimpCirc_Coord_p05.txt", sep="/"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Human_Coord, file=paste(path, "HumanCirc_Coord_p05.txt", sep="/"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Baboon_Coord, file=paste(path, "BaboonCirc_Coord_p05.txt", sep="/"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Lemur_Coord, file=paste(path, "LemurCirc_Coord_p05.txt", sep="/"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Marmoset_Coord, file=paste(path, "MarmosetCirc_Coord_p05.txt", sep="/"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(SquiMonkey_Coord, file=paste(path, "SquiMonkeyCirc_Coord_p05.txt", sep="/"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(Macaca_Coord, file=paste(path, "MacacaCirc_Coord_p05.txt", sep="/"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
