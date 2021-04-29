

#To make coordinates file for LiftOver

#We are going to use the files <Primate>ThatAreOrth_PSI.txt that have the Whippet PSI output of 
#all genes that are shared between primates, we want to filter the genes according to coordinates to
#use them agains LiftOVer

path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/"

Human_Shared = read.delim(file=paste(path,"HumanThatAreOrth_PSI.txt", sep=""), header = TRUE, as.is = TRUE )
Chimp_Shared = read.delim(file=paste(path,"ChimpThatAreOrth_PSI.txt", sep=""), header = TRUE, as.is = TRUE)
Lemur_Shared = read.delim(file=paste(path,"LemurThatAreOrth_PSI.txt", sep=""), header = TRUE, as.is = TRUE)
Macaque_Shared = read.delim(file=paste(path,"MacaqueThatAreOrth_PSI.txt", sep=""), header = TRUE, as.is = TRUE)
Macaca_Shared = read.delim(file=paste(path,"MacacaThatAreOrth_PSI.txt", sep=""), header = TRUE, as.is = TRUE)
Marmoset_Shared = read.delim(file=paste(path,"MarmosetThatAreOrth_PSI.txt", sep=""), header = TRUE, as.is = TRUE)
SquirrelMonkey_Shared = read.delim(file=paste(path,"SquirrelMonkeyThatAreOrth_PSI.txt", sep=""), header = TRUE, as.is = TRUE)
Baboon_Shared = read.delim(file=paste(path,"BaboonThatAreOrth_PSI.txt", sep=""), header = TRUE, as.is = TRUE)

####This was added on September 22nd of 2020:

#To get conserved circRNAs of primates we are going to get the coordinates of CE that had a 10 or 5 % PSI as filter
#The script to make such files is: PSI_AS_Shared.R

#Files with a PSI_10
path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_10/"


Human_10 = read.delim(file=paste(path,"HumanThatAreOrth_PSI_10.txt", sep=""), header = TRUE, as.is = TRUE )
Chimp_10 = read.delim(file=paste(path,"ChimpThatAreOrth_PSI_10.txt", sep=""), header = TRUE, as.is = TRUE)
Lemur_10 = read.delim(file=paste(path,"LemurThatAreOrth_PSI_10.txt", sep=""), header = TRUE, as.is = TRUE)
Macaque_10 = read.delim(file=paste(path,"MacaqueThatAreOrth_PSI_10.txt", sep=""), header = TRUE, as.is = TRUE)
Macaca_10 = read.delim(file=paste(path,"MacacaThatAreOrth_PSI_10.txt", sep=""), header = TRUE, as.is = TRUE)
Marmoset_10 = read.delim(file=paste(path,"MarmosetThatAreOrth_PSI_10.txt", sep=""), header = TRUE, as.is = TRUE)
SquirrelMonkey_10 = read.delim(file=paste(path,"SquirrelMonkeyThatAreOrth_PSI_10.txt", sep=""), header = TRUE, as.is = TRUE)
Baboon_10 = read.delim(file=paste(path,"BaboonThatAreOrth_PSI_10.txt", sep=""), header = TRUE, as.is = TRUE)



#Files with a PSI 5
path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_5/"


Human_5 = read.delim(file=paste(path,"HumanThatAreOrth_PSI_5.txt", sep=""), header = TRUE, as.is = TRUE )
Chimp_5 = read.delim(file=paste(path,"ChimpThatAreOrth_PSI_5.txt", sep=""), header = TRUE, as.is = TRUE)
Lemur_5 = read.delim(file=paste(path,"LemurThatAreOrth_PSI_5.txt", sep=""), header = TRUE, as.is = TRUE)
Macaque_5 = read.delim(file=paste(path,"MacaqueThatAreOrth_PSI_5.txt", sep=""), header = TRUE, as.is = TRUE)
Macaca_5 = read.delim(file=paste(path,"MacacaThatAreOrth_PSI_5.txt", sep=""), header = TRUE, as.is = TRUE)
Marmoset_5 = read.delim(file=paste(path,"MarmosetThatAreOrth_PSI_5.txt", sep=""), header = TRUE, as.is = TRUE)
SquirrelMonkey_5 = read.delim(file=paste(path,"SquirrelMonkeyThatAreOrth_PSI_5.txt", sep=""), header = TRUE, as.is = TRUE)
Baboon_5 = read.delim(file=paste(path,"BaboonThatAreOrth_PSI_5.txt", sep=""), header = TRUE, as.is = TRUE)


#Make a function that subsets Coord and Gene columns
#Splits Coords in three columns (Chr, Start, End) and make and ID as Gene_Coord

Get_Coords <- function(Org_Shared) {
  
  Coord = Org_Shared[,c("Coord", "Gene")]
  
  Coord$Chr = gsub(":.*", "", Coord$Coord)
  Coord$Start = gsub("-.*","",gsub(".*:", "", Coord$Coord))
  Coord$End = gsub(".*-", "", Coord$Coord)
  
  Coord$ID = paste(Coord$Gene, Coord$Coord, sep="_")
  
  Coord$Coord = NULL
  Coord$Gene = NULL
  
  return(Coord)
  
}

Hum_Coord = Get_Coords(Human_Shared)
Chimp_Coord = Get_Coords(Chimp_Shared)
Lemur_Coord = Get_Coords(Lemur_Shared)
Macaque_Coord = Get_Coords(Macaque_Shared)
Macaca_Coord = Get_Coords(Macaca_Shared)
Marmoset_Coord = Get_Coords(Marmoset_Shared)
SquirrelMonkey_Coord = Get_Coords(SquirrelMonkey_Shared)
Baboon_Coord = Get_Coords(Baboon_Shared)

###############Get Coords for PSI >= 10

Hum_Coord_10 = Get_Coords(Human_10)

Chimp_Coord_10 = Get_Coords(Chimp_10)

Lemur_Coord_10 = Get_Coords(Lemur_10)
Macaque_Coord_10 = Get_Coords(Macaque_10)
Macaca_Coord_10 = Get_Coords(Macaca_10)
Marmoset_Coord_10 = Get_Coords(Marmoset_10)
SquirrelMonkey_Coord_10 = Get_Coords(SquirrelMonkey_10)
Baboon_Coord_10 = Get_Coords(Baboon_10)

################Get Coords for PSI >= 5

Hum_Coord_5 = Get_Coords(Human_5)

Chimp_Coord_5 = Get_Coords(Chimp_5)

Lemur_Coord_5 = Get_Coords(Lemur_5)
Macaque_Coord_5 = Get_Coords(Macaque_5)
Macaca_Coord_5 = Get_Coords(Macaca_5)
Marmoset_Coord_5 = Get_Coords(Marmoset_5)
SquirrelMonkey_Coord_5 = Get_Coords(SquirrelMonkey_5)
Baboon_Coord_5 = Get_Coords(Baboon_5)


#Save Coords (WITHOUT HEADER)
path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/"


write.table(Hum_Coord, file=paste(path, "Coord_Hum", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Chimp_Coord, file=paste(path, "Coord_Chimp", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Lemur_Coord, file=paste(path, "Coord_Lemur", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Macaque_Coord, file=paste(path, "Coord_Macaque", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Macaca_Coord, file=paste(path, "Coord_Macaca", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(SquirrelMonkey_Coord, file=paste(path, "Coord_SquirrelMonkey", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Marmoset_Coord, file=paste(path, "Coord_Marmoset", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Baboon_Coord, file=paste(path, "Coord_Baboon", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")


###Save Coords for Primates analysis (PSI filters of 10 and 5)
path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_10/"


write.table(Hum_Coord_10, file=paste(path, "Coord_Hum10", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Chimp_Coord_10, file=paste(path, "Coord_Chimp10", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Lemur_Coord_10, file=paste(path, "Coord_Lemur10", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Macaque_Coord_10, file=paste(path, "Coord_Macaque10", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Macaca_Coord_10, file=paste(path, "Coord_Macaca10", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(SquirrelMonkey_Coord_10, file=paste(path, "Coord_SquirrelMonkey10", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Marmoset_Coord_10, file=paste(path, "Coord_Marmoset10", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Baboon_Coord_10, file=paste(path, "Coord_Baboon10", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")


path="/home/gaby/lab_Garvan/Primates/Whippet_CircRNAs/AS_CE/PSI_5/"


write.table(Hum_Coord_5, file=paste(path, "Coord_Hum5", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Chimp_Coord_5, file=paste(path, "Coord_Chimp5", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Lemur_Coord_5, file=paste(path, "Coord_Lemur5", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Macaque_Coord_5, file=paste(path, "Coord_Macaque5", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Macaca_Coord_5, file=paste(path, "Coord_Macaca5", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(SquirrelMonkey_Coord_5, file=paste(path, "Coord_SquirrelMonkey5", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Marmoset_Coord_5, file=paste(path, "Coord_Marmoset5", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

write.table(Baboon_Coord_5, file=paste(path, "Coord_Baboon5", sep=""), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

