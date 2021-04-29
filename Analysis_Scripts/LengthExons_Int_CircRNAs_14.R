
#Script to calculate length of exons of circRNAs of primates

#Files from bedintersect between exon coordinates and circRNA coordinates from Whippet (Primates data)

###NOT:

#Baboon = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Baboon_ExonsCircRNA.bed",
#                    header = FALSE, as.is = TRUE)
#Human = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Human_ExonsCircRNA.bed",
#                   header = FALSE, as.is = TRUE)
#Lemur = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Lemur_ExonsCircRNA.bed",
#                   header = FALSE, as.is = TRUE)
#Macaca = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Macaca_ExonsCircRNA.bed",
#                    header = FALSE, as.is = TRUE)
#Macaque = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Macaque_ExonsCircRNA.bed",
#                     header = FALSE, as.is = TRUE)
#Marmoset = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Marmoset_ExonsCircRNA.bed",
#                      header = FALSE, as.is = TRUE)
#SquiMonkey = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/SquiMonkey_ExonsCircRNA.bed",
#                        header = FALSE, as.is = TRUE)


#This was modified on September 29th:
#Read <Primate>_ExonsCircRNA_<PSIfilt>.bed files to calculate the length of the exons that frm circRNAs of primates and have certain PSI value

################NOT

###PSI 10

#Baboon10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Baboon_ExonsCircRNA_10.bed",
#                    header = FALSE, as.is = TRUE)
#Human10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Human_ExonsCircRNA_10.bed",
#                   header = FALSE, as.is = TRUE)
#Lemur10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Lemur_ExonsCircRNA_10.bed",
#                   header = FALSE, as.is = TRUE)
#Macaca10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Macaca_ExonsCircRNA_10.bed",
#                    header = FALSE, as.is = TRUE)
#Macaque10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Macaque_ExonsCircRNA_10.bed",
#                     header = FALSE, as.is = TRUE)
#Marmoset10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Marmoset_ExonsCircRNA_10.bed",
#                      header = FALSE, as.is = TRUE)
#SquiMonkey10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/SquiMonkey_ExonsCircRNA_10.bed",
#                        header = FALSE, as.is = TRUE)


###PSI 5

#Baboon5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Baboon_ExonsCircRNA_5.bed",
#                    header = FALSE, as.is = TRUE)
#Human5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Human_ExonsCircRNA_5.bed",
#                   header = FALSE, as.is = TRUE)
#Lemur5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Lemur_ExonsCircRNA_5.bed",
#                   header = FALSE, as.is = TRUE)
#Macaca5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Macaca_ExonsCircRNA_5.bed",
#                    header = FALSE, as.is = TRUE)
#Macaque5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Macaque_ExonsCircRNA_5.bed",
#                     header = FALSE, as.is = TRUE)
#Marmoset5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Marmoset_ExonsCircRNA_5.bed",
#                      header = FALSE, as.is = TRUE)
#SquiMonkey5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/SquiMonkey_ExonsCircRNA_5.bed",
#                        header = FALSE, as.is = TRUE)

######

##This script was modified on 7th of October of 2020. The coordinates we are using come from Whippet CE MERGED and GTF (intersected)

#PSI 10

Baboon10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Baboon_ExonsCircRNA_GTF-CE10.bed",
                    header = FALSE, as.is = TRUE)
Human10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Human_ExonsCircRNA_GTF-CE10.bed",
                   header = FALSE, as.is = TRUE)
Lemur10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Lemur_ExonsCircRNA_GTF-CE10.bed",
                   header = FALSE, as.is = TRUE)
Macaca10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Macaca_ExonsCircRNA_GTF-CE10.bed",
                    header = FALSE, as.is = TRUE)
Macaque10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Macaque_ExonsCircRNA_GTF-CE10.bed",
                     header = FALSE, as.is = TRUE)
Marmoset10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Marmoset_ExonsCircRNA_GTF-CE10.bed",
                      header = FALSE, as.is = TRUE)
SquiMonkey10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/SquiMonkey_ExonsCircRNA_GTF-CE10.bed",
                        header = FALSE, as.is = TRUE)

Chimp10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Chimp_ExonsCircRNA_GTF-CE10.bed",
                     header = FALSE, as.is = TRUE)

#PSI 5

#Baboon5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Baboon_ExonsCircRNA_GTF-CE5.bed",
#                    header = FALSE, as.is = TRUE)
#Human5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Human_ExonsCircRNA_GTF-CE5.bed",
#                   header = FALSE, as.is = TRUE)
#Lemur5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Lemur_ExonsCircRNA_GTF-CE5.bed",
#                   header = FALSE, as.is = TRUE)
#Macaca5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Macaca_ExonsCircRNA_GTF-CE5.bed",
#                    header = FALSE, as.is = TRUE)
#Macaque5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Macaque_ExonsCircRNA_GTF-CE5.bed",
#                     header = FALSE, as.is = TRUE)
#Marmoset5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Marmoset_ExonsCircRNA_GTF-CE5.bed",
#                      header = FALSE, as.is = TRUE)
#SquiMonkey5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/SquiMonkey_ExonsCircRNA_GTF-CE5.bed",
#                        header = FALSE, as.is = TRUE)

#Chimp5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Chimp_ExonsCircRNA_GTF-CE5.bed",
#                     header = FALSE, as.is = TRUE)
###Make a function that:
#Filter coordinates that match with circRNA coords
#Add colnames
#Calculate length of exons
cols = c("Chr_Exon", "Start_Exon", "End_Exon", "ID_Exon", "Chr_Circ", "Start_Circ", "End_Circ", "ID_Circ")

Filter_ExonCirc <- function(Primate, columns) {
  
  #Filter coordinates that match with circRNA coords
  PrimateFilt = Primate[grep("\\.", Primate$V8, invert=TRUE),]
  
  #Add colnames
  colnames(PrimateFilt) = columns
  
  #order
  PrimateFilt= PrimateFilt[order(PrimateFilt$ID_Exon),]
  
  #remove duplicated
  PrimateFilt = PrimateFilt[!duplicated(PrimateFilt),]
  
  #Calculate length of exons
  PrimateFilt$Length_Exon = PrimateFilt$End_Exon - PrimateFilt$Start_Exon
  
  return(PrimateFilt)
  
}

#Previous data
#BaboonFilt = Filter_ExonCirc(Baboon, cols) #199,220
#HumanFilt = Filter_ExonCirc(Human, cols) #1,182,977
#LemurFilt = Filter_ExonCirc(Lemur, cols) # 114,620
#MacacaFilt = Filter_ExonCirc(Macaca, cols) # 201,943
#MacaqueFilt = Filter_ExonCirc(Macaque, cols) # 239,041 
#MarmosetFilt = Filter_ExonCirc(Marmoset, cols) #187,704 
#SquiMonkeyFilt = Filter_ExonCirc(SquiMonkey, cols) #194,718


#PSI_10 (since Oct 7 2020 is GTF-CE)
BaboonFilt_10 = Filter_ExonCirc(Baboon10, cols) #59,271 (PSI) #62,959
HumanFilt_10 = Filter_ExonCirc(Human10, cols) #161,523 (PSI) #238,795
LemurFilt_10 = Filter_ExonCirc(Lemur10, cols) #40,804 (PSI) #43,258
MacacaFilt_10 = Filter_ExonCirc(Macaca10, cols) #57,828 (PSI) #61,910
MacaqueFilt_10 = Filter_ExonCirc(Macaque10, cols) #70,105 (PSI)  #75,514
MarmosetFilt_10 = Filter_ExonCirc(Marmoset10, cols) #58,414 (PSI)  #65,068
SquiMonkeyFilt_10 = Filter_ExonCirc(SquiMonkey10, cols) #58,953 (PSI) #62,005

ChimpFilt_10 = Filter_ExonCirc(Chimp10, cols) # (PSI) #67,844

#PSI_5
#BaboonFilt_5 = Filter_ExonCirc(Baboon5, cols) #59,294  (PSI) #62,982
#HumanFilt_5 = Filter_ExonCirc(Human5, cols) #161,817 (PSI) #239,129
#LemurFilt_5 = Filter_ExonCirc(Lemur5, cols) #40,826 (PSI) #43,280
#MacacaFilt_5 = Filter_ExonCirc(Macaca5, cols) #57,875 (PSI) #61,957
#MacaqueFilt_5 = Filter_ExonCirc(Macaque5, cols) #70,186 (PSI) #75,587
#MarmosetFilt_5 = Filter_ExonCirc(Marmoset5, cols) #58,508 (PSI) #65,162
#SquiMonkeyFilt_5 = Filter_ExonCirc(SquiMonkey5, cols) #58,995 (PSI) #62,047

#ChimpFilt_5 = Filter_ExonCirc(Chimp5, cols) # (PSI) #

#Save info

###Exons from GTF

#write.table(HumanFilt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Human_Int_ExonCirc.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(BaboonFilt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Baboon_Int_ExonCirc.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(LemurFilt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Lemur_Int_ExonCirc.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MacacaFilt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Macaca_Int_ExonCirc.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MacaqueFilt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Macaque_Int_ExonCirc.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MarmosetFilt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Marmoset_Int_ExonCirc.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(SquiMonkeyFilt, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/SquiMonkey_Int_ExonCirc.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


###Exons from Whippet CE info

#PSI_10
#write.table(HumanFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Human_Int_ExonCirc_10.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(BaboonFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Baboon_Int_ExonCirc_10.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(LemurFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Lemur_Int_ExonCirc_10.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MacacaFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Macaca_Int_ExonCirc_10.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MacaqueFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Macaque_Int_ExonCirc_10.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MarmosetFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/Marmoset_Int_ExonCirc_10.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(SquiMonkeyFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/SquiMonkey_Int_ExonCirc_10.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


#PSI_5
#write.table(HumanFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Human_Int_ExonCirc_5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(BaboonFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Baboon_Int_ExonCirc_5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(LemurFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Lemur_Int_ExonCirc_5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MacacaFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Macaca_Int_ExonCirc_5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MacaqueFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Macaque_Int_ExonCirc_5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MarmosetFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/Marmoset_Int_ExonCirc_5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(SquiMonkeyFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/SquiMonkey_Int_ExonCirc_5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")



####PSI GTF-CE

#PSI_10
write.table(HumanFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Human_Int_ExonCirc_GTF-CE10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(BaboonFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Baboon_Int_ExonCirc_GTF-CE10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(LemurFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Lemur_Int_ExonCirc_GTF-CE10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(MacacaFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Macaca_Int_ExonCirc_GTF-CE10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(MacaqueFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Macaque_Int_ExonCirc_GTF-CE10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(MarmosetFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Marmoset_Int_ExonCirc_GTF-CE10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(SquiMonkeyFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/SquiMonkey_Int_ExonCirc_GTF-CE10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(ChimpFilt_10, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_10/GTF-CE/Chimp_Int_ExonCirc_GTF-CE10.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#PSI_5
#write.table(HumanFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Human_Int_ExonCirc_GTF-CE5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(BaboonFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Baboon_Int_ExonCirc_GTF-CE5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(LemurFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Lemur_Int_ExonCirc_GTF-CE5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MacacaFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Macaca_Int_ExonCirc_GTF-CE5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MacaqueFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Macaque_Int_ExonCirc_GTF-CE5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(MarmosetFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Marmoset_Int_ExonCirc_GTF-CE5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(SquiMonkeyFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/SquiMonkey_Int_ExonCirc_GTF-CE5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#write.table(ChimpFilt_5, file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/PSI_5/GTF-CE/Chimp_Int_ExonCirc_GTF-CE5.txt",
#            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

