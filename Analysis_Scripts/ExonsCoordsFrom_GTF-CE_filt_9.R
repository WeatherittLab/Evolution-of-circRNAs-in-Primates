

#Script to filter coordinates from whippet CE intersected with GTFs coordinates.
#the result of the bedintersect of these files from primates data could have duplicated coordinates
#as the IDs of the GTF file is like: <Ensemble_transcID>_Coord and there could be many transcripts with the same exon...


library(gread)

#https://rdrr.io/github/openanalytics/gread/
#https://rdrr.io/github/openanalytics/gread/man/extract.html


###GTFs info

Human = read_format("~/lab_Garvan/GTFs/Homo_sapiens.GRCh38.97.gtf")
Baboon = read_format("~/lab_Garvan/GTFs/Papio_anubis.Panu_3.0.95.gtf")

Lemur = read_format("~/lab_Garvan/GTFs/micMur2.ensGene.gtf")

Macaca = read_format("~/lab_Garvan/GTFs/Macaca_fascicularis.Macaca_fascicularis_5.0.99.gtf")

Macaque = read_format("~/lab_Garvan/GTFs/Macaca_mulatta.Mmul_8.0.1.97.gtf")


Marmoset = read_format("~/lab_Garvan/GTFs/calJac3.ensGene.gtf")
SquirrelMonkey = read_format("~/lab_Garvan/GTFs/Saimiri_boliviensis_boliviensis.SaiBol1.0.99.gtf")

Chimp = read_format("~/lab_Garvan/GTFs/Pan_troglodytes.Pan_tro_3.0.99.gtf")


#### Files from whippet CE intersected with GTFs coordinates

#PSI_10

Human_10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Human_Exons_GTF-CE10.bed",
                      header = FALSE, as.is = TRUE)

Baboon_10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Baboon_Exons_GTF-CE10.bed",
                       header = FALSE, as.is = TRUE)

Lemur_10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Lemur_Exons_GTF-CE10.bed",
                      header = FALSE, as.is = TRUE)

Macaca_10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Macaca_Exons_GTF-CE10.bed",
                       header = FALSE, as.is = TRUE)

Macaque_10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Macaque_Exons_GTF-CE10.bed",
                        header = FALSE, as.is = TRUE)

Marmoset_10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Marmoset_Exons_GTF-CE10.bed",
                         header = FALSE, as.is = TRUE)

SquirrelMonkey_10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/SquirrelMonkey_Exons_GTF-CE10.bed",
                               header = FALSE, as.is = TRUE)

Chimp_10 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Chimp_Exons_GTF-CE10.bed",
                      header = FALSE, as.is = TRUE)

        
#PSI_5

#Human_5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Human_Exons_GTF-CE5.bed",
#                     header = FALSE, as.is = TRUE)

#Baboon_5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Baboon_Exons_GTF-CE5.bed",
#                      header = FALSE, as.is = TRUE)

#Lemur_5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Lemur_Exons_GTF-CE5.bed",
#                     header = FALSE, as.is = TRUE)

#Macaca_5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Macaca_Exons_GTF-CE5.bed",
#                      header = FALSE, as.is = TRUE)

#Macaque_5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Macaque_Exons_GTF-CE5.bed",
#                       header = FALSE, as.is = TRUE)

#Marmoset_5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Marmoset_Exons_GTF-CE5.bed",
#                        header = FALSE, as.is = TRUE)

#SquirrelMonkey_5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/SquiMonkey_Exons_GTF-CE5.bed",
#                              header = FALSE, as.is = TRUE)

#Chimp_5 = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Chimp_Exons_GTF-CE5.bed",
#                      header = FALSE, as.is = TRUE)

###Add colnames to <Primate>_<PSI>
cols = c("Chr", "Start", "End", "Exon_ID")

colnames(Human_10) = cols
#colnames(Human_5) = cols
colnames(Baboon_10) = cols
#colnames(Baboon_5) = cols
colnames(Lemur_10) = cols
#colnames(Lemur_5) = cols
colnames(Macaca_10) = cols
#colnames(Macaca_5) = cols
colnames(Macaque_10) = cols
#colnames(Macaque_5) = cols
colnames(Marmoset_10) = cols
#colnames(Marmoset_5) = cols
colnames(SquirrelMonkey_10) = cols
#colnames(SquirrelMonkey_5) = cols


colnames(Chimp_10) = cols
#colnames(Chimp_5)= cols

###Chroms info
Human_chroms = c(as.character(1:22), "X", "Y", "MT")
Baboon_chroms = c(as.character(1:20), "X")
Lemur_chroms = levels(seqnames(Lemur))
Macaca_chroms = c(as.character(1:20), "X", "MT")
Macaque_chroms = c(as.character(1:20), "X", "Y", "MT")

Marmoset_chroms = levels(seqnames(Marmoset))

SquirrelMonkey_chroms = levels(seqnames(SquirrelMonkey))
SquirrelMonkey_chroms = SquirrelMonkey_chroms[grep("JH", SquirrelMonkey_chroms)]


Chimp_chroms = c(as.character(1:22),"X", "Y", "MT")


#function to get Exons coords from GTF file

Get_ExonsCoords <- function(GTF, chroms) {
  
  Exons = extract(GTF, feature = "exon", type = "default", 
                  transcript_id = "transcript_id", gene_id = "gene_id")
  
  #Save exons info of coordinates in chromosomes 1:22,sexual and MT
  Exons = Exons[seqnames(Exons) %in% chroms] #1,362,975
  
  #Retrieve info to get a bed file (chromosome, start, end, ID)
  
  #1 Get coordinates
  Ranges = ranges(Exons)
  
  #2 From Ranges object get the Start and End
  BedFile = data.frame("Chrom" = as.character(seqnames(Exons)),
                       "Start" = start(ranges(Exons)),
                       "End" = end(ranges(Exons)))
  
  Transcr_IDs = paste(BedFile$Chrom, paste(BedFile$Start, BedFile$End, sep = "-"), sep = ":")
  Gene_IDs = Transcr_IDs
  
  
  BedFile$Transcr_ID = paste(Exons$transcript_id, Transcr_IDs, sep="_")
  BedFile$Gene_ID = paste(Exons$gene_id, Gene_IDs, sep="_")
  
  
  return(BedFile)
  
}


Human_Exons = Get_ExonsCoords(Human, Human_chroms)
Baboon_Exons = Get_ExonsCoords(Baboon, Baboon_chroms)
Lemur_Exons = Get_ExonsCoords(Lemur, Lemur_chroms)
Macaca_Exons = Get_ExonsCoords(Macaca, Macaca_chroms)
Macaque_Exons = Get_ExonsCoords(Macaque, Macaque_chroms)
SquirrelMonkey_Exons = Get_ExonsCoords(SquirrelMonkey, SquirrelMonkey_chroms)
Marmoset_Exons = Get_ExonsCoords(Marmoset, Marmoset_chroms)

Chimp_Exons = Get_ExonsCoords(Chimp, Chimp_chroms)


#For SquirrelMonkey remove the ".#" of Chrom
SquirrelMonkey_Exons$Chrom = gsub("\\..*","", SquirrelMonkey_Exons$Chrom)


####Add Gene_ID from <Primate>_Exons to <Primate>_<PSI> 

Human_10$ExonGene_ID = Human_Exons[match(Human_10$Exon_ID, Human_Exons$Transcr_ID), "Gene_ID"]
#Human_5$ExonGene_ID = Human_Exons[match(Human_5$Exon_ID, Human_Exons$Transcr_ID), "Gene_ID"]

Baboon_10$ExonGene_ID = Baboon_Exons[match(Baboon_10$Exon_ID, Baboon_Exons$Transcr_ID), "Gene_ID"]
#Baboon_5$ExonGene_ID = Baboon_Exons[match(Baboon_5$Exon_ID, Baboon_Exons$Transcr_ID), "Gene_ID"]

Lemur_10$ExonGene_ID = Lemur_Exons[match(Lemur_10$Exon_ID, Lemur_Exons$Transcr_ID), "Gene_ID"]
#Lemur_5$ExonGene_ID = Lemur_Exons[match(Lemur_5$Exon_ID, Lemur_Exons$Transcr_ID), "Gene_ID"]

Macaca_10$ExonGene_ID = Macaca_Exons[match(Macaca_10$Exon_ID, Macaca_Exons$Transcr_ID), "Gene_ID"]
#Macaca_5$ExonGene_ID = Macaca_Exons[match(Macaca_5$Exon_ID, Macaca_Exons$Transcr_ID), "Gene_ID"]

Macaque_10$ExonGene_ID = Macaque_Exons[match(Macaque_10$Exon_ID, Macaque_Exons$Transcr_ID), "Gene_ID"]
#Macaque_5$ExonGene_ID = Macaque_Exons[match(Macaque_5$Exon_ID, Macaque_Exons$Transcr_ID), "Gene_ID"]

SquirrelMonkey_10$ExonGene_ID = SquirrelMonkey_Exons[match(SquirrelMonkey_10$Exon_ID, SquirrelMonkey_Exons$Transcr_ID), "Gene_ID"]
#SquirrelMonkey_5$ExonGene_ID = SquirrelMonkey_Exons[match(SquirrelMonkey_5$Exon_ID, SquirrelMonkey_Exons$Transcr_ID), "Gene_ID"]

Marmoset_10$ExonGene_ID = Marmoset_Exons[match(Marmoset_10$Exon_ID, Marmoset_Exons$Transcr_ID), "Gene_ID"]
#Marmoset_5$ExonGene_ID = Marmoset_Exons[match(Marmoset_5$Exon_ID, Marmoset_Exons$Transcr_ID), "Gene_ID"]

Chimp_10$ExonGene_ID = Chimp_Exons[match(Chimp_10$Exon_ID, Chimp_Exons$Transcr_ID), "Gene_ID"]
#Chimp_5$ExonGene_ID = Chimp_Exons[match(Chimp_5$Exon_ID, Chimp_Exons$Transcr_ID), "Gene_ID"]


#Filter duplicates
#Human_5 = Human_5[which(!duplicated(Human_5[,c("Chr", "Start", "End", "ExonGene_ID")])),]
Human_10 = Human_10[which(!duplicated(Human_10[,c("Chr", "Start", "End", "ExonGene_ID")])),]

#Baboon_5 = Baboon_5[which(!duplicated(Baboon_5[,c("Chr", "Start", "End", "ExonGene_ID")])),]
Baboon_10 = Baboon_10[which(!duplicated(Baboon_10[,c("Chr", "Start", "End", "ExonGene_ID")])),]
  
#Lemur_5 =  Lemur_5[which(!duplicated(Lemur_5[,c("Chr", "Start", "End", "ExonGene_ID")])),]
Lemur_10 = Lemur_10[which(!duplicated(Lemur_10[,c("Chr", "Start", "End", "ExonGene_ID")])),]
  
#Macaca_5 =  Macaca_5[which(!duplicated(Macaca_5[,c("Chr", "Start", "End", "ExonGene_ID")])),]
Macaca_10 = Macaca_10[which(!duplicated(Macaca_10[,c("Chr", "Start", "End", "ExonGene_ID")])),]
  
#Macaque_5 = Macaque_5[which(!duplicated(Macaque_5[,c("Chr", "Start", "End", "ExonGene_ID")])),]
Macaque_10 = Macaque_10[which(!duplicated(Macaque_10[,c("Chr", "Start", "End", "ExonGene_ID")])),]
  
#SquirrelMonkey_5 = SquirrelMonkey_5[which(!duplicated(SquirrelMonkey_5[,c("Chr", "Start", "End", "ExonGene_ID")])),]
SquirrelMonkey_10 = SquirrelMonkey_10[which(!duplicated(SquirrelMonkey_10[,c("Chr", "Start", "End", "ExonGene_ID")])),]
  
#Marmoset_5 = Marmoset_5[which(!duplicated(Marmoset_5[,c("Chr", "Start", "End", "ExonGene_ID")])),]
Marmoset_10 = Marmoset_10[which(!duplicated(Marmoset_10[,c("Chr", "Start", "End", "ExonGene_ID")])),]

#Chimp_5 = Chimp_5[which(!duplicated(Chimp_5[,c("Chr","Start","End","ExonGene_ID")])),]
Chimp_10 = Chimp_10[which(!duplicated(Chimp_10[,c("Chr","Start","End","ExonGene_ID")])),]


#Save <Primate>_<PSI> with Chr, Start, End, ExonGene_ID

#PSI_5

write.table(Human_5[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Human_Exons_GTF-CE5.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Baboon_5[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Baboon_Exons_GTF-CE5.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Lemur_5[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Lemur_Exons_GTF-CE5.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Macaca_5[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Macaca_Exons_GTF-CE5.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Macaque_5[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Macaque_Exons_GTF-CE5.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Marmoset_5[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Marmoset_Exons_GTF-CE5.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(SquirrelMonkey_5[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/SquiMonkey_Exons_GTF-CE5.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Chimp_5[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_5/Chimp_Exons_GTF-CE5.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")


  #PSI_10

write.table(Human_10[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Human_Exons_GTF-CE10.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Baboon_10[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Baboon_Exons_GTF-CE10.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Lemur_10[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Lemur_Exons_GTF-CE10.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Macaca_10[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Macaca_Exons_GTF-CE10.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Macaque_10[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Macaque_Exons_GTF-CE10.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(Marmoset_10[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Marmoset_Exons_GTF-CE10.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(SquirrelMonkey_10[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/SquiMonkey_Exons_GTF-CE10.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


write.table(Chimp_10[,c(1,2,3,5)], file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/Exons_GTF-CE/PSI_10/Chimp_Exons_GTF-CE10.filt.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
