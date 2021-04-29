
#Script to get exons info from GTF file


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

  
###Chroms info
Human_chroms = c(as.character(1:22), "X", "Y", "MT")
Baboon_chroms = c(as.character(1:20), "X")
Lemur_chroms = levels(seqnames(Lemur))
Macaca_chroms = c(as.character(1:20), "X", "MT")
Macaque_chroms = c(as.character(1:20), "X", "Y", "MT")

Marmoset_chroms = levels(seqnames(Marmoset))

SquirrelMonkey_chroms = levels(seqnames(SquirrelMonkey))
SquirrelMonkey_chroms = temp[grep("JH", temp)]

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

  IDs = paste(BedFile$Chrom, paste(BedFile$Start, BedFile$End, sep = "-"), sep = ":")

  BedFile$ID = paste(Exons$transcript_id, IDs, sep="_")
  
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

##Save data
  path="/home/gaby/lab_Garvan/GTFs/ExonsInfo"

write.table(Human_Exons, file = paste(path, "Human_Exons.txt", sep = "/"), 
            col.names = FALSE, quote = FALSE, sep ="\t", row.names = FALSE)

write.table(Baboon_Exons, paste(path, "Baboon_Exons.txt", sep = "/"), 
            col.names = FALSE, quote = FALSE, sep ="\t", row.names = FALSE)
  
write.table(Lemur_Exons, paste(path, "Lemur_Exons.txt", sep = "/"), 
            col.names = FALSE, quote = FALSE, sep ="\t", row.names = FALSE)

write.table(Macaca_Exons, paste(path, "Macaca_Exons.txt", sep = "/"), 
            col.names = FALSE, quote = FALSE, sep ="\t", row.names = FALSE)

write.table(Macaque_Exons, paste(path, "Macaque_Exons.txt", sep = "/"), 
            col.names = FALSE, quote = FALSE, sep ="\t", row.names = FALSE)

write.table(Marmoset_Exons, paste(path, "Marmoset_Exons.txt", sep = "/"),
            col.names = FALSE, quote = FALSE, sep ="\t", row.names = FALSE)

write.table(SquirrelMonkey_Exons, paste(path, "SquirrelMonkey_Exons.txt", sep = "/"), 
            col.names = FALSE, quote = FALSE, sep ="\t", row.names = FALSE)

write.table(Chimp_Exons, paste(path, "Chimp_Exons.txt", sep="/"),
            col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)



####Check that coords match to whippet circRNA bed file