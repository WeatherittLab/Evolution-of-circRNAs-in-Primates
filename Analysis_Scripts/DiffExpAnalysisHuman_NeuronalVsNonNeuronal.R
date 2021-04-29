

#Script to perform differential expression analysis of human samples (from circRNAs dataset)

#The aim of this analysis is to find genes of conserved circRNAs are in the differential expressed genes

library(edgeR)
library(ggplot2)

#Load info of Conserved circRNAs
InfoConservedCircRNAs = read.delim(file="~/lab_Garvan/Primates/Whippet_CircRNAs/Intersect_ExonsGTFs_CircRNAs/ForLiftOver/LiftOverResults/PSI_10/GTF-CE/Conserved_BSJ/ConservedCircRNAs_Bob1.txt",
                                   header = TRUE, as.is = TRUE) #11,974



##Get Counts Matrix of Human samples
Human_Reads =read.delim(file="~/lab_Garvan/Primates/HumanGeneExp_Reads.txt",
            header = TRUE, as.is = TRUE) #50,927

data_info = data.frame("Samples" = colnames(Human_Reads),
                       "Tissue" = gsub("^(?:[^_]+_){1}([^_]+).*", "\\1", colnames(Human_Reads)),
                       "Replicate"= gsub(".*_","", colnames(Human_Reads)) )

###Make DGE List

HumanTissues_DGEList = DGEList(as.matrix(Human_Reads), genes= rownames(Human_Reads))

##Add tissue info to the DGEList
HumanTissues_DGEList$samples$group = data_info$Tissue


#filter low expressed genes
quantile(HumanTissues_DGEList$samples$lib.size)
#50%: 21,657,073

#Check how many are with 3 cpm in at least 2 samples
keep_cc = rowSums(cpm(HumanTissues_DGEList) > 3) >=2 #20,095

HumanTissues_DGEList = HumanTissues_DGEList[keep_cc, ,keep.lib.sizes = FALSE]

#CalcNormFators
HumanTissues_DGEList = calcNormFactors(HumanTissues_DGEList)

col.tissue = rainbow(length(levels(data_info$Tissue)))[data_info$Tissue]
plotMDS(cpm(HumanTissues_DGEList$counts, prior.count = 2, log = TRUE), cex = 0.8, col=col.tissue)

HumanTissues_DGEList$samples$group = factor(HumanTissues_DGEList$samples$group)



###Make design matrices and estimate dispersion Add Replicate info as a covariate
#group = levels(factor(data_info$Tissue))
#batch_rep = factor(data_info$Replicates)

design = model.matrix(~0+group, data = HumanTissues_DGEList$samples)
colnames(design) = levels(HumanTissues_DGEList$samples$group)

HumanTissues_DGEList_filt = estimateDisp(HumanTissues_DGEList, design = design)

plotMDS(cpm(HumanTissues_DGEList_filt$counts, prior.count = 2, log = TRUE), cex=0.8, col= col.tissue)


#Apply the general linear model
HumanTissues_fit = glmQLFit(HumanTissues_DGEList_filt, design)

##Define contrast of Neuronal Samples vs Non_Neuronal
NeuronalVsNon = makeContrasts("NeuronalVsNon" = ((Brain + Cerebellum + FrontalCortex)/3) - ((Colon + Heart + Liver + Lung + SkeletalMuscle + Spleen)/6), levels = design)

QLF.Neuron = glmQLFTest(HumanTissues_fit, contrast = NeuronalVsNon[,"NeuronalVsNon"])


##Add FDR 
QLF.Neuron$table$FDR = p.adjust(QLF.Neuron$table$PValue, method = "BH")


#Number of Differentially expressed genes
summary(decideTests(QLF.Neuron, lfc = log2(1.5),
                    adjust.method = "BH", p.value = 0.05))

#down: 4,383
#NotSig: 11,278
#up: 4,434

DecideTest = decideTests(QLF.Neuron, lfc = log2(1.5),
                         adjust.method = "BH", p.value = 0.05)

UpDecideTest = DecideTest[DecideTest == 1,]
colnames(UpDecideTest) = "NeuronalVsNon_Up"

DownDecideTest = DecideTest[DecideTest == -1,]
colnames(DownDecideTest) = "NeuronalVsNon_Down"

length(intersect(rownames(UpDecideTest), InfoConservedCircRNAs[grep("^Conserved$", InfoConservedCircRNAs$Category), "Ensembl_ID"]))  #143


###for hypergeometric test
tot_size = nrow(QLF.Neuron$table)
#sample = nrow(UpDecideTest)
tot_consv = length(intersect(rownames(QLF.Neuron$table), InfoConservedCircRNAs[grep("^Conserved$", InfoConservedCircRNAs$Category), "Ensembl_ID"]))
#sample_consv = length(intersect(rownames(UpDecideTest), InfoConservedCircRNAs[grep("^Conserved$", InfoConservedCircRNAs$Category), "Ensembl_ID"]))

#phyper(sample_consv, sample, tot_size-sample, tot_consv)

#phyper(tot_consv-1, sample, tot_size-sample, tot_consv)


#Make hypergeometric test taking all diff exp genes
sample_all = sum(nrow(UpDecideTest), nrow(DownDecideTest))
sample_consv = length(intersect(c(rownames(UpDecideTest), rownames(DownDecideTest)), InfoConservedCircRNAs[grep("^Conserved$", InfoConservedCircRNAs$Category), "Ensembl_ID"]))


phyper(sample_consv,  sample_all, tot_size-sample_all, tot_consv, lower.tail = FALSE)


###To make volcano plot. In this volcanoplot we want to show colored the UpRegulated genes (4,434) and the genes from consv. circRNAs

###Add info to color to the QLF.Neuron
QLF.Neuron$table$diffexpressed <- "Not Up Regulated"
QLF.Neuron$table$diffexpressed[QLF.Neuron$table$logFC >= log2(1.5) & QLF.Neuron$table$FDR <= 0.05] <- "Up Regulated"

QLF.Neuron$table[rownames(QLF.Neuron$table) %in% InfoConservedCircRNAs[grep("^Conserved$", InfoConservedCircRNAs$Category), "Ensembl_ID"],"diffexpressed"] <- "Conserved"


#Make a vector of colors of interest
mycolors = c("darkgoldenrod3", "darkmagenta", "skyblue1")
names(mycolors) = c("Conserved", "Up Regulated", "Not Up Regulated")


volcano = ggplot(data=QLF.Neuron$table, aes(x= logFC, y = -log10(PValue), col = diffexpressed )) +
          geom_point(size=4) + 
          theme_bw() + 
          theme(legend.title = element_text(size=20), legend.text = element_text(size = 20)) +
          geom_vline(xintercept =c(-0.6, 0.6), col="red") +
          geom_hline(yintercept = -log10(0.05), col="red" ) +
          scale_color_manual(values = mycolors)










