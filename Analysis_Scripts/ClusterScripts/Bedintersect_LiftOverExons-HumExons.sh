# -S /bin/sh

#Intersect Lifted Over Exons with human exons (All exons)

#Exons GTF files path
path="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords/PSI_10/LiftOver_AllExons"

###Files Lifted Over
baboon="LiftOver_Baboon2Hum_ExonGTF-CE10.txt"
chimp="LiftOver_Chimp2Hum_ExonGTF-CE10.txt"
lemur="LiftOver_Lemur2Hum_ExonGTF-CE10.txt"
macaca="LiftOver_Macaca2Hum_ExonGTF-CE10.txt"
macaque="LiftOver_Macaque2Hum_ExonGTF-CE10.txt"
marmoset="LiftOver_Marmoset2Hum_ExonGTF-CE10.txt"
squimonkey="LiftOver_SquirrelMonkey2Hum_ExonGTF-CE10.txt"

###Human Exons
path_human="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords/PSI_10/"

human_exons="Human_Exons_GTF-CE10.filt.bed"


###Comand per primate

###BABOON

BabInt="bedtools intersect -loj -a $path/$baboon -b $path_human/$human_exons > $path/Baboon_HumanExons_GTF-CE10.txt"

echo "$BabInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_BaboonHumExons -q short.q -V


###CHIMP

ChimpInt="bedtools intersect -loj -a $path/$chimp -b $path_human/$human_exons > $path/Chimp_HumanExons_GTF-CE10.txt"

echo "$ChimpInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_ChimpHumExons -q short.q -V


###LEMUR

LemurInt="bedtools intersect -loj -a $path/$lemur -b $path_human/$human_exons > $path/Lemur_HumanExons_GTF-CE10.txt"

echo "$LemurInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_LemurHumExons -q short.q -V


####MACACA

MacaInt="bedtools intersect -loj -a $path/$macaca -b $path_human/$human_exons > $path/Macaca_HumanExons_GTF-CE10.txt"

echo "$MacaInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_MacacaHumExons -q short.q -V


####MACAQUE

MacqInt="bedtools intersect -loj -a $path/$macaque -b $path_human/$human_exons > $path/Macaque_HumanExons_GTF-CE10.txt"

echo "$MacqInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_MacaqueHumExons -q short.q -V


####MARMOSET

MarmoInt="bedtools intersect -loj -a $path/$marmoset -b $path_human/$human_exons > $path/Marmoset_HumanExons_GTF-CE10.txt"

echo "$MarmoInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_MarmosetHumExons -q short.q -V


#####SQUIRREL MONKEY

SquiMonkInt="bedtools intersect -loj -a $path/$squimonkey -b $path_human/$human_exons > $path/SquirrelMonkey_HumanExons_GTF-CE10.txt"

echo "$SquiMonkInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_SquirrelMonkeyHumExons -q short.q -V
