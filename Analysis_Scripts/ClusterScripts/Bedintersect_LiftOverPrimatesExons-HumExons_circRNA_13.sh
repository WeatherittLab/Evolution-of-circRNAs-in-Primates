# -S /bin/sh

#Intersect Lifted Over Exons with human exons (primates circRNAs)

#Exons GTF files path
path="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords/PSI_10/"

###Files Lifted Over
baboon="LiftOver_Baboon2Hum_ExonCirc_GTF-CE10.txt"
chimp="LiftOver_Chimp2Hum_ExonCirc_GTF-CE10.txt"
lemur="LiftOver_Lemur2Hum_ExonCirc_GTF-CE10.txt"
macaca="LiftOver_Macaca2Hum_ExonCirc_GTF-CE10.txt"
macaque="LiftOver_Macaque2Hum_ExonCirc_GTF-CE10.txt"
marmoset="LiftOver_Marmoset2Hum_ExonCirc_GTF-CE10.txt"
squimonkey="LiftOver_SquirrelMonkey2Hum_ExonCirc_GTF-CE10.txt"

###Human Exons
human_exons="Human_ExonCoords_4LiftOver_GTF-CE10.bed"


###Comand per primate

###BABOON

#BabInt="bedtools intersect -loj -a $path/$baboon -b $path/$human_exons > $path/Baboon_Human_ExonsCircs_GTF-CE10.txt"

#echo "$BabInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_BaboonHumExonsCircs -q short.q -V


###CHIMP

#ChimpInt="bedtools intersect -loj -a $path/$chimp -b $path/$human_exons > $path/Chimp_Human_ExonsCircs_GTF-CE10.txt"

#echo "$ChimpInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_ChimpHumExonsCircs -q short.q -V


###LEMUR

LemurInt="bedtools intersect -loj -a $path/$lemur -b $path/$human_exons > $path/Lemur_Human_ExonsCircs_GTF-CE10.txt"

echo "$LemurInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_LemurHumExonsCircs -q short.q -V


####MACACA

#MacaInt="bedtools intersect -loj -a $path/$macaca -b $path/$human_exons > $path/Macaca_Human_ExonsCircs_GTF-CE10.txt"

#echo "$MacaInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_MacacaHumExonsCircs -q short.q -V


####MACAQUE

#MacqInt="bedtools intersect -loj -a $path/$macaque -b $path/$human_exons > $path/Macaque_Human_ExonsCircs_GTF-CE10.txt"

#echo "$MacqInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_MacaqueHumExonsCircs -q short.q -V


####MARMOSET

#MarmoInt="bedtools intersect -loj -a $path/$marmoset -b $path/$human_exons > $path/Marmoset_Human_ExonsCircs_GTF-CE10.txt"

#echo "$MarmoInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_MarmosetHumExonsCircs -q short.q -V


#####SQUIRREL MONKEY

#SquiMonkInt="bedtools intersect -loj -a $path/$squimonkey -b $path/$human_exons > $path/SquirrelMonkey_Human_ExonsCircs_GTF-CE10.txt"

#echo "$SquiMonkInt" | qsub -l mem_requested=20G,tmp_requested=20G -N BedtoolsInt_SquirrelMonkeyHumExonsCircs -q short.q -V
