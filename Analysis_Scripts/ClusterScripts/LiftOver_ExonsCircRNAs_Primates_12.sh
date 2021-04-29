# -S /bin/sh

#LiftOver exon coordinates that are within circRNAs coordinates of Primates to Human Coordinates

#path exons in circRNAs
path_2lift="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords/PSI_10"

baboon="Baboon_ExonCoords_4LiftOver_GTF-CE10.bed"
chimp="Chimp_ExonCoords_4LiftOver_GTF-CE10.bed"
lemur="Lemur_ExonCoords_4LiftOver_GTF-CE10.bed"
macaca="Macaca_ExonCoords_4LiftOver_GTF-CE10.bed"
macaque="Macaque_ExonCoords_4LiftOver_GTF-CE10.bed"
marmo="Marmoset_ExonCoords_4LiftOver_GTF-CE10.bed"
sqmonk="SquirrelMonkey_ExonCoords_4LiftOver_GTF-CE10.bed"

#liftover chains
path_liftover="/share/ScratchGeneral/gabrod/ChainFiles"

baboon_chn="papAnu4ToHg38.over.chain.gz"
chimp_chn="panTro5ToHg38.over.chain.gz"
lemur_chn="micMur2ToHg38.over.chain.gz"
macaca_chn="macFas5ToHg38.over.chain.gz"
macq_chn="rheMac8ToHg38.over.chain.gz"
marmo_chn="calJac3ToHg38.over.chain.gz"
sqmonk_chn="saiBol1ToHg38.over.chain.gz"


#########Commands for each primate

#BABOON
#baboonLO="liftOver $path_2lift/$baboon $path_liftover/$baboon_chn  $path_2lift/LiftOver_Baboon2Hum_ExonCirc_GTF-CE10.txt $path_2lift/unmapped_Baboon2Hum_ExonCirc_GTF-CE10"

#echo "$baboonLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveBaboonCirc10 -q short.q -V

#CHIMP
#chimpLO="liftOver $path_2lift/$chimp $path_liftover/$chimp_chn  $path_2lift/LiftOver_Chimp2Hum_ExonCirc_GTF-CE10.txt $path_2lift/unmapped_Chimp2Hum_ExonCirc_GTF-CE10"

#echo "$chimpLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveChimpCirc10 -q short.q -V



#LEMUR
lemurLO="liftOver $path_2lift/$lemur $path_liftover/$lemur_chn  $path_2lift/LiftOver_Lemur2Hum_ExonCirc_GTF-CE10.txt $path_2lift/unmapped_Lemur2Hum_ExonCirc_GTF-CE10"

echo "$lemurLO"	| qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveLemurCirc10 -q short.q -V


#MACACA
#macaLO="liftOver $path_2lift/$macaca $path_liftover/$macaca_chn  $path_2lift/LiftOver_Macaca2Hum_ExonCirc_GTF-CE10.txt $path_2lift/unmapped_Macaca2Hum_ExonCirc_GTF-CE10"

#echo "$macaLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveMacacaCirc10 -q short.q -V


#MACAQUE
#macqLO="liftOver $path_2lift/$macaque $path_liftover/$macq_chn  $path_2lift/LiftOver_Macaque2Hum_ExonCirc_GTF-CE10.txt $path_2lift/unmapped_Macaque2Hum_ExonCirc_GTF-CE10"

#echo "$macqLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveMacaqueCirc10 -q short.q -V


#MARMOSET
#marmoLO="liftOver $path_2lift/$marmo $path_liftover/$marmo_chn  $path_2lift/LiftOver_Marmoset2Hum_ExonCirc_GTF-CE10.txt $path_2lift/unmapped_Marmoset2Hum_ExonCirc_GTF-CE10"

#echo "$marmoLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveMarmosetCirc10 -q short.q -V

#SQUIRREL_MONKEY
#squimonkLO="liftOver $path_2lift/$sqmonk $path_liftover/$sqmonk_chn  $path_2lift/LiftOver_SquirrelMonkey2Hum_ExonCirc_GTF-CE10.txt $path_2lift/unmapped_SquirrelMonkey2Hum_ExonCirc_GTF-CE10"

#echo "$squimonkLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveSquirrelMonkeyCirc10 -q short.q -V

