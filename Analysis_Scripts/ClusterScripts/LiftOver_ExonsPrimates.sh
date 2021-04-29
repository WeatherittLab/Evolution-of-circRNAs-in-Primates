# -S /bin/sh

#LiftOver exon coordinates of Primates to Human Coordinates

#path exons
path_2lift="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords/PSI_10/"

#out_path
out="/share/ClusterShare/biodata/contrib/gabrod/Primates/WhippetOut_circ/ExonsCoords/PSI_10/LiftOver_AllExons/"

baboon="Baboon_Exons_GTF-CE10.filt.bed"
chimp="Chimp_Exons_GTF-CE10.filt.bed"
lemur="Lemur_Exons_GTF-CE10.filt.bed"
macaca="Macaca_Exons_GTF-CE10.filt.bed"
macaque="Macaque_Exons_GTF-CE10.filt.bed"
marmo="Marmoset_Exons_GTF-CE10.filt.bed"
sqmonk="SquiMonkey_Exons_GTF-CE10.filt.bed"

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
baboonLO="liftOver $path_2lift/$baboon $path_liftover/$baboon_chn  $out/LiftOver_Baboon2Hum_ExonGTF-CE10.txt $out/unmapped_Baboon2Hum_ExonGTF-CE10"

echo "$baboonLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveBaboonExons -q short.q -V

#CHIMP
chimpLO="liftOver $path_2lift/$chimp $path_liftover/$chimp_chn  $out/LiftOver_Chimp2Hum_ExonGTF-CE10.txt $out/unmapped_Chimp2Hum_ExonGTF-CE10"

echo "$chimpLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveChimpExons -q short.q -V


#LEMUR
lemurLO="liftOver $path_2lift/$lemur $path_liftover/$lemur_chn  $out/LiftOver_Lemur2Hum_ExonGTF-CE10.txt $out/unmapped_Lemur2Hum_ExonGTF-CE10"

echo "$lemurLO"	| qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveLemurExons -q short.q -V


#MACACA
macaLO="liftOver $path_2lift/$macaca $path_liftover/$macaca_chn  $out/LiftOver_Macaca2Hum_ExonGTF-CE10.txt $out/unmapped_Macaca2Hum_ExonGTF-CE10"

echo "$macaLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveMacacaExons -q short.q -V


#MACAQUE
macqLO="liftOver $path_2lift/$macaque $path_liftover/$macq_chn  $out/LiftOver_Macaque2Hum_ExonGTF-CE10.txt $out/unmapped_Macaque2Hum_ExonGTF-CE10"

echo "$macqLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveMacaqueExons -q short.q -V


#MARMOSET
marmoLO="liftOver $path_2lift/$marmo $path_liftover/$marmo_chn  $out/LiftOver_Marmoset2Hum_ExonGTF-CE10.txt $out/unmapped_Marmoset2Hum_ExonGTF-CE10"

echo "$marmoLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveMarmosetExons -q short.q -V


#SQUIRREL_MONKEY
squimonkLO="liftOver $path_2lift/$sqmonk $path_liftover/$sqmonk_chn  $out/LiftOver_SquirrelMonkey2Hum_ExonGTF-CE10.txt $out/unmapped_SquirrelMonkey2Hum_ExonGTF-CE10"

echo "$squimonkLO" | qsub -l mem_requested=20G,tmp_requested=20G -N LiftOveSquirrelMonkeyExons -q short.q -V

